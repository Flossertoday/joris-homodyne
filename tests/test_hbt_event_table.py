"""Synthetic-data tests for the merged event-table pass in extract_hbt_traces.

The event table is written by the same single pass that accumulates the
Fig. 2(a) traces (merge of the retired extract_hbt_event_table.py). These
tests pin the table schema, the refuse-overwrite/atomic write, the detection
column-range ROI, the per-frame CFD50 pulse calibration (tau = y - t0), and
the table-vs-npz self-consistency validation on a 6-frame synthetic dataset.
Also covers the per-frame sample-overlay frame sampling and rendering.
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

SCRIPT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))
import extract_hbt_traces as eht  # noqa: E402

FRAME_SHAPE = (1016, 1344)
BACKGROUND = 105.0  # flat streak background, counts/px (matches measured 104.9)
THRESHOLD = 45.0
N_FRAMES = 6
FAILED_FRAME = 3  # this frame gets no reference pulse -> t0 failure
BAND_WINDOWS = np.array([[490, 510], [890, 910]], dtype=np.float64)


def make_frame(seed, with_pulse=True, t0_row=150.0):
    """Synthetic raw frame: flat background + reference pulse + event blobs.

    Events (all well above THRESHOLD):
      A: 3x3 block, centroid ~(401, 501)   -- inside col range, band 0
      B: 2x2 block, centroid ~(700.5, 900.5) -- inside col range, band 1
      hot: single pixel (500, 600) -- rejected by min_size=2
      out: 3x3 block at columns 200-202 -- outside the detection col range
    """
    from scipy.special import erf

    rng = np.random.default_rng(seed)
    lam = np.full(FRAME_SHAPE, BACKGROUND)
    if with_pulse:
        rows = np.arange(FRAME_SHAPE[0], dtype=np.float64)
        profile = 0.5 * (1.0 + erf((rows - t0_row) / (2.0 * np.sqrt(2.0))))
        lam[:, 170:215] += (20000.0 * profile)[:, None] / 45.0
    lam[400:403, 500:503] += 200.0
    lam[700:702, 900:902] += 250.0
    lam[500, 600] += 250.0
    lam[100:103, 200:203] += 300.0
    return rng.poisson(lam).astype(np.int16)


def write_synthetic_h5(path):
    """6 frames; frame FAILED_FRAME has no reference pulse (t0 failure)."""
    background = np.full(FRAME_SHAPE, BACKGROUND)
    frames = np.stack(
        [
            make_frame(seed=100 + i, with_pulse=(i != FAILED_FRAME),
                       t0_row=150.0 + 0.25 * i)
            for i in range(N_FRAMES)
        ]
    )
    with h5py.File(path, "w") as handle:
        handle.create_dataset(eht.FRAMES_DATASET, data=frames)
        handle.create_dataset(eht.BACKGROUND_DATASET, data=background)
    return background


def run_synthetic_extraction(tmp):
    """One merged accumulate pass over the synthetic dataset."""
    h5_path = Path(tmp) / "synthetic.h5"
    background = write_synthetic_h5(h5_path)
    indices = np.arange(N_FRAMES, dtype=np.int64)
    acc = eht.accumulate(
        h5_path, background, indices, THRESHOLD, collect_events=True
    )
    frame_t0 = np.full(N_FRAMES, np.nan)
    frame_t0[acc["frame_indices"]] = acc["t0_rows"]
    failed = np.setdiff1d(indices, acc["frame_indices"]).astype(np.int32)
    return h5_path, acc, frame_t0, failed


class TestAccumulateCollectEvents(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.h5_path, cls.acc, cls.frame_t0, cls.failed = run_synthetic_extraction(
            cls.tmp
        )

    def test_failed_frame_excluded(self):
        np.testing.assert_array_equal(self.failed, [FAILED_FRAME])
        self.assertEqual(self.acc["frame_indices"].size, N_FRAMES - 1)
        self.assertTrue(np.isnan(self.frame_t0[FAILED_FRAME]))

    def test_two_events_per_kept_frame(self):
        events = self.acc["events"]
        frames, counts = np.unique(events["frame"], return_counts=True)
        np.testing.assert_array_equal(
            frames, [i for i in range(N_FRAMES) if i != FAILED_FRAME]
        )
        np.testing.assert_array_equal(counts, 2)

    def test_events_inside_detection_col_range(self):
        # 横轴 ROI: detection is restricted to EVENT_COL_RANGE; the 3x3 block
        # at columns 200-202 must never appear, nor the reference pulse.
        x = self.acc["events"]["x"]
        lo, hi = eht.EVENT_COL_RANGE
        self.assertTrue(np.all(x >= lo) and np.all(x < hi))

    def test_event_centroids_near_planted_blobs(self):
        events = self.acc["events"]
        first = events["frame"] == 0
        rows = np.sort(events["y"][first])
        cols = np.sort(events["x"][first])
        self.assertAlmostEqual(rows[0], 401.0, delta=1.0)
        self.assertAlmostEqual(rows[1], 700.5, delta=1.0)
        self.assertAlmostEqual(cols[0], 501.0, delta=1.0)
        self.assertAlmostEqual(cols[1], 900.5, delta=1.0)
        self.assertTrue(np.all(events["amplitude"] > 500.0))


class TestWriteEventsH5(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.h5_path, cls.acc, cls.frame_t0, cls.failed = run_synthetic_extraction(
            cls.tmp
        )
        cls.out = Path(cls.tmp) / "events.h5"
        attributes = {
            "source_h5": str(cls.h5_path),
            "threshold": THRESHOLD,
            "min_size": 2,
            "col_range": list(eht.EVENT_COL_RANGE),
            "t0_method": "cfd50",
            "cfd_fraction": eht.CFD_FRACTION,
            "t0_roi": list(eht.T0_ROI),
            "band_windows": BAND_WINDOWS,
            "n_frames": N_FRAMES,
            "n_frames_kept": int(cls.acc["frame_indices"].size),
            "n_events": int(cls.acc["events"]["frame"].size),
            "script": "extract_hbt_traces.py",
        }
        eht.write_events_h5(
            cls.acc["events"], cls.frame_t0, cls.failed, cls.out, attributes
        )

    def test_schema_and_dtypes(self):
        with h5py.File(self.out, "r") as handle:
            self.assertEqual(handle["events/frame"].dtype, np.int32)
            for name in ("x", "y", "tau", "amplitude"):
                self.assertEqual(handle["events"][name].dtype, np.float64)
            self.assertEqual(handle["frame_t0"].dtype, np.float64)
            self.assertEqual(handle["failed_frames"].dtype, np.int32)
            np.testing.assert_array_equal(handle["failed_frames"][:], [FAILED_FRAME])
            self.assertEqual(handle.attrs["t0_method"], "cfd50")
            np.testing.assert_array_equal(handle.attrs["col_range"], [350, 1150])
            np.testing.assert_array_equal(
                handle.attrs["band_windows"], BAND_WINDOWS
            )

    def test_pulse_calibrated_tau_consistent_with_frame_t0(self):
        # 脉冲校准: every event's tau must equal y - frame_t0[frame] (CFD50),
        # and the failed frame contributes no events.
        with h5py.File(self.out, "r") as handle:
            frame = handle["events/frame"][:]
            y = handle["events/y"][:]
            tau = handle["events/tau"][:]
            frame_t0 = handle["frame_t0"][:]
        self.assertFalse(np.any(frame == FAILED_FRAME))
        np.testing.assert_allclose(tau, y - frame_t0[frame], atol=1e-9)

    def test_refuses_overwrite_and_cleans_partial(self):
        with self.assertRaises(FileExistsError):
            eht.write_events_h5(
                self.acc["events"], self.frame_t0, self.failed, self.out, {}
            )
        self.assertFalse(self.out.with_name(self.out.name + ".partial").exists())


class TestValidateEventsAgainstNpz(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.h5_path, cls.acc, cls.frame_t0, cls.failed = run_synthetic_extraction(
            cls.tmp
        )
        mean_photon, sem_photon = eht.mean_sem(
            cls.acc["photon_sum"], cls.acc["photon_sumsq"], cls.acc["photon_count"]
        )
        band_traces, _ = eht.aggregate_band_traces(
            mean_photon, sem_photon, BAND_WINDOWS.tolist()
        )
        cls.npz = Path(cls.tmp) / "traces.npz"
        np.savez(
            cls.npz,
            time_ps=eht.TAU_GRID * eht.PS_PER_PX,
            band_windows=BAND_WINDOWS,
            band_traces_photon=band_traces,
        )
        cls.max_diff, cls.totals, cls.kept = eht.validate_events_against_npz(
            cls.acc["events"], cls.frame_t0, cls.npz
        )

    def test_same_pass_rebuild_is_exact(self):
        # Table and npz come from the same detection pass: no bin-edge ties.
        np.testing.assert_array_equal(self.max_diff, [0, 0])
        self.assertEqual(self.kept, N_FRAMES - 1)

    def test_band_event_totals(self):
        # Each kept frame contributes one blob A (band 0) and one blob B (band 1).
        self.assertEqual(self.totals, [N_FRAMES - 1, N_FRAMES - 1])


class TestSampleOverlays(unittest.TestCase):
    """Per-frame sample overlays: seeded random sampling, one PNG per frame."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.h5_path = Path(cls.tmp) / "synthetic.h5"
        cls.background = write_synthetic_h5(cls.h5_path)
        cls.frame_t0 = np.full(N_FRAMES, np.nan)
        for i in range(N_FRAMES):
            if i != FAILED_FRAME:
                cls.frame_t0[i] = 150.0 + 0.25 * i
        cls.events_h5 = Path(cls.tmp) / "detected-peaks" / "synthetic_events.h5"

    def test_pick_overlay_frames_seeded_excludes_failed(self):
        picked = eht.pick_overlay_frames(self.frame_t0, count=3, seed=7)
        self.assertEqual(picked.size, 3)
        self.assertNotIn(FAILED_FRAME, picked)
        np.testing.assert_array_equal(
            picked, eht.pick_overlay_frames(self.frame_t0, count=3, seed=7)
        )
        # count above the kept population is clipped to all kept frames.
        self.assertEqual(
            eht.pick_overlay_frames(self.frame_t0, count=99, seed=7).size,
            N_FRAMES - 1,
        )

    def test_one_full_frame_png_per_sampled_frame(self):
        sample = [0, 4]
        paths = eht.save_event_sample_overlays(
            self.h5_path, self.background, self.frame_t0, BAND_WINDOWS,
            sample, self.events_h5, THRESHOLD,
        )
        self.assertEqual(len(paths), len(sample))
        for path, index in zip(paths, sample):
            self.assertTrue(path.exists())
            self.assertEqual(
                path.name, f"synthetic_events_sample_overlay_frame{index:05d}.png"
            )
            self.assertGreater(path.stat().st_size, 10000)


class TestHelpers(unittest.TestCase):
    def test_default_events_h5(self):
        path = eht.default_events_h5("DATA/x/26-08-04_hbt_time_reference.h5")
        self.assertEqual(path.parent.name, "detected-peaks")
        self.assertEqual(path.name, "26-08-04_hbt_time_reference_events.h5")

    def test_assign_band_indices(self):
        x = np.array([349.9, 495.0, 900.0, 1150.0, 700.0])
        np.testing.assert_array_equal(
            eht.assign_band_indices(x, BAND_WINDOWS), [-1, 0, 1, -1, -1]
        )


if __name__ == "__main__":
    unittest.main()
