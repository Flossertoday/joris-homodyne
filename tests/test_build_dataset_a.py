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
import build_dataset_a as bda  # noqa: E402
import extract_hbt_traces as eht  # noqa: E402

FRAME_SHAPE = (1016, 1344)
BACKGROUND = 105.0


def noise_frame(seed):
    """Low-noise synthetic frame: flat background with mild Gaussian pixel
    noise (sigma 3 << threshold 25, so random above-threshold clusters are
    negligible and planted blobs dominate detection)."""
    rng = np.random.default_rng(seed)
    return rng.normal(0.0, 3.0, FRAME_SHAPE)


def add_blob(frame, row, col, amplitude=3000.0, sigma=2.5):
    r0, r1 = max(0, int(row) - 12), min(FRAME_SHAPE[0], int(row) + 13)
    c0, c1 = max(0, int(col) - 12), min(FRAME_SHAPE[1], int(col) + 13)
    rows, cols = np.mgrid[r0:r1, c0:c1]
    frame[r0:r1, c0:c1] += amplitude * np.exp(
        -((rows - row) ** 2 + (cols - col) ** 2) / (2.0 * sigma ** 2)
    )
    return frame


def add_reference_pulse(frame, row0=155.0):
    """Asymmetric reference pulse in the T0_ROI columns (CFD50 anchor)."""
    from scipy.special import erf

    rows = np.arange(FRAME_SHAPE[0], dtype=np.float64)
    profile = 0.5 * (1.0 + erf((rows - row0) / (2.0 * np.sqrt(2.0))))
    decay = np.exp(-np.maximum(rows - (row0 + 6.0), 0.0) / 25.0)
    profile = profile * np.where(rows > row0 + 6.0, decay, 1.0)
    frame[:, 170:215] += (20000.0 * profile)[:, None] / 45.0
    return frame


class DetectionFeaturesTest(unittest.TestCase):
    def test_centroids_and_amplitudes_match_production(self):
        for seed in (1, 2, 3):
            frame = noise_frame(seed)
            add_blob(frame, 300.4, 500.2)
            add_blob(frame, 700.1, 900.7, amplitude=800.0, sigma=1.5)
            feats = bda.detect_events_features(frame, 25.0)
            rows, cols, amps = eht.detect_events(frame, 25.0, with_counts=True)
            order_new = np.argsort(feats["rows"] * 1e4 + feats["cols"])
            order_ref = np.argsort(rows * 1e4 + cols)
            self.assertEqual(feats["rows"].size, rows.size)
            self.assertTrue(np.allclose(feats["rows"][order_new], rows[order_ref], atol=1e-12))
            self.assertTrue(np.allclose(feats["cols"][order_new], cols[order_ref], atol=1e-12))
            self.assertTrue(np.allclose(feats["amplitude"][order_new], amps[order_ref], atol=1e-9))

    def test_morphology_self_consistency(self):
        frame = noise_frame(11)
        add_blob(frame, 400.0, 600.0, amplitude=4000.0, sigma=3.0)
        feats = bda.detect_events_features(frame, 25.0)
        self.assertGreater(feats["rows"].size, 0)
        self.assertTrue(np.all(feats["size"] >= 2))
        self.assertTrue(np.all(feats["peak"] > 25.0))
        self.assertTrue(np.all(feats["amplitude"] >= feats["peak"]))
        bbox_area = (feats["r1"] - feats["r0"] + 1) * (feats["c1"] - feats["c0"] + 1)
        self.assertTrue(np.all(feats["size"] <= bbox_area))
        self.assertTrue(np.all(feats["r1"] >= feats["r0"]))
        self.assertTrue(np.all(feats["c1"] >= feats["c0"]))


class CurateTest(unittest.TestCase):
    def feats(self, x, r0=10, r1=20, c0=None, c1=None):
        x = np.asarray(x, dtype=np.float64)
        n = x.size
        return {
            "rows": np.full(n, 500.0),
            "cols": x,
            "amplitude": np.full(n, 100.0),
            "size": np.full(n, 5.0),
            "peak": np.full(n, 50.0),
            "r0": np.full(n, r0, dtype=np.int32),
            "r1": np.full(n, r1, dtype=np.int32),
            "c0": np.floor(x).astype(np.int32) - 3 if c0 is None else np.full(n, c0, dtype=np.int32),
            "c1": np.floor(x).astype(np.int32) + 3 if c1 is None else np.full(n, c1, dtype=np.int32),
        }

    def test_band_membership(self):
        keep, band, reason = bda.curate_events(self.feats([400.0, 500.0, 625.0]))
        self.assertEqual(keep.tolist(), [True, False, True])
        self.assertEqual(band.tolist(), [0, -1, 1])
        self.assertEqual(reason.tolist(), [0, 1, 0])

    def test_defect_column(self):
        keep, band, reason = bda.curate_events(self.feats([676.3, 677.0]))
        self.assertEqual(keep.tolist(), [False, True])
        self.assertEqual(reason.tolist(), [2, 0])

    def test_edge_touch(self):
        keep, _, reason = bda.curate_events(self.feats([400.0], r0=0, r1=5))
        self.assertEqual(keep.tolist(), [False])
        self.assertEqual(reason.tolist(), [3])
        keep, _, reason = bda.curate_events(self.feats([400.0], r0=1000, r1=1015))
        self.assertEqual(reason.tolist(), [3])

    def test_detection_boundary_touch(self):
        keep, _, reason = bda.curate_events(
            self.feats([400.0], c0=350, c1=360))
        self.assertEqual(reason.tolist(), [3])


class FloorFitTest(unittest.TestCase):
    def test_recovers_flat_floor_under_bump(self):
        cols = np.arange(*bda.EVENT_COL_RANGE)
        rate = np.full(cols.size, 0.010)
        for lo, hi in bda.DEFAULT_BANDS:
            center = 0.5 * (lo + hi)
            rate += 0.02 * np.exp(-((cols - center) ** 2) / (2 * 15.0 ** 2))
        rate_out = rate.copy()
        member = np.zeros(cols.size, dtype=bool)
        for lo, hi in bda.DEFAULT_BANDS:
            member |= (cols >= lo) & (cols <= hi)
        rate_out[member] = 0.0
        floor_mean, variants, report = bda.fit_band_floors(rate_out)
        for k in range(4):
            self.assertAlmostEqual(floor_mean[k], 0.010, places=4)
        # B1 has no usable left flank (detection edge) -> one-sided right
        self.assertEqual(report[0]["variants"][1]["one_sided"], "right")


class BuildRoundtripTest(unittest.TestCase):
    def make_h5(self, path, n_frames=6):
        rng = np.random.default_rng(5)
        frames = np.clip(
            BACKGROUND + rng.normal(0.0, 3.0, (n_frames,) + FRAME_SHAPE), 0, 65535
        ).astype(np.uint16)
        for i in range(n_frames):
            frame = frames[i].astype(np.float64) - BACKGROUND
            add_reference_pulse(frame, 155.0 + 0.1 * i)
            add_blob(frame, 300.0, 415.0)     # inside B1 -> kept
            add_blob(frame, 500.0, 500.0)     # out-of-band -> rejected
            add_blob(frame, 700.0, 676.2)     # defect column -> rejected
            frames[i] = np.clip(frame + BACKGROUND, 0, 65535).astype(np.uint16)
        with h5py.File(path, "w") as handle:
            handle.create_dataset(bda.FRAMES_DATASET, data=frames)
            handle.create_dataset(
                "conditions/background_and_streak_pulses_delay_calibration/background",
                data=np.full(FRAME_SHAPE, BACKGROUND),
            )
        return path

    def test_roundtrip_schema_and_cuts(self):
        with tempfile.TemporaryDirectory() as tmp:
            h5 = self.make_h5(Path(tmp) / "synthetic.h5")
            out = Path(tmp) / "dataset_a.h5"
            qc_dir = Path(tmp) / "qc"
            bda.build(h5, out, qc_dir, max_frames=None)
            with h5py.File(out, "r") as handle:
                events = handle["events"]
                for name in ("frame", "x", "y", "tau", "amplitude", "size",
                             "peak", "r0", "r1", "c0", "c1", "band"):
                    self.assertIn(name, events)
                x = events["x"][:]
                self.assertTrue(np.all((x >= 380) & (x <= 450)))  # only B1 blobs
                self.assertTrue(np.all(events["band"][:] == 0))
                t0 = handle["frame_t0"][:]
                self.assertTrue(np.all(np.isfinite(t0)))
                self.assertTrue(np.allclose(
                    events["tau"][:], events["y"][:] - t0[events["frame"][:]],
                    atol=1e-9))
                self.assertIn("noise_model", handle)
                self.assertEqual(handle["noise_model"]["floor_per_band"].shape, (4,))
                self.assertEqual(handle.attrs["n_events"], x.size)
            # refuse overwrite
            with self.assertRaises(FileExistsError):
                bda.build(h5, out, qc_dir)

    def test_t0_failure_frame_excluded(self):
        with tempfile.TemporaryDirectory() as tmp:
            h5 = self.make_h5(Path(tmp) / "synthetic.h5")
            with h5py.File(h5, "a") as handle:  # blank the pulse in frame 2
                frames = handle[bda.FRAMES_DATASET]
                frame = frames[2].astype(np.float64)
                frame[:, 170:215] = BACKGROUND
                frames[2] = frame.astype(np.uint16)
            out = Path(tmp) / "dataset_a.h5"
            bda.build(h5, out, Path(tmp) / "qc")
            with h5py.File(out, "r") as handle:
                self.assertEqual(handle["failed_frames"][:].tolist(), [2])
                self.assertTrue(np.isnan(handle["frame_t0"][2]))
                self.assertTrue(np.all(handle["events/frame"][:] != 2))


if __name__ == "__main__":
    unittest.main()
