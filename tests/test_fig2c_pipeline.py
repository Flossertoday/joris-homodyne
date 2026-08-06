import subprocess
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
import compute_g2_fig2c as cg2  # noqa: E402
import fig2c_g2_estimator as est  # noqa: E402

G2_SCRIPT = SCRIPT_DIR / "compute_g2_fig2c.py"


def synthetic_table(path, frames=500, seed=7):
    """Small synthetic event table: Poisson events in 4 bands + noise."""
    rng = np.random.default_rng(seed)
    tau_lo, tau_hi = est.TAU_RANGE_PX
    t0 = 155.0 + rng.normal(0.0, 3.0, frames)
    t0[3] = np.nan  # one t0 failure
    t0[37] = np.nan
    ev_frame = []
    ev_x = []
    ev_tau = []
    for band, (lo, hi) in enumerate(est.BAND_WINDOWS):
        # ~5 signal events/frame/band across the tau range
        counts = rng.poisson(5.0, frames)
        total = int(counts.sum())
        ev_frame.append(np.repeat(np.arange(frames), counts))
        ev_x.append(rng.uniform(lo, hi, total))
        ev_tau.append(rng.uniform(tau_lo, tau_hi, total))
    # residual noise: ~6 events/frame spread over the column range
    counts = rng.poisson(6.0, frames)
    total = int(counts.sum())
    ev_frame.append(np.repeat(np.arange(frames), counts))
    ev_x.append(rng.uniform(*est.EVENT_COL_RANGE, total))
    ev_tau.append(rng.uniform(tau_lo, tau_hi, total))
    frame = np.concatenate(ev_frame).astype(np.int32)
    x = np.concatenate(ev_x)
    tau = np.concatenate(ev_tau)
    with h5py.File(path, "w") as handle:
        group = handle.create_group("events")
        group.create_dataset("frame", data=frame)
        group.create_dataset("x", data=x)
        group.create_dataset("y", data=tau + np.nan_to_num(t0)[frame])
        group.create_dataset("tau", data=tau)
        group.create_dataset("amplitude", data=np.full(x.size, 40.0))
        handle.create_dataset("frame_t0", data=t0)
        handle.create_dataset(
            "failed_frames", data=np.flatnonzero(np.isnan(t0)).astype(np.int32)
        )
    return t0


class TestBuildCube(unittest.TestCase):
    def test_failed_frames_excluded_and_remapped(self):
        with tempfile.TemporaryDirectory() as tmp:
            table_path = Path(tmp) / "table.h5"
            t0 = synthetic_table(table_path, frames=200, seed=13)
            events, t0_loaded, _ = cg2.load_table(table_path)
            cube, nbins, kept = cg2.build_cube(events, t0_loaded, 4)
            self.assertEqual(198, kept)  # two t0 failures
            self.assertEqual((198, 4, nbins), cube.shape)
            self.assertEqual(cube.dtype, np.int16)
            # Every in-band event from kept frames is in the cube.
            band = est.assign_bands(events["x"])
            kept_events = np.isfinite(t0_loaded[events["frame"]]) & (band >= 0)
            self.assertEqual(int(cube.sum()), int(kept_events.sum()))


class TestComputeCli(unittest.TestCase):
    def test_end_to_end_poisson_is_unity(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            table_path = tmp / "table.h5"
            synthetic_table(table_path, frames=500, seed=17)
            out_dir = tmp / "out"
            command = [
                sys.executable,
                str(G2_SCRIPT),
                "--table",
                str(table_path),
                "--out-dir",
                str(out_dir),
            ]
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(0, result.returncode, result.stderr)
            with np.load(out_dir / "g2_results_bin4px.npz") as results:
                for key in ("t_ps", "g2_L", "g2_R", "sem_plugin_L", "sem_plugin_R"):
                    self.assertIn(key, results)
            cube = np.load(out_dir / "count_cube_bin4px.npy")
            self.assertEqual(np.int16, cube.dtype)
            self.assertEqual((500 - 2, 4), cube.shape[:2])
            import json

            with open(out_dir / "g2_summary.json") as handle:
                summary = json.load(handle)
            for pair in ("L", "R"):
                entry = summary["per_width"]["4"]["pairs"][pair]
                # Poisson table: pooled g2 consistent with 1.
                self.assertLess(
                    abs(entry["pooled_g2"] - 1.0),
                    4.0 * entry["pooled_sem_jackknife"],
                )


if __name__ == "__main__":
    unittest.main()
