import json
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
import fig2c_g2_estimator as est  # noqa: E402

CROSS_SCRIPT = SCRIPT_DIR / "compute_g2_cross_mesa.py"


def synthetic_table(path, frames=500, seed=7):
    """Small synthetic event table: independent Poisson events in 4 bands."""
    rng = np.random.default_rng(seed)
    tau_lo, tau_hi = est.TAU_RANGE_PX
    t0 = 155.0 + rng.normal(0.0, 3.0, frames)
    t0[3] = np.nan  # one t0 failure
    t0[37] = np.nan
    ev_frame = []
    ev_x = []
    ev_tau = []
    for band, (lo, hi) in enumerate(est.BAND_WINDOWS):
        counts = rng.poisson(5.0, frames)
        total = int(counts.sum())
        ev_frame.append(np.repeat(np.arange(frames), counts))
        ev_x.append(rng.uniform(lo, hi, total))
        ev_tau.append(rng.uniform(tau_lo, tau_hi, total))
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


class TestCrossMesaWiring(unittest.TestCase):
    def test_shared_b1_b4_coincidences_appear_only_in_lr(self):
        # A shared B1 x B4 component is a cross-mesa correlation; it must be
        # recovered by the LR pair and be absent from the RL pair.
        frames, bins, mu, q = 30000, 8, 0.3, 0.05
        rng = np.random.default_rng(97)
        shared = (rng.random((frames, bins)) < q).astype(np.float64)
        n_b1 = rng.poisson(mu, (frames, bins)) + shared
        n_b4 = rng.poisson(mu, (frames, bins)) + shared
        n_b2 = rng.poisson(mu, (frames, bins)).astype(np.float64)
        n_b3 = rng.poisson(mu, (frames, bins)).astype(np.float64)
        cube = np.stack([n_b1, n_b2, n_b3, n_b4], axis=1)
        beta = np.zeros((4, bins))
        true_g2 = (mu**2 + 2 * mu * q + q) / (mu + q) ** 2
        self.assertGreater(true_g2, 1.3)  # the test injects real correlation

        pooled_lr, sem_lr, _ = est.pooled_zero_delay(
            *est.pair_counts(cube, beta, "LR")
        )
        pooled_rl, sem_rl, _ = est.pooled_zero_delay(
            *est.pair_counts(cube, beta, "RL")
        )
        self.assertLess(abs(pooled_lr - true_g2), 4 * sem_lr)
        self.assertLess(abs(pooled_rl - 1.0), 4 * sem_rl)


class TestCrossMesaCli(unittest.TestCase):
    def test_end_to_end_poisson_is_unity(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            table_path = tmp / "table.h5"
            synthetic_table(table_path, frames=500, seed=19)
            out_dir = tmp / "out"
            command = [
                sys.executable,
                str(CROSS_SCRIPT),
                "--table",
                str(table_path),
                "--out-dir",
                str(out_dir),
            ]
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(0, result.returncode, result.stderr)
            with np.load(out_dir / "g2_cross_results_bin4px.npz") as results:
                for key in (
                    "t_ps",
                    "g2_LR",
                    "g2_RL",
                    "sem_plugin_LR",
                    "sem_plugin_RL",
                ):
                    self.assertIn(key, results)
            with open(out_dir / "g2_cross_summary.json") as handle:
                summary = json.load(handle)
            for pair in ("LR", "RL"):
                entry = summary["per_width"]["4"]["pairs"][pair]
                # Poisson table: pooled cross-mesa g2 consistent with 1.
                self.assertLess(
                    abs(entry["pooled_g2"] - 1.0),
                    4.0 * entry["pooled_sem_jackknife"],
                )


if __name__ == "__main__":
    unittest.main()
