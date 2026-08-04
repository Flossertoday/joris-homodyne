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
    / "26-07-31-derived-single-sweep"
)
sys.path.insert(0, str(SCRIPT_DIR))
import compute_g2_matrix as g2mod  # noqa: E402

G2_SCRIPT = SCRIPT_DIR / "compute_g2_matrix.py"


def simulate_counts(frames, bins, mu, q, beta, seed):
    """N1, N2 with signal Poisson(mu) + shared Bernoulli(q) + background Poisson(beta).

    E[S1] = mu + q; E[S1 S2] = mu^2 + 2 mu q + q (signal moments).
    """
    rng = np.random.default_rng(seed)
    shared = (rng.random((frames, bins)) < q).astype(np.float64)
    N1 = rng.poisson(mu, (frames, bins)) + shared
    N2 = rng.poisson(mu, (frames, bins)) + shared
    N1 = N1 + rng.poisson(beta, (frames, bins))
    N2 = N2 + rng.poisson(beta, (frames, bins))
    return N1.astype(np.float64), N2.astype(np.float64)


def analytic_g2(mu, q, beta):
    signal_mean = mu + q
    signal_product = mu ** 2 + 2 * mu * q + q
    true_g2 = signal_product / signal_mean ** 2
    measured_mean = signal_mean + beta
    measured_product = signal_product + 2 * signal_mean * beta + beta ** 2
    diluted_g2 = measured_product / measured_mean ** 2
    return true_g2, diluted_g2


class TestAssignmentAndBinning(unittest.TestCase):
    def test_assign_classes_edges_and_outside(self):
        intervals = np.array([[10.0, 20.0], [30.0, 40.0]])
        x = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 45.0])
        cls = g2mod.assign_classes(x, intervals)
        np.testing.assert_array_equal(cls, [-1, 0, 0, 0, -1, 1, 1, -1])

    def test_time_bin_index_covers_height_and_clips(self):
        index, nbins = g2mod.time_bin_index(
            np.array([0.0, 3.9, 4.0, 31.9, 31.99]), 32, 4
        )
        self.assertEqual(8, nbins)
        np.testing.assert_array_equal(index, [0, 0, 1, 7, 7])
        # Non-divisible height: ceil coverage, last bin partial.
        index, nbins = g2mod.time_bin_index(np.array([33.9]), 34, 4)
        self.assertEqual(9, nbins)
        np.testing.assert_array_equal(index, [8])

    def test_count_cube_sums_per_frame_class_bin(self):
        frame = np.array([0, 0, 1, 1, 1])
        cls = np.array([0, 0, 2, 3, -1])
        bins = np.array([1, 1, 1, 2, 0])
        cube = g2mod.count_cube(frame, cls, bins, nbins=3, frame_count=2)
        self.assertEqual((2, 4, 3), cube.shape)
        self.assertEqual(2.0, cube[0, 0, 1])
        self.assertEqual(1.0, cube[1, 2, 1])
        self.assertEqual(1.0, cube[1, 3, 2])
        self.assertEqual(0.0, cube.sum() - 4.0)


class TestCrossCorrelation(unittest.TestCase):
    def test_cross_corr_all_matches_direct(self):
        rng = np.random.default_rng(7)
        a = rng.normal(size=(3, 8))
        b = rng.normal(size=(3, 8))
        corr = g2mod._cross_corr_all(a, b)
        self.assertEqual((3, 15), corr.shape)
        taus = np.arange(-7, 8)
        for f in range(3):
            for k, tau in enumerate(taus):
                direct = sum(
                    a[f, n + tau] * b[f, n]
                    for n in range(8)
                    if 0 <= n + tau < 8
                )
                self.assertAlmostEqual(direct, corr[f, k], places=10)


class TestEstimatorRecovery(unittest.TestCase):
    def test_poisson_map_is_unity(self):
        frames, bins = 20000, 8
        N1, N2 = simulate_counts(frames, bins, mu=0.4, q=0.0, beta=0.0, seed=11)
        g2, sem = g2mod.g2_map(N1, N2, np.zeros(bins), np.zeros(bins))
        self.assertTrue(np.all(np.isfinite(g2)))
        self.assertLess(abs(np.median(g2) - 1.0), 0.02)

    def test_bunched_pooled_zero_delay_recovery(self):
        frames, bins, mu, q = 30000, 8, 0.3, 0.05
        true_g2, _ = analytic_g2(mu, q, beta=0.0)
        N1, N2 = simulate_counts(frames, bins, mu, q, beta=0.0, seed=13)
        tau, g2_tau, sem_tau, contrib = g2mod.g2_tau_profile(
            N1, N2, np.zeros(bins), np.zeros(bins)
        )
        zero = bins - 1
        self.assertEqual(0, tau[zero])
        np.testing.assert_array_equal(contrib, bins - np.abs(tau))
        self.assertLess(abs(g2_tau[zero] - true_g2), 4 * sem_tau[zero])
        self.assertGreater(true_g2, 1.3)  # the test injects real bunching

    def test_background_correction_restores_true_g2(self):
        frames, bins, mu, q, beta = 30000, 8, 0.3, 0.05, 0.1
        true_g2, diluted_g2 = analytic_g2(mu, q, beta)
        N1, N2 = simulate_counts(frames, bins, mu, q, beta, seed=17)
        betas = np.full(bins, float(beta))

        _, g2_corrected, sem_corrected, _ = g2mod.g2_tau_profile(N1, N2, betas, betas)
        _, g2_raw, _, _ = g2mod.g2_tau_profile(
            N1, N2, np.zeros(bins), np.zeros(bins)
        )
        zero = bins - 1
        self.assertLess(abs(g2_corrected[zero] - true_g2), 4 * sem_corrected[zero])
        self.assertLess(abs(g2_raw[zero] - diluted_g2), 0.02)
        self.assertLess(diluted_g2, true_g2)  # background dilutes toward 1

    def test_map_diagonal_matches_tau_zero(self):
        frames, bins, mu, q = 20000, 8, 0.3, 0.05
        N1, N2 = simulate_counts(frames, bins, mu, q, beta=0.0, seed=19)
        betas = np.zeros(bins)
        g2, _ = g2mod.g2_map(N1, N2, betas, betas)
        _, g2_tau, _, _ = g2mod.g2_tau_profile(N1, N2, betas, betas)
        # Pooled tau=0 is the ratio of diagonal sums of the map moments.
        m1 = N1.mean(axis=0)
        m2 = N2.mean(axis=0)
        product = (N1.T @ N2) / frames
        pooled = product.diagonal().sum() / (m1 * m2).sum()
        self.assertAlmostEqual(pooled, g2_tau[bins - 1], places=12)
        self.assertEqual((bins, bins), g2.shape)

    def test_tau_profile_arm_swap_reverses(self):
        frames, bins = 5000, 8
        N1, N2 = simulate_counts(frames, bins, mu=0.4, q=0.05, beta=0.0, seed=23)
        _, forward, _, _ = g2mod.g2_tau_profile(N1, N2, np.zeros(bins), np.zeros(bins))
        _, swapped, _, _ = g2mod.g2_tau_profile(N2, N1, np.zeros(bins), np.zeros(bins))
        np.testing.assert_allclose(forward, swapped[::-1], rtol=1e-10)


class TestErrorModel(unittest.TestCase):
    def test_jackknife_sem_matches_monte_carlo(self):
        reps, frames, bins = 40, 1500, 8
        estimates = []
        sems = []
        map_diag_estimates = []
        map_diag_sems = []
        for seed in range(reps):
            N1, N2 = simulate_counts(frames, bins, mu=0.5, q=0.05, beta=0.0,
                                     seed=1000 + seed)
            zeros = np.zeros(bins)
            _, g2_tau, sem_tau, _ = g2mod.g2_tau_profile(N1, N2, zeros, zeros)
            estimates.append(g2_tau[bins - 1])
            sems.append(sem_tau[bins - 1])
            g2, sem = g2mod.g2_map(N1, N2, zeros, zeros)
            map_diag_estimates.append(g2[3, 3])
            map_diag_sems.append(sem[3, 3])
        spread = np.std(estimates, ddof=1)
        mean_sem = np.mean(sems)
        self.assertLess(spread, 1.8 * mean_sem)
        self.assertGreater(spread, 0.5 * mean_sem)
        map_spread = np.std(map_diag_estimates, ddof=1)
        map_mean_sem = np.mean(map_diag_sems)
        self.assertLess(map_spread, 1.8 * map_mean_sem)
        self.assertGreater(map_spread, 0.5 * map_mean_sem)


class TestBackgroundEstimate(unittest.TestCase):
    def test_background_beta_counts_and_flatness(self):
        rng = np.random.default_rng(29)
        frame_count, height, width = 100, 64, 200
        intervals = np.array([[10.0, 20.0], [30.0, 40.0], [50.0, 60.0], [70.0, 80.0]])
        outside_width = width - 40.0
        density = 0.001  # per frame per x-px per y-bin(1px)
        n_bg = int(density * frame_count * outside_width * height)
        full_x = rng.uniform(0, width, n_bg)
        outside = np.ones(n_bg, dtype=bool)
        for lo, hi in intervals:
            outside &= ~((full_x >= lo) & (full_x <= hi))
        full_x = full_x[outside]
        full_y = rng.uniform(0, height, len(full_x))
        beta, rho_y, chi2_ndf = g2mod.background_beta(
            full_x, full_y, intervals, width, frame_count, 1, height, height
        )
        self.assertEqual((4, height), beta.shape)
        expected_rho = len(full_x) / (frame_count * outside_width * height)
        np.testing.assert_allclose(rho_y.mean(), expected_rho, rtol=1e-12)
        np.testing.assert_allclose(beta[0].mean(), expected_rho * 10.0, rtol=1e-12)
        self.assertLess(chi2_ndf, 2.0)  # uniform synthetic background


def write_peaks(path, frame, y, x, attributes):
    with h5py.File(path, "w") as handle:
        handle.create_dataset("frame", data=np.asarray(frame, dtype=np.int32))
        handle.create_dataset("y", data=np.asarray(y, dtype=np.float64))
        handle.create_dataset("x", data=np.asarray(x, dtype=np.float64))
        handle.create_dataset("mass", data=np.ones(len(x)))
        for key, value in attributes.items():
            handle.attrs[key] = value


class TestCli(unittest.TestCase):
    def test_cli_end_to_end_and_refuses_overwrite(self):
        rng = np.random.default_rng(31)
        frame_count, height, width = 150, 32, 120
        intervals = np.array(
            [[10.0, 20.0], [30.0, 40.0], [50.0, 60.0], [70.0, 80.0]]
        )
        frames_list = []
        ys = []
        xs = []
        for index, (lo, hi) in enumerate(intervals):
            counts = rng.poisson(1.2, frame_count)
            total = int(counts.sum())
            frames_list.append(np.repeat(np.arange(frame_count), counts))
            ys.append(rng.uniform(0, height, total))
            xs.append(rng.uniform(lo, hi, total))
        signal_frame = np.concatenate(frames_list)
        signal_y = np.concatenate(ys)
        signal_x = np.concatenate(xs)
        n_bg = 400
        bg_frame = rng.integers(0, frame_count, n_bg)
        bg_y = rng.uniform(0, height, n_bg)
        bg_x = rng.uniform(0, width, n_bg)

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            filtered_path = tmp / "filtered.h5"
            full_path = tmp / "full.h5"
            write_peaks(
                filtered_path, signal_frame, signal_y, signal_x,
                {
                    "condition_name": "synthetic",
                    "central_dot_strength": 0.18,
                    "x_range_intervals": intervals,
                },
            )
            write_peaks(
                full_path,
                np.concatenate([signal_frame, bg_frame]),
                np.concatenate([signal_y, bg_y]),
                np.concatenate([signal_x, bg_x]),
                {
                    "condition_name": "synthetic",
                    "central_dot_strength": 0.18,
                    "frame_shape": np.array([height, width]),
                    "frame_count": frame_count,
                },
            )
            output_path = tmp / "g2.h5"
            figure_dir = tmp / "figures"
            command = [
                sys.executable,
                str(G2_SCRIPT),
                str(filtered_path),
                str(full_path),
                "--output-h5",
                str(output_path),
                "--figure-dir",
                str(figure_dir),
                "--time-bins-px",
                "2",
                "4",
            ]
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(0, result.returncode, result.stderr)

            with h5py.File(output_path, "r") as handle:
                self.assertEqual("synthetic", handle.attrs["condition_name"])
                self.assertTrue(bool(handle.attrs["background_correction"]))
                np.testing.assert_array_equal(
                    handle.attrs["time_bins_px"], [2, 4]
                )
                group = handle["bin_px_2"]
                bins = height // 2
                self.assertEqual(bins, int(group.attrs["bin_count"]))
                self.assertEqual((bins, bins), handle["bin_px_2/g2_AA_map"].shape)
                self.assertEqual((bins,), handle["bin_px_2/g2_BB_zero_delay"].shape)
                self.assertEqual((2 * bins - 1,), handle["bin_px_2/g2_AB_tau"].shape)
                self.assertEqual((2 * bins - 1,), handle["bin_px_2/g2_BA_tau"].shape)
                np.testing.assert_allclose(
                    handle["bin_px_2/g2_BA_tau"][:],
                    handle["bin_px_2/g2_AB_tau"][:][::-1],
                )
                for pair in ("AA", "BB", "AB"):
                    self.assertIn(f"g2_{pair}_pooled_zero_delay", group.attrs)
                self.assertIn("bin_px_4", handle)
            figures = sorted(figure_dir.glob("*.png"))
            self.assertEqual(8, len(figures))  # 4 figure types x 2 bin widths

            rerun = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(0, rerun.returncode)
            self.assertIn("Refusing to overwrite", rerun.stderr)


if __name__ == "__main__":
    unittest.main()
