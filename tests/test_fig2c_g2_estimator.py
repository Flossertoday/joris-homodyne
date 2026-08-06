import sys
import unittest
from pathlib import Path

import numpy as np

SCRIPT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(
    0,
    str(
        Path(__file__).resolve().parent.parent
        / "analysis"
        / "26-07-31-derived-single-sweep"
    ),
)
import compute_g2_matrix as g2mod  # noqa: E402
import fig2c_g2_estimator as est  # noqa: E402


def simulate_counts(frames, bins, mu, q, beta, seed):
    """N1, N2 with signal Poisson(mu) + shared Bernoulli(q) + background Poisson(beta)."""
    rng = np.random.default_rng(seed)
    shared = (rng.random((frames, bins)) < q).astype(np.float64)
    N1 = rng.poisson(mu, (frames, bins)) + shared
    N2 = rng.poisson(mu, (frames, bins)) + shared
    N1 = N1 + rng.poisson(beta, (frames, bins))
    N2 = N2 + rng.poisson(beta, (frames, bins))
    return N1.astype(np.float64), N2.astype(np.float64)


def analytic_g2(mu, q, beta):
    signal_mean = mu + q
    signal_product = mu**2 + 2 * mu * q + q
    true_g2 = signal_product / signal_mean**2
    measured_mean = signal_mean + beta
    measured_product = signal_product + 2 * signal_mean * beta + beta**2
    diluted_g2 = measured_product / measured_mean**2
    return true_g2, diluted_g2


class TestBandAssignmentAndBinning(unittest.TestCase):
    def test_assign_bands_edges_and_outside(self):
        x = np.array([349.0, 380.0, 415.0, 450.0, 500.0, 570.0, 680.0, 1010.0, 1110.0, 1200.0])
        band = est.assign_bands(x)
        np.testing.assert_array_equal(band, [-1, 0, 0, 0, -1, 1, 1, 3, 3, -1])

    def test_aligned_bin_index_out_of_range(self):
        index, nbins = est.aligned_bin_index(
            np.array([-64.0, -63.5, -0.1, 0.0, 1079.9, 1080.0, 2000.0, -100.0]), 2
        )
        self.assertEqual((1080 - (-64)) // 2, nbins)
        np.testing.assert_array_equal(index, [0, 0, 31, 32, 571, -1, -1, -1])

    def test_count_cube_matches_per_band_histogram(self):
        rng = np.random.default_rng(5)
        frame_count = 50
        # Synthetic event table: uniform in tau, only bands 0 and 2 populated.
        n_events = 4000
        frame = rng.integers(0, frame_count, n_events)
        tau = rng.uniform(-64.0, 1080.0, n_events)
        band = rng.choice([0, 2], n_events)
        cube, nbins = est.count_cube_aligned(frame, band, tau, frame_count, 2)
        self.assertEqual((frame_count, 4, nbins), cube.shape)
        self.assertEqual(int(cube.sum()), n_events)
        # Cube summed over frames equals the per-band tau histogram.
        for b in (0, 2):
            bin_index, _ = est.aligned_bin_index(tau[band == b], 2)
            hist = np.bincount(bin_index, minlength=nbins)
            np.testing.assert_array_equal(cube[:, b, :].sum(axis=0), hist)
        # Bands 1 and 3 are empty.
        self.assertEqual(0, int(cube[:, 1, :].sum()))
        self.assertEqual(0, int(cube[:, 3, :].sum()))


class TestBackgroundEstimate(unittest.TestCase):
    def test_beta_from_outside_band_events(self):
        rng = np.random.default_rng(29)
        frame_count = 200
        bin_width = 4
        _, nbins = est.aligned_bin_index(np.array([0.0]), bin_width)
        col_lo, col_hi = est.EVENT_COL_RANGE
        outside_width = (col_hi - col_lo) - float(
            np.sum(est.BAND_WINDOWS[:, 1] - est.BAND_WINDOWS[:, 0])
        )
        density = 0.02  # per frame per x-px per bin
        n_bg = int(density * frame_count * outside_width * nbins)
        x = rng.uniform(col_lo, col_hi, n_bg)
        outside = np.ones(n_bg, dtype=bool)
        for lo, hi in est.BAND_WINDOWS:
            outside &= ~((x >= lo) & (x <= hi))
        x = x[outside]
        tau = rng.uniform(-64.0, 1080.0, len(x))
        beta, rho, chi2_ndf = est.background_beta_aligned(
            tau, x, frame_count, bin_width
        )
        self.assertEqual((4, nbins), beta.shape)
        expected_rho = len(x) / (frame_count * outside_width * nbins)
        np.testing.assert_allclose(rho.mean(), expected_rho, rtol=1e-12)
        width0 = est.BAND_WINDOWS[0, 1] - est.BAND_WINDOWS[0, 0]
        np.testing.assert_allclose(beta[0].mean(), expected_rho * width0, rtol=1e-12)
        self.assertLess(chi2_ndf, 2.0)  # uniform synthetic background


class TestConsistencyWith0731Estimator(unittest.TestCase):
    def test_zero_delay_matches_g2_map_diagonal(self):
        frames, bins, mu, q, beta = 3000, 12, 0.4, 0.05, 0.1
        N1, N2 = simulate_counts(frames, bins, mu, q, beta, seed=41)
        betas = np.full(bins, float(beta))
        g2_map, sem_map = g2mod.g2_map(N1, N2, betas, betas)
        g2_diag, sem_diag = est.g2_zero_delay(N1, N2, betas, betas)
        np.testing.assert_allclose(g2_diag, np.diag(g2_map), rtol=1e-12)
        np.testing.assert_allclose(sem_diag, np.diag(sem_map), rtol=1e-12)

    def test_pooled_matches_tau_profile_zero(self):
        frames, bins, mu, q, beta = 3000, 12, 0.3, 0.05, 0.1
        N1, N2 = simulate_counts(frames, bins, mu, q, beta, seed=43)
        betas = np.full(bins, float(beta))
        tau, g2_tau, sem_tau, _ = g2mod.g2_tau_profile(N1, N2, betas, betas)
        zero = bins - 1
        pooled, sem, n_used = est.pooled_zero_delay(N1, N2, betas, betas)
        self.assertEqual(bins, n_used)
        self.assertAlmostEqual(float(g2_tau[zero]), pooled, places=12)
        # Same jackknife construction; agree to numerical precision.
        self.assertAlmostEqual(float(sem_tau[zero]), sem, places=10)


class TestEstimatorRecovery(unittest.TestCase):
    def test_poisson_zero_delay_is_unity(self):
        # Moderate occupancy: g2 is unbiased; the paper-style three-term
        # propagation is CONSERVATIVE there (the neglected covariances between
        # <n_i> and <n1*n2> reduce the true variance), so chi2/ndf < 1.
        frames, bins = 20000, 64
        N1, N2 = simulate_counts(frames, bins, mu=0.4, q=0.0, beta=0.0, seed=51)
        zeros = np.zeros(bins)
        g2, sem = est.g2_zero_delay(N1, N2, zeros, zeros)
        finite = np.isfinite(g2) & (sem > 0.0)
        self.assertGreater(finite.sum(), 0.9 * bins)
        chi2 = float(np.mean(((g2[finite] - 1.0) / sem[finite]) ** 2))
        self.assertLess(chi2, 1.0)
        self.assertGreater(chi2, 0.2)  # conservative but not absurdly so
        self.assertLess(abs(np.nanmean(g2) - 1.0), 0.02)  # unbiased

    def test_low_occupancy_sem_ensemble_matches_monte_carlo(self):
        # The 26-08-04 regime (mu ~ 0.004/bin): per-realization SEMs degenerate
        # (pinned in test_low_occupancy_degenerates_to_zero_sem), but the
        # ensemble rms of the propagated SEM must track the Monte-Carlo
        # scatter of g2 across replicas — that ensemble property is what the
        # G3 power gate relies on for per-bin statements.
        reps, frames, bins, mu = 40, 20000, 8, 0.004
        g2s = []
        sems = []
        for seed in range(reps):
            N1, N2 = simulate_counts(frames, bins, mu=mu, q=0.0, beta=0.0,
                                     seed=300 + seed)
            zeros = np.zeros(bins)
            g, s = est.g2_zero_delay(N1, N2, zeros, zeros)
            g2s.append(g)
            sems.append(s)
        g2s = np.array(g2s)
        sems = np.array(sems)
        mc_spread = g2s.std(axis=0, ddof=1)
        rms_sem = np.sqrt((sems**2).mean(axis=0))
        ratio = rms_sem / mc_spread
        self.assertTrue(np.all(ratio > 0.75), f"SEM underpredicts: {ratio}")
        self.assertTrue(np.all(ratio < 1.25), f"SEM overpredicts: {ratio}")

    def test_low_occupancy_pooled_jackknife_matches_monte_carlo(self):
        # Same regime: the window-pooled g2(0) jackknife SEM is honest per
        # realization (the pooled statistic is coincidence-count dominated).
        reps, frames, bins, mu = 40, 20000, 8, 0.004
        pooled = []
        sems = []
        for seed in range(reps):
            N1, N2 = simulate_counts(frames, bins, mu=mu, q=0.0, beta=0.0,
                                     seed=400 + seed)
            zeros = np.zeros(bins)
            p, s, _ = est.pooled_zero_delay(N1, N2, zeros, zeros)
            pooled.append(p)
            sems.append(s)
        spread = float(np.std(pooled, ddof=1))
        mean_sem = float(np.mean(sems))
        self.assertLess(abs(spread - mean_sem) / mean_sem, 0.15)

    def test_low_occupancy_degenerates_to_zero_sem(self):
        # 26-08-04 per-bin occupancy regime (mu ~ 0.004): bins that observe no
        # coincidence over all frames give g2 = 0 with propagated SEM = 0, and
        # the g2 ensemble scatter stays far above that degenerate SEM. Pin this
        # here so the power-gate analysis is forced to use Monte-Carlo/plug-in
        # errors for per-bin statements instead of the paper-style propagated
        # per-bin SEM (which was validated at mu ~ 0.4 on the 26-07-31 data).
        frames, bins = 20000, 64
        N1, N2 = simulate_counts(frames, bins, mu=0.004, q=0.0, beta=0.0, seed=61)
        zeros = np.zeros(bins)
        g2, sem = est.g2_zero_delay(N1, N2, zeros, zeros)
        coincidence = (N1 * N2).sum(axis=0)  # total observed coincidences per bin
        no_coinc = coincidence == 0
        self.assertGreater(no_coinc.sum(), 0.5 * bins)  # most bins, this regime
        np.testing.assert_array_equal(g2[no_coinc], 0.0)
        np.testing.assert_array_equal(sem[no_coinc], 0.0)
        # Ensemble scatter across bins is order-unity: per-bin g2 is pure noise.
        self.assertGreater(np.nanstd(g2), 0.5)

    def test_bunched_pooled_recovery(self):
        frames, bins, mu, q = 30000, 8, 0.3, 0.05
        true_g2, _ = analytic_g2(mu, q, beta=0.0)
        N1, N2 = simulate_counts(frames, bins, mu, q, beta=0.0, seed=53)
        zeros = np.zeros(bins)
        pooled, sem, _ = est.pooled_zero_delay(N1, N2, zeros, zeros)
        self.assertGreater(true_g2, 1.3)  # the test injects real bunching
        self.assertLess(abs(pooled - true_g2), 4 * sem)

    def test_background_correction_restores_true_g2(self):
        frames, bins, mu, q, beta = 30000, 8, 0.3, 0.05, 0.1
        true_g2, diluted_g2 = analytic_g2(mu, q, beta)
        N1, N2 = simulate_counts(frames, bins, mu, q, beta, seed=57)
        betas = np.full(bins, float(beta))
        zeros = np.zeros(bins)
        pooled_corr, sem_corr, _ = est.pooled_zero_delay(N1, N2, betas, betas)
        pooled_raw, _, _ = est.pooled_zero_delay(N1, N2, zeros, zeros)
        self.assertLess(abs(pooled_corr - true_g2), 4 * sem_corr)
        self.assertLess(abs(pooled_raw - diluted_g2), 0.02)
        self.assertLess(diluted_g2, true_g2)  # background dilutes toward 1

    def test_jackknife_sem_matches_monte_carlo(self):
        reps, frames, bins = 40, 1500, 8
        estimates = []
        sems = []
        diag_estimates = []
        diag_sems = []
        for seed in range(reps):
            N1, N2 = simulate_counts(frames, bins, mu=0.5, q=0.05, beta=0.0,
                                     seed=2000 + seed)
            zeros = np.zeros(bins)
            pooled, sem, _ = est.pooled_zero_delay(N1, N2, zeros, zeros)
            estimates.append(pooled)
            sems.append(sem)
            g2, dsem = est.g2_zero_delay(N1, N2, zeros, zeros)
            diag_estimates.append(g2[3])
            diag_sems.append(dsem[3])
        spread = np.std(estimates, ddof=1)
        mean_sem = np.mean(sems)
        self.assertLess(spread, 1.8 * mean_sem)
        self.assertGreater(spread, 0.5 * mean_sem)
        diag_spread = np.std(diag_estimates, ddof=1)
        diag_mean_sem = np.mean(diag_sems)
        self.assertLess(diag_spread, 1.8 * diag_mean_sem)
        self.assertGreater(diag_spread, 0.5 * diag_mean_sem)


class TestPairCounts(unittest.TestCase):
    def test_same_mesa_pairs(self):
        cube = np.arange(2 * 4 * 3, dtype=np.int32).reshape(2, 4, 3)
        beta = np.full((4, 3), 0.1)
        N1, N2, b1, b2 = est.pair_counts(cube, beta, "L")
        np.testing.assert_array_equal(N1, cube[:, 0, :])
        np.testing.assert_array_equal(N2, cube[:, 2, :])
        N1, N2, b1, b2 = est.pair_counts(cube, beta, "R")
        np.testing.assert_array_equal(N1, cube[:, 1, :])
        np.testing.assert_array_equal(N2, cube[:, 3, :])


if __name__ == "__main__":
    unittest.main()
