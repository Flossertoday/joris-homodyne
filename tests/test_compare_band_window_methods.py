import sys
import unittest
from pathlib import Path

import numpy as np

SCRIPT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-08-band-window-methods"
)
HBT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(HBT_DIR))
import compare_band_window_methods as cbwm  # noqa: E402
import extract_hbt_traces as eht  # noqa: E402

MU = [415.0, 625.0, 845.0, 1060.0]
MU_SUPP = [395.0, 610.0, 815.0, 1035.0]
SIGMA_BAND = 36.0


def make_hist(floor_map, amps, mus=MU, sigma=SIGMA_BAND, seed=None, truncate=True):
    """Noiseless or Poisson-sampled event column histogram: local floor map
    plus four Gaussian bands, zero outside the detection column range.
    """
    cols = np.arange(cbwm.N_COLS)
    profile = np.asarray(floor_map, dtype=float).copy()
    for a, mu in zip(amps, mus):
        profile = profile + a * np.exp(-0.5 * ((cols - mu) / sigma) ** 2)
    if truncate:
        profile[: cbwm.DET_RANGE[0]] = 0.0
        profile[cbwm.DET_RANGE[1]:] = 0.0
    if seed is None:
        return profile
    rng = np.random.default_rng(seed)
    return rng.poisson(np.clip(profile, 0, None)).astype(float)


def flat_floor(level=200.0):
    return np.full(cbwm.N_COLS, level)


def shelf_floor():
    """Engineered supplement-like floor: 400 events/col baseline with a 620
    shelf under the first band (cols 350-470) — the structured out-of-band
    glow that breaks global-threshold methods, with B1 additionally truncated
    by the detection edge.
    """
    floor = np.full(cbwm.N_COLS, 400.0)
    floor[350:471] = 620.0
    return floor


def smooth_structured_floor():
    """Smooth +-20% sinusoidal structure (floor-fidelity check)."""
    cols = np.arange(cbwm.N_COLS)
    return (
        500.0
        + 55.0 * np.sin((cols - 300.0) / 180.0)
        + 35.0 * np.cos((cols - 500.0) / 90.0)
    )


class TestIncumbentReimplementation(unittest.TestCase):
    def test_matches_extract_hbt_traces_detect_bands(self):
        hist = make_hist(smooth_structured_floor(), [400, 600, 400, 550], seed=7)
        mine = cbwm.method_subwindow(hist, k=5.0)["windows"]
        ref = eht.detect_bands(
            hist, np.zeros(cbwm.N_COLS), band_sigma=5.0
        )["band_windows"]
        self.assertEqual(mine, [tuple(w) for w in ref])

    def test_flat_baseline_recovers_four_bands(self):
        hist = make_hist(flat_floor(), [400, 600, 400, 550], seed=3)
        windows = cbwm.method_subwindow(hist, k=5.0)["windows"]
        self.assertEqual(len(windows), 4)


class TestMethodsOnFlatFloor(unittest.TestCase):
    def setUp(self):
        self.hist = make_hist(flat_floor(), [400, 600, 400, 550], seed=11)
        self.floor = cbwm.local_floor(self.hist)

    def _check_centers(self, windows, mus=MU, tol=25.0):
        self.assertEqual(len(windows), 4)
        centers = [0.5 * (x0 + x1 - 1) for (x0, x1) in windows]
        for c, mu in zip(centers, mus):
            self.assertLess(abs(c - mu), tol)

    def _check_range(self, windows):
        for (x0, x1) in windows:
            self.assertGreaterEqual(x0, cbwm.DET_RANGE[0])
            self.assertLessEqual(x1, cbwm.DET_RANGE[1])

    def test_fdr_snr(self):
        w = cbwm.method_fdr_snr(self.hist, self.floor)["windows"]
        self._check_centers(w)
        self._check_range(w)

    def test_snr_aperture(self):
        w = cbwm.method_snr_aperture(self.hist, self.floor)["windows"]
        self._check_centers(w)
        self._check_range(w)

    def test_gauss4(self):
        # mu +/- k sigma windows are not range-limited by design
        self._check_centers(
            cbwm.method_gauss4(self.hist, MU, rule="mu2sigma")["windows"])

    def test_otsu(self):
        self._check_centers(cbwm.method_otsu(self.hist)["windows"])

    def test_findpeaks(self):
        self._check_centers(cbwm.method_findpeaks(self.hist, self.floor)["windows"])

    def test_matched(self):
        self._check_centers(cbwm.method_matched(self.hist, self.floor)["windows"])


class TestStructuredFloor(unittest.TestCase):
    """Supplement regime: shelf floor + weak edge-truncated first band."""

    def setUp(self):
        self.hist = make_hist(
            shelf_floor(), [550, 1000, 550, 800], mus=MU_SUPP, seed=13,
        )
        self.floor = cbwm.local_floor(self.hist)

    def _check_four(self, windows):
        self.assertEqual(len(windows), 4)
        centers = [0.5 * (x0 + x1 - 1) for (x0, x1) in windows]
        for c, mu in zip(centers, MU_SUPP):
            self.assertLess(abs(c - mu), 30.0)

    def test_fdr_snr_finds_all_four(self):
        windows = cbwm.method_fdr_snr(self.hist, self.floor)["windows"]
        self._check_four(windows)
        # windows stay inside the detection range and exclude hot edge columns
        for (x0, x1) in windows:
            self.assertGreaterEqual(x0, cbwm.DET_RANGE[0] + 2)
            self.assertLessEqual(x1, cbwm.DET_RANGE[1] - 4)

    def test_snr_aperture_finds_all_four(self):
        self._check_four(cbwm.method_snr_aperture(self.hist, self.floor)["windows"])

    def test_incumbent_fails_as_documented(self):
        # B1's thresholded run touches the detection edge and is discarded by
        # the incumbent's truncation rejection: the documented supplement
        # failure mode (three partial peak-top bands)
        windows = cbwm.method_subwindow(self.hist, k=5.0)["windows"]
        self.assertEqual(len(windows), 3)
        centers = [0.5 * (x0 + x1 - 1) for (x0, x1) in windows]
        for c in centers:
            self.assertGreater(c, 500.0)  # B1 (center 395) is missing

    def test_local_floor_tracks_smooth_structure(self):
        hist = make_hist(smooth_structured_floor(), [400, 600, 400, 550], seed=5)
        floor = cbwm.local_floor(hist)
        true = smooth_structured_floor()
        for (x0, x1) in [(660, 730), (900, 950)]:
            err = (floor[x0:x1] - true[x0:x1]).mean()
            self.assertLess(abs(err), 0.1 * true[x0:x1].mean())


class TestSnrOptimalWindow(unittest.TestCase):
    def test_realistic_contrast_optimum_near_1p6_sigma(self):
        # realistic per-column contrast (peak ~1.5x floor): the S/N-optimal
        # half-width of a Gaussian band lands near the classic ~1.6 sigma
        cols = np.arange(cbwm.N_COLS)
        floor = flat_floor(200.0)
        hist = floor + 100.0 * np.exp(-0.5 * ((cols - 625.0) / SIGMA_BAND) ** 2)
        x0, x1 = cbwm._snr_optimal_window(hist, floor, 625.0)
        half = 0.5 * (x1 - x0 - 1)
        self.assertGreater(half, 1.2 * SIGMA_BAND)
        self.assertLess(half, 2.2 * SIGMA_BAND)

    def test_pure_noise_fdr_detects_nothing(self):
        hist = make_hist(flat_floor(), [], seed=17)
        floor = cbwm.local_floor(hist)
        windows = cbwm.method_fdr_snr(hist, floor, q=0.01)["windows"]
        self.assertEqual(windows, [])


class TestEmissionSupport(unittest.TestCase):
    """Morphological emission-support walk (render_roi_morphology)."""

    def setUp(self):
        self.rwm = sys.modules.get("render_roi_morphology")
        if self.rwm is None:
            import render_roi_morphology as rwm
            self.rwm = rwm

    def test_gaussian_support_recovers_expected_extent(self):
        # threshold = max(3*sqrt(200)~42, 1.5% of 400 = 6) -> 42; for a
        # Gaussian with amp 400 and sigma 36 the walk ends at
        # 36*sqrt(2*ln(400/42)) ~ 75 px from the center
        hist = make_hist(flat_floor(), [400.0], mus=[625.0], truncate=False)
        floor = flat_floor()
        from scipy import ndimage

        res = ndimage.gaussian_filter1d(hist - floor, 3.0)
        x0, x1, tl, tr, thr = self.rwm.emission_support(res, floor, 625.0)
        self.assertFalse(tl)
        self.assertFalse(tr)
        self.assertLess(abs((625 - x0) - 75), 10)
        self.assertLess(abs((x1 - 625) - 75), 10)

    def test_support_truncation_flag_at_detection_edge(self):
        # B1-like band whose flank reaches the detection edge: the walk is
        # capped and flagged
        hist = make_hist(flat_floor(), [550.0], mus=[MU_SUPP[0]], seed=19)
        floor = cbwm.local_floor(hist)
        from scipy import ndimage

        res = ndimage.gaussian_filter1d(hist - floor, 3.0)
        x0, x1, tl, tr, thr = self.rwm.emission_support(
            res, floor, MU_SUPP[0], det_edges=cbwm.DET_RANGE)
        self.assertTrue(tl)
        self.assertEqual(x0, cbwm.DET_RANGE[0])

    def test_metrics_coverage_and_floor_fraction(self):
        hist = make_hist(flat_floor(), [400.0], mus=[625.0], truncate=False)
        floor = flat_floor()
        support = (550, 701, False, False, 42.0)
        narrow = self.rwm.method_metrics(hist, floor, [(590, 661)], [support])
        wide = self.rwm.method_metrics(hist, floor, [(545, 706)], [support])
        self.assertGreater(wide["coverage"], narrow["coverage"])
        self.assertGreater(wide["floor_fraction"], narrow["floor_fraction"])
        self.assertLessEqual(wide["coverage"], 1.0)
        # narrow window is fully inside the support: coverage < 1
        self.assertLess(narrow["coverage"], 0.9)


if __name__ == "__main__":
    unittest.main()
