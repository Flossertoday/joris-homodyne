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
import extract_hbt_traces as eht  # noqa: E402

FRAME_SHAPE = (1016, 1344)
BACKGROUND = 105.0  # flat streak background, counts/px (matches measured 104.9)
SPOT_COL = 192.5  # column centre of SPOT_ROI = (145, 190, 170, 215)


def background_frame(rng):
    """Background-subtracted pure-noise frame (Poisson around a flat level)."""
    return rng.poisson(BACKGROUND, FRAME_SHAPE).astype(np.float64) - BACKGROUND


def add_gaussian_spot(frame, row0, col0=SPOT_COL, amplitude=20000.0, sigma=3.0):
    """Add a 2D Gaussian reference spot sampled at pixel centres."""
    r_min = max(0, int(round(row0)) - 15)
    r_max = min(FRAME_SHAPE[0], int(round(row0)) + 16)
    c_min = max(0, int(round(col0)) - 15)
    c_max = min(FRAME_SHAPE[1], int(round(col0)) + 16)
    rows, cols = np.mgrid[r_min:r_max, c_min:c_max]
    signal = amplitude * np.exp(
        -((rows - row0) ** 2 + (cols - col0) ** 2) / (2.0 * sigma ** 2)
    )
    frame[r_min:r_max, c_min:c_max] += signal
    return frame


def make_spot_frame(row0, seed, col0=SPOT_COL, amplitude=20000.0, sigma=3.0):
    """Poisson-noise frame with a Gaussian spot; returned background-subtracted."""
    rng = np.random.default_rng(seed)
    lam = np.full(FRAME_SHAPE, BACKGROUND)
    lam = add_gaussian_spot(lam, row0, col0, amplitude, sigma)
    return rng.poisson(lam).astype(np.float64) - BACKGROUND


def make_band_frame(seed, band_x0=600, band_x1=630, row_center=400.0,
                    row_sigma=8.0, col_amplitude=200.0):
    """Band spanning columns [band_x0, band_x1) with a Gaussian row profile.

    Returns (frame_bg, expected_trace) where expected_trace is the noiseless
    column sum over [band_x0, band_x1) of the added signal.
    """
    rng = np.random.default_rng(seed)
    rows = np.arange(FRAME_SHAPE[0])
    profile = col_amplitude * np.exp(
        -((rows - row_center) ** 2) / (2.0 * row_sigma ** 2)
    )
    lam = np.full(FRAME_SHAPE, BACKGROUND)
    lam[:, band_x0:band_x1] += profile[:, None]
    frame_bg = rng.poisson(lam).astype(np.float64) - BACKGROUND
    expected = (band_x1 - band_x0) * profile
    return frame_bg, expected


def make_event_frame(seed, with_events=True):
    """Background-noise frame with compact bright event clusters.

    Events (all well above threshold 45, none near the boundary):
      A: 3x3 block, centroid (301.0, 501.0)
      B: 2x2 block, centroid (700.5, 900.5)
      hot: single pixel (500, 600) -- rejected by min_size=2
      out: 3x3 block at columns 200-202 -- outside default col_range
    """
    rng = np.random.default_rng(seed)
    lam = np.full(FRAME_SHAPE, BACKGROUND)
    if with_events:
        lam[300:303, 500:503] += 100.0
        lam[700:702, 900:902] += 150.0
        lam[500, 600] += 250.0
        lam[100:103, 200:203] += 100.0
    return rng.poisson(lam).astype(np.float64) - BACKGROUND


class TestEstimateT0(unittest.TestCase):
    def test_recovers_subpixel_spot_row(self):
        # Half-max centroid of a sigma=3 px Gaussian: discretization bias and
        # Poisson jitter at amplitude 20000 are both far below 0.5 px.
        for row0 in (160.3, 167.7):
            with self.subTest(row0=row0):
                frame = make_spot_frame(row0, seed=101)
                t0 = eht.estimate_t0(frame)
                self.assertTrue(np.isfinite(t0))
                self.assertLess(abs(t0 - row0), 0.5)

    def test_cosmic_ray_inside_roi_does_not_pull_centroid(self):
        # Compact cosmic-ray dot inside the ROI, dimmer than the spot peak:
        # the half-max mask excludes it, so the spot still dominates.
        row0 = 160.3
        frame = make_spot_frame(row0, seed=103)
        frame[175:177, 205:207] += 8000.0
        t0 = eht.estimate_t0(frame)
        self.assertTrue(np.isfinite(t0))
        self.assertLess(abs(t0 - row0), 1.0)

    def test_pure_background_returns_nan(self):
        frame = background_frame(np.random.default_rng(105))
        t0 = eht.estimate_t0(frame)
        self.assertTrue(np.isnan(t0))


class TestBandTrace(unittest.TestCase):
    def setUp(self):
        self.frame, self.expected = make_band_frame(seed=107)

    def test_recovers_row_profile_inside_band(self):
        trace = eht.band_trace(self.frame, 600, 630)
        self.assertEqual((FRAME_SHAPE[0],), trace.shape)
        self.assertEqual(np.float64, trace.dtype)
        # Noise on a 30-column sum at the band peak: sqrt(30*(105+200)) = 96.
        peak_row = 400
        self.assertLess(abs(trace[peak_row] - self.expected[peak_row]), 600.0)
        core = slice(384, 417)
        self.assertLess(
            np.max(np.abs(trace[core] - self.expected[core])), 700.0
        )
        # Far from the Gaussian row profile the trace is pure noise:
        # sigma = sqrt(30*105) = 56, 7-sigma bound below.
        far = np.concatenate([trace[:340], trace[460:]])
        self.assertLess(np.max(np.abs(far)), 400.0)

    def test_window_outside_band_is_zero(self):
        trace = eht.band_trace(self.frame, 700, 730)
        self.assertLess(np.max(np.abs(trace)), 400.0)

    def test_window_is_half_open(self):
        # [x0, x1): shifting the window by one column must drop exactly the
        # first band column; the difference is that column's own trace.
        full = eht.band_trace(self.frame, 600, 630)
        shifted = eht.band_trace(self.frame, 599, 629)
        diff = full - shifted
        # Column 600 signal at row 400 is 200 counts; its noise sigma is
        # sqrt(105+200) = 17.5, so 200 +- 105 is a 6-sigma band.
        self.assertLess(abs(diff[400] - 200.0), 105.0)
        far = np.concatenate([diff[:340], diff[460:]])
        self.assertLess(np.max(np.abs(far)), 200.0)


class TestDetectEvents(unittest.TestCase):
    def test_finds_two_clusters_at_weighted_centroids(self):
        frame = make_event_frame(seed=109)
        rows, cols = eht.detect_events(frame, threshold=45.0)
        self.assertEqual(2, len(rows))
        self.assertEqual(2, len(cols))
        order = np.argsort(rows)
        rows, cols = rows[order], cols[order]
        # Centroid noise on these ~1000-count clusters is <0.1 px.
        self.assertLess(abs(rows[0] - 301.0), 0.5)
        self.assertLess(abs(cols[0] - 501.0), 0.5)
        self.assertLess(abs(rows[1] - 700.5), 0.5)
        self.assertLess(abs(cols[1] - 900.5), 0.5)

    def test_rejects_single_pixel_hot_event(self):
        frame = make_event_frame(seed=111)
        rows, cols = eht.detect_events(frame, threshold=45.0, min_size=2)
        for r, c in zip(rows, cols):
            self.assertFalse(abs(r - 500.0) < 2.0 and abs(c - 600.0) < 2.0)

    def test_respects_col_range(self):
        frame = make_event_frame(seed=113)
        rows, cols = eht.detect_events(
            frame, threshold=45.0, col_range=(350, 1150)
        )
        self.assertTrue(np.all(cols >= 350.0))
        self.assertTrue(np.all(cols < 1150.0))
        # Widening the range recovers the third cluster near column 201.
        rows_wide, cols_wide = eht.detect_events(
            frame, threshold=45.0, col_range=(0, FRAME_SHAPE[1])
        )
        self.assertEqual(3, len(rows_wide))

    def test_background_only_returns_empty(self):
        frame = background_frame(np.random.default_rng(115))
        rows, cols = eht.detect_events(frame, threshold=45.0)
        self.assertEqual(0, len(rows))
        self.assertEqual(0, len(cols))


class TestInterpRows(unittest.TestCase):
    def setUp(self):
        rows = np.arange(FRAME_SHAPE[0], dtype=np.float64)
        self.rows = rows
        self.trace = 1000.0 * np.exp(-((rows - 400.0) ** 2) / (2.0 * 60.0 ** 2))

    def test_fractional_shift_matches_direct_interpolation(self):
        t0 = 165.7
        tau_grid = np.linspace(-150.0, 150.0, 77)
        out = eht.interp_rows(self.trace, t0, tau_grid)
        expected = np.interp(tau_grid + t0, self.rows, self.trace)
        # Both are linear interpolation of the same samples; 1% of peak.
        np.testing.assert_allclose(out, expected, atol=10.0, rtol=0.0)

    def test_integer_t0_samples_trace_exactly(self):
        tau_grid = np.arange(0.0, 1016.0 - 165.0)
        out = eht.interp_rows(self.trace, 165.0, tau_grid)
        np.testing.assert_allclose(
            out, self.trace[165:1016], atol=1e-12, rtol=0.0
        )

    def test_nan_outside_valid_rows(self):
        t0 = 165.7
        tau_grid = np.array([-166.0, -165.7, 0.0, 849.3, 850.0])
        # Absolute rows: -0.3, 0.0, 165.7, 1015.0, 1015.7.
        out = eht.interp_rows(self.trace, t0, tau_grid)
        np.testing.assert_array_equal(
            np.isnan(out), [True, False, False, False, True]
        )
        expected = np.interp(tau_grid + t0, self.rows, self.trace)
        np.testing.assert_allclose(out[1:4], expected[1:4], atol=1e-10, rtol=0.0)


class TestPinnedConventions(unittest.TestCase):
    def test_ps_per_px(self):
        self.assertEqual(1000.0 / 1016.0, eht.PS_PER_PX)

    def test_spot_roi(self):
        self.assertEqual((145, 190, 170, 215), tuple(eht.SPOT_ROI))


if __name__ == "__main__":
    unittest.main()
