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
import build_dataset_a as bda  # noqa: E402
import dataset_a_frame_injection as inj  # noqa: E402


def gaussian_blob(amplitude=3000.0, sigma=2.5):
    rows, cols = np.mgrid[0:15, 0:15]
    return amplitude * np.exp(
        -((rows - 7) ** 2 + (cols - 7) ** 2) / (2 * sigma ** 2)
    )


class StampDetectTest(unittest.TestCase):
    def test_stamped_blob_detected_near_center(self):
        rng = np.random.default_rng(1)
        frame = rng.normal(0.0, 3.0, (1016, 1344))
        blob = gaussian_blob()
        stamp = inj._stamp(frame, blob, 500.0, 415.0)
        self.assertIsNotNone(stamp)
        feats = bda.detect_events_features(frame, 25.0)
        self.assertGreaterEqual(feats["rows"].size, 1)
        i = int(np.argmax(feats["amplitude"]))
        self.assertLess(abs(feats["rows"][i] - 500.0), 2.0)
        self.assertLess(abs(feats["cols"][i] - 415.0), 2.0)

    def test_retention_matching(self):
        injected_rows = np.array([100.0, 200.0, 300.0])
        injected_cols = np.array([400.0, 500.0, 600.0])
        det_rows = np.array([101.0, 299.0])
        det_cols = np.array([401.0, 601.0])
        self.assertEqual(
            inj.match_retention(injected_rows, injected_cols, det_rows, det_cols),
            2,
        )


class RateModelTest(unittest.TestCase):
    def test_z_positions_respect_window_and_modulation(self):
        rng = np.random.default_rng(2)
        t0 = 155.0
        rows, cols = inj.draw_positions(
            rng, 0, 4000, t0, 9.84, 0.5, ramp_sign=+1.0)
        taus = rows - t0
        self.assertTrue(np.all(rows > inj.BLOB_HALF + 1))
        self.assertTrue(np.all(rows < 1016 - inj.BLOB_HALF - 1))
        self.assertTrue(np.all(taus >= 0) and np.all(taus <= 833))
        self.assertTrue(np.all(cols >= 380) and np.all(cols <= 450))
        # modulation sign check: fit A*sin(2 pi t / 9.84) with known phase
        t_ps = taus * inj.PS_PER_PX
        # bins fine enough that bin-averaging attenuation sinc(dt/T) ~ 0.98
        hist, edges = np.histogram(t_ps, bins=400)
        centers = 0.5 * (edges[:-1] + edges[1:])
        b = np.sin(2 * np.pi * centers / 9.84)
        slope = np.sum((hist - hist.mean()) * b) / np.sum(b**2) / hist.mean()
        self.assertGreater(slope, 0.35)  # injected +0.5
        rows2, _ = inj.draw_positions(
            rng, 0, 4000, t0, 9.84, 0.5, ramp_sign=-1.0)
        t_ps2 = (rows2 - t0) * inj.PS_PER_PX
        hist2, _ = np.histogram(t_ps2, bins=edges)
        slope2 = np.sum((hist2 - hist2.mean()) * b) / np.sum(b**2) / hist2.mean()
        self.assertLess(slope2, -0.35)

    def test_bunching_pair_rate_normalization(self):
        rng = np.random.default_rng(4)
        blobs = [gaussian_blob()]
        signal_rates = np.array([0.7, 1.2, 0.6, 1.0])
        a_pk = 6.0
        counts = []
        for _ in range(200):
            frame = np.zeros((1016, 1344))
            injected = inj.inject_bunching(
                frame, 155.0, signal_rates, a_pk, 9.84, blobs, rng)
            pairs_l = sum(1 for _, _, band in injected if band == 0)
            counts.append(pairs_l)
        win = inj.TAU_WIN_PX[1] - inj.TAU_WIN_PX[0]
        expected = (
            signal_rates[0] * signal_rates[2] / win
            * np.mean(0.75 + 0.25 * np.linspace(0, 1, 8))  # env mean ~0.875
            * a_pk / 2
        )
        self.assertAlmostEqual(np.mean(counts), expected, delta=0.35 * expected + 0.02)


if __name__ == "__main__":
    unittest.main()
