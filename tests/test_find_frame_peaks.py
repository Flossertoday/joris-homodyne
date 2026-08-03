import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ANALYSIS_MODULE = "analysis.26-07-31-derived-single-sweep.find_frame_peaks"
ANALYSIS_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-07-31-derived-single-sweep"
    / "find_frame_peaks.py"
)

BACKGROUND = 105.0
NOISE_SIGMA = 3.0
SPOT_SIGMA = 2.55  # FWHM ≈ 6 px, as measured on real frames
# Single-pixel impulse noise as measured on 0.18 frames (~500 events/frame):
# counts over background, p95 ≈ 43, max ≈ 106.
IMPULSE_RANGE = (25.0, 110.0)


def load_find_peaks():
    try:
        module = importlib.import_module(ANALYSIS_MODULE)
    except ModuleNotFoundError as exc:
        raise AssertionError(f"{ANALYSIS_MODULE} is not implemented") from exc
    return module.find_peaks


def add_spot(frame, y, x, peak, half=8):
    """Add one Gaussian spot centered at subpixel (y, x), clipped to the frame."""
    iy, ix = int(round(y)), int(round(x))
    y0, y1 = max(0, iy - half), min(frame.shape[0], iy + half + 1)
    x0, x1 = max(0, ix - half), min(frame.shape[1], ix + half + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    frame[y0:y1, x0:x1] += peak * np.exp(
        -((yy - y) ** 2 + (xx - x) ** 2) / (2 * SPOT_SIGMA**2)
    )


def make_synthetic_frame(
    rng,
    shape=(256, 256),
    n_spots=20,
    n_impulses=500,
    spot_peaks=(500.0, 3000.0),
    min_separation=15.0,
):
    """Background + Gaussian noise + Gaussian spots + single-pixel impulses."""
    frame = rng.normal(BACKGROUND, NOISE_SIGMA, shape)
    positions = []
    margin = 10
    while len(positions) < n_spots:
        y = rng.uniform(margin, shape[0] - margin)
        x = rng.uniform(margin, shape[1] - margin)
        if all(
            (y - py) ** 2 + (x - px) ** 2 >= min_separation**2 for py, px in positions
        ):
            positions.append((y, x))
    for y, x in positions:
        add_spot(frame, y, x, rng.uniform(*spot_peaks))
    for _ in range(n_impulses):
        frame[rng.integers(0, shape[0]), rng.integers(0, shape[1])] += rng.uniform(
            *IMPULSE_RANGE
        )
    return frame, np.array(positions)


def assert_matches_truth(test_case, peaks, truth):
    test_case.assertEqual((len(truth), 3), peaks.shape)
    remaining = peaks[:, :2].copy()
    for y, x in truth:
        distances = np.hypot(remaining[:, 0] - y, remaining[:, 1] - x)
        nearest = int(distances.argmin())
        test_case.assertLess(distances[nearest], 1.0)
        remaining = np.delete(remaining, nearest, axis=0)


class FindFramePeaksTests(unittest.TestCase):
    def test_recovers_synthetic_spots_among_impulse_noise(self):
        find_peaks = load_find_peaks()
        rng = np.random.default_rng(7)
        frame, truth = make_synthetic_frame(rng)

        peaks = find_peaks(frame)

        assert_matches_truth(self, peaks, truth)

    def test_returns_empty_when_only_impulse_noise(self):
        find_peaks = load_find_peaks()
        rng = np.random.default_rng(11)
        frame, _ = make_synthetic_frame(rng, n_spots=0)

        peaks = find_peaks(frame)

        self.assertEqual((0, 3), peaks.shape)

    def test_detects_spots_near_the_border(self):
        find_peaks = load_find_peaks()
        rng = np.random.default_rng(13)
        frame = rng.normal(BACKGROUND, NOISE_SIGMA, (256, 256))
        truth = np.array([(5.0, 5.0), (250.0, 250.0)])
        for y, x in truth:
            add_spot(frame, y, x, 1500.0)

        peaks = find_peaks(frame)

        assert_matches_truth(self, peaks, truth)

    def test_recovers_faint_spots_above_the_impulse_floor(self):
        find_peaks = load_find_peaks()
        rng = np.random.default_rng(19)
        # Peaks only ~50-100 counts over background, like the faintest events
        # seen in real frames; impulses stay within the measured real range.
        frame, truth = make_synthetic_frame(
            rng, n_spots=5, spot_peaks=(150.0, 400.0)
        )

        peaks = find_peaks(frame)

        assert_matches_truth(self, peaks, truth)

    def test_cli_writes_csv_and_overlay_for_selected_frame(self):
        rng = np.random.default_rng(17)
        empty_frame = rng.normal(BACKGROUND, NOISE_SIGMA, (256, 256))
        spot_frame, truth = make_synthetic_frame(rng, n_spots=5)
        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.h5"
            output_csv = Path(tmp) / "peaks.csv"
            overlay_png = Path(tmp) / "overlay.png"
            with h5py.File(input_path, "w") as handle:
                low = handle.create_group("conditions/central_dot_0p05")
                low.attrs["central_dot_strength"] = 0.05
                low.create_dataset(
                    "frames", data=np.zeros((1, 256, 256), dtype=np.uint16)
                )
                high = handle.create_group("conditions/central_dot_0p18")
                high.attrs["central_dot_strength"] = 0.18
                high.create_dataset(
                    "frames",
                    data=np.stack([empty_frame, spot_frame]).astype(np.uint16),
                )

            result = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT),
                    str(input_path),
                    str(output_csv),
                    "--central-dot-strength",
                    "0.18",
                    "--frame-index",
                    "1",
                    "--overlay-png",
                    str(overlay_png),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertIn("5 peaks", result.stdout)
            lines = output_csv.read_text().strip().splitlines()
            self.assertEqual("y,x,mass", lines[0])
            self.assertEqual(len(truth) + 1, len(lines))
            self.assertEqual(b"\x89PNG\r\n\x1a\n", overlay_png.read_bytes()[:8])


if __name__ == "__main__":
    unittest.main()
