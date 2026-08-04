import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ANALYSIS_MODULE = "analysis.26-07-31-derived-single-sweep.localization_uncertainty"
ANALYSIS_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-07-31-derived-single-sweep"
    / "localization_uncertainty.py"
)
PEAKS_MODULE = "analysis.26-07-31-derived-single-sweep.find_frame_peaks"

BACKGROUND = 105.0
NOISE_SIGMA = 3.0


def load_module(name):
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as exc:
        raise AssertionError(f"{name} is not implemented") from exc


def add_spot(frame, y, x, peak, sigma=2.55, half=8):
    iy, ix = int(round(y)), int(round(x))
    y0, y1 = max(0, iy - half), min(frame.shape[0], iy + half + 1)
    x0, x1 = max(0, ix - half), min(frame.shape[1], ix + half + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    frame[y0:y1, x0:x1] += peak * np.exp(
        -((yy - y) ** 2 + (xx - x) ** 2) / (2 * sigma**2)
    )


def make_frames_h5(path, frames):
    with h5py.File(path, "w") as handle:
        condition = handle.create_group("conditions/central_dot_0p18")
        condition.attrs["central_dot_strength"] = 0.18
        condition.create_dataset("frames", data=frames.astype(np.uint16))


def make_peaks_h5(path, frame_index, peaks):
    with h5py.File(path, "w") as handle:
        handle.create_dataset("frame", data=np.full(len(peaks), frame_index, np.int32))
        handle.create_dataset("y", data=peaks[:, 0])
        handle.create_dataset("x", data=peaks[:, 1])
        handle.create_dataset("mass", data=peaks[:, 2])
        handle.attrs["condition_name"] = "central_dot_0p18"
        handle.attrs["peak_count"] = len(peaks)
        handle.attrs["source_h5"] = "synthetic.h5"
        handle.attrs["source_h5_sha256"] = "0" * 64


class BandpassNoiseTests(unittest.TestCase):
    def test_recovers_gaussian_noise_scale(self):
        module = load_module(ANALYSIS_MODULE)
        rng = np.random.default_rng(5)
        frame = rng.normal(BACKGROUND, NOISE_SIGMA, (512, 512))

        noise = module.frame_background_noise(frame)

        # White noise sigma=3 through gaussian(sigma=1) minus 7x7 boxcar:
        # expected std is of order 3/sqrt(4*pi) ~ 0.85, far from raw 3.
        self.assertGreater(noise, 0.5)
        self.assertLess(noise, 1.5)

    def test_bright_spots_do_not_bias_the_estimate(self):
        module = load_module(ANALYSIS_MODULE)
        rng = np.random.default_rng(7)
        clean = rng.normal(BACKGROUND, NOISE_SIGMA, (512, 512))
        spotted = clean.copy()
        add_spot(spotted, 100, 100, 5000.0)
        add_spot(spotted, 300, 400, 2000.0)

        self.assertAlmostEqual(
            module.frame_background_noise(clean),
            module.frame_background_noise(spotted),
            delta=0.05,
        )


class StaticEpTests(unittest.TestCase):
    def test_ep_is_positive_and_shrinks_with_mass(self):
        module = load_module(ANALYSIS_MODULE)
        ep = module.static_ep([300.0, 1000.0, 3000.0, 10000.0], noise=0.85)

        self.assertTrue(np.all(np.isfinite(ep)))
        self.assertTrue(np.all(ep > 0))
        self.assertTrue(np.all(np.diff(ep) < 0))
        # Savin-Doyle scaling: ep inversely proportional to mass.
        self.assertAlmostEqual(ep[0] / ep[3], 10000.0 / 300.0, places=6)


class CliTests(unittest.TestCase):
    def test_cli_writes_ep_table_and_figure(self):
        module = load_module(ANALYSIS_MODULE)
        peaks_module = load_module(PEAKS_MODULE)
        rng = np.random.default_rng(11)
        frame = rng.normal(BACKGROUND, NOISE_SIGMA, (256, 256))
        for y, x, peak in [(60, 70, 800.0), (150, 180, 3000.0), (200, 90, 500.0)]:
            add_spot(frame, y, x, peak)
        peaks = peaks_module.find_peaks(frame)
        self.assertEqual(3, len(peaks))

        with tempfile.TemporaryDirectory() as tmp:
            frames_path = Path(tmp) / "frames.h5"
            peaks_path = Path(tmp) / "peaks.h5"
            output_path = Path(tmp) / "peaks_ep.h5"
            figure_path = Path(tmp) / "ep.png"
            make_frames_h5(frames_path, frame[None])
            make_peaks_h5(peaks_path, 0, peaks)

            result = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT),
                    str(peaks_path),
                    str(frames_path),
                    str(output_path),
                    "--figure",
                    str(figure_path),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(0, result.returncode, result.stderr)
            with h5py.File(output_path, "r") as handle:
                self.assertEqual(
                    {"frame", "y", "x", "mass", "ep", "frame_noise"},
                    set(handle.keys()),
                )
                ep = np.asarray(handle["ep"])
                self.assertEqual(3, len(ep))
                self.assertTrue(np.all(np.isfinite(ep)))
                self.assertTrue(np.all(ep > 0))
                self.assertTrue(np.all(ep < 1.0))
                self.assertEqual(1, len(handle["frame_noise"]))
                self.assertGreater(float(handle["frame_noise"][0]), 0.0)
            self.assertEqual(b"\x89PNG\r\n\x1a\n", figure_path.read_bytes()[:8])

            rerun = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT),
                    str(peaks_path),
                    str(frames_path),
                    str(output_path),
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(0, rerun.returncode)
            self.assertIn("Refusing to overwrite", rerun.stderr)


if __name__ == "__main__":
    unittest.main()
