import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ANALYSIS_MODULE = "analysis.26-07-31-derived-single-sweep.artificial_spot_test"
ANALYSIS_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-07-31-derived-single-sweep"
    / "artificial_spot_test.py"
)

BACKGROUND = 105.0
NOISE_SIGMA = 3.0


def load_module(name):
    try:
        return importlib.import_module(name)
    except ModuleNotFoundError as exc:
        raise AssertionError(f"{name} is not implemented") from exc


def make_frames_h5(path, frames, strength=0.18):
    with h5py.File(path, "w") as handle:
        condition = handle.create_group("conditions/central_dot_0p18")
        condition.attrs["central_dot_strength"] = strength
        condition.create_dataset("frames", data=frames.astype(np.uint16))


class SpotModelTests(unittest.TestCase):
    def test_unit_bandpassed_mass_is_a_small_positive_fraction(self):
        module = load_module(ANALYSIS_MODULE)
        unit = module.unit_bandpassed_mass(2.55, 2.55)

        self.assertGreater(unit, 0.05)
        self.assertLess(unit, 0.5)

    def test_draw_injections_respects_separation_and_ranges(self):
        module = load_module(ANALYSIS_MODULE)
        rng = np.random.default_rng(3)
        injections = module.draw_injections(
            rng, 6, (256, 256), mass_range=(2000.0, 5000.0)
        )

        self.assertEqual(6, len(injections))
        for index, spot in enumerate(injections):
            self.assertGreaterEqual(spot["y"], module.MARGIN)
            self.assertLessEqual(spot["y"], 256 - module.MARGIN)
            self.assertGreaterEqual(spot["mass_target"], 2000.0)
            self.assertLessEqual(spot["mass_target"], 5000.0)
            self.assertGreater(spot["total_counts"], spot["mass_target"])
            for other in injections[:index]:
                distance = np.hypot(spot["y"] - other["y"], spot["x"] - other["x"])
                self.assertGreaterEqual(distance, module.MIN_SEPARATION)


class InjectionAndMatchTests(unittest.TestCase):
    def test_bright_injections_are_recovered_on_truth(self):
        module = load_module(ANALYSIS_MODULE)
        peaks_module = load_module(
            "analysis.26-07-31-derived-single-sweep.find_frame_peaks"
        )
        rng = np.random.default_rng(9)
        frame = rng.normal(BACKGROUND, NOISE_SIGMA, (256, 256))
        injections = module.draw_injections(
            rng, 4, frame.shape, mass_range=(4000.0, 8000.0)
        )

        injected = module.inject_spots(frame, injections, rng)
        peaks = peaks_module.find_peaks(injected)
        matched = module.match_to_truth(injections, peaks)

        self.assertTrue(all(row["matched"] for row in matched))
        for row in matched:
            self.assertLess(abs(row["dy"]), 0.5)
            self.assertLess(abs(row["dx"]), 0.5)
            self.assertGreater(row["mass_det"], 300.0)

    def test_dim_injections_below_minmass_are_missed(self):
        module = load_module(ANALYSIS_MODULE)
        peaks_module = load_module(
            "analysis.26-07-31-derived-single-sweep.find_frame_peaks"
        )
        rng = np.random.default_rng(13)
        frame = rng.normal(BACKGROUND, NOISE_SIGMA, (256, 256))
        injections = module.draw_injections(
            rng, 3, frame.shape, mass_range=(60.0, 120.0)
        )

        injected = module.inject_spots(frame, injections, rng)
        peaks = peaks_module.find_peaks(injected)
        matched = module.match_to_truth(injections, peaks)

        self.assertFalse(any(row["matched"] for row in matched))
        for row in matched:
            self.assertTrue(np.isnan(row["dy"]))
            self.assertTrue(np.isnan(row["mass_det"]))


class CliTests(unittest.TestCase):
    def test_cli_end_to_end_on_synthetic_frames(self):
        rng = np.random.default_rng(17)
        frames = np.stack(
            [rng.normal(BACKGROUND, NOISE_SIGMA, (256, 256)) for _ in range(3)]
        )
        with tempfile.TemporaryDirectory() as tmp:
            frames_path = Path(tmp) / "frames.h5"
            output_csv = Path(tmp) / "spots.csv"
            figure_path = Path(tmp) / "errors.png"
            make_frames_h5(frames_path, frames)

            result = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT),
                    str(frames_path),
                    str(output_csv),
                    "--central-dot-strength",
                    "0.18",
                    "--n-frames",
                    "3",
                    "--spots-per-frame",
                    "2",
                    "--mass-min",
                    "4000",
                    "--mass-max",
                    "8000",
                    "--figure",
                    str(figure_path),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertIn("6 injections, 6 matched", result.stdout)
            lines = output_csv.read_text().strip().splitlines()
            self.assertEqual(
                "frame,y_true,x_true,mass_target,sigma,ecc,angle_deg,"
                "total_counts,matched,y_det,x_det,dy,dx,mass_det,ep_det,"
                "frame_noise",
                lines[0],
            )
            self.assertEqual(7, len(lines))
            self.assertEqual(b"\x89PNG\r\n\x1a\n", figure_path.read_bytes()[:8])

            rerun = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT),
                    str(frames_path),
                    str(output_csv),
                    "--central-dot-strength",
                    "0.18",
                    "--n-frames",
                    "3",
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(0, rerun.returncode)
            self.assertIn("Refusing to overwrite", rerun.stderr)


if __name__ == "__main__":
    unittest.main()
