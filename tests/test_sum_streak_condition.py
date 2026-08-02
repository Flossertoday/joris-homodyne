import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ANALYSIS_MODULE = "analysis.26-07-31-derived-single-sweep.sum_streak_condition"
ANALYSIS_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-07-31-derived-single-sweep"
    / "sum_streak_condition.py"
)


def load_sum_api():
    try:
        module = importlib.import_module(ANALYSIS_MODULE)
    except ModuleNotFoundError as exc:
        raise AssertionError(f"{ANALYSIS_MODULE} is not implemented") from exc
    return module.sum_condition_frames


def load_render_api():
    try:
        module = importlib.import_module(ANALYSIS_MODULE)
        save_sum_figure = module.save_sum_figure
    except (AttributeError, ModuleNotFoundError) as exc:
        raise AssertionError("save_sum_figure is not implemented") from exc
    return save_sum_figure


class SumStreakConditionTests(unittest.TestCase):
    def test_sums_the_selected_central_dot_condition_in_uint64(self):
        sum_condition_frames = load_sum_api()
        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.h5"
            with h5py.File(input_path, "w") as handle:
                conditions = handle.create_group("conditions")
                low = conditions.create_group("arbitrary_low_name")
                low.attrs["central_dot_strength"] = 0.05
                low.create_dataset("frames", data=np.full((2, 2, 3), 7, dtype=np.uint16))
                high = conditions.create_group("arbitrary_high_name")
                high.attrs["central_dot_strength"] = 0.18
                high.create_dataset(
                    "frames",
                    data=np.array(
                        [
                            [[1, 2, 3], [4, 5, 6]],
                            [[10, 20, 30], [40, 50, 60]],
                            [[100, 200, 300], [400, 500, 600]],
                        ],
                        dtype=np.uint16,
                    ),
                )

            condition_name, frame_count, summed = sum_condition_frames(
                input_path,
                central_dot_strength=0.18,
                batch_size=2,
            )

            self.assertEqual("arbitrary_high_name", condition_name)
            self.assertEqual(3, frame_count)
            self.assertEqual(np.dtype(np.uint64), summed.dtype)
            np.testing.assert_array_equal(
                np.array([[111, 222, 333], [444, 555, 666]], dtype=np.uint64),
                summed,
            )

    def test_saves_the_sum_as_a_png_in_a_new_output_directory(self):
        save_sum_figure = load_render_api()
        summed = np.arange(48, dtype=np.uint64).reshape(6, 8)

        with tempfile.TemporaryDirectory() as tmp:
            output_path = Path(tmp) / "nested" / "sum.png"

            created = save_sum_figure(
                summed,
                central_dot_strength=0.18,
                frame_count=1000,
                output_path=output_path,
            )

            self.assertEqual(output_path, created)
            self.assertEqual(b"\x89PNG\r\n\x1a\n", output_path.read_bytes()[:8])

    def test_cli_selects_a_strength_and_writes_the_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.h5"
            output_path = Path(tmp) / "result.png"
            with h5py.File(input_path, "w") as handle:
                condition = handle.create_group("conditions/central_dot_0p18")
                condition.attrs["central_dot_strength"] = 0.18
                condition.create_dataset(
                    "frames", data=np.arange(96, dtype=np.uint16).reshape(2, 6, 8)
                )
            result = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT),
                    str(input_path),
                    str(output_path),
                    "--central-dot-strength",
                    "0.18",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertTrue(output_path.exists())
            self.assertIn("central_dot_0p18", result.stdout)
            self.assertIn("2 frames", result.stdout)


if __name__ == "__main__":
    unittest.main()
