import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ANALYSIS_MODULE = "analysis.26-07-31-derived-single-sweep.visualize_detected_peaks"
ANALYSIS_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-07-31-derived-single-sweep"
    / "visualize_detected_peaks.py"
)


def load_module():
    try:
        return importlib.import_module(ANALYSIS_MODULE)
    except ModuleNotFoundError as exc:
        raise AssertionError(f"{ANALYSIS_MODULE} is not implemented") from exc


def write_peak_table(path: Path, xs, ys, attrs: dict) -> None:
    xs = np.asarray(xs, dtype=np.float64)
    ys = np.asarray(ys, dtype=np.float64)
    with h5py.File(path, "w") as handle:
        handle.create_dataset("frame", data=np.zeros(len(xs), dtype=np.int32))
        handle.create_dataset("x", data=xs)
        handle.create_dataset("y", data=ys)
        handle.create_dataset("mass", data=np.ones(len(xs), dtype=np.float64))
        for key, value in attrs.items():
            handle.attrs[key] = value


class DensityMapTests(unittest.TestCase):
    def test_bins_known_coordinates_on_the_detector_grid(self):
        density_map = load_module().density_map
        counts, x_edges, y_edges = density_map(
            np.array([2.0, 7.0]),
            np.array([3.0, 8.0]),
            frame_shape=(10, 12),
            bin_size=5,
        )

        self.assertEqual((2, 3), counts.shape)
        self.assertEqual(1, counts[0, 0])
        self.assertEqual(1, counts[1, 1])
        self.assertEqual(2, counts.sum())
        self.assertGreaterEqual(y_edges[-1], 10)
        self.assertGreaterEqual(x_edges[-1], 12)

    def test_rejects_nonpositive_bin_size(self):
        density_map = load_module().density_map
        with self.assertRaises(ValueError):
            density_map(
                np.array([1.0]), np.array([1.0]), frame_shape=(10, 10), bin_size=0
            )


class LoadPeakTableTests(unittest.TestCase):
    def test_reads_coordinates_and_optional_provenance(self):
        load_peak_table = load_module().load_peak_table
        with tempfile.TemporaryDirectory() as tmp:
            table_path = Path(tmp) / "peaks.h5"
            write_peak_table(
                table_path,
                xs=[10.5, 20.5],
                ys=[1.5, 2.5],
                attrs={
                    "peak_count": 2,
                    "frame_shape": np.array([16, 24]),
                    "x_range_intervals": np.array([[4.0, 12.0]]),
                    "sigma_multiple_used": 2.0,
                    "source_peak_count": 5,
                },
            )

            table = load_peak_table(table_path)

            np.testing.assert_array_equal(np.array([10.5, 20.5]), table["x"])
            np.testing.assert_array_equal(np.array([1.5, 2.5]), table["y"])
            self.assertEqual(2, table["peak_count"])
            self.assertEqual((16, 24), table["frame_shape"])
            np.testing.assert_array_equal(
                np.array([[4.0, 12.0]]), table["x_range_intervals"]
            )
            self.assertEqual(2.0, table["sigma_multiple_used"])
            self.assertEqual(5, table["source_peak_count"])

    def test_optional_attributes_are_absent_for_the_full_table(self):
        load_peak_table = load_module().load_peak_table
        with tempfile.TemporaryDirectory() as tmp:
            table_path = Path(tmp) / "peaks.h5"
            write_peak_table(
                table_path,
                xs=[10.5],
                ys=[1.5],
                attrs={"peak_count": 1, "frame_shape": np.array([16, 24])},
            )

            table = load_peak_table(table_path)

            self.assertNotIn("x_range_intervals", table)
            self.assertNotIn("source_peak_count", table)


class SaveOverviewFigureTests(unittest.TestCase):
    def build_tables(self):
        full = {
            "x": np.array([4.0, 6.0, 18.0]),
            "y": np.array([3.0, 9.0, 12.0]),
            "peak_count": 3,
            "frame_shape": (16, 24),
        }
        filtered = {
            "x": np.array([4.0, 6.0]),
            "y": np.array([3.0, 9.0]),
            "peak_count": 2,
            "x_range_intervals": np.array([[2.0, 8.0]]),
            "sigma_multiple_used": 2.0,
            "source_peak_count": 3,
        }
        return full, filtered

    def test_writes_a_png_in_a_new_output_directory(self):
        save_overview_figure = load_module().save_overview_figure
        full, filtered = self.build_tables()
        with tempfile.TemporaryDirectory() as tmp:
            output_path = Path(tmp) / "nested" / "overview.png"

            created = save_overview_figure(full, filtered, output_path, bin_size=2)

            self.assertEqual(output_path, created)
            self.assertEqual(b"\x89PNG\r\n\x1a\n", output_path.read_bytes()[:8])

    def test_refuses_to_overwrite_an_existing_figure(self):
        save_overview_figure = load_module().save_overview_figure
        full, filtered = self.build_tables()
        with tempfile.TemporaryDirectory() as tmp:
            output_path = Path(tmp) / "overview.png"
            output_path.write_bytes(b"existing")

            with self.assertRaises(FileExistsError):
                save_overview_figure(full, filtered, output_path, bin_size=2)

            self.assertEqual(b"existing", output_path.read_bytes())


class CliTests(unittest.TestCase):
    def test_cli_writes_the_overview_png(self):
        with tempfile.TemporaryDirectory() as tmp:
            full_path = Path(tmp) / "full.h5"
            filtered_path = Path(tmp) / "filtered.h5"
            output_path = Path(tmp) / "overview.png"
            write_peak_table(
                full_path,
                xs=[4.0, 6.0, 18.0],
                ys=[3.0, 9.0, 12.0],
                attrs={"peak_count": 3, "frame_shape": np.array([16, 24])},
            )
            write_peak_table(
                filtered_path,
                xs=[4.0, 6.0],
                ys=[3.0, 9.0],
                attrs={
                    "peak_count": 2,
                    "x_range_intervals": np.array([[2.0, 8.0]]),
                    "sigma_multiple_used": 2.0,
                    "source_peak_count": 3,
                },
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(ANALYSIS_SCRIPT),
                    str(full_path),
                    str(filtered_path),
                    str(output_path),
                    "--bin-size",
                    "2",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertTrue(output_path.exists())
            self.assertIn("3 peaks full", result.stdout)
            self.assertIn("2 peaks filtered", result.stdout)


if __name__ == "__main__":
    unittest.main()
