import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from scripts.convert_matlab_v73_to_analysis_hdf5 import convert_dataset, discover_image_blocks


def write_mat(
    path: Path,
    frames: np.ndarray,
    central_dot_strength: float,
    cam2_mode: str = "Off",
):
    with h5py.File(path, "w") as handle:
        struct = handle.create_group("structVar")
        struct.create_dataset("StreakCamImage", data=frames)
        struct.create_dataset("ZeroOrderPowerRatio", data=[[central_dot_strength]])
        chars = np.array([ord(char) for char in cam2_mode], dtype=np.uint16).reshape(-1, 1)
        mode = struct.create_dataset("Cam2ContinuousCameraMode", data=chars)
        mode.attrs["MATLAB_class"] = np.bytes_("char")


class ConverterTests(unittest.TestCase):
    def test_cli_converts_a_dataset_and_reports_progress(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "source"
            output = Path(tmp) / "analysis.h5"
            block = root / "images"
            block.mkdir(parents=True)
            write_mat(
                block / "12_image_savedStructure.mat",
                np.zeros((1, 2, 3), dtype=np.uint16),
                0.05,
            )
            script = (
                Path(__file__).resolve().parents[1]
                / "scripts"
                / "convert_matlab_v73_to_analysis_hdf5.py"
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    str(root),
                    str(output),
                    "--expected-frames-per-condition",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(0, result.returncode, result.stderr)
            self.assertTrue(output.exists())
            self.assertIn("12_image_savedStructure.mat", result.stdout)

    def test_discovers_numeric_image_blocks_without_hard_coded_folder_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            parameter_dir = root / "parameters renamed"
            block_a = root / "first condition renamed"
            block_b = root / "second condition renamed"
            parameter_dir.mkdir()
            block_a.mkdir()
            block_b.mkdir()

            (parameter_dir / "experiment_savedStructure.mat").touch()
            for number in (13, 12):
                (block_a / f"{number}_image_savedStructure.mat").touch()
            for number in (15, 14):
                (block_b / f"{number}_image_savedStructure.mat").touch()

            blocks = discover_image_blocks(root)

            self.assertEqual([[12, 13], [14, 15]], [block.file_numbers for block in blocks])
            self.assertEqual([block_a, block_b], [block.directory for block in blocks])

    def test_refuses_to_overwrite_existing_analysis_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "source"
            output = Path(tmp) / "analysis.h5"
            block = root / "images"
            block.mkdir(parents=True)
            write_mat(
                block / "12_image_savedStructure.mat",
                np.zeros((1, 2, 3), dtype=np.uint16),
                0.05,
            )
            output.write_bytes(b"keep me")

            with self.assertRaises(FileExistsError):
                convert_dataset(root, output)

            self.assertEqual(b"keep me", output.read_bytes())

    def test_validation_failure_leaves_no_analysis_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "source"
            output = Path(tmp) / "derived" / "analysis.h5"
            block = root / "images"
            block.mkdir(parents=True)
            write_mat(
                block / "12_image_savedStructure.mat",
                np.zeros((1, 2, 3), dtype=np.uint16),
                0.05,
            )

            with self.assertRaises(ValueError):
                convert_dataset(root, output, expected_frames_per_condition=2)

            self.assertFalse(output.exists())
            self.assertEqual([], list(output.parent.glob("*.partial")))

    def test_converts_only_streak_payload_with_condition_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "source"
            output = Path(tmp) / "derived" / "analysis.h5"
            parameter_dir = root / "parameter state"
            block_a = root / "arbitrary first block"
            block_b = root / "arbitrary second block"
            parameter_dir.mkdir(parents=True)
            block_a.mkdir()
            block_b.mkdir()
            (root / "README.txt").write_text("streak-only test", encoding="utf-8")

            with h5py.File(parameter_dir / "parameters.mat", "w") as handle:
                struct = handle.create_group("structVar")
                struct.create_dataset("MainCamImage", data=np.ones((2, 3, 4)))
                struct.create_dataset("SecondCamImage", data=np.ones((2, 3, 4)))
                struct.create_dataset("StreakCamImage", data=np.ones((2, 4, 6), dtype=np.uint16))

            expected = {}
            for number, strength, folder in (
                (12, 0.05, block_a),
                (13, 0.05, block_a),
                (14, 0.18, block_b),
                (15, 0.18, block_b),
            ):
                frames = (np.arange(48, dtype=np.uint16).reshape(2, 4, 6) + number * 100)
                cam2_mode = "Off" if number % 2 == 0 else "On"
                write_mat(
                    folder / f"{number}_image_savedStructure.mat",
                    frames,
                    strength,
                    cam2_mode=cam2_mode,
                )
                expected[number] = frames

            progress_messages = []
            convert_dataset(
                root,
                output,
                expected_frames_per_condition=4,
                progress=progress_messages.append,
            )

            with h5py.File(output, "r") as converted:
                self.assertEqual("streak_only", converted.attrs["raw_payload_scope"])
                self.assertEqual("MainCamImage,SecondCamImage", converted.attrs["excluded_camera_payloads"])
                self.assertEqual("streak-only test", converted["metadata/readme"].asstr()[()])
                self.assertEqual(
                    {"central_dot_0p05", "central_dot_0p18"},
                    set(converted["conditions"].keys()),
                )
                first = converted["conditions/central_dot_0p05"]
                second = converted["conditions/central_dot_0p18"]
                np.testing.assert_array_equal(
                    np.concatenate([expected[12], expected[13]]), first["frames"][:]
                )
                np.testing.assert_array_equal(
                    np.concatenate([expected[14], expected[15]]), second["frames"][:]
                )
                np.testing.assert_array_equal(
                    [12, 12, 13, 13], first["frame_source_file_number"][:]
                )
                np.testing.assert_array_equal([0, 1, 0, 1], first["frame_in_source_file"][:])
                np.testing.assert_array_equal([12, 13], first["source_files/file_number"][:])
                self.assertTrue(all(len(value) == 64 for value in first["source_files/sha256"].asstr()[:]))
                self.assertEqual(0.05, first.attrs["central_dot_strength"])
                self.assertEqual("ZeroOrderPowerRatio", first.attrs["central_dot_strength_source_field"])
                self.assertNotIn("auxiliary_cam2_continuous_mode", first.attrs)
                self.assertEqual("gzip", first["frames"].compression)
                self.assertEqual(4, first["frames"].compression_opts)
                self.assertTrue(first["frames"].shuffle)
                self.assertTrue(first["frames"].fletcher32)
                self.assertNotIn("MainCamImage", converted)
                self.assertNotIn("SecondCamImage", converted)
                self.assertEqual(
                    ["parameter state/parameters.mat"],
                    converted["metadata/parameter_files/relative_path"].asstr()[:].tolist(),
                )
            self.assertTrue(any("12_image_savedStructure.mat" in message for message in progress_messages))


if __name__ == "__main__":
    unittest.main()
