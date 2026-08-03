import hashlib
import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


ANALYSIS_MODULE = "analysis.26-07-31-derived-single-sweep.extract_condition_peaks"
PEAKS_MODULE = "analysis.26-07-31-derived-single-sweep.find_frame_peaks"
ANALYSIS_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-07-31-derived-single-sweep"
    / "extract_condition_peaks.py"
)

BACKGROUND = 105.0
NOISE_SIGMA = 3.0
SPOT_SIGMA = 2.55  # FWHM ≈ 6 px, as measured on real frames


def add_spot(frame, y, x, peak, half=8):
    iy, ix = int(round(y)), int(round(x))
    y0, y1 = max(0, iy - half), min(frame.shape[0], iy + half + 1)
    x0, x1 = max(0, ix - half), min(frame.shape[1], ix + half + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    frame[y0:y1, x0:x1] += peak * np.exp(
        -((yy - y) ** 2 + (xx - x) ** 2) / (2 * SPOT_SIGMA**2)
    )


def make_spot_frame(rng, shape=(256, 256), spots=((40.0, 50.0), (200.0, 180.0))):
    frame = rng.normal(BACKGROUND, NOISE_SIGMA, shape)
    for y, x in spots:
        add_spot(frame, y, x, 2000.0)
    return frame


def write_two_condition_h5(path, frames_0p18):
    with h5py.File(path, "w") as handle:
        low = handle.create_group("conditions/central_dot_0p05")
        low.attrs["central_dot_strength"] = 0.05
        low.create_dataset("frames", data=np.zeros((1, 256, 256), dtype=np.uint16))
        high = handle.create_group("conditions/central_dot_0p18")
        high.attrs["central_dot_strength"] = 0.18
        high.create_dataset("frames", data=frames_0p18.astype(np.uint16))


def load_modules():
    try:
        extract = importlib.import_module(ANALYSIS_MODULE)
        peaks_module = importlib.import_module(PEAKS_MODULE)
    except ModuleNotFoundError as exc:
        raise AssertionError(f"{ANALYSIS_MODULE} is not implemented") from exc
    return extract, peaks_module


class ExtractConditionPeaksTests(unittest.TestCase):
    def test_cli_extracts_every_frame_matching_direct_find_peaks(self):
        extract, peaks_module = load_modules()
        rng = np.random.default_rng(23)
        quiet_frame = rng.normal(BACKGROUND, NOISE_SIGMA, (256, 256))
        spot_frame = make_spot_frame(rng)
        frames = np.stack([quiet_frame, spot_frame]).astype(np.uint16)

        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.h5"
            output_path = Path(tmp) / "peaks.h5"
            write_two_condition_h5(input_path, frames)

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

            with h5py.File(output_path, "r") as handle:
                self.assertEqual("central_dot_0p18", handle.attrs["condition_name"])
                self.assertEqual(0.18, handle.attrs["central_dot_strength"])
                self.assertEqual(2, handle.attrs["frame_count"])
                self.assertEqual((256, 256), tuple(handle.attrs["frame_shape"]))
                self.assertEqual(7, handle.attrs["diameter"])
                self.assertEqual(300.0, handle.attrs["minmass"])
                self.assertEqual(
                    hashlib.sha256(input_path.read_bytes()).hexdigest(),
                    handle.attrs["source_h5_sha256"],
                )
                out_frame = np.asarray(handle["frame"])
                out_y = np.asarray(handle["y"])
                out_x = np.asarray(handle["x"])
                out_mass = np.asarray(handle["mass"])
                peak_count = int(handle.attrs["peak_count"])
            self.assertEqual(len(out_frame), peak_count)

            for index in range(2):
                direct = peaks_module.find_peaks(frames[index].astype(np.float64))
                selected = out_frame == index
                self.assertEqual(
                    len(direct),
                    int(selected.sum()),
                    f"frame {index} peak count mismatch",
                )
                if len(direct):
                    # find_peaks returns rows ordered like tp.locate; the
                    # extraction preserves that order within each frame.
                    np.testing.assert_allclose(out_y[selected], direct[:, 0])
                    np.testing.assert_allclose(out_x[selected], direct[:, 1])
                    np.testing.assert_allclose(out_mass[selected], direct[:, 2])

    def test_cli_refuses_to_overwrite_existing_output(self):
        rng = np.random.default_rng(29)
        frames = make_spot_frame(rng)[np.newaxis, :, :].astype(np.uint16)

        with tempfile.TemporaryDirectory() as tmp:
            input_path = Path(tmp) / "input.h5"
            output_path = Path(tmp) / "peaks.h5"
            write_two_condition_h5(input_path, frames)

            command = [
                sys.executable,
                str(ANALYSIS_SCRIPT),
                str(input_path),
                str(output_path),
                "--central-dot-strength",
                "0.18",
            ]
            first = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(0, first.returncode, first.stderr)
            second = subprocess.run(command, text=True, capture_output=True)
            self.assertNotEqual(0, second.returncode)
            self.assertIn("Refusing to overwrite", second.stderr)


if __name__ == "__main__":
    unittest.main()
