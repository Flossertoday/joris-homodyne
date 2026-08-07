"""Convert the supplementary HBT acquisition (7200 MAT files, 72000 frames) into a
normalized analysis HDF5 that mirrors the schema-2 layout of
``26-08-04_hbt_time_reference.h5`` closely enough for the fig2a/fig2c pipelines.

The supplement lives in ``DATA/26-08-04-hbt-time-reference/New folder (10)/New folder (10)``
and was verified (2026-08-06) to carry acquisition metadata identical to the original
``01_HBT_condensate_Josephson_regime`` condition except ``SweepNbStepsVar1`` (7200 vs 2000,
i.e. the sweep length itself). The measured single-sweep background
(``structVar/StreakCamBackG``) exists only in the original batch's condition 09; it was
verified against the supplement's off-signal pixels (50-frame median differs by
+0.02 counts on a ~105 counts/px level) and is copied into this file from the original
analysis HDF5 so the pipelines can subtract it unchanged.

Usage:
    python convert_hbt_supplement.py SOURCE_ROOT ORIGINAL_H5 OUTPUT_H5
"""

import argparse
from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import re

import h5py
import numpy as np


IMAGE_FILE_PATTERN = re.compile(r"(\d+)_image_savedStructure\.mat$")
CONDITION_NAME = "hbt_condensate_josephson_regime"
BACKGROUND_CONDITION = "background_and_streak_pulses_delay_calibration"
BACKGROUND_DATASET = f"conditions/{BACKGROUND_CONDITION}/background"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


@contextmanager
def _atomic_hdf5_output(output_path: Path):
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing analysis file: {output_path}")
    partial_path = output_path.with_name(f"{output_path.name}.partial")
    if partial_path.exists():
        raise FileExistsError(f"Interrupted conversion artifact already exists: {partial_path}")
    created = False
    try:
        with h5py.File(partial_path, "x") as output:
            created = True
            yield output
        if output_path.exists():
            raise FileExistsError(f"Analysis file appeared during conversion: {output_path}")
        partial_path.rename(output_path)
        created = False
    finally:
        if created and partial_path.exists():
            partial_path.unlink()


def _write_string_dataset(group: h5py.Group, name: str, values):
    group.create_dataset(name, data=list(values), dtype=h5py.string_dtype("utf-8"))


def discover_files(source_root: Path):
    matches = []
    for path in source_root.rglob("*.mat"):
        match = IMAGE_FILE_PATTERN.fullmatch(path.name)
        if match:
            matches.append((int(match.group(1)), path))
    if not matches:
        raise ValueError(f"No numbered image MAT files found under {source_root}")
    matches.sort()
    numbers = [number for number, _ in matches]
    if numbers != list(range(numbers[0], numbers[-1] + 1)):
        raise ValueError("Non-contiguous file numbers in supplement")
    directories = {path.parent.resolve() for _, path in matches}
    if len(directories) != 1:
        raise ValueError(f"Supplement files spread over multiple directories: {directories}")
    return numbers, [path for _, path in matches]


def convert_supplement(source_root: Path, original_h5: Path, output_path: Path,
                       verify_stride: int = 29, progress=print):
    source_root = source_root.resolve()
    original_h5 = original_h5.resolve()
    output_path = output_path.resolve()
    numbers, files = discover_files(source_root)

    file_shapes = []
    strengths = []
    for path in files:
        with h5py.File(path, "r") as source:
            struct = source["structVar"]
            frames = struct["StreakCamImage"]
            if frames.ndim != 3 or frames.dtype != np.dtype("uint16"):
                raise ValueError(f"Expected uint16 3-D StreakCamImage in {path}")
            file_shapes.append(frames.shape)
            strengths.append(float(struct["ZeroOrderPowerRatio"][0, 0]))
    detector_shapes = {(shape[1], shape[2]) for shape in file_shapes}
    if len(detector_shapes) != 1:
        raise ValueError(f"Inconsistent detector shapes: {detector_shapes}")
    if len(set(strengths)) != 1:
        raise ValueError("Inconsistent ZeroOrderPowerRatio in supplement")

    total_frames = sum(shape[0] for shape in file_shapes)
    height, width = next(iter(detector_shapes))
    first_output_frame = np.cumsum([0] + [shape[0] for shape in file_shapes[:-1]])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with _atomic_hdf5_output(output_path) as output, h5py.File(original_h5, "r") as orig:
        output.attrs["schema_name"] = "joris_homodyne_streak_analysis"
        output.attrs["schema_version"] = 2
        output.attrs["source_format"] = "MATLAB 7.3 / HDF5"
        output.attrs["source_root"] = str(source_root)
        output.attrs["raw_payload_scope"] = "streak_only"
        output.attrs["payload_dataset"] = "structVar/StreakCamImage"
        output.attrs["excluded_camera_payloads"] = "MainCamImage,SecondCamImage"
        output.attrs["axis_semantics_status"] = "uncalibrated"
        output.attrs["condition_mode"] = "directory"
        output.attrs["supplement_of"] = str(original_h5)
        output.attrs["background_source"] = (
            f"copied from {original_h5.name}:{BACKGROUND_DATASET}; verified against "
            "supplement off-signal pixels 2026-08-06 (50-frame median diff +0.02 counts)")

        conditions = output.create_group("conditions")
        condition = conditions.create_group(CONDITION_NAME)
        condition.attrs["central_dot_strength"] = strengths[0]
        condition.attrs["central_dot_strength_source_field"] = "ZeroOrderPowerRatio"
        condition.attrs["source_block_directory"] = files[0].parent.relative_to(source_root).as_posix()
        condition.attrs["measurement_scope"] = "streak_only"
        condition.attrs["supplement"] = True

        frames_out = condition.create_dataset(
            "frames",
            shape=(total_frames, height, width),
            dtype=np.uint16,
            chunks=(1, min(254, height), min(224, width)),
            compression="gzip",
            compression_opts=4,
            shuffle=True,
            fletcher32=True,
        )
        frames_out.attrs["axis_0"] = "acquisition_frame"
        frames_out.attrs["axis_1"] = "detector_axis_1"
        frames_out.attrs["axis_2"] = "detector_axis_2"
        frames_out.attrs["axis_semantics_status"] = "uncalibrated"
        frame_source = condition.create_dataset(
            "frame_source_file_number", shape=(total_frames,), dtype=np.int32)
        frame_in_source = condition.create_dataset(
            "frame_in_source_file", shape=(total_frames,), dtype=np.int16)

        source_files = condition.create_group("source_files")
        source_files.create_dataset("file_number", data=numbers, dtype=np.int32)
        _write_string_dataset(
            source_files, "relative_path",
            [path.relative_to(source_root).as_posix() for path in files])
        _write_string_dataset(source_files, "sha256", [_sha256(path) for path in files])
        source_files.create_dataset(
            "byte_size", data=[path.stat().st_size for path in files], dtype=np.int64)
        source_files.create_dataset("first_output_frame", data=first_output_frame, dtype=np.int32)
        source_files.create_dataset(
            "frame_count", data=[shape[0] for shape in file_shapes], dtype=np.int16)

        cursor = 0
        for number, path, shape in zip(numbers, files, file_shapes):
            count = shape[0]
            with h5py.File(path, "r") as source:
                frames_out[cursor : cursor + count] = source["structVar/StreakCamImage"][:]
            frame_source[cursor : cursor + count] = number
            frame_in_source[cursor : cursor + count] = np.arange(count, dtype=np.int16)
            cursor += count
            if progress is not None:
                progress(f"{CONDITION_NAME}: {path.name} -> frames {cursor - count}:{cursor}")

        bg_condition = conditions.create_group(BACKGROUND_CONDITION)
        for key, value in orig[f"conditions/{BACKGROUND_CONDITION}"].attrs.items():
            bg_condition.attrs[key] = value
        bg_condition.attrs["copied_from"] = str(original_h5)
        bg_condition.create_dataset("background", data=orig[BACKGROUND_DATASET][:])
        bg_condition["background"].attrs["source_note"] = (
            "measured StreakCamBackG of the original 26-08-04 batch; verified applicable "
            "to this supplement (see root attr background_source)")

    # Stratified exact verification: first/last frame of every verify_stride-th file.
    checked = 0
    mismatches = []
    with h5py.File(output_path, "r") as output:
        frames_out = output[f"conditions/{CONDITION_NAME}/frames"]
        for index in range(0, len(files), verify_stride):
            path = files[index]
            start = int(first_output_frame[index])
            count = int(file_shapes[index][0])
            with h5py.File(path, "r") as source:
                block = source["structVar/StreakCamImage"]
                for local in (0, count - 1):
                    if not np.array_equal(frames_out[start + local], block[local]):
                        mismatches.append(f"{path.name} frame {local}")
                    checked += 1
    report = {
        "files": len(files),
        "total_frames": total_frames,
        "frames_checked_exact": checked,
        "mismatches": mismatches,
        "output": str(output_path),
    }
    print(json.dumps(report, indent=2))
    if mismatches:
        raise ValueError(f"Verification failed: {mismatches[:5]}")
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument("original_h5", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("--verify-stride", type=int, default=29)
    arguments = parser.parse_args()
    convert_supplement(
        arguments.source_root, arguments.original_h5, arguments.output_path,
        verify_stride=arguments.verify_stride,
        progress=lambda message: print(message, flush=True),
    )
