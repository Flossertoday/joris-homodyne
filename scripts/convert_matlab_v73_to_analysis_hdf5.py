import argparse
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
from typing import Callable

import h5py
import numpy as np


IMAGE_FILE_PATTERN = re.compile(r"(\d+)_image_savedStructure\.mat$")


@dataclass(frozen=True)
class ImageBlock:
    directory: Path
    file_numbers: list[int]
    files: list[Path]


def discover_image_blocks(source_root: Path) -> list[ImageBlock]:
    blocks = []
    for directory in sorted(path for path in source_root.rglob("*") if path.is_dir()):
        numbered_files = []
        for path in directory.glob("*.mat"):
            match = IMAGE_FILE_PATTERN.fullmatch(path.name)
            if match:
                numbered_files.append((int(match.group(1)), path))
        if not numbered_files:
            continue
        numbered_files.sort()
        blocks.append(
            ImageBlock(
                directory=directory,
                file_numbers=[number for number, _ in numbered_files],
                files=[path for _, path in numbered_files],
            )
        )
    return sorted(blocks, key=lambda block: block.file_numbers[0])


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _condition_name(central_dot_strength: float) -> str:
    value = format(central_dot_strength, ".12g").replace("-", "m").replace(".", "p")
    return f"central_dot_{value}"


def _write_string_dataset(group: h5py.Group, name: str, values):
    group.create_dataset(name, data=list(values), dtype=h5py.string_dtype("utf-8"))


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


def convert_dataset(
    source_root: Path,
    output_path: Path,
    expected_frames_per_condition: int | None = None,
    progress: Callable[[str], None] | None = None,
):
    source_root = source_root.resolve()
    output_path = output_path.resolve()
    blocks = discover_image_blocks(source_root)
    if not blocks:
        raise ValueError(f"No numbered image MAT files found under {source_root}")

    image_files = {path.resolve() for block in blocks for path in block.files}
    parameter_files = sorted(
        path for path in source_root.rglob("*.mat") if path.resolve() not in image_files
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with _atomic_hdf5_output(output_path) as output:
        output.attrs["schema_name"] = "joris_homodyne_streak_analysis"
        output.attrs["schema_version"] = 1
        output.attrs["source_format"] = "MATLAB 7.3 / HDF5"
        output.attrs["source_root"] = str(source_root)
        output.attrs["raw_payload_scope"] = "streak_only"
        output.attrs["payload_dataset"] = "structVar/StreakCamImage"
        output.attrs["excluded_camera_payloads"] = "MainCamImage,SecondCamImage"
        output.attrs["axis_semantics_status"] = "uncalibrated"

        metadata = output.create_group("metadata")
        readme_path = source_root / "README.txt"
        readme = readme_path.read_text(encoding="utf-8") if readme_path.exists() else ""
        metadata.create_dataset("readme", data=readme, dtype=h5py.string_dtype("utf-8"))

        parameter_group = metadata.create_group("parameter_files")
        _write_string_dataset(
            parameter_group,
            "relative_path",
            [path.relative_to(source_root).as_posix() for path in parameter_files],
        )
        _write_string_dataset(
            parameter_group,
            "sha256",
            [_sha256(path) for path in parameter_files],
        )
        parameter_group.create_dataset(
            "byte_size", data=[path.stat().st_size for path in parameter_files], dtype=np.int64
        )

        conditions = output.create_group("conditions")
        used_condition_names = set()
        for block in blocks:
            expected_numbers = list(range(block.file_numbers[0], block.file_numbers[-1] + 1))
            if block.file_numbers != expected_numbers:
                raise ValueError(f"Non-contiguous file numbers in {block.directory}")

            file_shapes = []
            strengths = []
            for path in block.files:
                with h5py.File(path, "r") as source:
                    struct = source["structVar"]
                    frames = struct["StreakCamImage"]
                    if frames.ndim != 3 or frames.dtype != np.dtype("uint16"):
                        raise ValueError(
                            f"Expected uint16 3-D StreakCamImage in {path}, got {frames.shape} {frames.dtype}"
                        )
                    file_shapes.append(frames.shape)
                    strengths.append(float(struct["ZeroOrderPowerRatio"][0, 0]))

            detector_shapes = {(shape[1], shape[2]) for shape in file_shapes}
            if len(detector_shapes) != 1:
                raise ValueError(f"Inconsistent detector shapes in {block.directory}: {detector_shapes}")
            if len(set(strengths)) != 1:
                raise ValueError(f"Inconsistent ZeroOrderPowerRatio in {block.directory}")

            total_frames = sum(shape[0] for shape in file_shapes)
            if expected_frames_per_condition is not None and total_frames != expected_frames_per_condition:
                raise ValueError(
                    f"Expected {expected_frames_per_condition} frames in {block.directory}, got {total_frames}"
                )

            central_dot_strength = strengths[0]
            condition_name = _condition_name(central_dot_strength)
            if condition_name in used_condition_names:
                raise ValueError(f"Duplicate central-dot condition {condition_name}")
            used_condition_names.add(condition_name)

            height, width = next(iter(detector_shapes))
            condition = conditions.create_group(condition_name)
            condition.attrs["central_dot_strength"] = central_dot_strength
            condition.attrs["central_dot_strength_source_field"] = "ZeroOrderPowerRatio"
            condition.attrs["source_block_directory"] = block.directory.relative_to(source_root).as_posix()
            condition.attrs["measurement_scope"] = "streak_only"

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
                "frame_source_file_number", shape=(total_frames,), dtype=np.int32
            )
            frame_in_source = condition.create_dataset(
                "frame_in_source_file", shape=(total_frames,), dtype=np.int16
            )

            source_files = condition.create_group("source_files")
            source_files.create_dataset("file_number", data=block.file_numbers, dtype=np.int32)
            _write_string_dataset(
                source_files,
                "relative_path",
                [path.relative_to(source_root).as_posix() for path in block.files],
            )
            _write_string_dataset(source_files, "sha256", [_sha256(path) for path in block.files])
            source_files.create_dataset(
                "byte_size", data=[path.stat().st_size for path in block.files], dtype=np.int64
            )
            source_files.create_dataset(
                "first_output_frame",
                data=np.cumsum([0] + [shape[0] for shape in file_shapes[:-1]]),
                dtype=np.int32,
            )
            source_files.create_dataset(
                "frame_count", data=[shape[0] for shape in file_shapes], dtype=np.int16
            )

            cursor = 0
            for number, path, shape in zip(block.file_numbers, block.files, file_shapes):
                count = shape[0]
                with h5py.File(path, "r") as source:
                    frames_out[cursor : cursor + count] = source["structVar/StreakCamImage"][:]
                frame_source[cursor : cursor + count] = number
                frame_in_source[cursor : cursor + count] = np.arange(count, dtype=np.int16)
                cursor += count
                if progress is not None:
                    progress(f"{condition_name}: {path.name} -> frames {cursor - count}:{cursor}")

    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert MATLAB 7.3 streak-camera files into a normalized analysis HDF5."
    )
    parser.add_argument("source_root", type=Path)
    parser.add_argument("output_path", type=Path)
    parser.add_argument("--expected-frames-per-condition", type=int, default=None)
    arguments = parser.parse_args()
    created = convert_dataset(
        arguments.source_root,
        arguments.output_path,
        expected_frames_per_condition=arguments.expected_frames_per_condition,
        progress=lambda message: print(message, flush=True),
    )
    print(f"Created {created}", flush=True)
