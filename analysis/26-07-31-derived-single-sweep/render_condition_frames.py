import argparse
from pathlib import Path

import h5py
import imageio.v3 as iio
import numpy as np


def find_condition_name(input_path: Path, central_dot_strength: float) -> str:
    """Return the single condition matching the requested central-dot strength."""
    with h5py.File(input_path, "r") as handle:
        conditions = handle["conditions"]
        matches = [
            name
            for name, condition in conditions.items()
            if np.isclose(
                float(condition.attrs["central_dot_strength"]),
                central_dot_strength,
                rtol=0.0,
                atol=1e-12,
            )
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one condition with central_dot_strength={central_dot_strength}, "
                f"found {len(matches)}"
            )
        return matches[0]


def even_frame_indices(frame_count: int, image_count: int) -> np.ndarray:
    """Return `image_count` evenly spaced frame indices covering [0, frame_count)."""
    if not 0 < image_count <= frame_count:
        raise ValueError(
            f"image_count must be within [1, {frame_count}], got {image_count}"
        )
    return np.arange(image_count) * frame_count // image_count


def condition_intensity_range(frames: h5py.Dataset, batch_size: int = 100) -> tuple[int, int]:
    """Return the global (min, max) pixel value across every frame."""
    minimum, maximum = None, None
    for start in range(0, frames.shape[0], batch_size):
        block = frames[start : start + batch_size]
        low, high = int(block.min()), int(block.max())
        minimum = low if minimum is None else min(minimum, low)
        maximum = high if maximum is None else max(maximum, high)
    return minimum, maximum


def build_sqrt_lut(minimum: int, maximum: int) -> np.ndarray:
    """Return a strictly increasing LUT mapping [minimum, maximum] counts to 16-bit.

    LUT[v] = round(65535 * sqrt((v - minimum) / (maximum - minimum))), so the
    background noise and the sparse bright events are both visible. Strict
    monotonicity keeps intensity in one-to-one correspondence with counts.
    """
    if not minimum < maximum:
        raise ValueError(f"Expected minimum < maximum, got {minimum} >= {maximum}")
    values = np.arange(minimum, maximum + 1, dtype=np.float64)
    lut = np.round(np.sqrt((values - minimum) / (maximum - minimum)) * 65535).astype(
        np.uint16
    )
    if np.any(np.diff(lut.astype(np.int64)) <= 0):
        raise ValueError("sqrt LUT is not strictly increasing; mapping is not one-to-one")
    return lut


def render_condition_frames(
    input_path: Path,
    output_dir: Path,
    central_dot_strength: float,
    image_count: int = 100,
) -> tuple[str, list[Path], int, int]:
    """Write evenly spaced frames as 16-bit grayscale PNGs through a sqrt LUT.

    The LUT is defined once from the condition's global intensity range, so all
    images share one mapping; it is strictly increasing, keeping pixel
    intensity in one-to-one correspondence with detector counts. PNG
    row/column (y, x) corresponds to frame array element
    [detector_axis_1, detector_axis_2] with no flipping.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    condition_name = find_condition_name(input_path, central_dot_strength)

    written: list[Path] = []
    with h5py.File(input_path, "r") as handle:
        frames = handle["conditions"][condition_name]["frames"]
        if frames.ndim != 3 or frames.dtype != np.uint16:
            raise ValueError(
                f"Expected a 3-D uint16 frames dataset, got {frames.dtype} {frames.shape}"
            )
        minimum, maximum = condition_intensity_range(frames)
        lut = build_sqrt_lut(minimum, maximum)
        for index in even_frame_indices(frames.shape[0], image_count):
            index = int(index)
            output_path = output_dir / f"frame_{index:04d}.png"
            iio.imwrite(output_path, lut[frames[index] - minimum])
            written.append(output_path)

    return condition_name, written, minimum, maximum


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Render evenly spaced frames of one central-dot-strength condition as "
            "16-bit grayscale PNGs through a strictly increasing sqrt LUT defined "
            "by the condition's global intensity range."
        )
    )
    parser.add_argument("input_h5", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--central-dot-strength", type=float, required=True)
    parser.add_argument("--image-count", type=int, default=100)
    arguments = parser.parse_args()

    selected_name, created_paths, minimum, maximum = render_condition_frames(
        arguments.input_h5,
        arguments.output_dir,
        central_dot_strength=arguments.central_dot_strength,
        image_count=arguments.image_count,
    )
    print(
        f"Created {len(created_paths)} PNGs from {selected_name} "
        f"in {arguments.output_dir}\n"
        f"LUT: round(65535 * sqrt((counts - {minimum}) / {maximum - minimum})), "
        f"strictly increasing on [{minimum}, {maximum}]",
        flush=True,
    )
