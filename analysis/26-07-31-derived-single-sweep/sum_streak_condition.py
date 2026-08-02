import argparse
from pathlib import Path

import h5py
import matplotlib
import numpy as np


matplotlib.use("Agg")
from matplotlib import pyplot as plt


def sum_condition_frames(
    input_path: Path,
    central_dot_strength: float,
    batch_size: int = 25,
) -> tuple[str, int, np.ndarray]:
    """Return the uint64 pixel sum for one central-dot condition."""
    if batch_size <= 0:
        raise ValueError("batch_size must be positive")

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

        condition_name = matches[0]
        frames = conditions[condition_name]["frames"]
        if frames.ndim != 3:
            raise ValueError(f"Expected a 3-D frames dataset, got shape {frames.shape}")

        summed = np.zeros(frames.shape[1:], dtype=np.uint64)
        for start in range(0, frames.shape[0], batch_size):
            summed += frames[start : start + batch_size].sum(axis=0, dtype=np.uint64)

        return condition_name, int(frames.shape[0]), summed


def save_sum_figure(
    summed: np.ndarray,
    central_dot_strength: float,
    frame_count: int,
    output_path: Path,
) -> Path:
    """Render a linearly scaled PNG of an element-wise frame sum."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lower, upper = np.percentile(summed, [1.0, 99.9])

    figure, axis = plt.subplots(figsize=(9.6, 7.2), constrained_layout=True)
    image = axis.imshow(
        summed,
        origin="lower",
        aspect="auto",
        cmap="magma",
        interpolation="nearest",
        vmin=lower,
        vmax=upper,
    )
    axis.set_title(
        f"Central-dot intensity {central_dot_strength:g} — sum of {frame_count:,} streak frames"
    )
    axis.set_xlabel("Detector axis 2 [pixel]")
    axis.set_ylabel("Detector axis 1 [pixel]")
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label("Summed camera counts [a.u.]\n(linear, clipped to p1–p99.9)")
    figure.savefig(
        output_path,
        dpi=200,
        metadata={
            "Description": (
                f"Element-wise uint64 sum of {frame_count} streak frames at "
                f"central_dot_strength={central_dot_strength:g}."
            )
        },
    )
    plt.close(figure)
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sum every streak frame for one central-dot-strength condition."
    )
    parser.add_argument("input_h5", type=Path)
    parser.add_argument("output_png", type=Path)
    parser.add_argument("--central-dot-strength", type=float, required=True)
    parser.add_argument("--batch-size", type=int, default=25)
    arguments = parser.parse_args()

    selected_name, selected_count, selected_sum = sum_condition_frames(
        arguments.input_h5,
        central_dot_strength=arguments.central_dot_strength,
        batch_size=arguments.batch_size,
    )
    created_path = save_sum_figure(
        selected_sum,
        central_dot_strength=arguments.central_dot_strength,
        frame_count=selected_count,
        output_path=arguments.output_png,
    )
    print(f"Created {created_path} from {selected_name}: {selected_count} frames", flush=True)
