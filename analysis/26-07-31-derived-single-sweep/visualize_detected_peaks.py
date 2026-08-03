import argparse
import sys
from pathlib import Path

import h5py
import matplotlib
import numpy as np


matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import PowerNorm

# The x-acceptance rule is defined by the production range fit, imported so
# this figure cannot drift away from it.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_spot_x_ranges import acceptance_mask  # noqa: E402


def load_peak_table(peaks_path: Path) -> dict:
    """Read one detected-peaks HDF5 file: per-peak coordinates plus provenance."""
    with h5py.File(peaks_path, "r") as handle:
        table = {
            "frame": handle["frame"][...],
            "x": handle["x"][...],
            "y": handle["y"][...],
            "peak_count": int(handle.attrs["peak_count"]),
        }
        if "frame_count" in handle.attrs:
            table["frame_count"] = int(handle.attrs["frame_count"])
        if "frame_shape" in handle.attrs:
            table["frame_shape"] = tuple(
                int(value) for value in handle.attrs["frame_shape"]
            )
        if "x_range_intervals" in handle.attrs:
            table["x_range_intervals"] = np.asarray(
                handle.attrs["x_range_intervals"], dtype=np.float64
            )
        if "sigma_multiple_used" in handle.attrs:
            table["sigma_multiple_used"] = float(handle.attrs["sigma_multiple_used"])
        if "source_peak_count" in handle.attrs:
            table["source_peak_count"] = int(handle.attrs["source_peak_count"])
        return table


def density_map(
    y: np.ndarray,
    x: np.ndarray,
    frame_shape: tuple[int, int],
    bin_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Histogram peak coordinates onto the detector grid in bin_size pixel bins."""
    if bin_size <= 0:
        raise ValueError("bin_size must be positive")
    rows, cols = (int(frame_shape[0]), int(frame_shape[1]))
    y_edges = np.arange(0.0, rows + bin_size, bin_size)
    x_edges = np.arange(0.0, cols + bin_size, bin_size)
    counts, _, _ = np.histogram2d(y, x, bins=[y_edges, x_edges])
    return counts, x_edges, y_edges


def save_overview_figure(
    full: dict,
    filtered: dict,
    output_path: Path,
    bin_size: int = 4,
) -> Path:
    """Render full vs x-range-filtered peak density maps side by side."""
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing figure: {output_path}"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frame_shape = full.get("frame_shape")
    if frame_shape is None:
        raise ValueError("The full peak table must carry a frame_shape attribute")

    counts_full, x_edges, y_edges = density_map(
        full["y"], full["x"], frame_shape, bin_size
    )
    counts_filtered, _, _ = density_map(
        filtered["y"], filtered["x"], frame_shape, bin_size
    )
    extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
    norm = PowerNorm(
        gamma=0.5, vmin=0, vmax=max(counts_full.max(), counts_filtered.max(), 1)
    )

    figure, axes = plt.subplots(
        1, 2, figsize=(12.8, 6.0), constrained_layout=True, sharey=True
    )
    panels = [
        (
            counts_full,
            f"26-07-31_central_dot_0p18_peaks.h5\n{full['peak_count']:,} peaks (all detections)",
        ),
        (
            counts_filtered,
            "26-07-31_central_dot_0p18_peaks_xrange2p0sigma.h5\n"
            f"{filtered['peak_count']:,} peaks kept after the x-range filter",
        ),
    ]
    for axis, (counts, title) in zip(axes, panels):
        image = axis.imshow(
            counts,
            origin="lower",
            aspect="auto",
            cmap="magma",
            interpolation="nearest",
            norm=norm,
            extent=extent,
        )
        axis.set_title(title, fontsize=10)
        axis.set_xlabel("Detector axis 2 [pixel]")
    axes[0].set_ylabel("Detector axis 1 [pixel]")

    intervals = filtered.get("x_range_intervals")
    if intervals is not None:
        multiple = filtered.get("sigma_multiple_used")
        label = (
            f"kept x ranges (μ ± {multiple:g}σ)" if multiple is not None else "kept x ranges"
        )
        for index, (lower, upper) in enumerate(intervals):
            axes[0].axvspan(
                lower,
                upper,
                color="deepskyblue",
                alpha=0.13,
                label=label if index == 0 else None,
            )
        axes[0].legend(loc="upper right", fontsize=9)

    figure.suptitle(
        "central_dot_0p18 — detected peaks on the detector grid (1000 frames)"
    )
    colorbar = figure.colorbar(image, ax=axes, pad=0.02, shrink=0.95)
    colorbar.set_label(
        f"Peaks per {bin_size}×{bin_size} px bin\n(sqrt color scale, shared)"
    )
    figure.savefig(
        output_path,
        dpi=200,
        metadata={
            "Description": (
                "Density of detected sparse peaks on the detector grid: full "
                "peak table vs the x-range-filtered subset."
            )
        },
    )
    plt.close(figure)
    return output_path


def default_frames(frame_count: int, panels: int = 6) -> list[int]:
    """Uniformly spaced frame indices across the acquisition."""
    return np.linspace(0, frame_count - 1, panels).round().astype(int).tolist()


def save_frames_figure(
    full: dict,
    filtered: dict,
    output_path: Path,
    frames: list[int] | None = None,
) -> Path:
    """Scatter the detected peaks of individual frames, split by the x filter."""
    output_path = Path(output_path)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing figure: {output_path}"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frame_shape = full.get("frame_shape")
    if frame_shape is None:
        raise ValueError("The full peak table must carry a frame_shape attribute")
    if frames is None:
        frame_count = full.get("frame_count")
        if frame_count is None:
            raise ValueError(
                "Pass explicit frames: the full peak table lacks frame_count"
            )
        frames = default_frames(frame_count)
    frames = [int(frame_index) for frame_index in frames]

    intervals = filtered.get("x_range_intervals")
    rows = int(np.ceil(len(frames) / 3))
    figure, axes = plt.subplots(
        rows,
        3,
        figsize=(12.8, 3.1 * rows + 0.9),
        constrained_layout=True,
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    for panel, frame_index in enumerate(frames):
        axis = axes[panel // 3, panel % 3]
        in_frame = full["frame"] == frame_index
        xs = full["x"][in_frame]
        ys = full["y"][in_frame]
        if intervals is not None:
            kept = acceptance_mask(xs, intervals)
        else:
            kept = np.ones(xs.shape, dtype=bool)
        axis.scatter(
            xs[~kept],
            ys[~kept],
            s=10,
            c="0.55",
            alpha=0.7,
            linewidths=0,
            label="rejected by x-range filter",
        )
        axis.scatter(
            xs[kept],
            ys[kept],
            s=18,
            c="tab:red",
            linewidths=0,
            label="kept",
        )
        axis.set_title(
            f"frame {frame_index} — {len(xs)} peaks, {int(kept.sum())} kept",
            fontsize=10,
        )
        axis.set_xlim(0, frame_shape[1])
        axis.set_ylim(0, frame_shape[0])
        if panel % 3 == 0:
            axis.set_ylabel("Detector axis 1 [pixel]")
        if panel // 3 == rows - 1:
            axis.set_xlabel("Detector axis 2 [pixel]")
        if panel == 0:
            axis.legend(loc="upper right", fontsize=9)
    for panel in range(len(frames), rows * 3):
        axes[panel // 3, panel % 3].axis("off")

    figure.suptitle("central_dot_0p18 — detected peaks in individual frames")
    figure.savefig(
        output_path,
        dpi=200,
        metadata={
            "Description": (
                "Detected sparse peaks in individual frames: kept peaks vs "
                "peaks rejected by the x-range filter."
            )
        },
    )
    plt.close(figure)
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Render a side-by-side density overview of the full detected-peaks "
            "table and its x-range-filtered subset on the detector grid."
        )
    )
    parser.add_argument("full_peaks_h5", type=Path)
    parser.add_argument("filtered_peaks_h5", type=Path)
    parser.add_argument("output_png", type=Path)
    parser.add_argument("--bin-size", type=int, default=4)
    parser.add_argument(
        "--frames-png",
        type=Path,
        default=None,
        help="Also write a per-frame scatter overview to this path.",
    )
    parser.add_argument(
        "--frames",
        type=int,
        nargs="+",
        default=None,
        help="Frame indices for --frames-png (default: 6 uniformly spaced frames).",
    )
    arguments = parser.parse_args()

    full_table = load_peak_table(arguments.full_peaks_h5)
    filtered_table = load_peak_table(arguments.filtered_peaks_h5)
    created_path = save_overview_figure(
        full_table,
        filtered_table,
        output_path=arguments.output_png,
        bin_size=arguments.bin_size,
    )
    print(
        f"Created {created_path}: {full_table['peak_count']} peaks full, "
        f"{filtered_table['peak_count']} peaks filtered",
        flush=True,
    )
    if arguments.frames_png is not None:
        frames_path = save_frames_figure(
            full_table,
            filtered_table,
            arguments.frames_png,
            frames=arguments.frames,
        )
        print(f"Created {frames_path}: per-frame panels", flush=True)
