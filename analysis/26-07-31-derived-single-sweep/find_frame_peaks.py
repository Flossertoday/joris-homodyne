import argparse
import warnings
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import trackpy as tp


matplotlib.use("Agg")
from matplotlib import pyplot as plt


# Measured on central_dot_0p18 frames: flat ~105-count background, ~3-count
# read noise, spot FWHM ~6-7 px, and ~500 single-pixel impulse events per
# frame (max ~110 counts over background). diameter must be odd and slightly
# larger than the spot. trackpy's default percentile=64 cut would discard the
# faintest real spots on these sparse frames, so it is disabled; the bandpass
# plus the minmass cut (bandpassed integrated brightness) reject the impulse
# noise with a wide margin instead.
DEFAULT_DIAMETER = 7
DEFAULT_MINMASS = 300.0
DEFAULT_PERCENTILE = 0.0


def find_peaks(
    frame: np.ndarray,
    diameter: int = DEFAULT_DIAMETER,
    minmass: float = DEFAULT_MINMASS,
    percentile: float = DEFAULT_PERCENTILE,
) -> np.ndarray:
    """Return one row per Gaussian-like peak: subpixel y, x, and bandpassed mass."""
    if diameter % 2 == 0:
        raise ValueError("diameter must be odd")
    # A frame legitimately containing no peaks is an expected outcome, not a
    # warning condition; trackpy would otherwise warn once per empty frame.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="No maxima survived.*", category=UserWarning
        )
        features = tp.locate(
            np.asarray(frame, dtype=np.float64),
            diameter=diameter,
            minmass=minmass,
            percentile=percentile,
        )
    if features.empty:
        return np.empty((0, 3), dtype=np.float64)
    return features[["y", "x", "mass"]].to_numpy(dtype=np.float64)


def load_condition_frame(
    input_path: Path,
    central_dot_strength: float,
    frame_index: int,
) -> tuple[str, np.ndarray]:
    """Return the condition name and one frame selected by central-dot strength."""
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
        if not 0 <= frame_index < frames.shape[0]:
            raise ValueError(
                f"frame_index {frame_index} out of range for {frames.shape[0]} frames"
            )
        return condition_name, np.asarray(frames[frame_index])


def save_overlay_figure(
    frame: np.ndarray,
    peaks: np.ndarray,
    condition_name: str,
    frame_index: int,
    output_path: Path,
) -> Path:
    """Render the frame with detected peaks circled, for manual review."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lower, upper = np.percentile(frame, [1.0, 99.9])

    figure, axis = plt.subplots(figsize=(9.6, 7.2), constrained_layout=True)
    axis.imshow(
        frame,
        origin="upper",
        aspect="equal",
        cmap="magma",
        interpolation="nearest",
        vmin=lower,
        vmax=upper,
    )
    if len(peaks):
        axis.plot(
            peaks[:, 1],
            peaks[:, 0],
            "o",
            markerfacecolor="none",
            markeredgecolor="cyan",
            markeredgewidth=1.0,
            markersize=12,
        )
    axis.set_title(f"{condition_name} frame {frame_index} — {len(peaks)} peaks")
    axis.set_xlabel("Detector axis 2 [pixel]")
    axis.set_ylabel("Detector axis 1 [pixel]")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Locate Gaussian-like peaks in one streak frame (trackpy Crocker-Grier)."
    )
    parser.add_argument("input_h5", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--central-dot-strength", type=float, required=True)
    parser.add_argument("--frame-index", type=int, default=0)
    parser.add_argument("--diameter", type=int, default=DEFAULT_DIAMETER)
    parser.add_argument("--minmass", type=float, default=DEFAULT_MINMASS)
    parser.add_argument("--overlay-png", type=Path, default=None)
    arguments = parser.parse_args()

    condition_name, frame = load_condition_frame(
        arguments.input_h5,
        central_dot_strength=arguments.central_dot_strength,
        frame_index=arguments.frame_index,
    )
    peaks = find_peaks(
        frame, diameter=arguments.diameter, minmass=arguments.minmass
    )

    arguments.output_csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        arguments.output_csv,
        peaks,
        delimiter=",",
        header="y,x,mass",
        comments="",
        fmt="%.4f",
    )
    if arguments.overlay_png is not None:
        save_overlay_figure(
            frame,
            peaks,
            condition_name,
            arguments.frame_index,
            arguments.overlay_png,
        )
    print(
        f"Detected {len(peaks)} peaks in {condition_name} "
        f"frame {arguments.frame_index} -> {arguments.output_csv}",
        flush=True,
    )
