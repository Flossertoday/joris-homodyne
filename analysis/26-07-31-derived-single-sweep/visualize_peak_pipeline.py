"""Visualize each step of the Crocker-Grier peak-finding pipeline on one frame.

Audit/review companion to ``find_frame_peaks.py``: it reproduces the exact
operations ``tp.locate`` performs (bandpass via ``trackpy.preprocessing``,
candidate maxima via ``tp.locate`` with minmass=0 vs the production minmass,
and the iterative centroid refinement with the same circular mask, shift
threshold and iteration cap as ``trackpy.refine.center_of_mass``) and renders
two figures:

- ``--steps-png``: full-frame panels (a) raw, (b) gaussian lowpass, (c) boxcar
  background, (d) bandpass result, (e) candidate maxima with the minmass cut,
  (f) refined subpixel centroids overlaid on the raw frame.
- ``--zoom-png``: per-pixel view of the iterative centroid for two
  illustrative peaks (brightest, and largest centroid-vs-argmax offset).

One deliberate simplification: ``tp.locate`` internally refines the centroid
on a uint8-normalized copy of the bandpassed image (performance optimization);
this script visualizes the refinement on the float bandpassed image. The
quantization shifts centroids by <0.02 px and does not change any conclusion.

Reproduce the documented frame-230 figures with:

    python visualize_peak_pipeline.py <h5> --central-dot-strength 0.18 \\
        --frame-index 230 --steps-png <steps.png> --zoom-png <zoom.png>
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Circle
import trackpy as tp
from trackpy.preprocessing import bandpass, lowpass, boxcar
from trackpy.masks import binary_mask

# Same detector-calibrated defaults as the production peak finder.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_frame_peaks import (  # noqa: E402
    DEFAULT_DIAMETER,
    DEFAULT_MINMASS,
    DEFAULT_PERCENTILE,
    load_condition_frame,
)

# trackpy.refine.center_of_mass defaults: iterate at most 10 times and stop
# once the centroid is within 0.6 px of the window centre.
SHIFT_THRESH = 0.6
MAX_ITERATIONS = 10


def iterate_centroid(image: np.ndarray, start: np.ndarray, radius: int):
    """Replicate the center_of_mass refinement loop on a bandpassed image.

    Returns the final subpixel centroid and the list of integer window
    centres visited (length-1 = number of window shifts).
    """
    mask = binary_mask(radius, image.ndim)
    ogrid = [g.astype(float)
             for g in np.ogrid[[slice(0, i) for i in mask.shape]]]
    coord = np.asarray(start, dtype=int).copy()
    trajectory = [tuple(coord)]
    for _ in range(MAX_ITERATIONS):
        rect = tuple(slice(c - radius, c + radius + 1) for c in coord)
        neighborhood = mask * image[rect]
        normalizer = neighborhood.sum()
        cm_n = np.array([(neighborhood * g).sum() / normalizer for g in ogrid])
        cm_i = cm_n - radius + coord
        off = cm_n - radius
        if np.all(np.abs(off) < SHIFT_THRESH):
            break
        coord[off > SHIFT_THRESH] += 1
        coord[off < -SHIFT_THRESH] -= 1
        trajectory.append(tuple(coord))
    return cm_i, trajectory


def local_max_pixel(image: np.ndarray, y: float, x: float,
                    radius: int) -> np.ndarray:
    """Integer pixel of the bandpass local maximum nearest to (y, x)."""
    iy, ix = int(round(y)), int(round(x))
    win = image[iy - radius:iy + radius + 1, ix - radius:ix + radius + 1]
    ay, ax = np.unravel_index(np.argmax(win), win.shape)
    return np.array([iy - radius + ay, ix - radius + ax])


def save_steps_figure(frame: np.ndarray, low: np.ndarray, bg: np.ndarray,
                      bp: np.ndarray, all_cand, kept, rejected,
                      condition_name: str, frame_index: int, minmass: float,
                      output_path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16.5, 9.6), constrained_layout=True)

    ax = axes[0, 0]
    lo, hi = np.percentile(frame, [1, 99.9])
    ax.imshow(frame, origin="upper", cmap="magma", interpolation="nearest",
              vmin=lo, vmax=hi)
    ax.set_title("(a) raw frame\n(flat ~105-count background)")

    ax = axes[0, 1]
    lo, hi = np.percentile(low, [1, 99.9])
    ax.imshow(low, origin="upper", cmap="magma", interpolation="nearest",
              vmin=lo, vmax=hi)
    ax.set_title("(b) lowpass: gaussian $\\sigma$=1 px\n"
                 "(suppresses 1-px impulse noise)")

    ax = axes[0, 2]
    lo, hi = np.percentile(bg, [1, 99.9])
    im = ax.imshow(bg, origin="upper", cmap="magma", interpolation="nearest",
                   vmin=lo, vmax=hi)
    fig.colorbar(im, ax=ax, shrink=0.85, label="counts")
    ax.set_title(f"(c) background: 7$\\times$7 rolling mean\n"
                 f"(mean {bg.mean():.1f} counts, spot ghosts visible)")

    ax = axes[1, 0]
    lo, hi = np.percentile(bp[bp > 0], [1, 99.9]) if np.any(bp > 0) else (0, 1)
    ax.imshow(bp, origin="upper", cmap="magma", interpolation="nearest",
              vmin=0, vmax=hi)
    ax.set_title("(d) bandpass = (b) - (c), clipped $\\geq$0\n"
                 "(spots on zero background)")

    ax = axes[1, 1]
    ax.imshow(bp, origin="upper", cmap="magma", interpolation="nearest",
              vmin=0, vmax=hi)
    ax.plot(rejected["x"], rejected["y"], "x", color="red", ms=3, mew=0.6,
            alpha=0.6, label=f"rejected by minmass ({len(rejected)})")
    ax.plot(kept["x"], kept["y"], "o", mfc="none", mec="cyan", mew=1.0, ms=11,
            label=f"kept ({len(kept)})")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(f"(e) local maxima on (d): {len(all_cand)} candidates\n"
                 f"minmass>{minmass:.0f} keeps the real spots")

    ax = axes[1, 2]
    lo, hi = np.percentile(frame, [1, 99.9])
    ax.imshow(frame, origin="upper", cmap="magma", interpolation="nearest",
              vmin=lo, vmax=hi)
    ax.plot(kept["x"], kept["y"], "o", mfc="none", mec="cyan", mew=1.0, ms=12)
    ax.set_title("(f) refined subpixel centroids on raw frame\n"
                 f"({len(kept)} peaks; circle size is display-only)")

    for ax in axes.flat:
        ax.set_xlabel("Detector axis 2 [pixel]")
        ax.set_ylabel("Detector axis 1 [pixel]")
    fig.suptitle(f"{condition_name} frame {frame_index} — Crocker-Grier "
                 "pipeline steps (trackpy v0.7)")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print("saved", output_path)


def save_zoom_figure(bp: np.ndarray, kept, radius: int, condition_name: str,
                     frame_index: int, output_path: Path) -> None:
    offsets = []
    for _, p in kept.iterrows():
        lm = local_max_pixel(bp, p["y"], p["x"], radius)
        offsets.append(np.hypot(lm[0] - p["y"], lm[1] - p["x"]))
    offsets = np.array(offsets)
    chosen = [
        ("brightest peak", kept.iloc[int(np.argmax(kept["mass"].to_numpy()))]),
        ("largest centroid-vs-max offset", kept.iloc[int(np.argmax(offsets))]),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 6.8), constrained_layout=True)
    for ax, (tag, p) in zip(axes, chosen):
        start = local_max_pixel(bp, p["y"], p["x"], radius)
        cm, trajectory = iterate_centroid(bp, start, radius)
        r = 7
        y0, y1 = int(cm[0]) - r, int(cm[0]) + r + 1
        x0, x1 = int(cm[1]) - r, int(cm[1]) + r + 1
        ax.imshow(bp[y0:y1, x0:x1], origin="upper", cmap="magma",
                  interpolation="nearest",
                  extent=[x0 - 0.5, x1 - 0.5, y1 - 0.5, y0 - 0.5])
        ax.set_xticks(np.arange(x0, x1))
        ax.set_yticks(np.arange(y0, y1))
        ax.grid(color="gray", lw=0.3, alpha=0.6)
        for yy in range(y0, y1):
            for xx in range(x0, x1):
                ax.text(xx, yy, f"{bp[yy, xx]:.0f}", color="w", fontsize=5.5,
                        ha="center", va="center", alpha=0.85)
        ax.add_patch(Circle((start[1], start[0]), radius + 0.5, fill=False,
                            color="gray", ls="--", lw=1.0,
                            label=f"mask window (r={radius}) at start"))
        for (py0, px0), (py1, px1) in zip(trajectory[:-1], trajectory[1:]):
            ax.annotate("", xy=(px1, py1), xytext=(px0, py0),
                        arrowprops=dict(arrowstyle="->", color="lime", lw=1.2))
        ax.plot(start[1], start[0], "rx", ms=9, mew=1.5,
                label=f"bandpass local max ({start[0]}, {start[1]})")
        ax.plot(cm[1], cm[0], "c+", ms=11, mew=1.6,
                label=f"converged centroid ({cm[0]:.2f}, {cm[1]:.2f})")
        ax.plot(p["x"], p["y"], "o", mfc="none", mec="cyan", mew=1.0, ms=12,
                alpha=0.8, label="overlay circle in panel (f)")
        ax.legend(loc="upper right", fontsize=7)
        ax.set_title(f"{tag}: mass={p['mass']:.0f}, "
                     f"{len(trajectory) - 1} window shift(s)")
        ax.set_xlabel("Detector axis 2 [pixel]")
        ax.set_ylabel("Detector axis 1 [pixel]")

    fig.suptitle(f"{condition_name} frame {frame_index} — centroid refinement "
                 "detail (numbers = bandpassed counts)")
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print("saved", output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Render the Crocker-Grier pipeline steps of one frame."
    )
    parser.add_argument("input_h5", type=Path)
    parser.add_argument("--central-dot-strength", type=float, required=True)
    parser.add_argument("--frame-index", type=int, default=230,
                        help="frame of the documented figures (default 230)")
    parser.add_argument("--diameter", type=int, default=DEFAULT_DIAMETER)
    parser.add_argument("--minmass", type=float, default=DEFAULT_MINMASS)
    parser.add_argument("--percentile", type=float, default=DEFAULT_PERCENTILE)
    parser.add_argument("--steps-png", type=Path, default=None)
    parser.add_argument("--zoom-png", type=Path, default=None)
    arguments = parser.parse_args()

    radius = arguments.diameter // 2
    condition_name, frame = load_condition_frame(
        arguments.input_h5,
        arguments.central_dot_strength,
        arguments.frame_index,
    )
    frame = np.asarray(frame, dtype=np.float64)

    # Step 1: the exact bandpass components tp.locate computes
    # (noise_size=1 gaussian lowpass, smoothing_size=diameter boxcar).
    low = lowpass(frame, 1)
    bg = boxcar(frame, arguments.diameter)
    bp = bandpass(frame, 1, arguments.diameter)

    # Step 2: candidate maxima without and with the production minmass cut.
    all_cand = tp.locate(frame, diameter=arguments.diameter, minmass=0,
                         percentile=arguments.percentile)
    kept = tp.locate(frame, diameter=arguments.diameter,
                     minmass=arguments.minmass,
                     percentile=arguments.percentile)
    rejected = all_cand[all_cand["mass"] <= arguments.minmass]
    print(f"{condition_name} frame {arguments.frame_index}: "
          f"{len(all_cand)} candidates, {len(rejected)} rejected by minmass, "
          f"{len(kept)} kept")

    if arguments.steps_png is not None:
        save_steps_figure(frame, low, bg, bp, all_cand, kept, rejected,
                          condition_name, arguments.frame_index,
                          arguments.minmass, arguments.steps_png)
    if arguments.zoom_png is not None:
        save_zoom_figure(bp, kept, radius, condition_name,
                         arguments.frame_index, arguments.zoom_png)
