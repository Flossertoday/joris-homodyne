"""Sample frame-overlay preview of the 26-08-04 HBT streak condition.

One annotated PNG per sampled frame of the 20000-frame
``hbt_condensate_josephson_regime`` condition:

- raw streak frame, sqrt display with the color minimum pinned at the
  sampled noise floor (~105-count camera pedestal), NOT at 0 counts — with
  vmin at 0 the pedestal would sit mid-scale and dim events would vanish;
- detected photon events, split into kept (centroid column inside one of
  the four band x-ROIs) and rejected (outside);
- the four band x-ROIs B1-B4 used by the Fig. 2(a) trace extraction
  (L = {B1, B3}, R = {B2, B4}), plus the event-detection column range.

Band windows, event threshold and detection column range are the final
values of the CFD50 Fig. 2(a) run (drafts/26-08-06_fig2a_cfd50/
bands_report.json, docs/26-08-04-fig2a.md). Event detection and the
measured background come from extract_hbt_traces.py so the overlay cannot
drift away from the reviewed pipeline. Figure style follows the 26-07-31
group-overview overlays (drafts/26-07-31_central_dot_0p18_group_overview/
visualize_group_overview.py).

Frames are sampled evenly spaced over the whole run (default 20), so the
preview also covers the synchroscan phase wander. The vertical axis is the
streak time axis (1 ns window over 1016 px, metadata-conventions.toml;
row 0 rendered as the window start at 0 ps, in-data direction calibration
pending).
"""

import argparse
import sys
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from extract_hbt_traces import (  # noqa: E402
    DEFAULT_H5,
    EVENT_COL_RANGE,
    FRAMES_DATASET,
    N_COLS,
    PS_PER_PX,
    detect_events,
    load_background,
)

DEFAULT_OUT_DIR = SCRIPT_DIR.parents[1] / "hbt-time-reference-0804-preview" / "frame_overlays"

# Final band x-ROIs and event threshold of the CFD50 Fig. 2(a) run
# (drafts/26-08-06_fig2a_cfd50/bands_report.json; docs/26-08-04-fig2a.md).
# L/R is the mesa pairing: L = {B1, B3}, R = {B2, B4}.
BAND_WINDOWS = ((380, 450), (570, 680), (810, 880), (1010, 1110))
BAND_LABELS = ("B1 (L)", "B2 (R)", "B3 (L)", "B4 (R)")
EVENT_THRESHOLD = 25.0  # counts above the measured background, >= 2 px blobs

N_SAMPLES = 20
WINDOW_PS = 1000.0
# Color minimum: this percentile of the pooled sampled-frame pixels anchors
# the noise floor (pedestal ~105 counts), so the display never starts at 0.
VMIN_PERCENTILE = 1.0

KEPT_COLOR = "deepskyblue"
REJECTED_COLOR = "0.65"
BAND_COLOR = "deepskyblue"


def sample_frame_indices(n_frames: int, n_samples: int) -> np.ndarray:
    """Evenly spaced frame indices covering the whole run."""
    if n_samples < 1:
        raise ValueError(f"Expected n_samples >= 1, got {n_samples}")
    indices = np.unique(np.rint(np.linspace(0, n_frames - 1, n_samples)).astype(np.int64))
    return indices


def sample_display_range(frames: h5py.Dataset, indices: np.ndarray) -> tuple[float, float]:
    """(vmin, vmax) for the shared display: noise-floor percentile and max.

    vmin is the VMIN_PERCENTILE percentile of the pooled sampled-frame
    pixels, which lands on the low tail of the ~105-count camera pedestal
    (~97 counts) — deliberately not 0. vmax is the sampled global maximum.
    """
    pooled = np.concatenate([np.asarray(frames[int(i)]).ravel() for i in indices])
    vmin = float(np.percentile(pooled, VMIN_PERCENTILE))
    vmax = float(pooled.max())
    if not vmin < vmax:
        raise ValueError(f"Expected vmin < vmax, got {vmin} >= {vmax}")
    return vmin, vmax


def band_membership(cols: np.ndarray) -> np.ndarray:
    """True where the event centroid column falls inside any band x-ROI."""
    member = np.zeros(cols.shape, dtype=bool)
    for x0, x1 in BAND_WINDOWS:
        member |= (cols >= x0) & (cols < x1)
    return member


def draw_bands(axis) -> None:
    """Shade the four band x-ROIs and label them with the L/R pairing."""
    for (x0, x1), label in zip(BAND_WINDOWS, BAND_LABELS):
        axis.axvspan(x0, x1, color=BAND_COLOR, alpha=0.10)
        axis.text(
            0.5 * (x0 + x1),
            0.015 * WINDOW_PS,
            label,
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=8,
            color="white",
        )


def render_frame_overlay(
    frame: np.ndarray,
    frame_index: int,
    rows: np.ndarray,
    cols: np.ndarray,
    norm: PowerNorm,
    vmin: float,
    vmax: float,
    output_path: Path,
    run_label: str = "hbt_condensate_josephson_regime",
) -> None:
    """One full-resolution overlay PNG: frame + events + band x-ROIs."""
    kept = band_membership(cols)
    figure, axis = plt.subplots(figsize=(9.4, 7.0), constrained_layout=True)
    image = axis.imshow(
        frame,
        origin="lower",
        aspect="auto",
        cmap="inferno",
        interpolation="nearest",
        extent=(0.0, float(N_COLS), 0.0, WINDOW_PS),
        norm=norm,
    )
    axis.scatter(
        cols[~kept],
        rows[~kept] * PS_PER_PX,
        s=25,
        marker="x",
        c=REJECTED_COLOR,
        linewidths=0.7,
    )
    axis.scatter(
        cols[kept],
        rows[kept] * PS_PER_PX,
        s=60,
        facecolors="none",
        edgecolors=KEPT_COLOR,
        linewidths=1.0,
    )
    draw_bands(axis)
    for boundary in EVENT_COL_RANGE:
        axis.axvline(boundary, color="0.55", lw=0.7, ls="--")
    axis.set_xlim(0, N_COLS)
    axis.set_ylim(0, WINDOW_PS)
    axis.set_title(
        f"{run_label} frame {frame_index:05d} — "
        f"{len(rows)} events (>{EVENT_THRESHOLD:g} above bg), "
        f"{int(kept.sum())} in band ROIs\n"
        f"(sqrt display over [{vmin:.0f}, {vmax:.0f}] counts, "
        f"vmin = p{VMIN_PERCENTILE:g} noise floor, not 0)",
        fontsize=10,
    )
    axis.set_xlabel("Detector axis 2 [pixel]")
    axis.set_ylabel("Time [ps]")
    legend_handles = [
        Line2D([], [], marker="x", color=REJECTED_COLOR, ls="none",
               label=f"rejected ({int((~kept).sum())})"),
        Line2D([], [], marker="o", markerfacecolor="none",
               markeredgecolor=KEPT_COLOR, ls="none",
               label=f"kept ({int(kept.sum())})"),
        Patch(facecolor=BAND_COLOR, alpha=0.25,
              label="band x-ROIs (L: B1, B3; R: B2, B4)"),
        Line2D([], [], color="0.55", lw=0.7, ls="--",
               label=f"event-detection range {EVENT_COL_RANGE[0]}–{EVENT_COL_RANGE[1]}"),
    ]
    axis.legend(handles=legend_handles, loc="upper right", fontsize=8)
    # Explicit colorbar ticks: the default locator labels the bar bottom "0"
    # even with a nonzero vmin, contradicting the noise-floor mapping note.
    ticks = [vmin] + list(np.arange(500.0, vmax, 500.0))
    colorbar = figure.colorbar(image, ax=axis, pad=0.02, shrink=0.9, ticks=ticks)
    colorbar.ax.set_yticklabels([f"{v:.0f}" for v in ticks])
    colorbar.set_label(
        "Camera counts [a.u.]\n"
        f"(inferno, sqrt color scale, [{vmin:.0f}, {vmax:.0f}];\n"
        f"vmin = p{VMIN_PERCENTILE:g} of sampled frames, not 0)"
    )
    figure.savefig(
        output_path,
        dpi=200,
        metadata={
            "Description": (
                f"Streak frame {frame_index} of the 26-08-04 "
                "hbt_condensate_josephson_regime condition with detected "
                "photon events (kept vs rejected by the band x-ROIs) and "
                "the four band x-ROIs B1-B4 of the Fig. 2(a) extraction."
            )
        },
    )
    plt.close(figure)


def write_readme(out_dir: Path, indices: np.ndarray, n_frames: int, vmin: float, vmax: float) -> None:
    """Record the exact sampling, intensity mapping and selection rules."""
    lines = [
        "# 26-08-04 frame-overlay sample preview",
        "",
        "One annotated PNG per sampled frame of the 20000-frame",
        "`hbt_condensate_josephson_regime` condition, regenerated by",
        "`analysis/26-08-04-derived-hbt-time-reference/render_frame_overlays.py`.",
        "Style follows the 26-07-31 group-overview overlays",
        "(`drafts/26-07-31_central_dot_0p18_group_overview/`).",
        "",
        "## Sampling",
        "",
        f"{indices.size} frames, evenly spaced over the run of {n_frames:,}:",
        ", ".join(str(int(i)) for i in indices) + ".",
        "",
        "## Intensity mapping",
        "",
        "Raw camera counts, inferno, sqrt color scale (gamma 0.5), shared over the",
        f"sample: vmin = {vmin:.1f} = p{VMIN_PERCENTILE:g} of the pooled sampled-frame pixels,",
        f"vmax = {vmax:.0f} = sampled global maximum. vmin is pinned at the noise",
        "floor, NOT at 0: the frames carry a ~105-count camera pedestal, and a",
        "zero-based scale would lift the pedestal to mid-scale and hide dim events.",
        "Each figure carries the colorbar with this mapping. The vertical axis is",
        "the streak time axis (1 ns window over 1016 px, metadata-conventions.toml;",
        "row 0 rendered as the window start at 0 ps; in-data direction calibration",
        "pending).",
        "",
        "## Event selection",
        "",
        "Events are detected by `extract_hbt_traces.detect_events` on the",
        "background-subtracted frame (measured `StreakCamBackG` of the",
        "`background_and_streak_pulses_delay_calibration` condition): 8-connected",
        f"components of >= 2 pixels above {EVENT_THRESHOLD:g} counts over background,",
        f"count-weighted subpixel centroids, columns {EVENT_COL_RANGE[0]}–{EVENT_COL_RANGE[1]}",
        "(dashed lines). This is the final event selection of the CFD50 Fig. 2(a)",
        "run (`drafts/26-08-06_fig2a_cfd50/bands_report.json`, `docs/26-08-04-fig2a.md`).",
        "",
        "## Band x-ROIs",
        "",
        "Shaded vertical bands are the four Fig. 2(a) extraction ROIs (2 mesas x 2",
        "HBT arms): B1 380–450, B2 570–680, B3 810–880, B4 1010–1110 detector",
        "pixels, with the mesa pairing L = {B1, B3}, R = {B2, B4}. An event is",
        "kept (open circle) when its centroid column falls inside a band, and",
        "rejected (x) otherwise.",
        "",
        "Regenerate (from the repository root):",
        "",
        "    .venv/Scripts/python.exe analysis/26-08-04-derived-hbt-time-reference/render_frame_overlays.py",
        "",
    ]
    (out_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def run(args) -> None:
    out_dir = Path(args.out_dir)
    if out_dir.exists() and any(out_dir.glob("frame_*.png")):
        raise FileExistsError(f"Refusing to write into a non-empty overlays directory: {out_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    background = load_background(args.h5)
    with h5py.File(args.h5, "r") as handle:
        frames = handle[FRAMES_DATASET]
        n_frames = int(frames.shape[0])
        indices = sample_frame_indices(n_frames, args.n_samples)
        vmin, vmax = sample_display_range(frames, indices)
        print(
            f"Rendering {indices.size} of {n_frames} frames; "
            f"display range [{vmin:.1f}, {vmax:.0f}] counts",
            flush=True,
        )
        norm = PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax)
        for out_number, index in enumerate(indices):
            raw_frame = np.asarray(frames[int(index)])
            frame_bg = raw_frame.astype(np.float64) - background
            rows, cols = detect_events(frame_bg, args.threshold)
            render_frame_overlay(
                raw_frame,
                int(index),
                rows,
                cols,
                norm,
                vmin,
                vmax,
                out_dir / f"frame_{int(index):05d}.png",
            )
            print(f"  {out_number + 1}/{indices.size} frame {int(index):05d}", flush=True)
    write_readme(out_dir, indices, n_frames, vmin, vmax)
    print(f"Wrote overlays and README to {out_dir}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Sample frame-overlay preview of the 26-08-04 HBT condition: raw "
            "frame (sqrt display, vmin at the noise floor, not 0), detected "
            "events (kept/rejected by the band x-ROIs), and the four band "
            "x-ROIs B1-B4."
        )
    )
    parser.add_argument("--h5", type=Path, default=DEFAULT_H5)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--n-samples", type=int, default=N_SAMPLES,
                        help=f"evenly spaced frames to render (default {N_SAMPLES})")
    parser.add_argument("--threshold", type=float, default=EVENT_THRESHOLD,
                        help=f"event threshold in counts above background "
                             f"(default {EVENT_THRESHOLD:g}, the CFD50 Fig. 2(a) value)")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
