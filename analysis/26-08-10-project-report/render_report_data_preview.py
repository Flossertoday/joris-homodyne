"""Render Figure: data preview — raw streak frame, t0-aligned event accumulation,
and column profile for the 26-08-04 HBT time-reference batch.

Output: analysis/26-08-10-project-report/figures/fig_data_preview.png
"""

import sys
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import patches as mpatches
from matplotlib import pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.patheffects import withStroke

# Make project analysis modules importable from the worktree.
WORKTREE_ROOT = Path(__file__).resolve().parent.parent.parent
ANALYSIS_DIR = WORKTREE_ROOT / "analysis" / "26-08-04-derived-hbt-time-reference"
sys.path.insert(0, str(ANALYSIS_DIR))
import extract_hbt_traces as eht  # noqa: E402

DATA_H5 = Path(
    "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/"
    "26-08-04_hbt_time_reference.h5"
)
EVENTS_H5 = ANALYSIS_DIR / "detected-peaks" / "26-08-04_hbt_time_reference_events.h5"
OUT_PATH = Path(__file__).resolve().parent / "figures" / "fig_data_preview.png"

FRAME_INDEX = 2414
PS_PER_PX = 0.984251968503937  # 1 ns / 1016 px, agreed convention

BANDS = {
    "B1": (380, 450),
    "B2": (570, 680),
    "B3": (810, 880),
    "B4": (1010, 1110),
}
BAND_CENTERS = {"B1": 415, "B2": 625, "B3": 845, "B4": 1060}
# Project pairing: Left = B1 + B3, Right = B2 + B4.
LR_OF_BAND = {"B1": "L", "B2": "R", "B3": "L", "B4": "R"}
LR_COLORS = {"L": "tab:blue", "R": "tab:red"}
DETECTION_COLS = (350, 1150)
REFERENCE_ROI = (110, 210, 170, 215)


def load_frame_and_background():
    with h5py.File(DATA_H5, "r") as h:
        background = np.asarray(h[eht.BACKGROUND_DATASET], dtype=np.float64)
        frame = np.asarray(h[eht.FRAMES_DATASET][FRAME_INDEX], dtype=np.float64)
    return frame - background, background


def load_events():
    with h5py.File(EVENTS_H5, "r") as h:
        g = h["events"]
        x = np.asarray(g["x"], dtype=np.float64)
        tau = np.asarray(g["tau"], dtype=np.float64)
    return x, tau


def stroke_text(axis, x, y, s, **kwargs):
    """Text with a thin black outline for contrast on variable backgrounds."""
    defaults = {"fontsize": 10, "fontweight": "bold", "color": "white",
                "ha": "center", "va": "top"}
    defaults.update(kwargs)
    return axis.text(
        x, y, s,
        path_effects=[withStroke(linewidth=2, foreground="black")],
        **defaults,
    )


def add_band_shades(axis, y_frac=0.94, alpha=0.18):
    """ translucent vertical bands + top labels for B1-B4, L/R color scheme. """
    for name, (x0, x1) in BANDS.items():
        color = LR_COLORS[LR_OF_BAND[name]]
        axis.axvspan(x0, x1, color=color, alpha=alpha, zorder=1)
        stroke_text(
            axis, 0.5 * (x0 + x1), y_frac, name,
            transform=axis.get_xaxis_transform(),
            fontsize=10, fontweight="bold", color=color, zorder=5,
        )


def render(out_path=OUT_PATH):
    frame_bg, background = load_frame_and_background()
    x_ev, tau_ev = load_events()

    fig = plt.figure(figsize=(16, 10.67), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], width_ratios=[1.35, 1.0])
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # ---- Panel (a): single realization ----
    det_region = frame_bg[:, DETECTION_COLS[0]:DETECTION_COLS[1]]
    vmax_a = float(np.percentile(det_region[det_region > 0], 99.9))
    im_a = ax_a.imshow(
        frame_bg,
        origin="upper",
        aspect="auto",
        cmap="magma",
        interpolation="nearest",
        extent=(0, frame_bg.shape[1], frame_bg.shape[0], 0),
        norm=PowerNorm(gamma=0.5, vmin=0, vmax=vmax_a),
    )
    cbar_a = fig.colorbar(im_a, ax=ax_a, pad=0.01, fraction=0.046)
    cbar_a.set_label(
        f"sqrt counts (magma, 0..p99.9 = {vmax_a:.0f})"
    )

    # Reference-pulse ROI rectangle (dashed)
    r0, r1, c0, c1 = REFERENCE_ROI
    rect = mpatches.Rectangle(
        (c0, r0), c1 - c0, r1 - r0,
        linewidth=1.6, edgecolor="cyan", facecolor="none", linestyle="--", zorder=4,
    )
    ax_a.add_patch(rect)
    ax_a.annotate(
        "laser reference pulse\n(per-frame time origin)",
        xy=(0.5 * (c0 + c1), r0), xycoords="data",
        xytext=(c1 + 60, r0 + 30), textcoords="data",
        arrowprops=dict(arrowstyle="->", color="cyan", lw=1.2),
        color="cyan", fontsize=9, va="top",
    )

    # Detection column boundaries (dotted)
    for xc in DETECTION_COLS:
        ax_a.axvline(xc, color="white", linestyle=":", lw=1.2, zorder=3)
    ax_a.text(
        DETECTION_COLS[0] - 10, 0.02, "detection\ncolumns",
        transform=ax_a.get_xaxis_transform(),
        ha="right", va="bottom", fontsize=8, color="white",
    )

    add_band_shades(ax_a, y_frac=0.07, alpha=0.22)
    ax_a.set_xlim(0, frame_bg.shape[1])
    ax_a.set_ylim(frame_bg.shape[0], 0)
    ax_a.set_xlabel("detector column [px]")
    ax_a.set_ylabel("streak row (time increasing downward)")
    ax_a.set_title(f"(a) one realization — frame {FRAME_INDEX}, background subtracted")

    note_a = (
        "sqrt scale, 0..p99.9 of detection columns; reference pulse saturates by design"
    )
    ax_a.text(
        0.97, 0.02, note_a,
        transform=ax_a.transAxes,
        ha="right", va="bottom", fontsize=8, color="white",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="black", alpha=0.5, edgecolor="none"),
    )

    # ---- Panel (b): t0-aligned accumulation ----
    tau_ps_ev = tau_ev * PS_PER_PX
    ylim_b = (-63, 860)
    mask_b = (tau_ps_ev >= ylim_b[0]) & (tau_ps_ev <= ylim_b[1])

    bins_x = np.arange(0, frame_bg.shape[1] + 1, 4)
    bins_y = np.linspace(ylim_b[0], ylim_b[1], 421)
    h2d, _, _ = np.histogram2d(x_ev[mask_b], tau_ps_ev[mask_b], bins=[bins_x, bins_y])
    vmax_b = float(np.percentile(h2d[h2d > 0], 99.5) if h2d.max() > 0 else 1)

    im_b = ax_b.imshow(
        h2d.T,
        origin="lower",
        aspect="auto",
        cmap="inferno",
        interpolation="nearest",
        extent=(bins_x[0], bins_x[-1], bins_y[0], bins_y[-1]),
        norm=PowerNorm(gamma=0.5, vmin=0, vmax=max(vmax_b, 1)),
    )
    cbar_b = fig.colorbar(im_b, ax=ax_b, pad=0.01, fraction=0.046)
    cbar_b.set_label(f"events per bin (sqrt, 0..p99.5 = {vmax_b:.0f})")

    # Band rectangles (white dashed), B1-B4 labels, and L/R pairing labels
    for name, (x0, x1) in BANDS.items():
        ax_b.axvspan(x0, x1, color="white", alpha=0.08, zorder=2)
        rect_b = mpatches.Rectangle(
            (x0, ylim_b[0] + 5), x1 - x0, ylim_b[1] - ylim_b[0] - 10,
            linewidth=1.4, edgecolor="white", facecolor="none", linestyle="--", zorder=3,
        )
        ax_b.add_patch(rect_b)
        ax_b.text(
            0.5 * (x0 + x1), ylim_b[1] - 35, name,
            ha="center", va="top", fontsize=10, fontweight="bold",
            color="white", zorder=4,
        )
        stroke_text(
            ax_b, 0.5 * (x0 + x1), ylim_b[1] - 58, LR_OF_BAND[name],
            fontsize=11, fontweight="bold", color="white", zorder=4,
        )

    for xc in DETECTION_COLS:
        ax_b.axvline(xc, color="white", linestyle=":", lw=1.0, zorder=3)

    # Mesa brackets
    def bracket(ax, x0, x1, y, label):
        ax.plot([x0, x1], [y, y], color="white", lw=1.5, zorder=4)
        ax.plot([x0, x0], [y - 4, y], color="white", lw=1.5, zorder=4)
        ax.plot([x1, x1], [y - 4, y], color="white", lw=1.5, zorder=4)
        ax.text(0.5 * (x0 + x1), y + 6, label, ha="center", va="bottom",
                color="white", fontsize=9, fontweight="bold", zorder=4)

    y_brack = ylim_b[1] - 12
    bracket(ax_b, BANDS["B1"][0], BANDS["B3"][1], y_brack, "Left = B1 + B3")
    bracket(ax_b, BANDS["B2"][0], BANDS["B4"][1], y_brack, "Right = B2 + B4")

    # Reference-pulse arrow annotation
    ax_b.annotate(
        "reference pulse\n(not in detection columns)",
        xy=(0.5 * (REFERENCE_ROI[2] + REFERENCE_ROI[3]), ylim_b[0] + 30),
        xytext=(REFERENCE_ROI[3] + 90, ylim_b[0] + 180),
        arrowprops=dict(arrowstyle="->", color="cyan", lw=1.2),
        color="cyan", fontsize=9, va="top",
    )

    ax_b.set_xlim(0, frame_bg.shape[1])
    ax_b.set_ylim(ylim_b)
    ax_b.set_xlabel("detector column [px]")
    ax_b.set_ylabel("time after reference pulse [ps]")
    ax_b.set_title("(b) t0-aligned accumulation — all 275 218 photon events")

    # ---- Panel (c): column profile ----
    mask_c = (tau_ev >= 0) & (tau_ev <= 820 / PS_PER_PX)
    counts_c, edges_c = np.histogram(x_ev[mask_c], bins=np.arange(0, frame_bg.shape[1] + 1))
    centers_c = 0.5 * (edges_c[:-1] + edges_c[1:])

    ax_c.fill_between(centers_c, counts_c, color="black", step="mid", alpha=0.85)
    for name, (x0, x1) in BANDS.items():
        color = LR_COLORS[LR_OF_BAND[name]]
        ax_c.axvspan(x0, x1, color=color, alpha=0.18, zorder=1)
        stroke_text(
            ax_c, 0.5 * (x0 + x1), ax_c.get_ylim()[1] * 0.97, name,
            fontsize=10, fontweight="bold", color=color,
        )
    # L/R pairing legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color=LR_COLORS["L"], lw=6, alpha=0.5, label="Left = B1 + B3"),
        Line2D([0], [0], color=LR_COLORS["R"], lw=6, alpha=0.5, label="Right = B2 + B4"),
    ]
    ax_c.legend(handles=legend_elements, loc="upper right", fontsize=9)

    # Annotate spacings
    ax_c.annotate(
        "",
        xy=(BAND_CENTERS["B2"], ax_c.get_ylim()[1] * 0.72),
        xytext=(BAND_CENTERS["B1"], ax_c.get_ylim()[1] * 0.72),
        arrowprops=dict(arrowstyle="<->", color="black", lw=1.2),
    )
    ax_c.text(
        0.5 * (BAND_CENTERS["B1"] + BAND_CENTERS["B2"]),
        ax_c.get_ylim()[1] * 0.78,
        "210 px ≈ mesa spacing\n(204 px in real space)",
        ha="center", va="bottom", fontsize=8,
    )

    ax_c.annotate(
        "",
        xy=(BAND_CENTERS["B3"], ax_c.get_ylim()[1] * 0.45),
        xytext=(BAND_CENTERS["B1"], ax_c.get_ylim()[1] * 0.45),
        arrowprops=dict(arrowstyle="<->", color="black", lw=1.2),
    )
    ax_c.text(
        0.5 * (BAND_CENTERS["B1"] + BAND_CENTERS["B3"]),
        ax_c.get_ylim()[1] * 0.52,
        "430 px = HBT arm separation",
        ha="center", va="bottom", fontsize=8,
    )

    ax_c.set_xlim(300, 1200)
    ax_c.set_xlabel("detector column [px]")
    ax_c.set_ylabel("events (0–820 ps)")
    ax_c.set_title("(c) column profile — events per column")

    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    render()
