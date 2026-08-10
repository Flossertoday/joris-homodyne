"""Paper-style rendering of the cross-mesa g2(0)(t) (G9/§8.3 follow-up).

Purple = LR (B1 x B4), green = RL (B2 x B3) — two independent cross-arm
measurements of the same condensate-1 x condensate-2 correlation; shaded
bands are +/-1 plug-in Poisson SEM around the g2 = 1 reference (the
experiment's resolution; see render_fig2c.py for the statistics-limited
display policy). Layout mirrors the Fig. 2(c) replica panel.
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2c_g2_estimator as est  # noqa: E402
from render_fig2c import WINDOW_PS, Y_LIMITS  # noqa: E402

DEFAULT_OUT_DIR = Path("drafts/26-08-11_g2_cross_mesa")


def draw_cross_panel(axis, t_ps, g2_lr, sem_lr, g2_rl, sem_rl,
                     window=WINDOW_PS, y_limits=Y_LIMITS, pooled=None,
                     title_suffix=""):
    """Cross-mesa g2(0)(t) panel; same display policy as render_fig2c."""
    mask = (t_ps >= window[0]) & (t_ps <= window[1])
    t = t_ps[mask]

    def clip(values):
        return np.clip(values[mask], y_limits[0], y_limits[1])

    for g2, sem, color, label in (
        (g2_lr, sem_lr, "tab:purple", "LR = B1\u00d7B4"),
        (g2_rl, sem_rl, "tab:green", "RL = B2\u00d7B3"),
    ):
        band_lo = np.maximum(1.0 - sem[mask], y_limits[0])
        band_hi = np.minimum(1.0 + sem[mask], y_limits[1])
        axis.fill_between(t, band_lo, band_hi, color=color, alpha=0.18, lw=0)
        axis.plot(t, clip(g2), ".", color=color, ms=3.5, alpha=0.75,
                  label=f"$g^{{(2)}}(0)$ {label}")
    axis.axhline(1.0, color="black", lw=0.9, ls="--")
    axis.set_ylim(*y_limits)
    axis.set_xlim(window)
    axis.set_xlabel("time relative to reference pulse [ps]")
    axis.set_ylabel("cross-mesa $g^{(2)}(0)$")
    sem_med_lr = float(np.nanmedian(sem_lr[mask]))
    sem_med_rl = float(np.nanmedian(sem_rl[mask]))
    note = (
        f"shaded: plug-in ±1 SEM around the Poisson reference "
        f"(median {sem_med_lr:.1f} LR / {sem_med_rl:.1f} RL); per-bin points "
        f"clipped to {list(y_limits)};\nall bins statistics-limited "
        f"(SEM > {est.MAX_DISPLAY_SEM:g})"
    )
    if pooled is not None:
        note += (
            f"\nwindow-pooled g2(0): LR = {pooled['LR'][0]:.2f} ± "
            f"{pooled['LR'][1]:.2f}, RL = {pooled['RL'][0]:.2f} ± "
            f"{pooled['RL'][1]:.2f} (jackknife SEM)"
        )
    axis.text(0.02, 0.02, note, transform=axis.transAxes, fontsize=7.5,
              va="bottom", ha="left",
              bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"})
    axis.legend(fontsize=8, loc="upper right")
    if title_suffix:
        axis.set_title(title_suffix, fontsize=9)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--width", type=int, default=4)
    parser.add_argument("--title", default="26-08-04 HBT")
    args = parser.parse_args()

    import json

    data = np.load(args.out_dir / f"g2_cross_results_bin{args.width}px.npz")
    t_ps = np.asarray(data["t_ps"], dtype=np.float64)
    g2_lr = np.asarray(data["g2_LR"], dtype=np.float64)
    g2_rl = np.asarray(data["g2_RL"], dtype=np.float64)
    sem_lr = np.asarray(data["sem_plugin_LR"], dtype=np.float64)
    sem_rl = np.asarray(data["sem_plugin_RL"], dtype=np.float64)

    pooled = None
    summary_path = args.out_dir / "g2_cross_summary.json"
    if summary_path.exists():
        with open(summary_path) as handle:
            pairs = json.load(handle)["per_width"][str(args.width)]["pairs"]
        pooled = {
            pair: (pairs[pair]["pooled_g2"], pairs[pair]["pooled_sem_jackknife"])
            for pair in ("LR", "RL")
        }

    figure, axis = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    draw_cross_panel(
        axis, t_ps, g2_lr, sem_lr, g2_rl, sem_rl, pooled=pooled,
        title_suffix=(
            f"{args.title} cross-mesa g2(0)(t), {args.width} px bins "
            f"({args.width * est.PS_PER_PX:.2f} ps), plug-in ±1 SEM"
        ),
    )
    out = args.out_dir / f"fig2c_cross_mesa_bin{args.width}px.png"
    figure.savefig(out, dpi=200)
    plt.close(figure)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
