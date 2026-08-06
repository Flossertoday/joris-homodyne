"""G5: paper-style rendering of the 26-08-04 Fig. 2(c) g2(0)(t) replica.

Blue = left mesa (B1 x B3), red = right mesa (B2 x B4); shaded bands are
+/-1 plug-in Poisson SEM (validated against Monte Carlo in the G3 gate; the
paper-style sample-moment SEM degenerates at this occupancy). Linear y axis,
g2 = 1 reference line. Bins whose plug-in SEM exceeds MAX_DISPLAY_SEM are
marked statistics-limited (the analogue of the paper's SNR<20 gray band);
at the 26-08-04 occupancies that is essentially the whole window, which the
panel states explicitly instead of pretending a per-bin measurement.
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

DEFAULT_OUT_DIR = Path("drafts/26-08-06_fig2c")
WINDOW_PS = (0.0, 820.0)
Y_LIMITS = (0.0, 4.0)


def draw_fig2c_panel(axis, t_ps, g2_l, sem_l, g2_r, sem_r, window=WINDOW_PS,
                     y_limits=Y_LIMITS, pooled=None, title_suffix=""):
    """Paper-style g2(0)(t) panel: L blue, R red, +/-1 SEM bands, g2=1 line."""
    mask = (t_ps >= window[0]) & (t_ps <= window[1])
    t = t_ps[mask]

    def clip(values):
        return np.clip(values[mask], y_limits[0], y_limits[1])

    for g2, sem, color, label in (
        (g2_l, sem_l, "tab:blue", "L (B1$\\times$B3)"),
        (g2_r, sem_r, "tab:red", "R (B2$\\times$B4)"),
    ):
        sem_w = sem[mask]
        limited = sem_w > est.MAX_DISPLAY_SEM
        axis.plot(t, clip(g2), color=color, lw=0.9, label=f"$g^{{(2)}}_{label}(0)$")
        axis.fill_between(
            t, clip(g2 - sem), clip(g2 + sem), color=color, alpha=0.22, lw=0,
        )
        if np.any(limited):
            axis.fill_between(
                t[limited], y_limits[0], y_limits[0] + 0.06 * np.ptp(y_limits),
                color="gray", alpha=0.5, lw=0,
            )
    axis.axhline(1.0, color="black", lw=0.9, ls="--")
    axis.set_ylim(*y_limits)
    axis.set_xlim(window)
    axis.set_xlabel("time relative to reference pulse [ps]")
    axis.set_ylabel("$g^{(2)}(0)$")
    frac_limited = float(np.mean(sem_l[mask] > est.MAX_DISPLAY_SEM))
    note = (
        f"gray strip: plug-in SEM > {est.MAX_DISPLAY_SEM:g} (statistics-limited): "
        f"{frac_limited:.0%} of bins"
    )
    if pooled is not None:
        note += (
            f"\npooled g2(0): L = {pooled['L'][0]:.2f} ± {pooled['L'][1]:.2f}, "
            f"R = {pooled['R'][0]:.2f} ± {pooled['R'][1]:.2f} (±1 jackknife SEM)"
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
    args = parser.parse_args()

    data = np.load(args.out_dir / f"g2_results_bin{args.width}px.npz")
    t_ps = np.asarray(data["t_ps"], dtype=np.float64)
    g2_l = np.asarray(data["g2_L"], dtype=np.float64)
    g2_r = np.asarray(data["g2_R"], dtype=np.float64)
    sem_l = np.asarray(data["sem_plugin_L"], dtype=np.float64)
    sem_r = np.asarray(data["sem_plugin_R"], dtype=np.float64)

    figure, axis = plt.subplots(figsize=(9.5, 5.2), constrained_layout=True)
    draw_fig2c_panel(
        axis, t_ps, g2_l, sem_l, g2_r, sem_r,
        title_suffix=(
            f"26-08-04 HBT g2(0)(t), {args.width} px bins "
            f"({args.width * est.PS_PER_PX:.2f} ps), plug-in ±1 SEM"
        ),
    )
    out = args.out_dir / f"fig2c_replica_bin{args.width}px.png"
    figure.savefig(out, dpi=200)
    plt.close(figure)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
