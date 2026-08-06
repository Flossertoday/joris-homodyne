"""Additive-offset arm comparison figure for the 26-08-04 HBT Fig. 2(a) data.

Anti-phase diagnostic complement to render_fig2a.py: shifts the right-arm
trace by a constant c = mean(L) - mean(R) over the plotted window so both
condensates share the same window mean, then renders both arms on a zoomed
LINEAR intensity axis — full window plus a tighter zoom on the stationary
part (t >= --zoom-min) — with ±1 SEM bands and thin connecting lines, for
direct inspection of whether the two arms' up/down fluctuations are
anti-phase. Unlike multiplicative normalization (--normalize matched), the
additive shift preserves each arm's noise amplitude. Also prints the lag-0
Pearson correlation of the demeaned stationary-window fluctuations with a
two-sided permutation p-value.

Default output overwrites the tracked `fig2a_offset_zoom.png`/`.pdf` next to
this script. Reuses render_fig2a loading/rebin primitives.
"""

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import render_fig2a as rf  # noqa: E402

DEFAULT_NPZ = Path("drafts/26-08-06_fig2a_cfd50/hbt_traces.npz")
DEFAULT_OUT = SCRIPT_DIR / "fig2a_offset_zoom.png"


def permutation_pvalue(x: np.ndarray, y: np.ndarray, n_perm: int = 20000, seed: int = 7):
    """Two-sided permutation p-value for the lag-0 Pearson correlation."""
    r_obs = float(np.corrcoef(x, y)[0, 1])
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        r = float(np.corrcoef(x, rng.permutation(y))[0, 1])
        if abs(r) >= abs(r_obs):
            count += 1
    return r_obs, (count + 1) / (n_perm + 1)


def draw_panel(axis, time_ps, trace_l, sem_l, trace_r, sem_r, t_lo, t_hi, common_mean,
               show_legend: bool):
    """Draw both arms (R already offset) on a linear axis zoomed to [t_lo, t_hi]."""
    mask = (time_ps >= t_lo) & (time_ps <= t_hi)
    t, l, sl = time_ps[mask], trace_l[mask], sem_l[mask]
    r, sr = trace_r[mask], sem_r[mask]
    for values, errors, color, label in (
        (l, sl, rf.LEFT_COLOR, "Left (B1+B3)"),
        (r, sr, rf.RIGHT_COLOR, "Right (B2+B4) + const"),
    ):
        axis.fill_between(t, values - errors, values + errors, color=color, alpha=0.12,
                          linewidth=0, zorder=1)
        axis.plot(t, values, color=color, alpha=0.55, linewidth=0.8, zorder=2)
        axis.scatter(t, values, s=10, c=color, edgecolors="none", zorder=3, label=label)
    axis.axhline(common_mean, color="0.4", linewidth=0.7, linestyle="--", zorder=0)
    v_min = min(float(l.min()), float(r.min()))
    v_max = max(float(l.max()), float(r.max()))
    pad = 0.12 * (v_max - v_min)
    axis.set_ylim(v_min - pad, v_max + pad)
    axis.set_xlim(t[0], t[-1])
    axis.set_ylabel("Intensity (a.u.)", fontsize=8)
    axis.tick_params(labelsize=7)
    for spine in axis.spines.values():
        spine.set_linewidth(0.8)
    if show_legend:
        axis.legend(loc="lower right", fontsize=7, frameon=False, scatterpoints=1,
                    markerscale=2.0, handletextpad=0.4)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Additive-offset (constant-shift) arm comparison on a zoomed "
            "linear intensity axis, for visual anti-phase inspection."
        )
    )
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--source", choices=("photon", "analog"), default="photon")
    parser.add_argument("--bands-l", default="0,2")
    parser.add_argument("--bands-r", default="1,3")
    parser.add_argument("--rebin", type=int, default=4)
    parser.add_argument("--t-min", type=float, default=0.0)
    parser.add_argument("--t-max", type=float, default=820.0)
    parser.add_argument(
        "--zoom-min", type=float, default=150.0,
        help="Start of the stationary-part zoom panel (ps)",
    )
    arguments = parser.parse_args(argv)

    bands_l = rf.parse_band_selection(arguments.bands_l)
    bands_r = rf.parse_band_selection(arguments.bands_r)
    data = rf.load_band_traces(arguments.npz, arguments.source)
    trace_l, sem_l = rf.combine_bands(data["traces"], data["sems"], bands_l)
    trace_r, sem_r = rf.combine_bands(data["traces"], data["sems"], bands_r)
    time_win, trace_l, sem_l = rf.select_time_window(
        data["time_ps"], trace_l, sem_l, arguments.t_min, arguments.t_max
    )
    _, trace_r, sem_r = rf.select_time_window(
        data["time_ps"], trace_r, sem_r, arguments.t_min, arguments.t_max
    )
    time_ps, trace_l, sem_l = rf.rebin_trace(time_win, trace_l, sem_l, arguments.rebin)
    _, trace_r, sem_r = rf.rebin_trace(time_win, trace_r, sem_r, arguments.rebin)

    offset = float(trace_l.mean() - trace_r.mean())
    trace_r = trace_r + offset  # additive shift; sem_r unchanged
    common_mean = float(trace_l.mean())

    figure, (ax_full, ax_zoom) = plt.subplots(
        2, 1, figsize=(7.2, 6.0), constrained_layout=True
    )
    draw_panel(ax_full, time_ps, trace_l, sem_l, trace_r, sem_r,
               time_ps[0], time_ps[-1], common_mean, show_legend=True)
    ax_full.set_title(
        f"Additive offset: Right shifted by {offset:+.4g} so both arms share "
        f"the window mean; ±1 SEM bands",
        fontsize=8.5, loc="left",
    )
    draw_panel(ax_zoom, time_ps, trace_l, sem_l, trace_r, sem_r,
               arguments.zoom_min, time_ps[-1], common_mean, show_legend=False)
    ax_zoom.set_title(f"Stationary part, t ≥ {arguments.zoom_min:g} ps (tighter y)",
                      fontsize=8.5, loc="left")
    ax_zoom.set_xlabel("Time (ps)", fontsize=8)

    arguments.out.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(arguments.out, dpi=300)
    figure.savefig(arguments.out.with_suffix(".pdf"))
    plt.close(figure)

    stat = time_ps >= arguments.zoom_min
    fluct_l = trace_l[stat] - trace_l[stat].mean()
    fluct_r = trace_r[stat] - trace_r[stat].mean()
    r_obs, p_val = permutation_pvalue(fluct_l, fluct_r)

    bin_width = (time_ps[-1] - time_ps[0]) / max(len(time_ps) - 1, 1)
    print(
        f"source={arguments.source} rebin={arguments.rebin} ({bin_width:.2f} ps/bin), "
        f"window {time_ps[0]:.1f}–{time_ps[-1]:.1f} ps",
        flush=True,
    )
    print(
        f"means: L={trace_l.mean():.5g}, R(raw)={trace_r.mean() - offset:.5g}; "
        f"offset c={offset:+.5g} ({offset / common_mean:+.2%} of the common mean)",
        flush=True,
    )
    print(
        f"stationary (t>={arguments.zoom_min:g} ps, {stat.sum()} bins) demeaned "
        f"fluctuation correlation: r={r_obs:+.3f}, permutation p={p_val:.4f}",
        flush=True,
    )
    print(f"Wrote {arguments.out} and {arguments.out.with_suffix('.pdf')}", flush=True)


if __name__ == "__main__":
    main()
