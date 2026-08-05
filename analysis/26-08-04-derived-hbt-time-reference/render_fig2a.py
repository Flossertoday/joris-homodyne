"""Render a Fig. 2(a) replica (Adiyatullin et al. 2017) from extracted HBT band traces.

Reads the trace archive written by extract_hbt_traces.py
(drafts/26-08-04_fig2a/hbt_traces.npz), sums the selected detector bands into
left/right mesa intensities I_L(t) and I_R(t), and renders them in the paper's
panel-a style: small unconnected dots (left blue, right red) on a log10
intensity axis peak-normalized to a 10^0 top, with gray horizontal shading
over the intensity region where the per-bin signal-to-noise ratio
(trace / SEM) falls below the SNR threshold.

Rebinning merges `k` adjacent time bins as their count-weighted mean
(trace_merged = sum(trace) / k) with independent-bin SEM propagation
(sem_merged = sqrt(sum(sem^2)) / k), matching the paper's 2-streak-pixel time
bins; the merged bin's SNR is invariant to the sum-vs-mean convention.
Non-positive intensities are masked before log plotting.
"""

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.ticker import NullLocator

DEFAULT_NPZ = Path("drafts/26-08-04_fig2a/hbt_traces.npz")
DEFAULT_OUT = Path("drafts/26-08-04_fig2a/fig2a_replica.png")

LEFT_COLOR = "#1a5fb4"  # paper: left mesa in blue
RIGHT_COLOR = "#d62728"  # paper: right mesa in red
SHADE_COLOR = "0.85"
# Axis limits are fractions of the larger trace peak, so peak normalization
# (peak = 1) reproduces the paper's 10^0 top with 10^-2 near the bottom.
Y_TOP_HEADROOM = 1.4
Y_BOTTOM_MAX_FRAC = 8e-3
Y_BOTTOM_FLOOR_FRAC = 3e-3


def parse_band_selection(text: str) -> list[int]:
    """Parse a comma-separated band index list like '0,2'."""
    try:
        indices = [int(part.strip()) for part in text.split(",") if part.strip()]
    except ValueError as exc:
        raise ValueError(f"Invalid band list {text!r}: {exc}") from exc
    if not indices:
        raise ValueError(f"Band list {text!r} selects no bands")
    if len(set(indices)) != len(indices):
        raise ValueError(f"Band list {text!r} contains duplicate indices")
    return indices


def _scalar_text(value: np.ndarray) -> str:
    """Return a stored JSON string from a 0-d or bytes npz scalar."""
    if isinstance(value, np.ndarray) and value.shape == ():
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    return str(value)


def load_band_traces(npz_path: Path, source: str) -> dict:
    """Load time axis, band traces, and SEMs for `source` ('photon' or 'analog')."""
    with np.load(npz_path) as archive:
        time_ps = archive["time_ps"].astype(np.float64)
        traces = archive[f"band_traces_{source}"].astype(np.float64)
        sems = archive[f"band_sem_{source}"].astype(np.float64)
        meta = json.loads(_scalar_text(archive["meta_json"])) if "meta_json" in archive else {}
        n_frames_total = (
            int(archive["n_frames_total"]) if "n_frames_total" in archive else None
        )
    if traces.ndim != 2 or sems.shape != traces.shape:
        raise ValueError(
            f"Expected 2-D band_traces_{source} with matching SEMs, "
            f"got traces {traces.shape}, sems {sems.shape}"
        )
    if time_ps.ndim != 1 or time_ps.shape[0] != traces.shape[1]:
        raise ValueError(
            f"time_ps shape {time_ps.shape} does not match trace length {traces.shape[1]}"
        )
    return {
        "time_ps": time_ps,
        "traces": traces,
        "sems": sems,
        "meta": meta,
        "n_frames_total": n_frames_total,
    }


def combine_bands(
    traces: np.ndarray, sems: np.ndarray, indices: list[int]
) -> tuple[np.ndarray, np.ndarray]:
    """Sum the selected bands into one trace; independent SEMs add in quadrature."""
    trace = traces[indices].sum(axis=0)
    sem = np.sqrt((sems[indices] ** 2).sum(axis=0))
    return trace, sem


def select_time_window(
    time_ps: np.ndarray,
    trace: np.ndarray,
    sem: np.ndarray,
    t_min: float,
    t_max: float | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Restrict all arrays to bins with t_min <= time <= t_max (None = full range)."""
    upper = float(time_ps[-1]) if t_max is None else t_max
    mask = (time_ps >= t_min) & (time_ps <= upper)
    if not np.any(mask):
        raise ValueError(f"Time window [{t_min}, {upper}] ps selects no bins")
    return time_ps[mask], trace[mask], sem[mask]


def rebin_trace(
    time_ps: np.ndarray, trace: np.ndarray, sem: np.ndarray, factor: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Merge `factor` adjacent bins: count-weighted mean with quadrature SEM / k.

    trace_merged = sum(trace) / k and sem_merged = sqrt(sum(sem^2)) / k for a
    group of k independent bins (the last group keeps its own smaller k), so
    the merged SNR equals sum(trace) / sqrt(sum(sem^2)) either way.
    """
    if factor < 1:
        raise ValueError(f"rebin factor must be >= 1, got {factor}")
    if factor == 1:
        return time_ps, trace, sem
    group_count = (len(time_ps) + factor - 1) // factor
    out_time = np.empty(group_count)
    out_trace = np.empty(group_count)
    out_sem = np.empty(group_count)
    for group in range(group_count):
        window = slice(group * factor, min((group + 1) * factor, len(time_ps)))
        k = window.stop - window.start
        out_time[group] = time_ps[window].mean()
        out_trace[group] = trace[window].sum() / k
        out_sem[group] = np.sqrt((sem[window] ** 2).sum()) / k
    return out_time, out_trace, out_sem


def normalize_peak(
    trace_l: np.ndarray, trace_r: np.ndarray, sem_l: np.ndarray, sem_r: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Scale traces and SEMs so the largest peak of the two traces is 1.0 (10^0)."""
    peak = max(float(trace_l.max()), float(trace_r.max()))
    if peak <= 0:
        raise ValueError("Cannot peak-normalize: largest trace value is <= 0")
    return trace_l / peak, trace_r / peak, sem_l / peak, sem_r / peak


def compute_snr(trace: np.ndarray, sem: np.ndarray) -> np.ndarray:
    """Per-bin SNR = trace / sem; bins with sem <= 0 get SNR 0."""
    snr = np.zeros_like(trace)
    np.divide(trace, sem, out=snr, where=sem > 0)
    return snr


def snr_shade_ceiling(
    traces: list[np.ndarray], snrs: list[np.ndarray], threshold: float
) -> float | None:
    """Highest positive intensity whose bin SNR is below `threshold`, over all traces.

    Defines the top edge of the gray SNR band; None when every positive bin
    clears the threshold.
    """
    low_levels = [
        trace[(snr < threshold) & (trace > 0)] for trace, snr in zip(traces, snrs)
    ]
    low_levels = [levels for levels in low_levels if levels.size]
    if not low_levels:
        return None
    return float(max(levels.max() for levels in low_levels))


def _place_curve_label(
    axis,
    time_ps: np.ndarray,
    trace: np.ndarray,
    x_fraction: float,
    text: str,
    color: str,
    y_offset: float = 1.5,
) -> None:
    """Put a colored curve name just above the trace's local maximum.

    The anchor bin is the one nearest `x_fraction` of the displayed range; the
    text sits above the highest positive bin within +/-5% of the range around
    it, so oscillation peaks cannot overlap the label.
    """
    valid = np.flatnonzero(trace > 0)
    if not valid.size:
        return
    span = time_ps[-1] - time_ps[0]
    target = time_ps[0] + x_fraction * span
    index = valid[np.argmin(np.abs(time_ps[valid] - target))]
    window = np.abs(time_ps - time_ps[index]) <= 0.05 * span
    window &= trace > 0
    local_peak = float(trace[window].max()) if np.any(window) else float(trace[index])
    axis.text(
        time_ps[index],
        local_peak * y_offset,
        text,
        color=color,
        fontsize=8,
        ha="center",
        va="bottom",
    )


def render_fig2a(
    time_ps: np.ndarray,
    trace_l: np.ndarray,
    trace_r: np.ndarray,
    sem_l: np.ndarray,
    sem_r: np.ndarray,
    snr_threshold: float,
    output_path: Path,
) -> tuple[Path, Path, float | None]:
    """Render the Fig. 2(a)-style panel; write PNG (dpi 300) and PDF alongside.

    Returns (png_path, pdf_path, snr_shade_ceiling).
    """
    snr_l = compute_snr(trace_l, sem_l)
    snr_r = compute_snr(trace_r, sem_r)
    ceiling = snr_shade_ceiling([trace_l, trace_r], [snr_l, snr_r], snr_threshold)

    positive = np.concatenate([trace_l[trace_l > 0], trace_r[trace_r > 0]])
    if not positive.size:
        raise ValueError("Both traces are non-positive everywhere; nothing to plot")
    peak = max(float(trace_l.max()), float(trace_r.max()))
    y_top = Y_TOP_HEADROOM * peak
    y_bottom = min(
        Y_BOTTOM_MAX_FRAC * peak,
        max(Y_BOTTOM_FLOOR_FRAC * peak, 0.6 * float(positive.min())),
    )

    figure, axis = plt.subplots(figsize=(3.5, 2.7), constrained_layout=True)
    if ceiling is not None and ceiling > y_bottom:
        shade_top = min(ceiling, y_top)
        axis.axhspan(y_bottom, shade_top, color=SHADE_COLOR, zorder=0)
        axis.text(
            time_ps[0] + 0.82 * (time_ps[-1] - time_ps[0]),
            shade_top * 0.92,
            f"SNR<{snr_threshold:g}",
            color="0.35",
            fontsize=6.5,
            ha="center",
            va="top",
            zorder=4,
            bbox=dict(facecolor=SHADE_COLOR, edgecolor="none", pad=1.2),
        )

    mask_l = trace_l > 0
    mask_r = trace_r > 0
    axis.scatter(
        time_ps[mask_l], trace_l[mask_l], s=7, c=LEFT_COLOR, edgecolors="none", zorder=3
    )
    axis.scatter(
        time_ps[mask_r], trace_r[mask_r], s=7, c=RIGHT_COLOR, edgecolors="none", zorder=3
    )

    axis.set_yscale("log")
    axis.set_ylim(y_bottom, y_top)
    half_bin = 0.5 * float(np.median(np.diff(time_ps))) if len(time_ps) > 1 else 0.5
    axis.set_xlim(time_ps[0] - half_bin, time_ps[-1] + half_bin)
    # Whole-decade ticks inside the limits (10^0/10^-1/10^-2 when normalized).
    low_decade = int(np.ceil(np.log10(y_bottom)))
    high_decade = int(np.floor(np.log10(y_top)))
    axis.set_yticks([10.0 ** decade for decade in range(low_decade, high_decade + 1)])
    axis.yaxis.set_minor_locator(NullLocator())
    axis.set_xlabel("Time (ps)", fontsize=8)
    axis.set_ylabel("Intensity (a.u.)", fontsize=8)
    axis.tick_params(labelsize=7)
    for spine in axis.spines.values():
        spine.set_linewidth(0.8)

    _place_curve_label(axis, time_ps, trace_l, 0.26, "Left", LEFT_COLOR)
    _place_curve_label(axis, time_ps, trace_r, 0.36, "Right", RIGHT_COLOR, y_offset=1.25)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = output_path.with_suffix(".pdf")
    figure.savefig(output_path, dpi=300)
    figure.savefig(pdf_path)
    plt.close(figure)
    return output_path, pdf_path, ceiling


def build_lr_traces(
    data: dict,
    bands_l: list[int],
    bands_r: list[int],
    t_min: float,
    t_max: float | None,
    rebin: int,
    normalize: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Combine bands, apply the time window and rebinning, and normalize.

    Returns (time_ps, trace_l, trace_r, sem_l, sem_r).
    """
    band_count = data["traces"].shape[0]
    for index in bands_l + bands_r:
        if not 0 <= index < band_count:
            raise ValueError(f"Band index {index} out of range for {band_count} bands")
    if set(bands_l) & set(bands_r):
        raise ValueError("--bands-l and --bands-r must not share a band index")

    trace_l, sem_l = combine_bands(data["traces"], data["sems"], bands_l)
    trace_r, sem_r = combine_bands(data["traces"], data["sems"], bands_r)
    time_win, trace_l, sem_l = select_time_window(
        data["time_ps"], trace_l, sem_l, t_min, t_max
    )
    _, trace_r, sem_r = select_time_window(data["time_ps"], trace_r, sem_r, t_min, t_max)
    time_ps, trace_l, sem_l = rebin_trace(time_win, trace_l, sem_l, rebin)
    _, trace_r, sem_r = rebin_trace(time_win, trace_r, sem_r, rebin)
    if normalize == "peak":
        trace_l, trace_r, sem_l, sem_r = normalize_peak(trace_l, trace_r, sem_l, sem_r)
    return time_ps, trace_l, trace_r, sem_l, sem_r


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Render a Fig. 2(a) replica (Adiyatullin et al. 2017) from the "
            "extracted HBT band traces: left/right mesa intensities as blue/red "
            "dots on a log intensity axis with gray SNR shading."
        )
    )
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--source", choices=("photon", "analog"), default="photon")
    parser.add_argument(
        "--bands-l",
        required=True,
        help="Comma-separated band indices summed into the left mesa trace",
    )
    parser.add_argument(
        "--bands-r",
        required=True,
        help="Comma-separated band indices summed into the right mesa trace",
    )
    parser.add_argument("--rebin", type=int, default=2)
    parser.add_argument("--t-min", type=float, default=0.0)
    parser.add_argument("--t-max", type=float, default=None)
    parser.add_argument("--snr-threshold", type=float, default=20.0)
    parser.add_argument("--normalize", choices=("peak", "none"), default="peak")
    arguments = parser.parse_args(argv)

    bands_l = parse_band_selection(arguments.bands_l)
    bands_r = parse_band_selection(arguments.bands_r)
    data = load_band_traces(arguments.npz, arguments.source)
    time_ps, trace_l, trace_r, sem_l, sem_r = build_lr_traces(
        data,
        bands_l,
        bands_r,
        t_min=arguments.t_min,
        t_max=arguments.t_max,
        rebin=arguments.rebin,
        normalize=arguments.normalize,
    )
    png_path, pdf_path, ceiling = render_fig2a(
        time_ps, trace_l, trace_r, sem_l, sem_r, arguments.snr_threshold, arguments.out
    )

    bin_width = (time_ps[-1] - time_ps[0]) / max(len(time_ps) - 1, 1)
    print(
        f"source={arguments.source} bands L={bands_l} R={bands_r} rebin={arguments.rebin} "
        f"({bin_width:.3f} ps/bin), bins={len(time_ps)} over "
        f"{time_ps[0]:.1f}–{time_ps[-1]:.1f} ps, normalize={arguments.normalize}",
        flush=True,
    )
    print(
        f"peaks: L={trace_l.max():.4g}, R={trace_r.max():.4g}; "
        f"SNR<{arguments.snr_threshold:g} shade top: "
        + (f"{ceiling:.4g}" if ceiling is not None else "none (all bins clear threshold)"),
        flush=True,
    )
    if data["n_frames_total"] is not None:
        print(f"n_frames_total={data['n_frames_total']}", flush=True)
    print(f"Wrote {png_path} and {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
