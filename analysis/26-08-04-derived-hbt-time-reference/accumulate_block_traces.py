"""Block-resolved accumulation of the 26-08-04 HBT bands (no t0 alignment).

Motivation: the t0-aligned averaged traces (extract_hbt_traces.py) show no
anti-phase intensity oscillation, while the data collector reports an
intensity oscillation in this dataset. Two ways the aligned average could
hide it: (a) the oscillation phase drifts between frame blocks (averaging
cancels it), (b) the per-frame t0 alignment itself smears or mis-tracks the
true optical timing. This script makes one pass over all 20000 frames and
produces, WITHOUT any t0 alignment:

- the background-subtracted total sum image (for the preview figure);
- band row traces per 10-frame block (for block power spectra: averaging
  |FFT|^2 over blocks keeps an oscillation whose phase drifts between blocks);
- a two-panel preview figure `_overview_sum_bg_subtracted.png` for the
  preview tree (full frame + spot-masked band contrast);
- a block power-spectrum diagnostic figure (drafts, git-ignored).

Bands (from the full-statistics extraction, see docs/26-08-04-fig2a.md):
B1 380-450, B2 570-680, B3 810-880, B4 1010-1110; L={B1,B3}, R={B2,B4}.
"""

import argparse
import json
import time
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from streak_time_axis import ps_per_px

# Acquisition metadata (metadata-conventions.toml [streak_camera.time_range_0804]):
# the HBT condition is single-sweep, 1 ns window over 1016 rows.
PS_PER_PX = ps_per_px("01_HBT_condensate_Josephson_regime")
BANDS = [(380, 450), (570, 680), (810, 880), (1010, 1110)]
SPOT_MASK_COLS = (150, 270)
BLOCK = 10  # frames per block
BATCH = 100  # frames per H5 read
ROW_RANGE = (50, 950)  # rows used for FFT diagnostics


def accumulate(h5_path, max_frames=None):
    with h5py.File(h5_path, "r") as h:
        background = h["conditions/background_and_streak_pulses_delay_calibration/background"][:]
        frames = h["conditions/hbt_condensate_josephson_regime/frames"]
        n_total = frames.shape[0] if max_frames is None else min(max_frames, frames.shape[0])
        n_total -= n_total % BLOCK
        total_sum = np.zeros(frames.shape[1:], np.float64)
        blocks = []
        t_start = time.time()
        for start in range(0, n_total, BATCH):
            stop = min(start + BATCH, n_total)
            batch = frames[start:stop].astype(np.float64)
            total_sum += batch.sum(axis=0)
            block_sums = batch.reshape(-1, BLOCK, *frames.shape[1:]).sum(axis=1)
            n_blocks = block_sums.shape[0]
            block_sums -= BLOCK * background  # background-subtracted block sums
            traces = np.stack([block_sums[:, :, x0:x1].sum(axis=2) for x0, x1 in BANDS], axis=1)
            blocks.append(traces)  # [n_blocks, 4, 1016]
            print(f"{stop}/{n_total} frames, {time.time()-t_start:.0f}s", flush=True)
        block_traces = np.concatenate(blocks, axis=0)
    return total_sum, block_traces, background, n_total


def save_preview_figure(total_sum, background, n_frames, out_path):
    bg_sub = total_sum - n_frames * background
    t_ps = np.arange(bg_sub.shape[0]) * PS_PER_PX
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.4), constrained_layout=True, sharey=True)

    lo, hi = np.percentile(bg_sub, [1.0, 99.9])
    im = axes[0].imshow(
        bg_sub, origin="lower", aspect="auto", cmap="magma", interpolation="nearest",
        extent=(0, bg_sub.shape[1], 0, t_ps[-1]), vmin=lo, vmax=hi,
    )
    axes[0].set_title("background-subtracted, full columns")
    cb = fig.colorbar(im, ax=axes[0], pad=0.02)
    cb.set_label(f"counts - {n_frames}xStreakCamBackG\n(linear, magma, p1-p99.9 = {lo:,.0f}..{hi:,.0f})")

    masked = bg_sub.copy()
    masked[:, SPOT_MASK_COLS[0] : SPOT_MASK_COLS[1]] = np.nan
    lo2, hi2 = np.nanpercentile(masked, [50.0, 99.9])
    im2 = axes[1].imshow(
        masked, origin="lower", aspect="auto", cmap="magma", interpolation="nearest",
        extent=(0, masked.shape[1], 0, t_ps[-1]), vmin=lo2, vmax=hi2,
    )
    axes[1].set_title(f"spot columns {SPOT_MASK_COLS[0]}-{SPOT_MASK_COLS[1]} masked")
    cb2 = fig.colorbar(im2, ax=axes[1], pad=0.02)
    cb2.set_label(f"counts - {n_frames}xStreakCamBackG\n(linear, magma, p50-p99.9 = {lo2:,.0f}..{hi2:,.0f})")

    for ax in axes:
        ax.set_xlabel("Detector axis 2 [pixel]")
        for x0, x1 in BANDS:
            ax.axvspan(x0, x1, color="cyan", alpha=0.12, lw=0)
    axes[0].set_ylabel("Time [ps] (row 0 = window start, direction uncalibrated)")
    fig.suptitle(f"26-08-04 01_HBT_condensate_Josephson_regime — bg-subtracted sum of {n_frames:,} frames")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def power_spectrum_diagnostics(block_traces, out_path):
    r0, r1 = ROW_RANGE
    rows = np.arange(r0, r1)
    t = rows * PS_PER_PX
    # band traces per block -> L, R, S, z per block
    L = block_traces[:, [0, 2], r0:r1].sum(axis=1)
    R = block_traces[:, [1, 3], r0:r1].sum(axis=1)
    S = L + R

    def detrend(y):
        c = np.polyfit(t, y, 1)
        return y - np.polyval(c, t)

    def avg_power(series):
        # series: [n_blocks, n_rows]; mean |FFT|^2 of detrended blocks
        p = np.zeros((series.shape[1] // 2) + 1)
        win = np.hanning(series.shape[1])
        for block in series:
            p += np.abs(np.fft.rfft(detrend(block) * win)) ** 2
        return p / series.shape[0]

    pS = avg_power(S)
    pz = avg_power((L - R) / np.maximum(S, 1e-9))
    # coherent reference: |FFT|^2 of the grand-average traces
    gc_S = np.abs(np.fft.rfft(detrend(S.mean(axis=0)) * np.hanning(len(t)))) ** 2
    gc_z = np.abs(np.fft.rfft(detrend(((L - R) / np.maximum(S, 1e-9)).mean(axis=0)) * np.hanning(len(t)))) ** 2

    freq = np.fft.rfftfreq(len(t), d=PS_PER_PX)
    period = 1.0 / freq[1:]

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    for col, (p, gc, name) in enumerate([(pS, gc_S, "S = L+R (common mode)"), (pz, gc_z, "z = (L-R)/(L+R)")]):
        axes[0, col].semilogy(period, p[1:], lw=0.8, label="mean block |FFT|$^2$ (phase-drift-proof)")
        axes[0, col].semilogy(period, gc[1:], lw=0.8, label="|FFT|$^2$ of grand average")
        axes[0, col].set_xlim(1.5, 120)
        axes[0, col].set_xlabel("period [ps]")
        axes[0, col].set_title(name)
        axes[0, col].legend(fontsize=8)
        # top peaks of the block-averaged power in 2-100 ps
        sel = (period >= 2) & (period <= 100)
        med = np.median(p[1:][sel])
        top = np.argsort(p[1:][sel])[-3:][::-1]
        txt = "  ".join(f"P={period[sel][i]:.1f}ps x{p[1:][sel][i]/med:.0f}med" for i in top)
        axes[0, col].text(0.02, 0.95, txt, transform=axes[0, col].transAxes, fontsize=8, va="top")
        print(f"{name}: top block-power peaks {txt}")
        # grand-average traces
        axes[1, col].plot(t, S.mean(axis=0) if col == 0 else ((L - R) / np.maximum(S, 1e-9)).mean(axis=0), lw=0.8)
        axes[1, col].set_xlabel("t [ps] (unaligned, row 0 = window start)")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", required=True)
    parser.add_argument("--out-preview", required=True, help="preview PNG path (tracked preview tree)")
    parser.add_argument("--out-dir", required=True, help="git-ignored dir for npy/diagnostics")
    parser.add_argument("--max-frames", type=int, default=None)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    total_sum, block_traces, background, n_frames = accumulate(args.h5, args.max_frames)
    np.save(out_dir / "unaligned_total_sum.npy", total_sum)
    np.save(out_dir / "unaligned_block_traces.npy", block_traces)
    (out_dir / "unaligned_meta.json").write_text(json.dumps({
        "h5": args.h5, "n_frames": n_frames, "block": BLOCK, "bands": BANDS,
        "row_range_fft": ROW_RANGE, "ps_per_px": PS_PER_PX,
    }, indent=2) + "\n")
    save_preview_figure(total_sum, background, n_frames, Path(args.out_preview))
    power_spectrum_diagnostics(block_traces, out_dir / "unaligned_block_power.png")
    print(f"wrote {args.out_preview} and diagnostics to {out_dir}")


if __name__ == "__main__":
    main()
