"""Measure the intensity oscillation in the synchroscan k-space condition.

The 26-08-04 batch DOES contain an intensity oscillation — directly visible in
`synchroscan_kspace_image` (accumulated in-camera, so per-frame timing jitter
does not apply): tilted two-condensate interference fringes (relative phase
winding in time). At a fixed k column the intensity oscillates with a dominant
period of ~74-84 ps over rows ~60-400 (2 ns synchroscan window per
metadata-conventions.toml; measured 74.4 ps fixed-k / 74.4 and 83.7 ps 2D-FFT
peaks). The no-pulse k-space calibration shows
no fringes, ruling out a detector artifact.

This script renders the summary figure and prints the quantitative
measurements used in docs/26-08-04-fig2a.md.
"""

import argparse
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from streak_time_axis import ps_per_px

# Both analyzed conditions (02 synchroscan_kspace_image and 05
# kspace_no_pulse_intensity_calibration) are synchroscan-mode: 2 ns window
# (metadata-conventions.toml [streak_camera.time_range_0804]).
PS_PER_PX = ps_per_px("02_Synchroscan_Kspace_image")
FRINGE_ROWS = (60, 400)
FRINGE_COLS = (400, 1000)
K_BINS = [(c, c + 10) for c in range(450, 950, 50)]  # fixed-k intensity traces


def load_sums(h5_path):
    with h5py.File(h5_path, "r") as h:
        pulse = h["conditions/synchroscan_kspace_image/frames"][:].astype(np.float64)
        calib = h["conditions/kspace_no_pulse_intensity_calibration/frames"][:].astype(np.float64)
    return pulse.sum(axis=0), calib.sum(axis=0), pulse


def fixed_k_periods(region):
    """Dominant oscillation period of fixed-k intensity traces (rows axis)."""
    r0, r1 = FRINGE_ROWS
    ny = r1 - r0
    rows = np.arange(ny)
    freq = np.fft.rfftfreq(ny, d=PS_PER_PX)
    specs, periods = [], []
    for c0, c1 in K_BINS:
        v = region[:, c0 - FRINGE_COLS[0] : c1 - FRINGE_COLS[0]].mean(axis=1)
        v = v - np.polyval(np.polyfit(rows, v, 2), rows)
        sp = np.abs(np.fft.rfft(v * np.hanning(ny))) ** 2
        specs.append(sp)
        periods.append(1.0 / freq[1 + np.argmax(sp[1:])])
    return freq, np.mean(specs, axis=0), periods


def fft2_peaks(region, n_peaks=6):
    ny, nx = region.shape
    F = np.abs(np.fft.fft2(region - region.mean()))
    F[:3, :] = F[-3:, :] = 0
    F[:, :3] = 0
    F[:, -3:] = 0
    out = []
    for i in np.argsort(F.ravel())[::-1][:n_peaks]:
        a, b = np.unravel_index(i, F.shape)
        ft = a if a <= ny // 2 else ny - a
        fx = b if b <= nx // 2 else nx - b
        out.append((ny / ft * PS_PER_PX, nx / fx, F[a, b]))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", required=True)
    parser.add_argument("--out", required=True, help="output figure path (drafts)")
    args = parser.parse_args()

    pulse_sum, calib_sum, pulse_frames = load_sums(args.h5)
    r0, r1 = FRINGE_ROWS
    c0, c1 = FRINGE_COLS
    region = pulse_sum[r0:r1, c0:c1]

    freq, mean_spec, periods = fixed_k_periods(region)
    peaks = fft2_peaks(region)
    per_frame = []
    for n in range(pulse_frames.shape[0]):
        fr = pulse_frames[n][r0:r1, c0:c1]
        F = np.abs(np.fft.fft2(fr - fr.mean()))
        F[:5, :] = F[-5:, :] = 0
        F[:, :5] = 0
        F[:, -5:] = 0
        per_frame.append(F.max() / np.median(F))

    print(f"fixed-k oscillation periods per k-bin: {[f'{p:.1f}' for p in periods]} ps")
    print("top 2D-FFT peaks (temporal period ps, spatial period px, power):")
    for tp, sp, pw in peaks:
        print(f"  T={tp:.1f} ps  K={sp:.0f} px  |F|={pw:.3e}")
    print(f"per-frame 2D-FFT peak/median: {[f'{v:.0f}' for v in per_frame]}")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), constrained_layout=True)
    extent = (c0, c1, r0 * PS_PER_PX, r1 * PS_PER_PX)
    lo, hi = np.percentile(region, [1, 99.9])
    im = axes[0].imshow(region, origin="lower", aspect="auto", cmap="magma",
                        interpolation="nearest", extent=extent, vmin=lo, vmax=hi)
    axes[0].set_title("synchroscan k-space (sum of 10 frames)")
    axes[0].set_xlabel("Detector axis 2 [pixel]")
    axes[0].set_ylabel("Time [ps] (row 0 = window start, direction uncalibrated)")
    cb = fig.colorbar(im, ax=axes[0], pad=0.02)
    cb.set_label(f"summed counts (linear, magma, p1-p99.9 = {lo:,.0f}..{hi:,.0f})")

    calib_region = calib_sum[r0:r1, c0:c1]
    lo2, hi2 = np.percentile(calib_region, [1, 99.9])
    im2 = axes[1].imshow(calib_region, origin="lower", aspect="auto", cmap="magma",
                         interpolation="nearest", extent=extent, vmin=lo2, vmax=hi2)
    axes[1].set_title("no-pulse k-space calibration (sum of 10 frames)")
    axes[1].set_xlabel("Detector axis 2 [pixel]")
    cb2 = fig.colorbar(im2, ax=axes[1], pad=0.02)
    cb2.set_label(f"summed counts (linear, magma, p1-p99.9 = {lo2:,.0f}..{hi2:,.0f})")

    period_axis = 1.0 / freq[1:]
    axes[2].plot(period_axis, mean_spec[1:], lw=1)
    axes[2].set_xlim(10, 120)
    axes[2].set_xlabel("period [ps]")
    axes[2].set_ylabel("mean |FFT|$^2$ (detrended, Hann)")
    axes[2].set_title(f"fixed-k intensity oscillation, peak at {np.median(periods):.1f} ps")
    fig.savefig(args.out, dpi=200)
    plt.close(fig)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
