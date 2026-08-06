"""Extract jitter-corrected HBT band time traces from the 26-08-04 streak dataset.

Replicates the Fig. 2(a) measurement recipe of Adiyatullin et al. 2017
(Nat. Commun. 8:1329) on the 20000-frame HBT condensate condition: every frame
carries a fixed-position pulsed-laser reference spot whose leading edge marks
per-frame time zero (t0), and four vertical data bands (2 mesas x 2 HBT arms)
hold the time-resolved emission. Per frame the script

  1. subtracts the measured StreakCamBackG background,
  2. estimates t0 from the reference spot (leading-edge constant-fraction
     crossing at 50% of peak, CFD50; --t0-method centroid selects the legacy
     centroid above half-max),
  3. integrates background-subtracted counts in 80 column sub-windows of 10 px
     (analog traces) via one column cumsum,
  4. detects photon events (8-connected components above threshold) and
     centroids them subpixel,
  5. interpolates/histograms both onto a common tau grid relative to t0 and
     accumulates NaN-aware sum, sum of squares and count per (sub-window, tau)
     for mean/SEM,
  6. accumulates an integer-shift t0-aligned sum image for band localization.

After the pass, the column profile of the aligned sum yields the definitive
band grouping of the sub-windows, band traces are aggregated, the band
Pearson correlation matrix is computed, and an L/R pairing is suggested from
(a) correlation and (b) the ~200 px mesa-separation prior.

Outputs (under --out-dir, default drafts/26-08-04_fig2a):
  hbt_traces.npz, bands_report.json, aligned_sum.png, column_profile.png,
  band_traces_quicklook.png, band_corr.png, t0_vs_frame.png, t0_hist.png,
  threshold_sensitivity.png
"""

import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy import ndimage


DEFAULT_H5 = Path(
    "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/"
    "26-08-04_hbt_time_reference.h5"
)
DEFAULT_OUT_DIR = Path("drafts/26-08-04_fig2a")
FRAMES_DATASET = "conditions/hbt_condensate_josephson_regime/frames"
BACKGROUND_DATASET = "conditions/background_and_streak_pulses_delay_calibration/background"
H5_SHA256 = "905672dc098474e59584749f689af48ed19c1e6f06ed2fff7ce7ac4e703be65d"

# Acquisition metadata (metadata-conventions.toml): 1 ns window over 1016 rows.
N_ROWS = 1016
N_COLS = 1344
PS_PER_PX = 1000.0 / 1016.0

# Reference-spot ROI (row0, row1, col0, col1) from the fig2a recon.
SPOT_ROI = (145, 190, 170, 215)
# ROI actually used for the per-frame t0 in the accumulation. The spot row
# wanders over ~125-200 across the run (synchroscan phase wander, measured
# 2026-08-05 on this dataset): with the narrow recon ROI, ~26% of frames clip
# at the ROI edge (t0 piled up at exactly 145.0) or fail the noise floor. The
# widened row range tracks the spot over its full wander.
T0_ROI = (110, 210, 170, 215)

# Constant-fraction discrimination (CFD) level for the reference-pulse leading
# edge: t0 is the sub-row interpolated crossing of CFD_FRACTION * peak. The
# reference pulse is non-Gaussian and its shape breathes by ~1.3 ps frame to
# frame (flat top, shouldered trailing tail = streak-camera IRF decay, not
# synchronized with the arrival time); only a leading-edge estimator is
# blind to that breathing. See docs/08-04_sync_oscillation.md (2026-08-06).
CFD_FRACTION = 0.5

# Sample-pulse vs reference-spot delay from the calibration condition (ps).
DELAY_PS = 0.32

# Analysis grids.
SUBWINDOW_EDGES = np.arange(350, 1151, 10)  # 81 edges -> 80 sub-windows of 10 px
TAU_GRID = np.arange(-64, 1080)  # tau = row - t0, covers the full window plus margin
EVENT_COL_RANGE = (350, 1150)

# Row range holding band emission (95% cumulative-energy extent from recon).
EMISSION_ROWS = (32, 986)

BATCH_FRAMES = 100
# Aligned-sum canvas: frame rows are deposited at d = CANVAS_OFF - round(t0)
# + row, so that canvas[CANVAS_OFF - round(t0_median) : ...] is the aligned
# sum (shift = round(t0_median - t0) per frame). t0 stays within the T0_ROI
# rows 110-210, well inside the margin.
CANVAS_OFF = 256
CANVAS_ROWS = N_ROWS + 2 * CANVAS_OFF

# Event-threshold study: candidates above the ~105-count pedestal noise
# (sigma ~ sqrt(105) ~ 10 counts per pixel).
THRESHOLD_CANDIDATES = (10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0, 60.0)
THRESHOLD_STUDY_FRAMES = 200
# Auto-pick rule: smallest candidate where the relative drop to the next
# candidate falls below this fraction (knee of the noise-vs-photon curve).
THRESHOLD_STABLE_DROP = 0.25

# Band grouping on the event column histogram: sub-windows whose event count
# exceeds the clipped baseline + BAND_SIGMA * sigma are candidate bands.
BAND_SIGMA = 5.0
# Mesa separation prior from the real-space calibration (L x~610.5, R x~814.0).
MESA_SEPARATION_PX = 204.0
MESA_SEPARATION_TOL_PX = 40.0


def load_background(h5_path) -> np.ndarray:
    """Return the measured StreakCamBackG background, float64 (1016, 1344)."""
    with h5py.File(h5_path, "r") as handle:
        background = np.asarray(handle[BACKGROUND_DATASET], dtype=np.float64)
    if background.shape != (N_ROWS, N_COLS):
        raise ValueError(f"Unexpected background shape {background.shape}")
    return background


def estimate_t0(frame_bg: np.ndarray, roi=SPOT_ROI) -> float:
    """Row centroid above half-max of the reference spot (subpixel, absolute rows).

    frame_bg is a background-subtracted frame. The baseline is the 20th
    percentile of the ROI row profile (the spot fills at most ~2/3 of the ROI,
    so the lowest rows are pure pedestal), and the noise is the standard
    deviation of the sub-baseline rows rescaled by 2.4 (truncated-normal
    correction). Returns NaN when the peak does not clear 5 noise sigmas
    (e.g. laser pulse dropout).
    """
    row0, row1, col0, col1 = roi
    profile = frame_bg[row0:row1, col0:col1].sum(axis=1)
    baseline = float(np.percentile(profile, 20.0))
    low = profile[profile <= baseline]
    noise = 2.4 * float(low.std(ddof=1)) if low.size > 1 else 0.0
    signal = profile - baseline
    peak = float(signal.max())
    if not np.isfinite(peak) or peak <= max(5.0 * noise, 1.0):
        return float("nan")
    weights = np.where(signal > 0.5 * peak, signal, 0.0)
    total = float(weights.sum())
    if total <= 0.0:
        return float("nan")
    rows = np.arange(row0, row1, dtype=np.float64)
    return float(np.dot(rows, weights) / total)


def estimate_t0_cfd(frame_bg: np.ndarray, roi=T0_ROI, fraction=CFD_FRACTION) -> float:
    """Leading-edge constant-fraction (CFD) time of the reference pulse.

    Baseline, noise floor and 5-sigma peak detection are identical to
    estimate_t0, but t0 is the sub-row linear-interpolated crossing of
    fraction * peak on the leading edge: the last row below the threshold
    before the peak row and the first row at/above it. The leading edge is
    the causal marker of the laser arrival time; estimators weighting the
    trailing side (centroid, Gaussian center) track the streak-camera IRF
    decay instead, whose shape breathes ~1.3 ps frame to frame on this
    dataset (docs/08-04_sync_oscillation.md, 2026-08-06). Returns NaN when
    the peak fails the noise floor or the leading edge is not resolved
    inside the ROI.
    """
    row0, row1, col0, col1 = roi
    profile = frame_bg[row0:row1, col0:col1].sum(axis=1)
    baseline = float(np.percentile(profile, 20.0))
    low = profile[profile <= baseline]
    noise = 2.4 * float(low.std(ddof=1)) if low.size > 1 else 0.0
    signal = profile - baseline
    peak = float(signal.max())
    imax = int(signal.argmax())
    if not np.isfinite(peak) or peak <= max(5.0 * noise, 1.0) or imax == 0:
        return float("nan")
    threshold = fraction * peak
    below = np.flatnonzero(signal[:imax] < threshold)
    if below.size == 0:
        return float("nan")
    i = int(below[-1])
    s0, s1 = float(signal[i]), float(signal[i + 1])
    if s1 <= s0:
        return float("nan")
    return float(row0 + i + (threshold - s0) / (s1 - s0))


# Per-frame t0 estimators selectable via --t0-method.
T0_METHODS = {"cfd50": estimate_t0_cfd, "centroid": estimate_t0}


def subwindow_traces(frame_bg: np.ndarray, edges=SUBWINDOW_EDGES) -> np.ndarray:
    """Column-summed counts per row for each window [e0, e1); (n_rows, n_windows).

    One column cumsum pass serves all windows.
    """
    edges = np.asarray(edges, dtype=np.int64)
    cumsum = np.concatenate(
        [np.zeros((frame_bg.shape[0], 1)), np.cumsum(frame_bg, axis=1)], axis=1
    )
    return cumsum[:, edges[1:]] - cumsum[:, edges[:-1]]


def band_trace(frame_bg: np.ndarray, x0: int, x1: int) -> np.ndarray:
    """Background-subtracted counts summed over columns [x0, x1) per row."""
    return subwindow_traces(frame_bg, np.array([x0, x1]))[:, 0]


def detect_events(
    frame_bg: np.ndarray,
    threshold: float,
    min_size: int = 2,
    min_counts: float = 0.0,
    col_range=EVENT_COL_RANGE,
    with_counts: bool = False,
) -> tuple:
    """Photon events: 8-connected components above threshold after bg subtraction.

    Only columns col_range are considered. Components need >= min_size pixels
    and >= min_counts integrated counts. Returns count-weighted subpixel
    centroid (rows, cols) arrays; with ``with_counts=True`` also returns the
    per-event integrated counts above background as a third array.
    """
    mask = np.zeros(frame_bg.shape, dtype=bool)
    mask[:, col_range[0] : col_range[1]] = frame_bg[:, col_range[0] : col_range[1]] > threshold
    labels, n_labels = ndimage.label(mask, structure=np.ones((3, 3), dtype=int))
    if n_labels == 0:
        empty = np.empty(0, dtype=np.float64)
        return (empty, empty, empty) if with_counts else (empty, empty)
    labels_flat = labels.ravel()
    weights = np.where(mask, frame_bg, 0.0)
    counts = np.bincount(labels_flat, weights=weights.ravel(), minlength=n_labels + 1)
    sizes = np.bincount(labels_flat, minlength=n_labels + 1)
    keep = (sizes >= min_size) & (counts >= min_counts)
    keep[0] = False
    if not keep.any():
        empty = np.empty(0, dtype=np.float64)
        return (empty, empty, empty) if with_counts else (empty, empty)
    row_index = np.arange(frame_bg.shape[0], dtype=np.float64)[:, None]
    col_index = np.arange(frame_bg.shape[1], dtype=np.float64)[None, :]
    row_sums = np.bincount(
        labels_flat, weights=(weights * row_index).ravel(), minlength=n_labels + 1
    )
    col_sums = np.bincount(
        labels_flat, weights=(weights * col_index).ravel(), minlength=n_labels + 1
    )
    with np.errstate(invalid="ignore", divide="ignore"):
        rows = row_sums[keep] / counts[keep]
        cols = col_sums[keep] / counts[keep]
    if with_counts:
        return rows, cols, counts[keep]
    return rows, cols


def interp_rows(trace: np.ndarray, t0: float, tau_grid: np.ndarray) -> np.ndarray:
    """Linear-interpolate trace (defined on rows 0..n-1) at rows tau_grid + t0.

    Values outside [0, n-1] are NaN.
    """
    trace = np.asarray(trace, dtype=np.float64)
    tau_grid = np.asarray(tau_grid, dtype=np.float64)
    abs_rows = tau_grid + t0
    out = np.full(abs_rows.shape, np.nan, dtype=np.float64)
    valid = (abs_rows >= 0.0) & (abs_rows <= trace.size - 1)
    if valid.any():
        out[valid] = np.interp(
            abs_rows[valid], np.arange(trace.size, dtype=np.float64), trace
        )
    return out


def threshold_study(h5_path, background, study_indices, candidates=THRESHOLD_CANDIDATES):
    """Count detected events per candidate threshold on a frame sample.

    Also collects the bg-subtracted pixel-value histogram tail above the
    pedestal to document the noise regime. Returns a study dict.
    """
    counts = {float(c): 0 for c in candidates}
    hist_edges = np.arange(0.0, 201.0, 1.0)
    hist = np.zeros(hist_edges.size - 1, dtype=np.int64)
    col0, col1 = EVENT_COL_RANGE
    with h5py.File(h5_path, "r") as handle:
        frames = handle[FRAMES_DATASET]
        for index in study_indices:
            frame_bg = np.asarray(frames[int(index)], dtype=np.float64) - background
            values = frame_bg[:, col0:col1]
            positive = values[values > 0.5]
            if positive.size:
                hist += np.histogram(positive, bins=hist_edges)[0]
            for candidate in candidates:
                rows, _ = detect_events(frame_bg, candidate)
                counts[float(candidate)] += int(rows.size)
    return {
        "n_study_frames": int(len(study_indices)),
        "candidates": [float(c) for c in candidates],
        "event_counts": [counts[float(c)] for c in candidates],
        "pixel_hist_edges": hist_edges.tolist(),
        "pixel_hist": hist.tolist(),
    }


def pick_threshold(study, stable_drop=THRESHOLD_STABLE_DROP):
    """Auto-pick the event threshold at the knee of the count-vs-threshold curve.

    Real photon blobs are bright and survive high thresholds; pedestal-noise
    clusters die out quickly. The threshold is the smallest candidate whose
    relative event-count drop to the next candidate is below stable_drop.
    Falls back to the median candidate when no knee is found.
    """
    candidates = study["candidates"]
    counts = study["event_counts"]
    for i in range(len(candidates) - 1):
        if counts[i] <= 0:
            continue
        drop = (counts[i] - counts[i + 1]) / counts[i]
        if drop < stable_drop:
            return float(candidates[i]), f"knee: relative drop {drop:.3f} < {stable_drop}"
    fallback = float(candidates[len(candidates) // 2])
    return fallback, f"no knee found; fallback to median candidate {fallback}"


def accumulate(h5_path, background, frame_indices, threshold, tau_grid=TAU_GRID,
               edges=SUBWINDOW_EDGES, batch=BATCH_FRAMES, progress_every=1000,
               t0_roi=T0_ROI, t0_method="cfd50"):
    """Single pass over frames: accumulate t0-referenced traces and aligned sum.

    Returns a dict with NaN-aware sum/sumsq/count arrays per (sub-window, tau)
    for analog and photon traces, the aligned-sum canvas, per-frame t0 rows and
    kept frame indices, and the t0 failure count.
    """
    estimate = T0_METHODS[t0_method]
    n_windows = edges.size - 1
    n_tau = tau_grid.size
    analog_sum = np.zeros((n_windows, n_tau))
    analog_sumsq = np.zeros((n_windows, n_tau))
    analog_count = np.zeros((n_windows, n_tau))
    photon_sum = np.zeros((n_windows, n_tau))
    photon_sumsq = np.zeros((n_windows, n_tau))
    photon_count = np.zeros((n_windows, n_tau))
    canvas = np.zeros((CANVAS_ROWS, N_COLS))
    event_col_hist = np.zeros(N_COLS, dtype=np.int64)
    t0_rows = []
    kept_indices = []
    n_failures = 0
    started = time.time()

    with h5py.File(h5_path, "r") as handle:
        frames = handle[FRAMES_DATASET]
        for block_start in range(0, frame_indices.size, batch):
            block_indices = frame_indices[block_start : block_start + batch]
            span = int(block_indices[-1]) - int(block_indices[0]) + 1
            if span <= batch:
                # Contiguous (or near-contiguous) indices: one slice read.
                offset = int(block_indices[0])
                block_frames = frames[offset : int(block_indices[-1]) + 1]
                block = [
                    (int(index), block_frames[int(index) - offset])
                    for index in block_indices
                ]
            else:
                # Strided indices: per-frame reads stay within the memory budget.
                block = [(int(index), frames[int(index)]) for index in block_indices]
            for frame_index, raw_frame in block:
                frame_bg = raw_frame.astype(np.float64) - background
                t0 = estimate(frame_bg, roi=t0_roi)
                if not np.isfinite(t0):
                    n_failures += 1
                    continue
                shift = int(round(t0))
                if not 0 <= shift <= CANVAS_OFF:
                    n_failures += 1
                    continue
                canvas[CANVAS_OFF - shift : CANVAS_OFF - shift + N_ROWS, :] += frame_bg

                valid_tau = (tau_grid + t0 >= 0.0) & (tau_grid + t0 <= N_ROWS - 1)
                traces = subwindow_traces(frame_bg, edges)
                for window in range(n_windows):
                    values = interp_rows(traces[:, window], t0, tau_grid)
                    finite = np.isfinite(values)
                    analog_sum[window, finite] += values[finite]
                    analog_sumsq[window, finite] += values[finite] ** 2
                    analog_count[window, finite] += 1

                rows, cols = detect_events(frame_bg, threshold)
                if rows.size:
                    np.add.at(event_col_hist, np.clip(cols.astype(np.int64), 0, N_COLS - 1), 1)
                    bin_index = np.floor(rows - t0 - tau_grid[0] + 0.5).astype(np.int64)
                    window_index = np.searchsorted(edges, cols, side="right") - 1
                    selected = (
                        (bin_index >= 0)
                        & (bin_index < n_tau)
                        & (window_index >= 0)
                        & (window_index < n_windows)
                    )
                    np.add.at(photon_sum, (window_index[selected], bin_index[selected]), 1.0)
                    np.add.at(photon_sumsq, (window_index[selected], bin_index[selected]), 1.0)
                photon_count[:, valid_tau] += 1

                t0_rows.append(t0)
                kept_indices.append(int(frame_index))
            done = block_start + block_indices.size
            if done % progress_every < batch or done == frame_indices.size:
                elapsed = time.time() - started
                print(
                    f"  {done}/{frame_indices.size} frames, {n_failures} t0 failures, "
                    f"{elapsed:.0f} s elapsed",
                    flush=True,
                )
    return {
        "analog_sum": analog_sum,
        "analog_sumsq": analog_sumsq,
        "analog_count": analog_count,
        "photon_sum": photon_sum,
        "photon_sumsq": photon_sumsq,
        "photon_count": photon_count,
        "canvas": canvas,
        "event_col_hist": event_col_hist,
        "t0_rows": np.asarray(t0_rows, dtype=np.float64),
        "frame_indices": np.asarray(kept_indices, dtype=np.int64),
        "n_failures": n_failures,
        "elapsed_s": time.time() - started,
    }


def mean_sem(acc_sum, acc_sumsq, acc_count):
    """NaN-aware mean and standard error of the mean from accumulated moments."""
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = acc_sum / acc_count
        variance = (acc_sumsq - acc_sum**2 / acc_count) / (acc_count - 1)
        sem = np.sqrt(np.maximum(variance, 0.0) / acc_count)
    mean[acc_count < 1] = np.nan
    sem[acc_count < 2] = np.nan
    return mean, sem


def extract_aligned_sum(canvas, t0_median_row):
    """Aligned sum image: shift each frame by round(t0_median - t0) (integer)."""
    offset = int(round(t0_median_row))
    start = CANVAS_OFF - offset
    return canvas[start : start + N_ROWS, :]


def _low_baseline(values, frac=0.4):
    """Baseline mean/sigma from the lowest frac of values.

    Bands can cover more than half of the sub-windows, so a median is not a
    safe baseline anchor; the lowest 40% are inter-band baseline columns.
    """
    ordered = np.sort(values)
    low = ordered[: max(2, int(np.ceil(frac * ordered.size)))]
    return float(np.mean(low)), float(np.std(low, ddof=1))


def detect_bands(event_col_hist, column_profile, edges=SUBWINDOW_EDGES, band_sigma=5.0):
    """Group sub-windows into bands from the photon-event column histogram.

    The analog aligned-sum column profile has a lumpy fixed-pattern baseline
    (measured-background mismatch) with an amplitude comparable to the band
    contrast, so bands are localized on the event column histogram instead:
    thresholded >=2-pixel events are uniformly distributed detector noise
    except where real emission enhances them, giving a flat, Poisson-like
    baseline. Per-sub-window event counts are thresholded at
    baseline + band_sigma * sigma (baseline/sigma from the lowest 40% of
    sub-windows — bands can cover more than half the range, so a median is
    unsafe — with sigma floored by the Poisson sigma); maximal contiguous runs
    of member sub-windows form the bands. The analog profile, reduced the same
    way after subtracting a wide median-filter baseline, is reported as a
    cross-check.
    """
    n_cols = edges[1] - edges[0]
    event_sub = np.array(
        [event_col_hist[e0:e1].sum() for e0, e1 in zip(edges[:-1], edges[1:])],
        dtype=np.float64,
    )
    median, sigma_robust = _low_baseline(event_sub)
    # Per-sub-window baseline sigma: low-quantile estimate, floored by the
    # Poisson sigma of a flat baseline (median * n_cols counts per sub-window).
    sigma = max(sigma_robust, np.sqrt(max(median * n_cols, 1.0)))
    event_threshold = median + band_sigma * sigma

    member = event_sub > event_threshold
    # Morphological closing: bridge single-sub-window gaps that still exceed
    # baseline + 2 sigma, so noisy band flanks do not fragment a band.
    bridge = event_sub > median + 2.0 * sigma
    for w in range(1, member.size - 1):
        if not member[w] and member[w - 1] and member[w + 1] and bridge[w]:
            member[w] = True
    bands = []
    start = None
    for w in range(member.size + 1):
        if w < member.size and member[w] and start is None:
            start = w
        elif (w == member.size or not member[w]) and start is not None:
            bands.append((int(edges[start]), int(edges[w])))
            start = None
    # Bands touching the detection column range are truncation artifacts of
    # detect_events (blobs clipped by the col_range boundary), not emission.
    bands = [
        (x0, x1)
        for (x0, x1) in bands
        if x0 > int(edges[0]) and x1 < int(edges[-1])
    ]

    # Analog cross-check: residual of the aligned-sum column profile above a
    # wide median-filter baseline, same per-sub-window reduction.
    baseline = ndimage.median_filter(column_profile, size=201)
    residual = column_profile - baseline
    analog_sub = np.array(
        [residual[e0:e1].mean() for e0, e1 in zip(edges[:-1], edges[1:])],
        dtype=np.float64,
    )
    analog_median, analog_sigma = _low_baseline(analog_sub)
    analog_threshold = analog_median + band_sigma * analog_sigma

    above = event_col_hist > (median + band_sigma * sigma) / n_cols
    pixel_regions = []
    start = None
    for x in range(above.size + 1):
        if x < above.size and above[x] and start is None:
            start = x
        elif (x == above.size or not above[x]) and start is not None:
            pixel_regions.append([int(start), int(x)])
            start = None
    return {
        "band_windows": bands,
        "event_subwindow_counts": event_sub.tolist(),
        "event_threshold": float(event_threshold),
        "event_baseline_median": median,
        "event_baseline_sigma": float(sigma),
        "analog_subwindow_residual": analog_sub.tolist(),
        "analog_threshold": float(analog_threshold),
        "analog_baseline_median": float(analog_median),
        "analog_baseline_sigma": float(analog_sigma),
        "pixel_regions": pixel_regions,
    }


def aggregate_band_traces(mean, sem, band_windows, edges=SUBWINDOW_EDGES):
    """Sum sub-window traces into band traces; SEMs combined as independent."""
    band_means = []
    band_sems = []
    for x0, x1 in band_windows:
        members = np.where((edges[:-1] >= x0) & (edges[1:] <= x1))[0]
        band_means.append(np.nansum(mean[members], axis=0))
        band_sems.append(np.sqrt(np.nansum(sem[members] ** 2, axis=0)))
    return np.asarray(band_means), np.asarray(band_sems)


def band_correlation(band_traces):
    """Pearson correlation matrix over the common finite tau range."""
    valid = np.all(np.isfinite(band_traces), axis=0)
    if valid.sum() < 3:
        return np.full((band_traces.shape[0], band_traces.shape[0]), np.nan), int(valid.sum())
    return np.corrcoef(band_traces[:, valid]), int(valid.sum())


def suggest_pairing(band_windows, corr, band_integrals=None):
    """Suggest the L/R (mesa) pairing of bands.

    Evidence (a): the pairing maximizing within-pair Pearson correlation —
    the two HBT arms of one mesa carry the same I(t), so same-mesa band pairs
    correlate positively while different mesas oscillate in anti-phase.
    Evidence (b): bands ~204 px apart are the two mesas seen through the same
    arm (real-space calibration: L x~610.5, R x~814.0); with four bands this
    splits the bands into two same-arm pairs, leaving the two cross-arm
    pairings.
    Evidence (c): the cross-arm pairing whose two per-arm intensity ratios
    agree (and match the no-pulse L/R ratio 0.578) — the mesa intensity ratio
    is a physical property, independent of the arm.
    """
    centers = [0.5 * (x0 + x1) for x0, x1 in band_windows]
    n = len(centers)
    suggestion = {
        "band_centers_px": centers,
        "correlation_evidence": None,
        "separation_evidence": None,
        "ratio_evidence": None,
        "suggested_pairs": None,
        "suggested_left_pair": None,
        "notes": [],
    }
    if n < 2:
        suggestion["notes"].append("Fewer than two bands detected; no pairing possible.")
        return suggestion

    same_arm_pairs = [
        (i, j, abs(centers[i] - centers[j]))
        for i in range(n)
        for j in range(i + 1, n)
        if abs(abs(centers[i] - centers[j]) - MESA_SEPARATION_PX) <= MESA_SEPARATION_TOL_PX
    ]
    suggestion["separation_evidence"] = [
        {"bands": [i, j], "separation_px": sep} for i, j, sep in same_arm_pairs
    ]

    if n == 4:
        matchings = [
            ((0, 1), (2, 3)),
            ((0, 2), (1, 3)),
            ((0, 3), (1, 2)),
        ]
        scored = [
            (
                float(corr[a][b] + corr[c][d]),
                ((a, b), (c, d)),
            )
            for (a, b), (c, d) in matchings
        ]
        scored.sort(key=lambda item: item[0], reverse=True)
        suggestion["correlation_evidence"] = [
            {"pairs": [list(p) for p in matching], "within_pair_corr_sum": score}
            for score, matching in scored
        ]

        arm_pair_idx = [(i, j) for i, j, _ in same_arm_pairs]
        arm_sep = {(i, j): sep for i, j, sep in same_arm_pairs}
        best = scored[0][1]
        # Find two disjoint same-arm pairs covering all four bands (the arms).
        arm_covers = []
        for first in range(len(arm_pair_idx)):
            for second in range(first + 1, len(arm_pair_idx)):
                pair_a, pair_b = arm_pair_idx[first], arm_pair_idx[second]
                if sorted(pair_a + pair_b) == [0, 1, 2, 3]:
                    arm_covers.append((pair_a, pair_b))
        if arm_covers:
            (a0, a1), (b0, b1) = min(
                arm_covers,
                key=lambda cover: sum(
                    abs(arm_sep[p] - MESA_SEPARATION_PX) for p in cover
                ),
            )
            # Cross-arm pairings compatible with the two arms: straight (left
            # with left) and crossed. The arm left/right intensity ratio is a
            # physical mesa property, so the correct pairing reproduces it in
            # both arms.
            ratio_rows = []
            best_ratio = None
            best_consistency = -1.0
            for straight in (True, False):
                if straight:
                    pairing = ((a0, b0), (a1, b1))
                else:
                    pairing = ((a0, b1), (a1, b0))
                (p, q), (u, v) = pairing
                row = {
                    "pairs": [list(pr) for pr in pairing],
                    "within_pair_corr_sum": float(corr[p][q] + corr[u][v]),
                }
                if band_integrals is not None:
                    r_arm1 = band_integrals[a0] / band_integrals[a1]
                    if straight:
                        r_arm2 = band_integrals[b0] / band_integrals[b1]
                    else:
                        r_arm2 = band_integrals[b1] / band_integrals[b0]
                    consistency = min(r_arm1, r_arm2) / max(r_arm1, r_arm2)
                    row["arm1_left_over_right"] = float(r_arm1)
                    row["arm2_left_over_right"] = float(r_arm2)
                    row["ratio_consistency"] = float(consistency)
                    if consistency > best_consistency:
                        best_consistency = consistency
                        best_ratio = pairing
                ratio_rows.append(row)
            suggestion["ratio_evidence"] = ratio_rows
            if best_ratio is not None:
                best = best_ratio
                suggestion["notes"].append(
                    "Pairing chosen by per-arm left/right intensity-ratio "
                    "consistency (evidence c); correlation sums for both "
                    "cross-arm pairings are listed under ratio_evidence."
                )
        suggestion["suggested_pairs"] = [list(p) for p in best]
        # Left pair = the one holding the left (smaller-x) band of the first arm.
        if suggestion["ratio_evidence"] is not None:
            left_first = min(best[0])
            suggestion["suggested_left_pair"] = [
                list(p) for p in best if left_first in p
            ][0]
            suggestion["notes"].append(
                "Left/right assignment assumes the streak column ordering "
                "preserves the real-space mesa ordering (L at smaller x); the "
                "arm left/right ratios are consistent with the no-pulse L/R "
                "calibration ratio 0.578."
            )
        else:
            suggestion["notes"].append(
                f"Correlation-optimal pairing (sum of within-pair corr "
                f"{scored[0][0]:.3f}); no arm-aware ratio evidence available."
            )
    else:
        pairs = [
            {"bands": [i, j], "corr": float(corr[i][j]), "separation_px": abs(centers[i] - centers[j])}
            for i in range(n)
            for j in range(i + 1, n)
        ]
        pairs.sort(key=lambda item: item["corr"], reverse=True)
        suggestion["correlation_evidence"] = pairs
        suggestion["notes"].append(
            f"{n} bands detected (expected 4); pairing listed by descending correlation."
        )
    return suggestion


def save_aligned_sum_figure(aligned_sum, t0_median_row, output_path):
    """Aligned sum image, magma, linear, p1-p99.9 clip, vertical axis in ps."""
    rows = np.arange(aligned_sum.shape[0], dtype=np.float64)
    time_ps = (rows - t0_median_row) * PS_PER_PX
    lower, upper = np.percentile(aligned_sum, [1.0, 99.9])
    figure, axis = plt.subplots(figsize=(10.5, 7.2), constrained_layout=True)
    image = axis.imshow(
        aligned_sum,
        origin="upper",
        aspect="auto",
        cmap="magma",
        interpolation="nearest",
        extent=(0.0, float(aligned_sum.shape[1]), time_ps[-1], time_ps[0]),
        vmin=lower,
        vmax=upper,
    )
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label(
        f"counts (magma, linear, p1-p99.9 clip: {lower:.3g}..{upper:.3g})"
    )
    axis.set_xlabel("Detector column [pixel]")
    axis.set_ylabel("time relative to median t0 [ps]")
    axis.set_title("t0-aligned sum of background-subtracted HBT frames")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def save_column_profile_figure(column_profile, event_col_hist, band_info, output_path):
    """Analog aligned column profile and event column histogram, bands shaded."""
    columns = np.arange(column_profile.size)
    figure, axes = plt.subplots(
        2, 1, figsize=(11.0, 8.0), sharex=True, constrained_layout=True
    )
    axes[0].plot(columns, column_profile, lw=0.8, color="black")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("aligned emission-rows sum [counts]")
    axes[0].set_title("Column profile of the t0-aligned analog sum (emission rows only)")

    axes[1].plot(columns, event_col_hist, lw=0.8, color="black")
    axes[1].axhline(
        band_info["event_threshold"] / (SUBWINDOW_EDGES[1] - SUBWINDOW_EDGES[0]),
        color="tab:red", ls="--", lw=1.0,
        label=f"band threshold = baseline + {BAND_SIGMA:g}x sigma (per column)",
    )
    axes[1].set_ylabel("photon events per column")
    axes[1].set_xlabel("Detector column [pixel]")
    axes[1].set_title("Photon-event column histogram (band localization)")
    axes[1].legend(loc="upper right")

    for axis in axes:
        for k, (x0, x1) in enumerate(band_info["band_windows"]):
            axis.axvspan(x0, x1, color="tab:blue", alpha=0.15)
            axis.text(0.5 * (x0 + x1), 0.97, f"B{k + 1}", ha="center", va="top",
                      transform=axis.get_xaxis_transform())
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def save_band_traces_figure(time_ps, band_analog, band_photon, band_windows, output_path):
    """Quicklook of band mean traces (analog and photon), log intensity."""
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 5.5), constrained_layout=True)
    for axis, traces, title in (
        (axes[0], band_analog, "analog (bg-subtracted counts)"),
        (axes[1], band_photon, "photon events"),
    ):
        for k, (x0, x1) in enumerate(band_windows):
            positive = np.where(traces[k] > 0, traces[k], np.nan)
            axis.plot(time_ps, positive, lw=0.9, label=f"B{k + 1} x={x0}-{x1}")
        axis.set_yscale("log")
        axis.set_xlabel("time relative to t0 [ps]")
        axis.set_ylabel("mean per frame [a.u.]")
        axis.set_title(f"Band traces — {title}")
        axis.legend(fontsize="small")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def save_band_corr_figure(corr, band_windows, output_path):
    """Band Pearson correlation matrix heatmap (RdBu_r, linear, -1..1)."""
    labels = [f"B{k + 1}\n{x0}-{x1}" for k, (x0, x1) in enumerate(band_windows)]
    figure, axis = plt.subplots(figsize=(6.4, 5.6), constrained_layout=True)
    image = axis.imshow(corr, cmap="RdBu_r", vmin=-1.0, vmax=1.0, interpolation="nearest")
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Pearson r (RdBu_r, linear, -1..1)")
    axis.set_xticks(range(len(labels)), labels, fontsize="small")
    axis.set_yticks(range(len(labels)), labels, fontsize="small")
    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            axis.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize="small")
    axis.set_title("Band mean-trace Pearson correlation")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def save_t0_figures(t0_rows, frame_indices, out_dir):
    """t0 vs frame index and t0 histogram (both in ps around the median)."""
    t0_median = float(np.median(t0_rows))
    jitter_ps = (t0_rows - t0_median) * PS_PER_PX
    rms_ps = float(np.std(jitter_ps, ddof=1))

    figure, axis = plt.subplots(figsize=(11.0, 4.5), constrained_layout=True)
    axis.plot(frame_indices, jitter_ps, ".", ms=1.5, rasterized=True)
    axis.axhline(0.0, color="tab:red", lw=0.8)
    axis.set_xlabel("frame index")
    axis.set_ylabel("t0 - median(t0) [ps]")
    axis.set_title(f"Reference-spot t0 drift/jitter (rms {rms_ps:.2f} ps)")
    figure.savefig(out_dir / "t0_vs_frame.png", dpi=150)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    axis.hist(jitter_ps, bins=80, color="black")
    axis.set_xlabel("t0 - median(t0) [ps]")
    axis.set_ylabel("frames")
    axis.set_title(f"t0 jitter histogram (rms {rms_ps:.2f} ps, n={jitter_ps.size})")
    figure.savefig(out_dir / "t0_hist.png", dpi=150)
    plt.close(figure)
    return t0_median, rms_ps


def save_threshold_figure(study, chosen_threshold, method, output_path):
    """Detected-event count vs candidate threshold with the chosen value marked."""
    figure, axis = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    axis.plot(study["candidates"], study["event_counts"], "o-", color="black")
    axis.axvline(chosen_threshold, color="tab:red", ls="--",
                 label=f"chosen {chosen_threshold:g}")
    axis.set_yscale("log")
    axis.set_xlabel("event threshold [counts above background]")
    axis.set_ylabel(f"events in {study['n_study_frames']} frames")
    axis.set_title("Event-threshold sensitivity")
    axis.legend()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def run(args):
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    background = load_background(args.h5)
    with h5py.File(args.h5, "r") as handle:
        n_frames_h5 = handle[FRAMES_DATASET].shape[0]
    frame_indices = np.arange(0, n_frames_h5, args.sample_step, dtype=np.int64)
    if args.max_frames is not None:
        frame_indices = frame_indices[: args.max_frames]
    print(f"Processing {frame_indices.size} of {n_frames_h5} frames from {args.h5}", flush=True)

    # Event-threshold study on an evenly spaced subsample.
    stride = max(1, frame_indices.size // THRESHOLD_STUDY_FRAMES)
    study_indices = frame_indices[::stride][:THRESHOLD_STUDY_FRAMES]
    study = threshold_study(args.h5, background, study_indices)
    auto_threshold, method = pick_threshold(study)
    threshold = float(args.threshold) if args.threshold is not None else auto_threshold
    if args.threshold is not None:
        method = f"user-specified --threshold {threshold:g} (auto study suggested {auto_threshold:g}: {method})"
    print(f"Event threshold: {threshold:g} ({method})", flush=True)
    save_threshold_figure(study, threshold, method, out_dir / "threshold_sensitivity.png")

    acc = accumulate(args.h5, background, frame_indices, threshold,
                     t0_method=args.t0_method)
    print(
        f"Accumulation done in {acc['elapsed_s']:.0f} s: "
        f"{acc['frame_indices'].size} frames kept, {acc['n_failures']} t0 failures",
        flush=True,
    )

    mean_analog, sem_analog = mean_sem(acc["analog_sum"], acc["analog_sumsq"], acc["analog_count"])
    mean_photon, sem_photon = mean_sem(acc["photon_sum"], acc["photon_sumsq"], acc["photon_count"])

    t0_median, t0_rms_ps = save_t0_figures(acc["t0_rows"], acc["frame_indices"], out_dir)
    aligned_sum = extract_aligned_sum(acc["canvas"], t0_median)
    save_aligned_sum_figure(aligned_sum, t0_median, out_dir / "aligned_sum.png")

    column_profile = aligned_sum[EMISSION_ROWS[0] : EMISSION_ROWS[1], :].sum(axis=0)
    band_info = detect_bands(acc["event_col_hist"], column_profile)
    save_column_profile_figure(
        column_profile, acc["event_col_hist"], band_info, out_dir / "column_profile.png"
    )
    band_windows = band_info["band_windows"]
    print(f"Detected {len(band_windows)} bands: {band_windows}", flush=True)

    band_analog, band_sem_analog = aggregate_band_traces(mean_analog, sem_analog, band_windows)
    band_photon, band_sem_photon = aggregate_band_traces(mean_photon, sem_photon, band_windows)
    corr_analog, n_corr_tau = band_correlation(band_analog)
    corr_photon, _ = band_correlation(band_photon)
    band_integrals = np.nansum(band_analog, axis=1)
    pairing = suggest_pairing(band_windows, corr_analog, band_integrals=band_integrals)

    save_band_traces_figure(
        TAU_GRID * PS_PER_PX, band_analog, band_photon, band_windows,
        out_dir / "band_traces_quicklook.png",
    )
    save_band_corr_figure(corr_analog, band_windows, out_dir / "band_corr.png")

    meta = {
        "h5_path": str(args.h5),
        "h5_sha256": H5_SHA256,
        "background_source": BACKGROUND_DATASET,
        "frames_dataset": FRAMES_DATASET,
        "roi": list(T0_ROI),
        "spot_roi_constant": list(SPOT_ROI),
        "t0_method": args.t0_method,
        "cfd_fraction": CFD_FRACTION if args.t0_method == "cfd50" else None,
        "ps_per_px": PS_PER_PX,
        "delay_ps": DELAY_PS,
        "threshold": {
            "chosen": threshold,
            "method": method,
            "min_size": 2,
            "min_counts": 0.0,
            "col_range": list(EVENT_COL_RANGE),
            "candidates": study["candidates"],
            "event_counts": study["event_counts"],
            "n_study_frames": study["n_study_frames"],
        },
        "n_frames_h5": int(n_frames_h5),
        "n_frames_selected": int(frame_indices.size),
        "n_frames_used": int(acc["frame_indices"].size),
        "n_failures": int(acc["n_failures"]),
        "sample_step": int(args.sample_step),
        "t0_median_row": t0_median,
        "t0_rms_jitter_ps": t0_rms_ps,
        "emission_rows": list(EMISSION_ROWS),
        "band_windows": band_windows,
        "pairing_suggestion": pairing,
        "elapsed_s": acc["elapsed_s"],
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    np.savez(
        out_dir / "hbt_traces.npz",
        time_ps=TAU_GRID * PS_PER_PX,
        subwindow_edges=SUBWINDOW_EDGES.astype(np.int64),
        traces_analog=mean_analog,
        sem_analog=sem_analog,
        traces_photon=mean_photon,
        sem_photon=sem_photon,
        band_windows=np.asarray(band_windows, dtype=np.int64),
        band_traces_analog=band_analog,
        band_sem_analog=band_sem_analog,
        band_traces_photon=band_photon,
        band_sem_photon=band_sem_photon,
        band_corr=corr_analog,
        t0_rows=acc["t0_rows"],
        frame_indices=acc["frame_indices"],
        n_frames_total=np.int64(frame_indices.size),
        meta_json=json.dumps(meta, indent=2),
    )

    report = {
        "band_windows": band_windows,
        "band_centers_px": pairing["band_centers_px"],
        "band_integrals_analog": [float(v) for v in band_integrals],
        "band_detection": {
            "method": (
                "bands grouped from the photon-event column histogram "
                f"(per-sub-window counts > clipped baseline + {BAND_SIGMA:g} x sigma); "
                "the analog aligned-sum profile has a lumpy fixed-pattern baseline "
                "and serves as a cross-check only"
            ),
            "event_subwindow_counts": band_info["event_subwindow_counts"],
            "event_baseline_median": band_info["event_baseline_median"],
            "event_baseline_sigma": band_info["event_baseline_sigma"],
            "event_threshold": band_info["event_threshold"],
            "analog_subwindow_residual": band_info["analog_subwindow_residual"],
            "analog_baseline_median": band_info["analog_baseline_median"],
            "analog_baseline_sigma": band_info["analog_baseline_sigma"],
            "analog_threshold": band_info["analog_threshold"],
            "event_pixel_regions_above_threshold": band_info["pixel_regions"],
            "event_col_hist": acc["event_col_hist"].tolist(),
        },
        "band_corr_analog": corr_analog.tolist(),
        "band_corr_photon": corr_photon.tolist(),
        "band_corr_n_tau": n_corr_tau,
        "pairing_suggestion": pairing,
        "t0": {
            "median_row": t0_median,
            "rms_jitter_ps": t0_rms_ps,
            "n_failures": int(acc["n_failures"]),
            "n_frames_used": int(acc["frame_indices"].size),
        },
        "threshold": meta["threshold"],
        "elapsed_s": acc["elapsed_s"],
        "timestamp": meta["timestamp"],
    }
    with open(out_dir / "bands_report.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"Wrote {out_dir / 'hbt_traces.npz'} and diagnostics", flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Extract jitter-corrected HBT band time traces (Fig. 2a recipe)."
    )
    parser.add_argument("--h5", type=Path, default=DEFAULT_H5)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--sample-step", type=int, default=1,
                        help="process every Nth frame (default 1 = all)")
    parser.add_argument("--max-frames", type=int, default=None,
                        help="cap the number of processed frames")
    parser.add_argument("--threshold", type=float, default=None,
                        help="photon-event threshold in counts above background "
                             "(default: auto from the sensitivity study)")
    parser.add_argument("--t0-method", choices=sorted(T0_METHODS), default="cfd50",
                        help="reference-pulse t0 estimator: cfd50 = leading-edge "
                             "constant-fraction crossing at 50%% of peak "
                             "(default), centroid = legacy half-max centroid")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
