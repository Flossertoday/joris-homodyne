"""Build dataset A: the curated per-photon event set of the 26-08-04 HBT condition.

Purpose (docs/26-08-08-dataset-a.md): answer "is the Fig. 2(a)/2(c) null a
pipeline artifact?" by constructing an event set in which every event can
physically be condensate signal, and running the downstream analyses on that
set only. The per-frame chain is IDENTICAL to the production extraction
(extract_hbt_traces.py: measured-background subtraction, CFD50 t0 from the
reference pulse, 8-connected detection at 25 counts above background,
min 2 px, columns 350-1150); this script adds (a) per-event morphology the
old table lacks and (b) the curation cuts decided from the 26-08-08
investigation:

  1. horizontal ROI: the event centroid x must fall inside one of the four
     definitive band windows (a window-resize scan showed the current windows
     already maximize S/sqrt(S+N) to <0.6%, so they are kept as-is);
  2. defect columns: 1-px-wide hot columns are detector defects (uniform in
     time, amplitude spectrum identical to neighbors); the in-band one
     (col 676 inside B2, z=14.5) is excluded, out-of-band spikes are already
     removed by the ROIs;
  3. frame-edge touching events (bbox touches row 0/1015 or the detection
     column boundary): clipped blobs with unreliable position;
  4. NO morphology cut: in-band and out-of-band events are the same physical
     objects (MCP clouds); the best measured discriminator (event size)
     rejects noise only 1.7x faster than signal, below the ~3x break-even
     for any sensitivity gain — every morphology cut would LOSE sensitivity.

Cuts 1-3 remove only events that cannot be (usable) signal. What remains is
dataset A: signal + an irreducible, per-event indistinguishable dark-count
floor. That floor is quantified per band during the same pass from the
out-of-band columns (local linear fit with 30-px guards that exclude the
measured ~25-px signal tails) and stored in the file's noise_model group;
downstream analyses read ONLY this file (events + noise model), never the
raw frames.

Outputs: the dataset-A HDF5 (refuses to overwrite, atomic via .partial),
plus QC figures/JSON in --qc-dir.
"""
import argparse
import hashlib
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

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_hbt_traces import (  # noqa: E402
    DEFAULT_H5,
    EVENT_COL_RANGE,
    FRAMES_DATASET,
    N_ROWS,
    PS_PER_PX,
    T0_ROI,
    estimate_t0_cfd,
    load_background,
)

SCHEMA_VERSION = 1
# Definitive band windows (docs/26-08-04-fig2a.md); S/N-optimal per the
# 26-08-08 window-resize scan (tmp/dataset_a_invest/part1b_summary.json).
DEFAULT_BANDS = ((380, 450), (570, 680), (810, 880), (1010, 1110))
# In-band hot column (z=14.5 narrow spike, time-uniform, amplitude spectrum
# identical to neighbors -> detector defect, not condensate emission).
DEFAULT_DEFECT_COLS = (676,)
# All narrow spike columns (incl. out-of-band and detection-range edges),
# excluded from the noise-floor fits. From the same investigation.
DEFAULT_FLOOR_EXCLUDE_COLS = (350, 351, 567, 676, 744, 970, 1147, 1148)
DETECTION_THRESHOLD = 25.0  # production value, knee of the sensitivity scan
TAU_GRID = np.arange(-64, 1080)  # 1-px aligned grid, as extract_hbt_traces
FLOOR_GUARD_PX = 30   # measured signal tails reach ~25-30 px past band edges
FLOOR_FIT_PX = 35
FLOOR_VARIANTS = ((25, 30), (30, 35), (35, 40))
EDGE_MARGIN_PX = 1    # bbox touching row 0/1015 or the col_range boundary


def detect_events_features(frame_bg, threshold=DETECTION_THRESHOLD, min_size=2,
                           col_range=EVENT_COL_RANGE):
    """Production detection + per-event morphology.

    Same mask/label/centroid arithmetic as extract_hbt_traces.detect_events
    (regression-checked against it in tests), plus per-event size [px], peak
    pixel [counts] and bounding box (r0, r1, c0, c1 inclusive).
    """
    mask = np.zeros(frame_bg.shape, dtype=bool)
    mask[:, col_range[0]:col_range[1]] = (
        frame_bg[:, col_range[0]:col_range[1]] > threshold
    )
    labels, n_labels = ndimage.label(mask, structure=np.ones((3, 3), dtype=int))
    empty = {
        key: np.empty(0, dtype=np.float64)
        for key in ("rows", "cols", "amplitude", "size", "peak")
    } | {key: np.empty(0, dtype=np.int32) for key in ("r0", "r1", "c0", "c1")}
    if n_labels == 0:
        return empty
    labels_flat = labels.ravel()
    weights = np.where(mask, frame_bg, 0.0)
    counts = np.bincount(labels_flat, weights=weights.ravel(), minlength=n_labels + 1)
    sizes = np.bincount(labels_flat, minlength=n_labels + 1)
    peaks = np.asarray(
        ndimage.maximum(frame_bg, labels, index=np.arange(n_labels + 1)),
        dtype=np.float64,
    )
    keep = (sizes >= min_size)
    keep[0] = False
    if not keep.any():
        return empty
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
    bbox = {"r0": [], "r1": [], "c0": [], "c1": []}
    kept_labels = np.flatnonzero(keep)
    objects = ndimage.find_objects(labels)
    for label in kept_labels:
        slices = objects[label - 1]
        bbox["r0"].append(slices[0].start)
        bbox["r1"].append(slices[0].stop - 1)
        bbox["c0"].append(slices[1].start)
        bbox["c1"].append(slices[1].stop - 1)
    return {
        "rows": rows,
        "cols": cols,
        "amplitude": counts[keep],
        "size": sizes[keep].astype(np.float64),
        "peak": peaks[keep],
        "r0": np.asarray(bbox["r0"], dtype=np.int32),
        "r1": np.asarray(bbox["r1"], dtype=np.int32),
        "c0": np.asarray(bbox["c0"], dtype=np.int32),
        "c1": np.asarray(bbox["c1"], dtype=np.int32),
    }


def curate_events(feats, bands=DEFAULT_BANDS, defect_cols=DEFAULT_DEFECT_COLS,
                  col_range=EVENT_COL_RANGE, n_rows=N_ROWS,
                  edge_margin=EDGE_MARGIN_PX):
    """Apply the dataset-A cuts. Returns (keep_mask, band_index, reason codes).

    reason codes per event: 0 kept, 1 out-of-band, 2 defect column,
    3 frame-edge/boundary touching (out-of-band takes precedence, then
    defect, then edge).
    """
    x = feats["cols"]
    n = x.size
    band = np.full(n, -1, dtype=np.int64)
    for index, (lo, hi) in enumerate(bands):
        band[(x >= lo) & (x <= hi)] = index
    reason = np.zeros(n, dtype=np.int64)
    reason[band < 0] = 1
    defect = np.isin(np.floor(x).astype(np.int64), np.asarray(defect_cols))
    reason[(reason == 0) & defect] = 2
    edge = (
        (feats["r0"] <= 0)
        | (feats["r1"] >= n_rows - 1)
        | (feats["c0"] <= col_range[0])
        | (feats["c1"] >= col_range[1] - 1)
    )
    reason[(reason == 0) & edge] = 3
    return reason == 0, band, reason


def fit_band_floors(col_hist_out_per_frame, bands=DEFAULT_BANDS,
                    col_range=EVENT_COL_RANGE,
                    exclude_cols=DEFAULT_FLOOR_EXCLUDE_COLS,
                    variants=FLOOR_VARIANTS):
    """Per-band residual-noise floor from the out-of-band column histogram.

    col_hist_out_per_frame: events/frame/col on the col_range grid, with
    out-of-band columns measured and in-band columns arbitrary (unused).
    For each band and each (guard, fit_width) variant, a line is fit to the
    two flank windows [edge-guard-fit, edge-guard) and (edge+guard,
    edge+guard+fit] (excluding defect/edge columns) and integrated under the
    band; if a flank has <10 usable columns (B1 sits 30 px from the detection
    edge) the other flank's mean is held constant (one-sided). Returns
    (floor_mean_per_band [events/frame/col], per-variant array, report dict).
    """
    cols = np.arange(col_range[0], col_range[1])
    rate = np.asarray(col_hist_out_per_frame, dtype=np.float64)
    excluded = np.isin(cols, np.asarray(exclude_cols))

    def region(c0, c1):
        cc = np.arange(max(c0, col_range[0]), min(c1, col_range[1]))
        cc = cc[~np.isin(cc, np.asarray(exclude_cols))]
        return cc

    floors = np.zeros((len(bands), len(variants)))
    report = []
    for k, (lo, hi) in enumerate(bands):
        sel = (cols >= lo) & (cols <= hi)
        width = float(sel.sum())  # inclusive membership -> hi-lo+1 columns
        variant_info = []
        for v, (guard, fit) in enumerate(variants):
            left = region(lo - guard - fit, lo - guard)
            right = region(hi + 1 + guard, hi + 1 + guard + fit)
            one_sided = None
            if left.size >= 10 and right.size >= 10:
                fx = np.concatenate([left, right])
                fy = rate[np.isin(cols, fx)]
                a, b = np.polyfit(fx, fy, 1)
                total = float((a * cols[sel] + b).sum())
            elif right.size >= 10:
                one_sided = "right"
                total = float(rate[np.isin(cols, right)].mean() * width)
            elif left.size >= 10:
                one_sided = "left"
                total = float(rate[np.isin(cols, left)].mean() * width)
            else:
                total = float("nan")
            floors[k, v] = total / width
            variant_info.append({
                "guard": guard, "fit_width": fit,
                "events_per_frame_in_band": total,
                "one_sided": one_sided,
                "left_cols": int(left.size), "right_cols": int(right.size),
            })
        report.append({
            "band": [int(lo), int(hi)],
            "floor_events_per_frame_per_col_mean": float(np.mean(floors[k])),
            "floor_halfspread": float((floors[k].max() - floors[k].min()) / 2),
            "variants": variant_info,
        })
    return floors.mean(axis=1), floors, report


def build(h5_path, out_path, qc_dir, bands=DEFAULT_BANDS,
          defect_cols=DEFAULT_DEFECT_COLS, threshold=DETECTION_THRESHOLD,
          sample_step=1, max_frames=None, progress_every=1000,
          floor_exclude_cols=DEFAULT_FLOOR_EXCLUDE_COLS):
    """Full curation pass. Writes the dataset-A HDF5 + QC artifacts."""
    out_path = Path(out_path)
    qc_dir = Path(qc_dir)
    qc_dir.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing file: {out_path}")
    partial = out_path.with_name(out_path.name + ".partial")
    if partial.exists():
        raise FileExistsError(f"Interrupted artifact already exists: {partial}")

    background = load_background(h5_path)
    with h5py.File(h5_path, "r") as handle:
        n_frames_total = handle[FRAMES_DATASET].shape[0]
    frame_indices = np.arange(0, n_frames_total, sample_step, dtype=np.int64)
    if max_frames is not None:
        frame_indices = frame_indices[:max_frames]
    print(f"Curating {frame_indices.size} of {n_frames_total} frames from {h5_path}",
          flush=True)

    n_tau = TAU_GRID.size
    n_cols = EVENT_COL_RANGE[1] - EVENT_COL_RANGE[0]
    noise_col_hist = np.zeros(n_cols, dtype=np.int64)  # out-of-band, all cols
    noise_tau_hist = np.zeros(n_tau, dtype=np.int64)   # out-of-band, 1-px bins
    band_tau_hist = np.zeros((len(bands), n_tau), dtype=np.int64)  # kept events
    detected_col_hist = np.zeros(n_cols, dtype=np.int64)  # all detected events
    cut_counts = {"total_detected": 0, "kept": 0, "out_of_band": 0,
                  "defect_column": 0, "edge_touch": 0}
    store = {key: [] for key in
             ("frame", "x", "y", "tau", "amplitude", "size", "peak", "r0", "r1",
              "c0", "c1", "band")}
    frame_t0 = np.full(frame_indices.size, np.nan)
    failed = []
    started = time.time()

    with h5py.File(h5_path, "r") as handle:
        frames = handle[FRAMES_DATASET]
        for done, index in enumerate(frame_indices):
            frame_bg = np.asarray(frames[int(index)], dtype=np.float64) - background
            t0 = estimate_t0_cfd(frame_bg, roi=T0_ROI)
            if not np.isfinite(t0):
                failed.append(int(index))
                continue
            frame_t0[done] = t0
            feats = detect_events_features(frame_bg, threshold)
            n_ev = feats["rows"].size
            if n_ev == 0:
                continue
            keep_mask, band, reason = curate_events(
                feats, bands=bands, defect_cols=defect_cols
            )
            cut_counts["total_detected"] += n_ev
            cut_counts["kept"] += int(keep_mask.sum())
            cut_counts["out_of_band"] += int((reason == 1).sum())
            cut_counts["defect_column"] += int((reason == 2).sum())
            cut_counts["edge_touch"] += int((reason == 3).sum())
            col_idx = np.clip(
                np.floor(feats["cols"]).astype(np.int64) - EVENT_COL_RANGE[0],
                0, n_cols - 1)
            np.add.at(detected_col_hist, col_idx, 1)
            out = band < 0
            if out.any():
                np.add.at(noise_col_hist, col_idx[out], 1)
                tau_out = feats["rows"][out] - t0
                bin_idx = np.floor(tau_out - TAU_GRID[0]).astype(np.int64)
                ok = (bin_idx >= 0) & (bin_idx < n_tau)
                np.add.at(noise_tau_hist, bin_idx[ok], 1)
            if keep_mask.any():
                tau = feats["rows"][keep_mask] - t0
                bin_idx = np.floor(tau - TAU_GRID[0]).astype(np.int64)
                ok = (bin_idx >= 0) & (bin_idx < n_tau)
                np.add.at(band_tau_hist, (band[keep_mask][ok], bin_idx[ok]), 1)
                store["frame"].append(
                    np.full(int(keep_mask.sum()), index, dtype=np.int32))
                store["x"].append(feats["cols"][keep_mask])
                store["y"].append(feats["rows"][keep_mask])
                store["tau"].append(tau)
                store["amplitude"].append(feats["amplitude"][keep_mask])
                store["size"].append(feats["size"][keep_mask].astype(np.int32))
                store["peak"].append(feats["peak"][keep_mask])
                for key in ("r0", "r1", "c0", "c1"):
                    store[key].append(feats[key][keep_mask])
                store["band"].append(band[keep_mask].astype(np.int8))
            if (done + 1) % progress_every == 0 or done + 1 == frame_indices.size:
                print(f"  {done + 1}/{frame_indices.size} frames, "
                      f"{cut_counts['kept']} events kept, "
                      f"{time.time() - started:.0f}s", flush=True)

    n_kept_frames = int(np.isfinite(frame_t0).sum())
    events = {key: (np.concatenate(value) if value else np.empty(0))
              for key, value in store.items()}

    # --- noise model ---------------------------------------------------------
    member = np.zeros(n_cols, dtype=bool)
    for lo, hi in bands:
        member |= (np.arange(n_cols) + EVENT_COL_RANGE[0] >= lo) & \
                  (np.arange(n_cols) + EVENT_COL_RANGE[0] <= hi)
    out_width = float((~member).sum())
    col_rate = noise_col_hist / max(n_kept_frames, 1)
    col_rate_for_fit = col_rate.copy()
    col_rate_for_fit[member] = 0.0  # unused by the fit
    floor_mean, floor_variants, floor_report = fit_band_floors(
        col_rate_for_fit, bands=bands, exclude_cols=floor_exclude_cols)
    rho_tau_global = noise_tau_hist / max(n_kept_frames, 1) / out_width
    # flatness of the out-of-band density along tau (window 0-820 ps)
    win = (TAU_GRID >= 0) & (TAU_GRID <= 820 / PS_PER_PX)
    expected = noise_tau_hist[win].mean()
    chi2_ndf = float(
        np.sum((noise_tau_hist[win] - expected) ** 2) / max(expected, 1e-12)
        / max(int(win.sum()) - 1, 1)
    )

    band_totals = [int((events["band"] == k).sum()) for k in range(len(bands))]
    noise_fraction = {
        f"B{k + 1}": {
            "events": band_totals[k],
            "noise_events_per_frame": float(
                floor_mean[k] * (bands[k][1] - bands[k][0] + 1)),
            "noise_fraction": float(
                floor_mean[k] * (bands[k][1] - bands[k][0] + 1) * n_kept_frames
                / max(band_totals[k], 1)
            ),
        }
        for k in range(len(bands))
    }

    attributes = {
        "schema_version": SCHEMA_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "source_h5": str(Path(h5_path).resolve()),
        "source_h5_sha256": _sha256(h5_path),
        "detection_threshold_counts": threshold,
        "detection_min_size_px": 2,
        "detection_col_range": list(EVENT_COL_RANGE),
        "t0_method": "cfd50",
        "t0_roi": list(T0_ROI),
        "ps_per_px": PS_PER_PX,
        "band_windows": np.asarray(bands, dtype=np.float64),
        "defect_columns": np.asarray(defect_cols, dtype=np.int64),
        "floor_exclude_columns": np.asarray(floor_exclude_cols, dtype=np.int64),
        "floor_guard_fit_variants": np.asarray(FLOOR_VARIANTS, dtype=np.int64),
        "n_frames_total": n_frames_total,
        "n_frames_processed": int(frame_indices.size),
        "n_frames_kept": n_kept_frames,
        "n_frames_failed": len(failed),
        "n_events": int(events["frame"].size),
        "sample_step": sample_step,
        "max_frames": -1 if max_frames is None else max_frames,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with h5py.File(partial, "x") as handle:
            for key, value in attributes.items():
                handle.attrs[key] = value
            group = handle.create_group("events")
            for name, data in events.items():
                group.create_dataset(
                    name, data=data, compression="gzip", shuffle=True)
            handle.create_dataset("frame_t0", data=frame_t0)
            handle.create_dataset("frame_indices", data=frame_indices)
            handle.create_dataset("failed_frames",
                                  data=np.asarray(failed, dtype=np.int32))
            nm = handle.create_group("noise_model")
            nm.create_dataset("floor_per_band", data=floor_mean)
            nm.create_dataset("floor_variants", data=floor_variants)
            nm.create_dataset("rho_tau_global", data=rho_tau_global)
            nm.create_dataset("tau_grid", data=TAU_GRID.astype(np.float64))
            nm.create_dataset("out_of_band_col_rate", data=col_rate)
            nm.create_dataset("detected_col_hist", data=detected_col_hist)
            nm.create_dataset("noise_tau_hist", data=noise_tau_hist)
            nm.create_dataset("band_tau_hist", data=band_tau_hist)
            nm.attrs["flatness_chi2_ndf_window_0_820ps"] = chi2_ndf
            nm.attrs["out_of_band_width_cols"] = out_width
            nm.attrs["floor_report_json"] = json.dumps(floor_report)
            for key, value in cut_counts.items():
                nm.attrs[f"cut_{key}"] = value
        if out_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing file: {out_path}")
        partial.rename(out_path)
    except BaseException:
        if partial.exists():
            partial.unlink()
        raise

    print(f"dataset A written: {out_path} ({events['frame'].size} events, "
          f"{n_kept_frames} frames, {len(failed)} t0 failures)", flush=True)
    qc = {
        "cut_counts": cut_counts,
        "noise_fraction_per_band": noise_fraction,
        "floor_report": floor_report,
        "noise_tau_flatness_chi2_ndf_window_0_820ps": chi2_ndf,
        "band_event_totals": band_totals,
        "n_frames_kept": n_kept_frames,
        "failed_frames": failed,
    }
    with open(qc_dir / "dataset_a_qc.json", "w") as stream:
        json.dump(qc, stream, indent=2)
    return out_path, qc


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        while chunk := stream.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def render_qc_figures(h5_path, dataset_a_path, qc_dir, bands=DEFAULT_BANDS,
                      defect_cols=DEFAULT_DEFECT_COLS,
                      sample_frames=(100, 5000, 12000),
                      threshold=DETECTION_THRESHOLD):
    """QC figures: x-profile with ROIs/floors and kept-vs-rejected overlays."""
    qc_dir = Path(qc_dir)
    with h5py.File(dataset_a_path, "r") as handle:
        nm = handle["noise_model"]
        detected = np.asarray(nm["detected_col_hist"])
        floor = np.asarray(nm["floor_per_band"])
        n_kept = int(handle.attrs["n_frames_kept"])
        band_tau = np.asarray(nm["band_tau_hist"])
        noise_tau = np.asarray(nm["noise_tau_hist"])
        out_width = float(nm.attrs["out_of_band_width_cols"])
        ev_x = np.asarray(handle["events/x"])
        ev_tau = np.asarray(handle["events/tau"])
        ev_band = np.asarray(handle["events/band"])

    cols = np.arange(EVENT_COL_RANGE[0], EVENT_COL_RANGE[1])
    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                             constrained_layout=True)
    axes[0].plot(cols, detected / n_kept, lw=0.6, color="black")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("detected events / frame / col")
    axes[1].plot(cols, detected / n_kept, lw=0.6, color="black")
    axes[1].set_ylabel("detected events / frame / col")
    axes[1].set_xlabel("detector column [px]")
    for ax in axes:
        for k, (lo, hi) in enumerate(bands):
            ax.axvspan(lo, hi, color="tab:green", alpha=0.15)
            ax.hlines(floor[k], lo, hi, color="tab:red", lw=1.2)
            axes[1].text(0.5 * (lo + hi), 0.96, f"A{k + 1}", ha="center",
                         va="top", transform=axes[1].get_xaxis_transform())
    axes[0].set_title(
        "Dataset-A horizontal ROIs (green) with per-band noise floors (red)")
    fig.savefig(qc_dir / "dataset_a_xprofile.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 4.5), constrained_layout=True)
    tau_ps = TAU_GRID * PS_PER_PX
    ax.plot(tau_ps, noise_tau / n_kept / out_width, lw=0.7, color="tab:red",
            label="out-of-band noise / col")
    for k in range(len(bands)):
        width = bands[k][1] - bands[k][0] + 1
        ax.plot(tau_ps, band_tau[k] / n_kept / width, lw=0.7,
                label=f"A{k + 1} kept / col")
        ax.hlines(floor[k], tau_ps[0], tau_ps[-1], ls="--", lw=0.7)
    ax.axvspan(0, 820, color="tab:green", alpha=0.06)
    ax.set_xlabel("tau [ps]")
    ax.set_ylabel("events / frame / col / px-tau")
    ax.legend(fontsize="small", ncol=2)
    ax.set_title("Per-band kept-event tau profiles vs noise floors "
                 "(green band = analysis window)")
    fig.savefig(qc_dir / "dataset_a_tau_floor.png", dpi=150)
    plt.close(fig)

    background = load_background(h5_path)
    fig, axes = plt.subplots(1, len(sample_frames),
                             figsize=(6.2 * len(sample_frames), 5.2),
                             constrained_layout=True)
    with h5py.File(h5_path, "r") as handle:
        frames = handle[FRAMES_DATASET]
        for ax, index in zip(axes, sample_frames):
            frame_bg = np.asarray(frames[int(index)], dtype=np.float64) - background
            t0 = estimate_t0_cfd(frame_bg, roi=T0_ROI)
            feats = detect_events_features(frame_bg, threshold)
            keep_mask, band, reason = curate_events(
                feats, bands=bands, defect_cols=defect_cols)
            crop = frame_bg[:, EVENT_COL_RANGE[0]:EVENT_COL_RANGE[1]]
            vmax = float(np.percentile(crop, 99.9))
            image = ax.imshow(crop, origin="upper", aspect="auto", cmap="magma",
                              vmin=0.0, vmax=vmax)
            fig.colorbar(image, ax=ax, shrink=0.85).set_label(
                f"counts above bg (magma, linear, 0..p99.9={vmax:.0f})")
            ax.plot(feats["cols"][keep_mask] - EVENT_COL_RANGE[0],
                    feats["rows"][keep_mask], "o", ms=5, mfc="none",
                    color="lime", label="dataset A")
            ax.plot(feats["cols"][~keep_mask] - EVENT_COL_RANGE[0],
                    feats["rows"][~keep_mask], "x", ms=4,
                    color="cyan", label="rejected")
            for lo, hi in bands:
                ax.axvline(lo - EVENT_COL_RANGE[0], color="white", lw=0.5,
                           alpha=0.4)
                ax.axvline(hi - EVENT_COL_RANGE[0], color="white", lw=0.5,
                           alpha=0.4)
            ax.axhline(t0, color="tab:red", lw=0.8)
            ax.set_title(f"frame {index}: kept {keep_mask.sum()}, "
                         f"rejected {(~keep_mask).sum()}", fontsize=9)
            ax.set_xlabel(f"column - {EVENT_COL_RANGE[0]} [px]")
            ax.set_ylabel("row [px]")
            ax.legend(fontsize=7, loc="upper right")
    fig.savefig(qc_dir / "dataset_a_overlay.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", type=Path, default=DEFAULT_H5)
    parser.add_argument("--out", type=Path, default=None,
                        help="default: detected-peaks/<h5 stem>_dataset_a.h5 "
                             "next to this script")
    parser.add_argument("--qc-dir", type=Path,
                        default=Path("drafts/26-08-08_dataset_a"))
    parser.add_argument("--bands", type=str, default=None,
                        help="comma list of 2N column edges; default: "
                             "definitive 26-08-04 windows")
    parser.add_argument("--defect-cols", type=str, default=None)
    parser.add_argument("--floor-exclude-cols", type=str, default=None,
                        help="comma list of spike/edge columns excluded from "
                             "the noise-floor fits; default: the 26-08-04 "
                             "original-batch list")
    parser.add_argument("--threshold", type=float, default=DETECTION_THRESHOLD)
    parser.add_argument("--sample-step", type=int, default=1)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args()

    bands = DEFAULT_BANDS
    if args.bands:
        values = [float(v) for v in args.bands.split(",") if v.strip()]
        if len(values) % 2 or not values:
            raise ValueError("--bands needs an even number of edges")
        bands = tuple((int(lo), int(hi)) for lo, hi in
                      np.asarray(values).reshape(-1, 2))
    defect_cols = DEFAULT_DEFECT_COLS
    if args.defect_cols:
        defect_cols = tuple(int(v) for v in args.defect_cols.split(",") if v.strip())
    floor_exclude_cols = DEFAULT_FLOOR_EXCLUDE_COLS
    if args.floor_exclude_cols:
        floor_exclude_cols = tuple(
            int(v) for v in args.floor_exclude_cols.split(",") if v.strip())
    out_path = args.out
    if out_path is None:
        out_path = (Path(__file__).resolve().parent / "detected-peaks"
                    / f"{Path(args.h5).stem}_dataset_a.h5")

    out_path, qc = build(args.h5, out_path, args.qc_dir, bands=bands,
                         defect_cols=defect_cols, threshold=args.threshold,
                         sample_step=args.sample_step,
                         max_frames=args.max_frames,
                         floor_exclude_cols=floor_exclude_cols)
    print(json.dumps(qc["noise_fraction_per_band"], indent=2), flush=True)
    if not args.skip_figures:
        render_qc_figures(args.h5, out_path, args.qc_dir, bands=bands,
                          defect_cols=defect_cols)
        print(f"QC figures in {args.qc_dir}", flush=True)


if __name__ == "__main__":
    main()
