"""Fig. 2(a)-style band traces from dataset A ONLY (no raw-frame access).

Reads the curated dataset-A event table (build_dataset_a.py) and bins the
kept events onto the production 1-px aligned tau grid, per band and per
10-px subwindow, producing an npz in the exact hbt_traces.npz schema of
extract_hbt_traces.py so the audited downstream consumers (render_fig2a.py,
hbt_z_targeted_fits.py, ...) run unchanged.

Differences vs the raw-frame extraction, by construction:
  * events are the curated set (band ROIs, defect columns and frame-edge
    events removed) — nothing out-of-band can enter;
  * the analog channel does not exist (dataset A stores no pixel data):
    traces_analog/sem_analog are NaN arrays and meta_json carries
    "analog_unavailable": true;
  * per-bin frame counts account for each frame's reachable tau range from
    its stored t0 (same validity rule as the production accumulation).

The residual in-band noise floor (noise_model group) is NOT subtracted in
the traces — they are the raw curated counts, directly comparable to the
old pipeline's photon traces; the floor enters the science interpretation
through the dilution factor quoted by dataset_a_fig2a_analysis.py.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_dataset_a import TAU_GRID  # noqa: E402
from extract_hbt_traces import (  # noqa: E402
    N_ROWS,
    PS_PER_PX,
    SUBWINDOW_EDGES,
    aggregate_band_traces,
    band_correlation,
    mean_sem,
)

DEFAULT_TABLE = (
    Path(__file__).resolve().parent
    / "detected-peaks/26-08-04_hbt_time_reference_dataset_a.h5"
)
DEFAULT_OUT = Path("drafts/26-08-08_dataset_a/hbt_traces_dataset_a.npz")


def build_traces(table_path):
    """Bin dataset-A events into (subwindow, tau) mean/SEM traces."""
    with h5py.File(table_path, "r") as handle:
        events = {
            name: np.asarray(handle["events"][name])
            for name in ("frame", "x", "tau", "band")
        }
        frame_t0 = np.asarray(handle["frame_t0"])
        frame_indices = np.asarray(handle["frame_indices"])
        band_windows = np.asarray(handle.attrs["band_windows"], dtype=np.float64)
        attrs = dict(handle.attrs)

    n_tau = TAU_GRID.size
    n_sub = SUBWINDOW_EDGES.size - 1
    kept_mask = np.isfinite(frame_t0)
    kept_t0 = frame_t0[kept_mask]
    n_kept = kept_t0.size

    # per-frame validity: bin tau is reachable when 0 <= tau + t0 <= N_ROWS-1
    valid = (TAU_GRID[None, :] + kept_t0[:, None] >= 0.0) & (
        TAU_GRID[None, :] + kept_t0[:, None] <= N_ROWS - 1
    )  # (kept frames, n_tau)

    kept_ordinal = np.full(frame_t0.size, -1, dtype=np.int64)
    kept_ordinal[kept_mask] = np.arange(n_kept)
    frame_ord = kept_ordinal[np.searchsorted(frame_indices, events["frame"])]
    sub = np.searchsorted(SUBWINDOW_EDGES, events["x"], side="right") - 1
    bin_index = np.floor(events["tau"] - TAU_GRID[0]).astype(np.int64)
    ok = (
        (frame_ord >= 0)
        & (sub >= 0)
        & (sub < n_sub)
        & (bin_index >= 0)
        & (bin_index < n_tau)
    )
    # guard: only count events in bins the frame can actually reach (a real
    # event's bin is always reachable; the mask costs nothing)
    ok[ok] &= valid[frame_ord[ok], bin_index[ok]]

    # moments per (subwindow, tau): per-frame counts are 0/1 (rarely 2) at
    # this occupancy, so square via the sparse per-(frame, cell) counts
    cell = frame_ord[ok] * (n_sub * n_tau) + (sub[ok] * n_tau + bin_index[ok])
    uniq, c = np.unique(cell, return_counts=True)
    cell_index = uniq % (n_sub * n_tau)
    acc_sum = np.bincount(cell_index, weights=c, minlength=n_sub * n_tau)
    acc_sumsq = np.bincount(cell_index, weights=c**2, minlength=n_sub * n_tau)
    acc_sum = acc_sum.reshape(n_sub, n_tau)
    acc_sumsq = acc_sumsq.reshape(n_sub, n_tau)
    count_per_tau = valid.astype(np.float64).sum(axis=0)
    count_full = np.broadcast_to(count_per_tau, (n_sub, n_tau)).copy()
    mean, sem = mean_sem(acc_sum, acc_sumsq, count_full)

    band_means, band_sems = aggregate_band_traces(
        mean, sem, band_windows, SUBWINDOW_EDGES
    )
    corr, n_overlap = band_correlation(
        np.where(np.isfinite(band_means), band_means, np.nan)
    )
    return {
        "time_ps": TAU_GRID.astype(np.float64) * PS_PER_PX,
        "subwindow_edges": SUBWINDOW_EDGES,
        "traces_photon": mean,
        "sem_photon": sem,
        "band_windows": band_windows.astype(np.int64),
        "band_traces_photon": band_means,
        "band_sem_photon": band_sems,
        "band_corr": corr,
        "t0_rows": kept_t0,
        "frame_indices": frame_indices[kept_mask],
        "n_frames_total": np.int64(frame_t0.size),
        "attrs": attrs,
        "n_overlap": n_overlap,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    result = build_traces(args.table)
    attrs = result.pop("attrs")
    n_overlap = result.pop("n_overlap")
    nan_traces = np.full_like(result["traces_photon"], np.nan)
    meta = {
        "source": "dataset_a_only",
        "dataset_a_table": str(args.table),
        "analog_unavailable": True,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "band_corr_overlap_bins": n_overlap,
        "dataset_a_attrs": {
            key: (value.tolist() if hasattr(value, "tolist") else value)
            for key, value in attrs.items()
        },
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        args.out,
        traces_analog=nan_traces,
        sem_analog=nan_traces.copy(),
        band_traces_analog=np.full_like(result["band_traces_photon"], np.nan),
        band_sem_analog=np.full_like(result["band_sem_photon"], np.nan),
        meta_json=json.dumps(meta),
        **result,
    )
    in_window = (result["time_ps"] >= 0) & (result["time_ps"] <= 820)
    print(f"wrote {args.out}")
    for k in range(result["band_traces_photon"].shape[0]):
        trace = result["band_traces_photon"][k]
        print(f"  B{k + 1}: mean events/frame/bin in 0-820 ps window: "
              f"{np.nanmean(trace[in_window]):.5f}")


if __name__ == "__main__":
    main()
