"""Build an hbt_traces-schema npz from a full curated event table under
explicit band windows (the 2026-08-10 FDR-SNR window switch of the project
report).

Same binning logic as dataset_a_traces.py (which reads the curated dataset-A
h5 with band windows baked in as an attribute); here the band windows are a
CLI argument, so any window set can be evaluated from the window-agnostic
curated event table without re-detecting events on the raw frames. The
analog channel does not exist for an event table -> NaN arrays (same
convention as dataset_a_traces.py), so downstream consumers
(render_fig2a.py, hbt_z_targeted_fits.py, ...) run unchanged with
--source photon.
"""

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ANALYSIS_0804 = HERE.parent / "26-08-04-derived-hbt-time-reference"
sys.path.insert(0, str(ANALYSIS_0804))
from build_dataset_a import TAU_GRID  # noqa: E402
from extract_hbt_traces import (  # noqa: E402
    N_ROWS,
    PS_PER_PX,
    SUBWINDOW_EDGES,
    aggregate_band_traces,
    band_correlation,
    mean_sem,
)


def build_traces(table_path, band_windows):
    """Bin full-table events into (subwindow, tau) mean/SEM traces.

    Mirrors dataset_a_traces.build_traces; differences: the full curated
    event table has no per-event band column and no band_windows attribute
    (bands are assigned downstream), and no frame_indices group (frames are
    the raw 0..N-1 indices).
    """
    with h5py.File(table_path, "r") as handle:
        events = {
            name: np.asarray(handle["events"][name])
            for name in ("frame", "x", "tau")
        }
        frame_t0 = np.asarray(handle["frame_t0"])
        attrs = dict(handle.attrs)
    band_windows = np.asarray(band_windows, dtype=np.float64)
    frame_indices = np.arange(frame_t0.size)

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
    frame_ord = kept_ordinal[events["frame"]]
    sub = np.searchsorted(SUBWINDOW_EDGES, events["x"], side="right") - 1
    bin_index = np.floor(events["tau"] - TAU_GRID[0]).astype(np.int64)
    ok = (
        (frame_ord >= 0)
        & (sub >= 0)
        & (sub < n_sub)
        & (bin_index >= 0)
        & (bin_index < n_tau)
    )
    ok[ok] &= valid[frame_ord[ok], bin_index[ok]]

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
    parser.add_argument("--table", type=Path, required=True)
    parser.add_argument(
        "--bands",
        type=str,
        required=True,
        help="comma list of 2N column edges, e.g. 368,471,560,685,798,903,999,1118",
    )
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    values = [float(v) for v in args.bands.split(",") if v.strip()]
    if len(values) % 2 or not values:
        raise ValueError("--bands needs an even number of edges")
    band_windows = np.asarray(values).reshape(-1, 2)

    result = build_traces(args.table, band_windows)
    attrs = result.pop("attrs")
    n_overlap = result.pop("n_overlap")
    nan_traces = np.full_like(result["traces_photon"], np.nan)
    meta = {
        "source": "full_event_table_rebinned",
        "event_table": str(args.table),
        "band_windows": band_windows.tolist(),
        "analog_unavailable": True,
        "n_frames_used": int(result["t0_rows"].size),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
        "band_corr_overlap_bins": n_overlap,
        "table_attrs": {
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
