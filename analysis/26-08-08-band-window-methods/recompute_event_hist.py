"""Recompute the photon-event column histogram with a shifted detection left edge.

Motivation: the archived bands_report.json histograms (both 26-08-04 HBT
batches) were accumulated with EVENT_COL_RANGE = (350, 1150), so events with
centroids left of column 350 were never recorded and components straddling
column 350 were truncated by the detection mask. The full-column analog
profile showed the supplement B1 emission support extends to ~column 337
(docs/26-08-08-band-window-methods.md, section on the 350/1150 provenance),
so the 350 left edge clips real signal. This script re-runs ONLY the event
histogram accumulation with a wider column range (default 330-1150), keeping
every other pipeline choice bit-identical:

- same H5 frames + measured background (extract_hbt_traces datasets);
- same per-frame t0 gating (cfd50 estimate on T0_ROI, failures dropped);
- same detect_events call (threshold 25.0, min_size 2, min_counts 0.0),
  only ``col_range`` widened;
- histogram accumulated over all 1344 columns exactly as in accumulate().

Because every frame is read anyway, the full-column analog profile
(rows 32-986 summed, bg-subtracted; the compute_full_column_profiles.py
product) is accumulated in the same pass at negligible extra cost.

QA cross-check against the archived 350-1150 report: outside the columns a
widened left edge can touch (>= ~370), the new histogram must equal the
archived one exactly, and the t0-kept frame count must match the archived
n_frames_used.

Outputs (--out-dir, default drafts/26-08-08_band_window_methods_col330):
- bands_report_col330_<batch>.json   band_detection.event_col_hist + provenance
- full_column_profile_<batch>.npy    full-column analog profile
- full_column_profiles_meta_col330.json

Run (from the repository/worktree root):

    python analysis/26-08-08-band-window-methods/recompute_event_hist.py \
        --batch main --col-range 330 1150
    python analysis/26-08-08-band-window-methods/recompute_event_hist.py \
        --batch supplement --col-range 330 1150
"""

import argparse
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

HBT_DIR = (
    Path(__file__).resolve().parent.parent / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(HBT_DIR))
import extract_hbt_traces as eht  # noqa: E402

MAIN_H5 = "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5"
SUPP_H5 = "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference_supplement.h5"
MAIN_REPORT_350 = "E:/Documents/Codex/joris-homodyne/drafts/26-08-06_fig2a_cfd50/bands_report.json"
SUPP_REPORT_350 = "E:/Documents/Codex/joris-homodyne/drafts/26-08-06_fig2a_supplement/bands_report.json"

BATCHES = {
    "main": {"h5": MAIN_H5, "report350": MAIN_REPORT_350},
    "supplement": {"h5": SUPP_H5, "report350": SUPP_REPORT_350},
}

THRESHOLD = 25.0  # threshold.chosen in both archived reports
BATCH = eht.BATCH_FRAMES
EMISSION_ROWS = eht.EMISSION_ROWS


def run_batch(h5_path, col_range, threshold, max_frames=None):
    """One pass: event column histogram + full-column analog profile."""
    estimate = eht.T0_METHODS["cfd50"]
    with h5py.File(h5_path, "r") as handle:
        background = np.asarray(handle[eht.BACKGROUND_DATASET], dtype=np.float64)
        frames = handle[eht.FRAMES_DATASET]
        n_total = frames.shape[0] if max_frames is None else min(max_frames, frames.shape[0])
        n_rows, n_cols = frames.shape[1], frames.shape[2]
        hist = np.zeros(n_cols, dtype=np.int64)
        total_raw = np.zeros((n_rows, n_cols), dtype=np.float64)
        n_used = n_failures = 0
        started = time.time()
        for start in range(0, n_total, BATCH):
            stop = min(start + BATCH, n_total)
            block = frames[start:stop].astype(np.float64)
            total_raw += block.sum(axis=0)
            for i in range(stop - start):
                frame_bg = block[i] - background
                t0 = estimate(frame_bg, roi=eht.T0_ROI)
                shift = int(round(t0)) if np.isfinite(t0) else -1
                if not 0 <= shift <= eht.CANVAS_OFF:
                    n_failures += 1
                    continue
                rows, cols = eht.detect_events(
                    frame_bg, threshold, col_range=col_range
                )
                if rows.size:
                    np.add.at(
                        hist, np.clip(cols.astype(np.int64), 0, n_cols - 1), 1
                    )
                n_used += 1
            if stop % 2000 < BATCH or stop == n_total:
                print(
                    f"  {stop}/{n_total} frames, {n_failures} t0 failures, "
                    f"{time.time() - started:.0f} s elapsed",
                    flush=True,
                )
    profile = (
        total_raw - n_total * background
    )[EMISSION_ROWS[0] : EMISSION_ROWS[1], :].sum(axis=0)
    return {
        "hist": hist,
        "profile": profile,
        "n_frames": n_total,
        "n_used": n_used,
        "n_failures": n_failures,
        "elapsed_s": time.time() - started,
    }


def qa_against_350(hist, report350_path, n_used, col_range):
    """Compare with the archived 350-1150 histogram where both must agree."""
    old = json.load(open(report350_path))
    old_hist = np.asarray(old["band_detection"]["event_col_hist"], dtype=np.int64)
    qa = {
        "archived_report": report350_path,
        "archived_n_frames_used": old["t0"]["n_frames_used"],
        "n_frames_used_match": bool(n_used == old["t0"]["n_frames_used"]),
    }
    # components touching the widened left edge can change centroid/extent;
    # right of col_range[0] + 40 the two histograms must be identical
    stable = slice(col_range[0] + 40, eht.EVENT_COL_RANGE[1])
    diff = np.nonzero(hist[stable] != old_hist[stable])[0]
    qa["stable_region"] = [col_range[0] + 40, eht.EVENT_COL_RANGE[1]]
    qa["stable_region_diff_columns"] = diff.size
    qa["stable_region_diff_abs_events"] = int(
        np.abs(hist[stable] - old_hist[stable]).sum()
    )
    qa["new_region_events"] = int(hist[col_range[0] : eht.EVENT_COL_RANGE[0]].sum())
    return qa


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch", choices=sorted(BATCHES), required=True)
    ap.add_argument("--col-range", type=int, nargs=2, default=[330, 1150])
    ap.add_argument("--threshold", type=float, default=THRESHOLD)
    ap.add_argument("--max-frames", type=int, default=None)
    ap.add_argument(
        "--out-dir", default="drafts/26-08-08_band_window_methods_col330"
    )
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    info = BATCHES[args.batch]
    col_range = tuple(args.col_range)
    print(f"== {args.batch}: {info['h5']}, col_range={col_range}")
    res = run_batch(info["h5"], col_range, args.threshold, args.max_frames)
    print(
        f"  kept {res['n_used']}/{res['n_frames']} frames "
        f"({res['n_failures']} t0 failures), {res['elapsed_s']:.0f} s"
    )

    qa = None
    if args.max_frames is None:
        qa = qa_against_350(res["hist"], info["report350"], res["n_used"], col_range)
        print("  QA vs archived 350 report:", json.dumps(qa, indent=2))

    tag = f"col{col_range[0]}"
    report = {
        "band_detection": {
            "method": (
                "photon-event column histogram, recomputed with widened "
                "detection column range; all other detection parameters "
                "identical to the archived 350-1150 extraction"
            ),
            "col_range": list(col_range),
            "threshold": args.threshold,
            "min_size": 2,
            "min_counts": 0.0,
            "event_col_hist": res["hist"].tolist(),
        },
        "t0": {
            "method": "cfd50",
            "n_frames": res["n_frames"],
            "n_frames_used": res["n_used"],
            "n_failures": res["n_failures"],
        },
        "provenance": {
            "h5": info["h5"],
            "frames_dataset": eht.FRAMES_DATASET,
            "background_dataset": eht.BACKGROUND_DATASET,
            "script": str(Path(__file__).resolve()),
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "elapsed_s": res["elapsed_s"],
        },
        "qa_vs_archived_350": qa,
    }
    report_path = out / f"bands_report_{tag}_{args.batch}.json"
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print("wrote", report_path)

    np.save(out / f"full_column_profile_{args.batch}.npy", res["profile"])
    meta_path = out / f"full_column_profiles_meta_{tag}.json"
    meta = {}
    if meta_path.exists():
        meta = json.load(open(meta_path))
    meta[args.batch] = {
        "source": info["h5"],
        "n_frames": res["n_frames"],
        "rows": list(EMISSION_ROWS),
        "note": "same pass as the event histogram recompute",
    }
    with open(meta_path, "w") as fh:
        json.dump(meta, fh, indent=2)
    print("wrote", out / f"full_column_profile_{args.batch}.npy")


if __name__ == "__main__":
    main()
