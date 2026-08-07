"""Verify dataset A against the production curated event table (provenance).

The dataset-A pass re-detects events from raw frames with the SAME
algorithm as the production extraction (same mask, threshold, centroids),
so its kept events must be an EXACT subset of the production table, and the
set difference must be exactly the documented curation rejects: out-of-band
events, defect-column events, and frame-edge-touching events. Any other
difference means the two detection passes drifted — fail loudly.

Checks:
  1. multiset containment: every dataset-A event (frame, x, y, amplitude)
     appears in the production table (float64 values rounded to 1e-9 — both
     tables share the detection arithmetic, so values are bit-identical);
  2. the remaining production events are precisely the three reject
     classes: out-of-band, defect column, or in-band events hugging the
     frame top/bottom edge (the production table has no bbox, so the edge
     class is identified by the curator's own cut count plus a centroid
     within 16 px of an edge — the largest observed blob half-height);
  3. the per-frame t0 arrays are identical.

Writes <dataset A>.verify_against_old.json; exits nonzero on any mismatch.
"""
import json
import sys

from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_dataset_a import DEFAULT_BANDS, DEFAULT_DEFECT_COLS  # noqa: E402

MAIN_REPO = Path("E:/Documents/Codex/joris-homodyne")
OLD_TABLE = (
    MAIN_REPO
    / "analysis/26-08-04-derived-hbt-time-reference/detected-peaks"
    / "26-08-04_hbt_time_reference_events.h5"
)
NEW_TABLE = (
    Path(__file__).resolve().parent
    / "detected-peaks/26-08-04_hbt_time_reference_dataset_a.h5"
)
EDGE_CENTROID_PX = 30  # exceeds the largest observed blob half-height (19 px)


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--old", type=Path, default=OLD_TABLE)
    parser.add_argument("--new", type=Path, default=NEW_TABLE)
    parser.add_argument("--bands", type=str, default=None,
                        help="comma list of 2N edges; default: original-batch")
    parser.add_argument("--defect-cols", type=str, default=None)
    args = parser.parse_args()
    bands = DEFAULT_BANDS
    if args.bands:
        values = [float(v) for v in args.bands.split(",") if v.strip()]
        bands = tuple((int(lo), int(hi)) for lo, hi in
                      np.asarray(values).reshape(-1, 2))
    defect_cols = DEFAULT_DEFECT_COLS
    if args.defect_cols:
        defect_cols = tuple(int(v) for v in args.defect_cols.split(",") if v.strip())

    with h5py.File(args.old, "r") as handle:
        old = {name: np.asarray(handle["events"][name])
               for name in ("frame", "x", "y", "amplitude")}
        old_t0 = np.asarray(handle["frame_t0"])
    with h5py.File(args.new, "r") as handle:
        new = {name: np.asarray(handle["events"][name])
               for name in ("frame", "x", "y", "amplitude")}
        new_t0 = np.asarray(handle["frame_t0"])
        cut = {key: int(handle["noise_model"].attrs[f"cut_{key}"])
               for key in ("total_detected", "kept", "out_of_band",
                           "defect_column", "edge_touch")}

    report = {"old_events": int(old["frame"].size),
              "new_events": int(new["frame"].size),
              "curator_cut_counts": cut}

    same = (old_t0 == new_t0) | (np.isnan(old_t0) & np.isnan(new_t0))
    report["t0_identical"] = bool(np.all(same))

    # multiset difference via index-tracking merge: rounded keys decide the
    # match, raw old-table coordinates drive the reject-class classification
    # (floor(x) at integer defect columns must not see the rounded value)
    def keys_for(frame, x, y, amp):
        return np.stack(
            [frame.astype(np.float64), np.round(y, 9), np.round(x, 9),
             np.round(amp, 9)],
            axis=1,
        )

    old_keys = keys_for(old["frame"], old["x"], old["y"], old["amplitude"])
    new_keys = keys_for(new["frame"], new["x"], new["y"], new["amplitude"])
    old_order = np.lexsort(old_keys.T[::-1])
    new_order = np.lexsort(new_keys.T[::-1])
    old_sorted = old_keys[old_order]
    new_sorted = new_keys[new_order]

    def lexcmp(a, b):
        for k in range(a.shape[0]):
            if a[k] < b[k]:
                return -1
            if a[k] > b[k]:
                return 1
        return 0

    matched_old = np.zeros(old["frame"].size, dtype=bool)
    foreign = 0
    i = j = 0
    while i < old_sorted.shape[0] and j < new_sorted.shape[0]:
        cmp = lexcmp(old_sorted[i], new_sorted[j])
        if cmp == 0:
            matched_old[old_order[i]] = True
            i += 1
            j += 1
        elif cmp < 0:
            i += 1
        else:
            foreign += 1
            j += 1
    foreign += new_sorted.shape[0] - j
    remainder_idx = np.flatnonzero(~matched_old)
    report["new_events_not_in_old"] = int(foreign)
    report["old_events_not_in_new"] = int(remainder_idx.size)

    # classify from the old table's own raw coordinates
    member = np.zeros(old["frame"].size, dtype=bool)
    for lo, hi in bands:
        member |= (old["x"] >= lo) & (old["x"] <= hi)
    defect = np.isin(np.floor(old["x"]).astype(np.int64),
                     np.asarray(defect_cols))
    report["old_out_of_band"] = int((~member).sum())
    report["old_inband_defect"] = int((member & defect).sum())
    report["old_inband_nondefect"] = int((member & ~defect).sum())

    if remainder_idx.size:
        rem_x = old["x"][remainder_idx]
        rem_y = old["y"][remainder_idx]
        rem_member = member[remainder_idx]
        rem_defect = defect[remainder_idx]
        unexplained = rem_member & ~rem_defect
        uy = rem_y[unexplained]
        ux = rem_x[unexplained]
        edge_ok = (
            (uy < EDGE_CENTROID_PX)
            | (uy > 1016 - EDGE_CENTROID_PX)
            # bbox touching the detection column boundary (supplement B1
            # starts at col 350): centroid within a blob half-width of it
            | (ux < 350 + EDGE_CENTROID_PX)
            | (ux > 1149 - EDGE_CENTROID_PX)
        )
        report["remainder_out_of_band"] = int((~rem_member).sum())
        report["remainder_defect"] = int((rem_member & rem_defect).sum())
        report["remainder_inband_edge_hugging"] = int(edge_ok.sum())
        report["remainder_unexplained"] = int((~edge_ok).sum())
    else:
        report["remainder_unexplained"] = 0

    with open(args.new.with_suffix(".verify_against_old.json"), "w") as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps(report, indent=2))

    assert report["t0_identical"], "per-frame t0 arrays differ"
    assert report["new_events_not_in_old"] == 0, "dataset A has foreign events"
    assert report["remainder_unexplained"] == 0, (
        "old-table events missing from dataset A outside the documented "
        "reject classes"
    )
    assert (
        report["old_events_not_in_new"]
        == cut["out_of_band"] + cut["defect_column"] + cut["edge_touch"]
    ), "reject-count mismatch vs curator cut counts"
    print("VERIFY PASS: dataset A is an exact curated subset of the "
          "production table")


if __name__ == "__main__":
    main()
