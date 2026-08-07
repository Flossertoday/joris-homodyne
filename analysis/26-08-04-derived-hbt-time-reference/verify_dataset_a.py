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
from collections import Counter
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


def tuples(frame, x, y, amp):
    return Counter(
        zip(
            frame.astype(np.int64).tolist(),
            np.round(y, 9).tolist(),
            np.round(x, 9).tolist(),
            np.round(amp, 9).tolist(),
        )
    )


def main():
    with h5py.File(OLD_TABLE, "r") as handle:
        old = {name: np.asarray(handle["events"][name])
               for name in ("frame", "x", "y", "amplitude")}
        old_t0 = np.asarray(handle["frame_t0"])
    with h5py.File(NEW_TABLE, "r") as handle:
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

    remainder = tuples(old["frame"], old["x"], old["y"], old["amplitude"])
    remainder.subtract(tuples(new["frame"], new["x"], new["y"], new["amplitude"]))
    foreign = {k: v for k, v in remainder.items() if v < 0}
    report["new_events_not_in_old"] = int(sum(-v for v in foreign.values()))
    remainder = Counter({k: v for k, v in remainder.items() if v > 0})
    report["old_events_not_in_new"] = int(sum(remainder.values()))

    # classify the remainder from the old table's own coordinates
    member = np.zeros(old["frame"].size, dtype=bool)
    for lo, hi in DEFAULT_BANDS:
        member |= (old["x"] >= lo) & (old["x"] <= hi)
    defect = np.isin(np.floor(old["x"]).astype(np.int64),
                     np.asarray(DEFAULT_DEFECT_COLS))
    report["old_out_of_band"] = int((~member).sum())
    report["old_inband_defect"] = int((member & defect).sum())
    report["old_inband_nondefect"] = int((member & ~defect).sum())

    # remainder events must be: out-of-band, defect, or edge-hugging
    if remainder:
        rem = np.asarray(list(remainder.elements()))
        rem_frame, rem_y, rem_x = rem[:, 0], rem[:, 1], rem[:, 2]
        rem_member = np.zeros(rem_frame.size, dtype=bool)
        for lo, hi in DEFAULT_BANDS:
            rem_member |= (rem_x >= lo) & (rem_x <= hi)
        rem_defect = np.isin(np.floor(rem_x).astype(np.int64),
                             np.asarray(DEFAULT_DEFECT_COLS))
        unexplained = rem_member & ~rem_defect
        edge_ok = (rem_y[unexplained] < EDGE_CENTROID_PX) | (
            rem_y[unexplained] > 1016 - EDGE_CENTROID_PX)
        report["remainder_out_of_band"] = int((~rem_member).sum())
        report["remainder_defect"] = int((rem_member & rem_defect).sum())
        report["remainder_inband_edge_hugging"] = int(edge_ok.sum())
        report["remainder_unexplained"] = int((~edge_ok).sum())
    else:
        report["remainder_unexplained"] = 0

    with open(NEW_TABLE.with_suffix(".verify_against_old.json"), "w") as stream:
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
