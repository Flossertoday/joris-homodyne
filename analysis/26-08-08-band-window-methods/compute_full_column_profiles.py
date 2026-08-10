"""Full-column analog emission profiles for both 26-08-04 HBT batches.

Motivation: the photon-event column histogram (the input of all band-window
methods) only exists inside the detection column range 350-1150, so it cannot
answer whether that range cuts into the B1/B4 emission. The analog
(background-subtracted, row-summed) profile covers all 1344 columns; the
column direction is unaffected by the per-frame t0 jitter, so an UNALIGNED
sum suffices for column morphology.

- main batch: reuses the verified unaligned total sum
  drafts/26-08-04_fig2a/block_full/unaligned_total_sum.npy (20000 frames,
  raw sum) and subtracts n_frames x measured background from the H5;
- supplement batch: accumulates the plain sum from the supplement H5
  (72000 frames; same measured background, copied into that H5).

Row range: EMISSION_ROWS (32, 986) from extract_hbt_traces.py, matching the
pipeline's analog column-profile convention.

Output: full_column_profile_{main,supplement}.npy + a json sidecar with
provenance, written to --out-dir.
"""

import argparse
import json
import time
from pathlib import Path

import h5py
import numpy as np

MAIN_H5 = "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5"
SUPP_H5 = "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference_supplement.h5"
MAIN_SUM_NPY = "E:/Documents/Codex/joris-homodyne/drafts/26-08-04_fig2a/block_full/unaligned_total_sum.npy"
MAIN_SUM_FRAMES = 20000

EMISSION_ROWS = (32, 986)  # EMISSION_ROWS in extract_hbt_traces.py
BATCH = 100
BG_PATH = "conditions/background_and_streak_pulses_delay_calibration/background"
FRAMES_PATH = "conditions/hbt_condensate_josephson_regime/frames"


def accumulate_sum(h5_path, max_frames=None):
    with h5py.File(h5_path, "r") as h:
        background = h[BG_PATH][:].astype(np.float64)
        frames = h[FRAMES_PATH]
        n_total = frames.shape[0] if max_frames is None else min(max_frames, frames.shape[0])
        total = np.zeros(frames.shape[1:], np.float64)
        t0 = time.time()
        for start in range(0, n_total, BATCH):
            stop = min(start + BATCH, n_total)
            total += frames[start:stop].astype(np.float64).sum(axis=0)
            if stop % 5000 == 0 or stop == n_total:
                print(f"{stop}/{n_total} frames, {time.time()-t0:.0f}s", flush=True)
    return total, background, n_total


def column_profile(total_sum, background, n_frames):
    bg_sub = total_sum - n_frames * background
    return bg_sub[EMISSION_ROWS[0] : EMISSION_ROWS[1], :].sum(axis=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default="drafts/26-08-08_band_window_methods")
    ap.add_argument("--main-h5", default=MAIN_H5)
    ap.add_argument("--supp-h5", default=SUPP_H5)
    ap.add_argument("--main-sum-npy", default=MAIN_SUM_NPY)
    ap.add_argument("--skip-supplement", action="store_true")
    args = ap.parse_args()
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # main batch from the verified archived sum
    total = np.load(args.main_sum_npy)
    with h5py.File(args.main_h5, "r") as h:
        background = h[BG_PATH][:].astype(np.float64)
    prof = column_profile(total, background, MAIN_SUM_FRAMES)
    np.save(out / "full_column_profile_main.npy", prof)
    meta = {
        "main": {
            "source": args.main_sum_npy,
            "n_frames": MAIN_SUM_FRAMES,
            "background": f"{args.main_h5}:{BG_PATH}",
            "rows": list(EMISSION_ROWS),
        }
    }
    print("main profile: total %.3e counts" % prof.sum())

    if not args.skip_supplement:
        total, background, n = accumulate_sum(args.supp_h5)
        prof = column_profile(total, background, n)
        np.save(out / "full_column_profile_supplement.npy", prof)
        meta["supplement"] = {
            "source": args.supp_h5,
            "n_frames": n,
            "background": f"{args.supp_h5}:{BG_PATH}",
            "rows": list(EMISSION_ROWS),
        }
        print("supplement profile: total %.3e counts" % prof.sum())

    with open(out / "full_column_profiles_meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)


if __name__ == "__main__":
    main()
