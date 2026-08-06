"""Extract the per-photon event table for the 26-08-04 Fig. 2(c) g2 replica.

G1 of docs/adiyatullin2017-fig2c.md: re-run the Fig. 2(a) detection chain
(measured-background subtraction, CFD50 t0, 8-connected event detection at
25 counts / min 2 px) on all 20000 HBT frames and store every detected
event: (frame index, subpixel x, subpixel y, integrated amplitude), plus the
per-frame t0. The Fig. 2(a) pipeline only kept accumulated traces; the g2
estimator needs per-frame coincidences, hence the event table.

Validation (pre-registered in the same doc): band-summed photon traces
rebuilt from this table must agree with the Fig. 2(a) accumulated photon
traces (drafts/26-08-06_fig2a_cfd50/hbt_traces.npz) within numerical
tolerance; t0 failures are excluded and recorded.

Output: HDF5 with datasets
  events/frame  int32    frame index of each event
  events/x      float32  subpixel column (detector px)
  events/y      float32  subpixel row (detector px, absolute)
  events/tau    float32  y - t0(frame)  (aligned coordinate, px)
  events/amplitude float32 integrated counts above background
  frame_t0      float64  per-frame CFD50 t0 row (NaN on failure)
  failed_frames int32    frames excluded (t0 failure)
plus parameter/provenance attributes.
"""
import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import extract_hbt_traces as eht  # noqa: E402
import fig2c_g2_estimator as est  # noqa: E402

DEFAULT_OUT = Path("drafts/26-08-06_fig2c/hbt_event_table.h5")
DEFAULT_NPZ = Path(
    "E:/Documents/Codex/joris-homodyne/drafts/26-08-06_fig2a_cfd50/hbt_traces.npz"
)
BATCH = 100
THRESHOLD = 25.0  # Fig. 2(a) final detection parameter (knee of the study)
TAU_GRID0 = eht.TAU_GRID[0]  # -64: bin origin of the Fig. 2(a) accumulation


def extract(h5_path, threshold=THRESHOLD, max_frames=None, progress_every=2000):
    """Single pass: returns event arrays, per-frame t0, failed frame list."""
    background = eht.load_background(h5_path)
    with h5py.File(h5_path, "r") as handle:
        n_frames = handle[eht.FRAMES_DATASET].shape[0]
    if max_frames is not None:
        n_frames = min(n_frames, max_frames)
    ev_frame = []
    ev_x = []
    ev_y = []
    ev_amp = []
    t0 = np.full(n_frames, np.nan)
    started = time.time()
    with h5py.File(h5_path, "r") as handle:
        frames = handle[eht.FRAMES_DATASET]
        for start in range(0, n_frames, BATCH):
            block = frames[start : min(start + BATCH, n_frames)]
            for offset, raw in enumerate(block):
                index = start + offset
                frame_bg = raw.astype(np.float64) - background
                t = eht.estimate_t0_cfd(frame_bg)
                if np.isfinite(t):
                    t0[index] = t
                    rows, cols, amps = eht.detect_events(
                        frame_bg, threshold, with_counts=True
                    )
                    if rows.size:
                        ev_frame.append(
                            np.full(rows.size, index, dtype=np.int32)
                        )
                        ev_x.append(cols)
                        ev_y.append(rows)
                        ev_amp.append(amps)
            done = min(start + BATCH, n_frames)
            if done % progress_every < BATCH or done == n_frames:
                print(f"  {done}/{n_frames} frames, "
                      f"{int(np.isnan(t0).sum())} t0 failures, "
                      f"{time.time() - started:.0f} s elapsed", flush=True)
    if ev_frame:
        events = {
            "frame": np.concatenate(ev_frame),
            "x": np.concatenate(ev_x),
            "y": np.concatenate(ev_y),
            "amplitude": np.concatenate(ev_amp),
        }
    else:
        events = {
            "frame": np.empty(0, np.int32),
            "x": np.empty(0),
            "y": np.empty(0),
            "amplitude": np.empty(0),
        }
    events["tau"] = events["y"] - t0[events["frame"]]
    return events, t0


def write_h5(events, t0, out_path, attributes):
    """Write the event table, refusing to overwrite; atomic via .partial."""
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing file: {out_path}")
    partial = out_path.with_name(out_path.name + ".partial")
    if partial.exists():
        raise FileExistsError(f"Interrupted artifact already exists: {partial}")
    failed = np.flatnonzero(np.isnan(t0)).astype(np.int32)
    try:
        with h5py.File(partial, "x") as handle:
            for key, value in attributes.items():
                handle.attrs[key] = value
            group = handle.create_group("events")
            for name in ("frame", "x", "y", "tau", "amplitude"):
                group.create_dataset(
                    name, data=events[name], compression="gzip", shuffle=True
                )
            handle.create_dataset("frame_t0", data=t0)
            handle.create_dataset("failed_frames", data=failed)
        if out_path.exists():
            raise FileExistsError(f"Refusing to overwrite existing file: {out_path}")
        partial.rename(out_path)
    except BaseException:
        if partial.exists():
            partial.unlink()
        raise
    return out_path


def validate_against_npz(events, t0, npz_path):
    """Rebuild band photon traces from the event table and compare to the
    Fig. 2(a) accumulated npz (same binning as extract_hbt_traces.accumulate:
    bin = floor(y - t0 - tau_grid0 + 0.5) over the 1-px tau grid).
    Returns per-band max abs difference over the 0-820 ps window.

    Tolerance note: isolated bins can differ by 1-2 counts (float64
    bin-edge ties; measured 32 of 165k in-window band events, 0.019%),
    so the pass criterion is 2/kept frames, not bit equality.
    """
    data = np.load(npz_path, allow_pickle=True)
    tau_ps = np.asarray(data["time_ps"])
    reference = np.nan_to_num(np.asarray(data["band_traces_photon"]), nan=0.0)
    n_tau = tau_ps.size
    kept = int(np.sum(np.isfinite(t0)))
    band = est.assign_bands(events["x"])
    bin_index = np.floor(events["tau"] - TAU_GRID0 + 0.5).astype(np.int64)
    valid = (band >= 0) & (bin_index >= 0) & (bin_index < n_tau)
    flat = band[valid] * n_tau + bin_index[valid]
    rebuilt = np.bincount(flat, minlength=4 * n_tau).reshape(4, n_tau) / kept
    window = (tau_ps >= 0.0) & (tau_ps <= 820.0)
    max_diff = np.abs(rebuilt - reference)[:, window].max(axis=1)
    totals = [int((band == k).sum()) for k in range(4)]
    return max_diff, totals, kept


def save_sample_overlay(h5_path, background, t0, frame_indices, out_png):
    """Overlay detected events on two sample frames (band crops)."""
    n = len(frame_indices)
    figure, axes = plt.subplots(
        1, n, figsize=(6.2 * n, 5.2), constrained_layout=True
    )
    if n == 1:
        axes = [axes]
    with h5py.File(h5_path, "r") as handle:
        frames = handle[eht.FRAMES_DATASET]
        for axis, index in zip(axes, frame_indices):
            frame_bg = np.asarray(frames[int(index)], dtype=np.float64) - background
            rows, cols, amps = eht.detect_events(
                frame_bg, THRESHOLD, with_counts=True
            )
            crop = frame_bg[:, est.EVENT_COL_RANGE[0] : est.EVENT_COL_RANGE[1]]
            vmax = float(np.percentile(crop, 99.9))
            image = axis.imshow(
                crop, origin="upper", aspect="auto", cmap="magma",
                vmin=0.0, vmax=vmax,
            )
            figure.colorbar(image, ax=axis, shrink=0.85).set_label(
                f"counts above bg (magma, linear, 0..p99.9={vmax:.0f})"
            )
            axis.plot(
                cols - est.EVENT_COL_RANGE[0], rows, ".", ms=4,
                color="cyan", mfc="none", label="events",
            )
            for k, (lo, hi) in enumerate(est.BAND_WINDOWS):
                axis.axvline(lo - est.EVENT_COL_RANGE[0], color="white",
                             lw=0.5, alpha=0.4)
                axis.axvline(hi - est.EVENT_COL_RANGE[0], color="white",
                             lw=0.5, alpha=0.4)
            axis.axhline(t0[int(index)], color="lime", lw=0.8,
                         label=f"t0 = {t0[int(index)]:.1f}")
            axis.set_title(f"frame {index}: {rows.size} events", fontsize=9)
            axis.set_xlabel("column - 350 [px]")
            axis.set_ylabel("row [px]")
            axis.legend(fontsize=7, loc="upper right")
    figure.savefig(out_png, dpi=150)
    plt.close(figure)
    return out_png


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", type=Path, default=eht.DEFAULT_H5)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ,
                        help="Fig. 2(a) traces for the rebuild validation")
    parser.add_argument("--threshold", type=float, default=THRESHOLD)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--sample-frames", type=int, nargs="+",
                        default=[100, 5000, 12000])
    parser.add_argument("--validate-only", action="store_true",
                        help="skip extraction; re-run the npz rebuild validation "
                             "on the existing --out table")
    args = parser.parse_args()

    if args.validate_only:
        with h5py.File(args.out, "r") as handle:
            events = {name: handle["events"][name][:]
                      for name in ("frame", "x", "tau")}
            t0 = handle["frame_t0"][:]
        max_diff, totals, kept = validate_against_npz(events, t0, args.npz)
        print(f"rebuild validation vs {args.npz.name} (0-820 ps window, "
              f"{kept} frames):")
        for k, name in enumerate(est.BAND_NAMES):
            print(f"  {name}: {totals[k]} events, max |rebuilt - npz| "
              f"= {max_diff[k]:.3e} events/frame/bin "
              f"({max_diff[k] * kept:.2f} counts)")
        passed = bool(np.all(max_diff * kept <= 2.0))
        print(f"validation {'PASS' if passed else 'FAIL'} (tolerance 2 counts/bin)")
        return

    print(f"Extracting event table from {args.h5} (threshold {args.threshold:g})",
          flush=True)
    events, t0 = extract(args.h5, args.threshold, args.max_frames)
    n_events = events["frame"].size
    print(f"{n_events} events from {np.isfinite(t0).sum()} frames "
          f"({int(np.isnan(t0).sum())} t0 failures)", flush=True)

    attributes = {
        "source_h5": str(args.h5),
        "source_h5_sha256": eht.H5_SHA256,
        "frames_dataset": eht.FRAMES_DATASET,
        "background_dataset": eht.BACKGROUND_DATASET,
        "threshold": float(args.threshold),
        "min_size": 2,
        "col_range": list(eht.EVENT_COL_RANGE),
        "t0_method": "cfd50",
        "cfd_fraction": eht.CFD_FRACTION,
        "t0_roi": list(eht.T0_ROI),
        "n_frames": int(t0.size),
        "n_frames_kept": int(np.isfinite(t0).sum()),
        "n_events": int(n_events),
        "band_windows": est.BAND_WINDOWS.tolist(),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
    }
    out = write_h5(events, t0, args.out, attributes)
    print(f"wrote {out}", flush=True)

    # G1 validation: rebuilt band traces vs the Fig. 2(a) npz. Only
    # meaningful on the full extraction (same frame set as the npz).
    if args.npz.exists() and args.max_frames is None:
        max_diff, totals, kept = validate_against_npz(events, t0, args.npz)
        print(f"rebuild validation vs {args.npz.name} (0-820 ps window, "
              f"{kept} frames):")
        for k, name in enumerate(est.BAND_NAMES):
            print(f"  {name}: {totals[k]} events, max |rebuilt - npz| "
              f"= {max_diff[k]:.3e} events/frame/bin")
        check = {
            "npz": str(args.npz),
            "kept_frames": kept,
            "band_event_totals": dict(zip(est.BAND_NAMES, totals)),
            "max_abs_diff": dict(zip(est.BAND_NAMES, [float(v) for v in max_diff])),
            # 2 counts per bin: float64 bin-edge ties between independent
            # summation orders (measured 32 of 165k in-window band events).
            "pass": bool(np.all(max_diff * kept <= 2.0)),
        }
        check_path = Path(args.out).with_name("event_table_validation.json")
        with open(check_path, "w") as handle:
            json.dump(check, handle, indent=2)
        print(f"validation {'PASS' if check['pass'] else 'FAIL'} "
              f"(tolerance 2 counts/bin); wrote {check_path}", flush=True)

    # Sample-frame overlay for visual inspection.
    background = eht.load_background(args.h5)
    overlay = Path(args.out).with_name("event_table_sample_overlay.png")
    sample = [i for i in args.sample_frames if i < t0.size]
    print(f"wrote {save_sample_overlay(args.h5, background, t0, sample, overlay)}",
          flush=True)


if __name__ == "__main__":
    main()
