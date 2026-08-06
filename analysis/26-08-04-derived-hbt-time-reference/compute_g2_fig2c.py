"""G2+G4: count cubes and g2(0)(t) for the 26-08-04 Fig. 2(c) replica.

Reads the G1 event table (drafts/26-08-06_fig2c/hbt_event_table.h5), bins
events into the jitter-corrected coordinate tau = y - t0 at 1/2/4 px, and
computes the same-mesa zero-delay g2 curves:

  L pair = B1 x B3 (the two HBT arms of the left mesa)
  R pair = B2 x B4

Per width and pair the script stores
  * g2(0)(t) with the paper-style three-term propagated SEM (sample
    moments; degenerates to 0 in coincidence-free bins at this occupancy),
  * the plug-in Poisson SEM from the measured rates (validated against
    Monte Carlo in the G3 gate; used as the display band),
  * the window-pooled g2(0) with leave-one-frame-out jackknife SEM,
  * the residual-background beta per band and bin (from out-of-band events),
and saves the per-frame count cube itself (int16) for later 2-D/cross-mesa
extensions (docs/adiyatullin2017-fig2c.md G9).

G2 validations printed: cube summed over frames equals the per-band trace;
out-of-band bins fluctuate around the beta expectation; t0-failed frames
excluded and counted.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2c_g2_estimator as est  # noqa: E402
from fig2c_power_gate import sem_plugin  # noqa: E402

DEFAULT_TABLE = Path("drafts/26-08-06_fig2c/hbt_event_table.h5")
DEFAULT_OUT_DIR = Path("drafts/26-08-06_fig2c")
BIN_WIDTHS = (1, 2, 4)
WINDOW_PS = (0.0, 820.0)


def load_table(path):
    with h5py.File(path, "r") as handle:
        events = {name: handle["events"][name][:] for name in
                  ("frame", "x", "tau", "amplitude")}
        t0 = handle["frame_t0"][:]
        attributes = dict(handle.attrs)
    return events, t0, attributes


def build_cube(events, t0, bin_width):
    """(kept_frames, 4, nbins) int16 cube; frames are remapped to kept ordinals."""
    kept_index = np.flatnonzero(np.isfinite(t0))
    remap = np.full(t0.size, -1, dtype=np.int64)
    remap[kept_index] = np.arange(kept_index.size)
    band = est.assign_bands(events["x"])
    frame_ord = remap[events["frame"]]
    valid = frame_ord >= 0
    cube, nbins = est.count_cube_aligned(
        frame_ord[valid], band[valid], events["tau"][valid],
        kept_index.size, bin_width,
    )
    return cube.astype(np.int16), nbins, kept_index.size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    events, t0, attributes = load_table(args.table)
    n_kept = int(np.isfinite(t0).sum())
    print(f"{args.table}: {events['frame'].size} events, "
          f"{n_kept} kept frames ({t0.size - n_kept} t0 failures)")

    summary = {
        "table": str(args.table),
        "n_frames_kept": n_kept,
        "window_ps": list(WINDOW_PS),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "per_width": {},
    }

    for width in BIN_WIDTHS:
        cube, nbins, frames = build_cube(events, t0, width)
        bin_ps = width * est.PS_PER_PX
        t_ps = (est.TAU_RANGE_PX[0] + (np.arange(nbins) + 0.5) * width) * est.PS_PER_PX
        window = (t_ps >= WINDOW_PS[0]) & (t_ps <= WINDOW_PS[1])

        beta, rho, chi2_ndf = est.background_beta_aligned(
            events["tau"], events["x"], frames, width
        )

        # G2 validation: cube summed over frames equals the band tau histogram
        # of the kept-frame events (failed frames excluded as in build_cube).
        band = est.assign_bands(events["x"])
        bin_index, _ = est.aligned_bin_index(events["tau"], width)
        kept = np.isfinite(t0)
        valid = (band >= 0) & (bin_index >= 0) & kept[events["frame"]]
        flat = band[valid] * nbins + bin_index[valid]
        direct = np.bincount(flat, minlength=4 * nbins).reshape(4, nbins)
        cube_sum = cube.sum(axis=0)
        assert np.array_equal(cube_sum, direct), "cube frame-sum mismatch"

        mean_counts = cube.astype(np.float64).mean(axis=0)
        width_entry = {
            "bin_width_px": width,
            "bin_width_ps": float(bin_ps),
            "nbins": int(nbins),
            "background_flatness_chi2_ndf": chi2_ndf,
            "beta_window_per_band": {
                name: float(beta[k][window].sum())
                for k, name in enumerate(est.BAND_NAMES)
            },
            "mean_counts_window_per_band": {
                name: float(mean_counts[k][window].sum())
                for k, name in enumerate(est.BAND_NAMES)
            },
            "pairs": {},
        }

        results = {
            "t_ps": t_ps.astype(np.float32),
            "mean_counts": mean_counts.astype(np.float32),
            "beta": beta.astype(np.float32),
        }
        for pair in est.SAME_MESA_PAIRS:
            N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
            N1 = N1.astype(np.float64)
            N2 = N2.astype(np.float64)
            g2, sem_prop = est.g2_zero_delay(N1, N2, b1, b2)
            m1 = N1.mean(axis=0)
            m2 = N2.mean(axis=0)
            sem_plug = sem_plugin(m1, m2, b1, b2, frames)
            pooled, pooled_sem, n_used = est.pooled_zero_delay(
                N1, N2, b1, b2, mask=window
            )
            results[f"g2_{pair}"] = g2.astype(np.float32)
            results[f"sem_prop_{pair}"] = sem_prop.astype(np.float32)
            results[f"sem_plugin_{pair}"] = sem_plug.astype(np.float32)
            width_entry["pairs"][pair] = {
                "pooled_g2": pooled,
                "pooled_sem_jackknife": pooled_sem,
                "pooled_bins_used": n_used,
                "nan_bins_window": int(np.isnan(g2[window]).sum()),
                "sem_plugin_median_window": float(np.nanmedian(sem_plug[window])),
            }
            print(f"{width} px, pair {pair}: pooled g2(0) = {pooled:.4f} +/- "
                  f"{pooled_sem:.4f} (jackknife, {n_used} bins); median plug-in "
                  f"SEM in window {np.nanmedian(sem_plug[window]):.2f}; "
                  f"NaN bins {np.isnan(g2[window]).sum()}")

        np.savez(out_dir / f"g2_results_bin{width}px.npz", **results)
        cube_path = out_dir / f"count_cube_bin{width}px.npy"
        np.save(cube_path, cube)
        width_entry["cube_file"] = cube_path.name
        width_entry["cube_shape"] = list(cube.shape)
        summary["per_width"][str(width)] = width_entry
        print(f"  saved {cube_path} {cube.shape} and g2_results_bin{width}px.npz")

    summary["source_table_attrs"] = {
        key: (value.tolist() if hasattr(value, "tolist") else value)
        for key, value in attributes.items()
    }
    with open(out_dir / "g2_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    print(f"wrote {out_dir / 'g2_summary.json'}")


if __name__ == "__main__":
    main()
