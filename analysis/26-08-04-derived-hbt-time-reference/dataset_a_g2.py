"""Fig. 2(c)-style g2(0)(t) from dataset A ONLY (no raw-frame access).

Mirrors compute_g2_fig2c.py but every input comes from the curated
dataset-A HDF5 (build_dataset_a.py):

  * per-frame count cubes are built from the curated events (stored band
    index; cross-checked against est.assign_bands on the stored x);
  * the residual-background beta per band is the file's OWN noise model
    (noise_model/floor_per_band, local-linear fit with 30-px signal-tail
    guards, tau-flat validated) times the band's member-column count and
    bin width — not a fresh estimate from out-of-band events, because
    dataset A deliberately contains none;
  * a second beta variant (the old pipeline's global out-of-band mean
    density, from noise_model/rho_tau_global) is carried through the
    pooled measurement to bracket the beta-choice systematic.

Outputs (same schema as compute_g2_fig2c.py, so fig2c_science_tests.py,
fig2c_synthetic_injection.py and render_fig2c.py run on them unchanged):
  count_cube_bin{1,2,4}px.npy, g2_results_bin{1,2,4}px.npz, g2_summary.json
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
from extract_hbt_traces import N_ROWS  # noqa: E402
from fig2c_power_gate import sem_plugin  # noqa: E402

DEFAULT_TABLE = (
    Path(__file__).resolve().parent
    / "detected-peaks/26-08-04_hbt_time_reference_dataset_a.h5"
)
DEFAULT_OUT_DIR = Path("drafts/26-08-08_dataset_a/g2")
BIN_WIDTHS = (1, 2, 4)
WINDOW_PS = (0.0, 820.0)


def load_dataset_a(path):
    with h5py.File(path, "r") as handle:
        events = {
            name: np.asarray(handle["events"][name])
            for name in ("frame", "x", "tau", "band")
        }
        frame_t0 = np.asarray(handle["frame_t0"])
        frame_indices = np.asarray(handle["frame_indices"])
        floor = np.asarray(handle["noise_model/floor_per_band"])
        rho_tau_global = np.asarray(handle["noise_model/rho_tau_global"])
        tau_grid = np.asarray(handle["noise_model/tau_grid"])
        band_windows = np.asarray(handle.attrs["band_windows"], dtype=np.float64)
        attrs = dict(handle.attrs)
    return {
        "events": events,
        "frame_t0": frame_t0,
        "frame_indices": frame_indices,
        "floor": floor,
        "rho_tau_global": rho_tau_global,
        "tau_grid": tau_grid,
        "band_windows": band_windows,
        "attrs": attrs,
    }


def build_cube(events, frame_t0, frame_indices, bin_width):
    """(kept_frames, 4, nbins) int16 cube from the curated events."""
    kept_mask = np.isfinite(frame_t0)
    kept_ordinal = np.full(frame_t0.size, -1, dtype=np.int64)
    kept_ordinal[kept_mask] = np.arange(int(kept_mask.sum()))
    frame_ord = kept_ordinal[np.searchsorted(frame_indices, events["frame"])]
    valid = frame_ord >= 0
    cube, nbins = est.count_cube_aligned(
        frame_ord[valid],
        events["band"][valid].astype(np.int64),
        events["tau"][valid],
        int(kept_mask.sum()),
        bin_width,
    )
    return cube.astype(np.int16), nbins, int(kept_mask.sum())


def beta_from_noise_model(floor, band_windows, nbins, bin_width):
    """beta[4, nbins] = per-row floor density x member columns x bin width.

    floor_per_band is events/frame/column integrated over the frame rows, so
    the per-row density is floor / N_ROWS; the tau-flat noise then gives
    floor/N_ROWS * width_cols * bin_width per bin.
    """
    widths = band_windows[:, 1] - band_windows[:, 0] + 1.0
    beta = floor[:, None] / N_ROWS * widths[:, None] * bin_width
    return np.broadcast_to(beta, (len(floor), nbins)).copy()


def beta_global_variant(rho_tau_global, tau_grid, band_windows, nbins, bin_width):
    """Old-pipeline-style beta: global out-of-band density per tau bin,
    rebinned to bin_width, times hi-lo window width."""
    rebinned = []
    for k in range(nbins):
        lo = est.TAU_RANGE_PX[0] + k * bin_width
        sel = (tau_grid >= lo) & (tau_grid < lo + bin_width)
        rebinned.append(float(rho_tau_global[sel].sum()))
    rho = np.asarray(rebinned)  # events/frame/col per bin
    widths = band_windows[:, 1] - band_windows[:, 0]
    return rho[None, :] * widths[:, None]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    data = load_dataset_a(args.table)
    events = data["events"]
    n_kept = int(np.isfinite(data["frame_t0"]).sum())

    # stored band index must match the estimator's x-based assignment
    est_band = est.assign_bands(events["x"], data["band_windows"])
    assert np.array_equal(est_band, events["band"].astype(np.int64)), (
        "stored band index disagrees with x-based assignment"
    )
    print(f"{args.table}: {events['frame'].size} curated events, "
          f"{n_kept} kept frames", flush=True)

    summary = {
        "table": str(args.table),
        "source": "dataset_a_only",
        "n_frames_kept": n_kept,
        "window_ps": list(WINDOW_PS),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "per_width": {},
    }

    for width in BIN_WIDTHS:
        cube, nbins, frames = build_cube(
            events, data["frame_t0"], data["frame_indices"], width)
        bin_ps = width * est.PS_PER_PX
        t_ps = (est.TAU_RANGE_PX[0] + (np.arange(nbins) + 0.5) * width) * est.PS_PER_PX
        window = (t_ps >= WINDOW_PS[0]) & (t_ps <= WINDOW_PS[1])

        beta = beta_from_noise_model(
            data["floor"], data["band_windows"], nbins, width)
        beta_old = beta_global_variant(
            data["rho_tau_global"], data["tau_grid"], data["band_windows"],
            nbins, width)

        # validation: cube frame-sum equals the direct band x bin histogram
        bin_index, _ = est.aligned_bin_index(events["tau"], width)
        kept = np.isfinite(data["frame_t0"])[
            np.searchsorted(data["frame_indices"], events["frame"])]
        valid = (events["band"] >= 0) & (bin_index >= 0) & kept
        direct = np.bincount(
            events["band"][valid].astype(np.int64) * nbins + bin_index[valid],
            minlength=4 * nbins,
        ).reshape(4, nbins)
        assert np.array_equal(cube.sum(axis=0), direct), "cube frame-sum mismatch"

        mean_counts = cube.astype(np.float64).mean(axis=0)
        width_entry = {
            "bin_width_px": width,
            "bin_width_ps": float(bin_ps),
            "nbins": int(nbins),
            "beta_source": "dataset_a noise_model/floor_per_band (local fit)",
            "beta_window_per_band": {
                name: float(beta[k][window].sum())
                for k, name in enumerate(est.BAND_NAMES)
            },
            "beta_window_per_band_global_variant": {
                name: float(beta_old[k][window].sum())
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
                N1, N2, b1, b2, mask=window)
            _, _, b1o, b2o = est.pair_counts(cube, beta_old, pair)
            pooled_old, pooled_sem_old, _ = est.pooled_zero_delay(
                N1, N2, b1o, b2o, mask=window)
            results[f"g2_{pair}"] = g2.astype(np.float32)
            results[f"sem_prop_{pair}"] = sem_prop.astype(np.float32)
            results[f"sem_plugin_{pair}"] = sem_plug.astype(np.float32)
            width_entry["pairs"][pair] = {
                "pooled_g2": pooled,
                "pooled_sem_jackknife": pooled_sem,
                "pooled_g2_global_beta_variant": pooled_old,
                "pooled_sem_jackknife_global_beta_variant": pooled_sem_old,
                "pooled_bins_used": n_used,
                "nan_bins_window": int(np.isnan(g2[window]).sum()),
                "sem_plugin_median_window": float(np.nanmedian(sem_plug[window])),
            }
            print(f"{width} px, pair {pair}: pooled g2(0) = {pooled:.4f} +/- "
                  f"{pooled_sem:.4f} (local-floor beta) | {pooled_old:.4f} +/- "
                  f"{pooled_sem_old:.4f} (global beta); median plug-in SEM "
                  f"{np.nanmedian(sem_plug[window]):.2f}", flush=True)

        np.savez(out_dir / f"g2_results_bin{width}px.npz", **results)
        cube_path = out_dir / f"count_cube_bin{width}px.npy"
        np.save(cube_path, cube)
        width_entry["cube_file"] = cube_path.name
        width_entry["cube_shape"] = list(cube.shape)
        summary["per_width"][str(width)] = width_entry

    summary["source_table_attrs"] = {
        key: (value.tolist() if hasattr(value, "tolist") else value)
        for key, value in data["attrs"].items()
    }
    with open(out_dir / "g2_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    print(f"wrote {out_dir / 'g2_summary.json'}", flush=True)


if __name__ == "__main__":
    main()
