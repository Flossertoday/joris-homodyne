"""Cross-mesa g2(0)(t) from dataset A ONLY (local-floor beta).

Companion to compute_g2_cross_mesa.py: same cross-arm pairs LR = B1 x B4 /
RL = B2 x B3, but every input comes from the curated dataset-A HDF5 and
the beta is the file's own noise model (local-linear floor per band), with
the old-pipeline global-beta variant carried through as a systematic
bracket — exactly the dataset-A conventions of dataset_a_g2.py. See that
module and docs/26-08-08-dataset-a-noise-model.md for the noise model.

Outputs: g2_cross_results_bin{1,2,4}px.npz + g2_cross_summary.json in
--out-dir (same key schema as compute_g2_cross_mesa.py).
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2c_g2_estimator as est  # noqa: E402
from dataset_a_g2 import (  # noqa: E402
    DEFAULT_TABLE,
    beta_from_noise_model,
    beta_global_variant,
    build_cube,
    load_dataset_a,
)
from fig2c_power_gate import sem_plugin  # noqa: E402

DEFAULT_OUT_DIR = Path("drafts/26-08-11_g2_cross_mesa/dataset_a")
BIN_WIDTHS = (1, 2, 4)
WINDOW_PS = (0.0, 820.0)


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
        "pairs": {pair: list(index) for pair, index in
                  est.CROSS_MESA_PAIRS.items()},
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
            "pairs": {},
        }

        results = {
            "t_ps": t_ps.astype(np.float32),
            "mean_counts": mean_counts.astype(np.float32),
            "beta": beta.astype(np.float32),
        }
        for pair in est.CROSS_MESA_PAIRS:
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

        np.savez(out_dir / f"g2_cross_results_bin{width}px.npz", **results)
        summary["per_width"][str(width)] = width_entry
        print(f"  saved g2_cross_results_bin{width}px.npz", flush=True)

    summary["source_table_attrs"] = {
        key: (value.tolist() if hasattr(value, "tolist") else value)
        for key, value in data["attrs"].items()
    }
    with open(out_dir / "g2_cross_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    print(f"wrote {out_dir / 'g2_cross_summary.json'}", flush=True)


if __name__ == "__main__":
    main()
