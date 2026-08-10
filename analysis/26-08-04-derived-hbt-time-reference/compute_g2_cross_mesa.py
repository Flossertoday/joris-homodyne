"""Cross-mesa g2(0)(t) for the 26-08-04 HBT dataset (G9/§8.3 follow-up).

The same-mesa pairs L = B1 x B3 and R = B2 x B4 of compute_g2_fig2c.py are
the auto-correlations (self-g2) of condensates 1 and 2. This script computes
the missing condensate-1 x condensate-2 cross-correlation from the same
event tables with the same estimator (docs/adiyatullin2017-fig2c.md G9
item 4; project report §8.3), using the cross-arm pairs

  LR = B1 x B4  (condensate 1 in arm A x condensate 2 in arm B)
  RL = B2 x B3  (condensate 2 in arm A x condensate 1 in arm B)

Both estimate the same g2_LR(0)(t); cross-arm pairs keep the two detection
channels independent (same-arm cross pairs B1 x B2 / B3 x B4 would share
one streak trace and its fixed-pattern noise), so LR and RL are two
independent measurements of the cross-correlation — an internal
consistency check. Under independent condensates g2_LR = 1; number-
conserving Josephson population dynamics anti-correlates the two
condensates (g2_LR < 1).

Estimator, residual-background correction beta, plug-in SEM and the pooled
jackknife are exactly those of fig2c_g2_estimator; only the pair
definition differs. Set FIG2C_BAND_WINDOWS for the supplement batch as in
compute_g2_fig2c.py. Per bin width the script stores g2_LR/RL with the
paper-style propagated SEM, the plug-in Poisson SEM, the window-pooled
g2(0) with jackknife SEM, and writes g2_cross_results_bin{width}px.npz
plus g2_cross_summary.json into --out-dir.

Validation printed: cube summed over frames equals the per-band tau
histogram of kept-frame events (same G2 check as compute_g2_fig2c).
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2c_g2_estimator as est  # noqa: E402
from compute_g2_fig2c import DEFAULT_TABLE, build_cube, load_table  # noqa: E402
from fig2c_power_gate import sem_plugin  # noqa: E402

DEFAULT_OUT_DIR = Path("drafts/26-08-11_g2_cross_mesa")
BIN_WIDTHS = (1, 2, 4)
WINDOW_PS = (0.0, 820.0)


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
        "pairs": {pair: list(index) for pair, index in
                  est.CROSS_MESA_PAIRS.items()},
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

        # Same G2 validation as compute_g2_fig2c: the cube summed over
        # frames equals the band tau histogram of the kept-frame events.
        band = est.assign_bands(events["x"])
        bin_index, _ = est.aligned_bin_index(events["tau"], width)
        kept = np.isfinite(t0)
        valid = (band >= 0) & (bin_index >= 0) & kept[events["frame"]]
        flat = band[valid] * nbins + bin_index[valid]
        direct = np.bincount(flat, minlength=4 * nbins).reshape(4, nbins)
        assert np.array_equal(cube.sum(axis=0), direct), "cube frame-sum mismatch"

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
        for pair in est.CROSS_MESA_PAIRS:
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

        np.savez(out_dir / f"g2_cross_results_bin{width}px.npz", **results)
        summary["per_width"][str(width)] = width_entry
        print(f"  saved g2_cross_results_bin{width}px.npz")

    summary["source_table_attrs"] = {
        key: (value.tolist() if hasattr(value, "tolist") else value)
        for key, value in attributes.items()
    }
    with open(out_dir / "g2_cross_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    print(f"wrote {out_dir / 'g2_cross_summary.json'}")


if __name__ == "__main__":
    main()
