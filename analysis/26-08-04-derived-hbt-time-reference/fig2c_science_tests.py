"""G6 science verdicts for the 26-08-04 Fig. 2(c) g2 replica.

Implements the four pre-registered tests of docs/26-08-04-fig2c.md (written
before the real-data g2 computation, 2026-08-06):

  (a) bunching existence: per-bin z = (g2-1)/SEM_plugin, statistic max|z|
      over the 0-820 ps window;
  (b) oscillation: two-quadrature weighted sinusoid fit, period scan
      3-200 ps (plus targeted 9.84/37.2/41.8 ps), statistic max |A|/sigma_A;
  (c) anti-phase: zero-lag Pearson r of the detrended g2_L-1 / g2_R-1
      fluctuations (one-sided: anti-phase = significantly negative);
  (d) envelope growth: pooled (g2-1) late (400-820 ps) minus early
      (0-400 ps).

Null distributions come from 200 surrogates that cyclically shift one arm's
per-frame counts against the other (breaks same-frame coincidences, keeps
per-arm rates and t0 structure); the whole statistic pipeline is recomputed
per surrogate. p = (1 + #{surr >= obs}) / (1 + N_surr). Default bin width
4 px (finest MC-validated width); 2 px as cross-check (R validated, L
qualitative only).

Headline: window-pooled g2(0) +/- jackknife SEM for L and R with z vs 1.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2c_g2_estimator as est  # noqa: E402
from fig2c_power_gate import fit_matrices  # noqa: E402

DEFAULT_OUT_DIR = Path("drafts/26-08-06_fig2c")
WINDOW_PS = (0.0, 820.0)
EARLY_PS = (0.0, 400.0)
LATE_PS = (400.0, 820.0)
PRIOR_PERIODS_PS = (9.84, 37.2, 41.8)
N_SURROGATES = 200
# Scan grid: over an 820 ps window the frequency resolution is ~1.2e-3 /ps;
# the fig2a-era 0.5 ps spacing below 20 ps undersamples it by ~4x (a 9.84 ps
# line would lose ~80% of its amplitude at the nearest grid point). The
# null calibration is unaffected (surrogates see the same grid).
SCAN_PERIODS_PS = np.concatenate(
    [np.arange(3.0, 20.0, 0.1), np.arange(20.0, 200.0, 0.5)]
)


def load_width(out_dir, width):
    results = np.load(out_dir / f"g2_results_bin{width}px.npz")
    cube = np.load(out_dir / f"count_cube_bin{width}px.npy")
    return results, cube


def window_mask(t_ps, window):
    return (t_ps >= window[0]) & (t_ps <= window[1])


def per_bin_z_stat(g2, sem, mask):
    z = (g2 - 1.0) / sem
    z[~np.isfinite(z)] = 0.0
    return float(np.max(np.abs(z[mask]))), z


def fit_complex(t_ps, g2, sem, bin_ps, mats, periods=SCAN_PERIODS_PS):
    """Complex two-quadrature amplitudes per period, NaN-safe.

    Fast path uses the precomputed matrices (requires all-finite g2 on the
    SEM mask); otherwise falls back to a per-period fit dropping non-finite
    bins (mirrors the gate's slow path).
    """
    if np.all(np.isfinite(g2)):
        out = np.full(len(mats), np.nan + 0j)
        for k, entry in enumerate(mats):
            if entry is None:
                continue
            cov, xt_w, finite = entry
            beta_hat = cov @ (xt_w @ g2[finite])
            out[k] = beta_hat[3] + 1j * beta_hat[4]
        return out
    finite = np.isfinite(g2) & np.isfinite(sem) & (sem > 0.0)
    t = t_ps[finite]
    y = g2[finite]
    w = 1.0 / sem[finite] ** 2
    out = np.full(len(mats), np.nan + 0j)
    for k, period in enumerate(periods):
        attenuation = np.sinc(bin_ps / period)
        if abs(attenuation) < 0.05:
            continue  # degenerate under bin averaging; no sensitivity
        basis = attenuation * np.exp(2j * np.pi * t / period)
        design = np.column_stack([np.ones_like(t), t, t**2, basis.real, basis.imag])
        xt_w = design.T * w
        try:
            cov = np.linalg.inv(xt_w @ design)
        except np.linalg.LinAlgError:
            continue
        beta_hat = cov @ (xt_w @ y)
        out[k] = beta_hat[3] + 1j * beta_hat[4]
    return out


def amp_scan_stat(t_ps, g2, sem, bin_ps, mats, sigmas):
    """max over periods of |A|/sigma_A and the amplitude vector."""
    amps = np.abs(fit_complex(t_ps, g2, sem, bin_ps, mats))
    ratio = amps / sigmas
    return float(np.nanmax(ratio)), amps


def combined_amp_scan(t_ps, curves, sems, bin_ps, mats_l, sig_l, mats_r, sig_r):
    """Anti-phase combined amplitude: z = (z_L/sL^2 - z_R/sR^2)/(1/sL^2+1/sR^2)."""
    z_l = fit_complex(t_ps, curves[0], sems[0], bin_ps, mats_l)
    z_r = fit_complex(t_ps, curves[1], sems[1], bin_ps, mats_r)
    w_l = 1.0 / sig_l**2
    w_r = 1.0 / sig_r**2
    z = (z_l * w_l - z_r * w_r) / (w_l + w_r)
    sigma = 1.0 / np.sqrt(w_l + w_r)
    ratio = np.abs(z) / sigma
    return float(np.nanmax(ratio)), np.abs(z), sigma


def detrend_weighted(t, y, sem):
    finite = np.isfinite(y) & np.isfinite(sem) & (sem > 0)
    w = 1.0 / sem[finite] ** 2
    design = np.column_stack(
        [np.ones(finite.sum()), t[finite], t[finite] ** 2]
    )
    xt_w = design.T * w
    cov = np.linalg.pinv(xt_w @ design)
    beta = cov @ (xt_w @ y[finite])
    trend = np.column_stack([np.ones_like(t), t, t**2]) @ beta
    return y - trend, finite


def antiphase_r(g2_l, g2_r, sem_l, sem_r, t, mask):
    y_l, fin_l = detrend_weighted(t, g2_l - 1.0, sem_l)
    y_r, fin_r = detrend_weighted(t, g2_r - 1.0, sem_r)
    both = fin_l & fin_r & mask
    if both.sum() < 10:
        return float("nan")
    return float(np.corrcoef(y_l[both], y_r[both])[0, 1])


def envelope_stat(cube, beta, t_ps, pair):
    N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
    N1 = N1.astype(np.float64)
    N2 = N2.astype(np.float64)
    g_early, _, _ = est.pooled_zero_delay(
        N1, N2, b1, b2, mask=window_mask(t_ps, EARLY_PS)
    )
    g_late, _, _ = est.pooled_zero_delay(
        N1, N2, b1, b2, mask=window_mask(t_ps, LATE_PS)
    )
    return g_late - g_early


def shift_surrogate(N, rng):
    """Cyclic frame shift of one arm against the other (breaks coincidences)."""
    shift = int(rng.integers(1, N.shape[0]))
    return np.roll(N, shift, axis=0)


def evaluate_width(out_dir, width, n_surr, seed):
    results, cube = load_width(out_dir, width)
    t_ps = np.asarray(results["t_ps"], dtype=np.float64)
    beta = np.asarray(results["beta"], dtype=np.float64)
    mask = window_mask(t_ps, WINDOW_PS)
    bin_ps = float(np.median(np.diff(t_ps)))
    rng = np.random.default_rng(seed)

    out = {"bin_width_px": width, "bin_ps": bin_ps, "n_bins_window": int(mask.sum())}
    per_pair = {}
    mats = {}
    sigmas = {}
    curves = {}
    for pair in est.SAME_MESA_PAIRS:
        g2 = np.asarray(results[f"g2_{pair}"], dtype=np.float64)
        sem = np.asarray(results[f"sem_plugin_{pair}"], dtype=np.float64)
        curves[pair] = (g2, sem)
        mats[pair] = fit_matrices(t_ps, sem, bin_ps, SCAN_PERIODS_PS)
        sig = np.full(SCAN_PERIODS_PS.size, np.nan)
        for k, entry in enumerate(mats[pair]):
            if entry is not None:
                cov = entry[0]
                sig[k] = np.sqrt(max(cov[3, 3], cov[4, 4]))
        sigmas[pair] = sig

    # --- observed statistics ---
    obs = {}
    for pair in est.SAME_MESA_PAIRS:
        g2, sem = curves[pair]
        obs[f"maxz_{pair}"], _ = per_bin_z_stat(g2, sem, mask)
        obs[f"amp_{pair}"], amps = amp_scan_stat(
            t_ps, g2, sem, bin_ps, mats[pair], sigmas[pair]
        )
        per_pair[pair] = {"amps": amps}
    obs["amp_comb"], comb_amps, comb_sig = combined_amp_scan(
        t_ps,
        (curves["L"][0], curves["R"][0]),
        (curves["L"][1], curves["R"][1]),
        bin_ps, mats["L"], sigmas["L"], mats["R"], sigmas["R"],
    )
    obs["r_antiphase"] = antiphase_r(
        curves["L"][0], curves["R"][0],
        curves["L"][1], curves["R"][1], t_ps, mask,
    )
    for pair in est.SAME_MESA_PAIRS:
        obs[f"envelope_{pair}"] = envelope_stat(cube, beta, t_ps, pair)

    # --- surrogate null distributions ---
    surr = {key: [] for key in obs}
    for _ in range(n_surr):
        cube_s = cube.copy()
        # Shift one arm of each same-mesa pair (bands 2 and 3) against the rest.
        for band in (2, 3):
            cube_s[:, band, :] = shift_surrogate(cube_s[:, band, :], rng)
        g2_s = {}
        for pair in est.SAME_MESA_PAIRS:
            N1, N2, b1, b2 = est.pair_counts(cube_s, beta, pair)
            g2_s[pair], _ = est.g2_zero_delay(
                N1.astype(np.float64), N2.astype(np.float64), b1, b2
            )
            surr[f"envelope_{pair}"].append(
                envelope_stat(cube_s, beta, t_ps, pair)
            )
            stat_z, _ = per_bin_z_stat(g2_s[pair], curves[pair][1], mask)
            surr[f"maxz_{pair}"].append(stat_z)
            stat_a, _ = amp_scan_stat(
                t_ps, g2_s[pair], curves[pair][1], bin_ps, mats[pair], sigmas[pair]
            )
            surr[f"amp_{pair}"].append(stat_a)
        stat_c, _, _ = combined_amp_scan(
            t_ps,
            (g2_s["L"], g2_s["R"]),
            (curves["L"][1], curves["R"][1]),
            bin_ps, mats["L"], sigmas["L"], mats["R"], sigmas["R"],
        )
        surr["amp_comb"].append(stat_c)
        surr["r_antiphase"].append(
            antiphase_r(
                g2_s["L"], g2_s["R"], curves["L"][1], curves["R"][1], t_ps, mask
            )
        )

    # --- p values ---
    verdicts = {}
    for key, observed in obs.items():
        distribution = np.asarray(surr[key])
        distribution = distribution[np.isfinite(distribution)]
        if key == "r_antiphase":
            # One-sided: anti-phase = significantly NEGATIVE r.
            p = (1 + np.sum(distribution <= observed)) / (1 + distribution.size)
        else:
            p = (1 + np.sum(distribution >= observed)) / (1 + distribution.size)
        verdicts[key] = {
            "observed": float(observed),
            "null_median": float(np.median(distribution)),
            "null_q95": float(np.quantile(distribution, 0.95)),
            "p_value": float(p),
        }
    out["tests"] = verdicts

    # --- oscillation upper limits (fig2a convention: |A| + 2 sigma) ---
    # Targeted fits at the exact prior periods (the scan grid point nearest
    # to 9.84 ps would lose ~80% of a line's amplitude over this window).
    prior = np.asarray(PRIOR_PERIODS_PS, dtype=np.float64)
    for pair in est.SAME_MESA_PAIRS:
        amps = per_pair[pair]["amps"]
        ul95 = amps + 2.0 * sigmas[pair]
        mats_t = fit_matrices(t_ps, curves[pair][1], bin_ps, prior)
        z_t = fit_complex(t_ps, curves[pair][0], curves[pair][1], bin_ps,
                          mats_t, prior)
        sig_t = np.array([
            np.sqrt(max(entry[0][3, 3], entry[0][4, 4])) if entry else np.nan
            for entry in mats_t
        ])
        entry = {
            "ul95_median": float(np.nanmedian(ul95)),
            "targeted": {
                f"{period}ps": {
                    "amplitude": float(abs(z_t[j])),
                    "sigma": float(sig_t[j]),
                    "ul95": float(abs(z_t[j]) + 2.0 * sig_t[j]),
                }
                for j, period in enumerate(PRIOR_PERIODS_PS)
            },
        }
        per_pair[pair]["ul95"] = entry
    ul95_comb = comb_amps + 2.0 * comb_sig
    # combined targeted: inverse-variance combination of the exact fits
    mats_tl = fit_matrices(t_ps, curves["L"][1], bin_ps, prior)
    mats_tr = fit_matrices(t_ps, curves["R"][1], bin_ps, prior)
    z_tl = fit_complex(t_ps, curves["L"][0], curves["L"][1], bin_ps, mats_tl, prior)
    z_tr = fit_complex(t_ps, curves["R"][0], curves["R"][1], bin_ps, mats_tr, prior)
    sig_tl = np.array([
        np.sqrt(max(e[0][3, 3], e[0][4, 4])) if e else np.nan for e in mats_tl
    ])
    sig_tr = np.array([
        np.sqrt(max(e[0][3, 3], e[0][4, 4])) if e else np.nan for e in mats_tr
    ])
    w_l = 1.0 / sig_tl**2
    w_r = 1.0 / sig_tr**2
    z_comb = (z_tl * w_l - z_tr * w_r) / (w_l + w_r)
    sig_comb = 1.0 / np.sqrt(w_l + w_r)
    combined_entry = {
        "ul95_median": float(np.nanmedian(ul95_comb)),
        "targeted": {
            f"{period}ps": {
                "amplitude": float(abs(z_comb[j])),
                "sigma": float(sig_comb[j]),
                "ul95": float(abs(z_comb[j]) + 2.0 * sig_comb[j]),
            }
            for j, period in enumerate(PRIOR_PERIODS_PS)
        },
    }
    out["oscillation"] = {
        "periods_ps": [float(v) for v in SCAN_PERIODS_PS],
        "pairs": {
            pair: {
                "amplitudes": [float(v) for v in per_pair[pair]["amps"]],
                "ul95": per_pair[pair]["ul95"],
            }
            for pair in est.SAME_MESA_PAIRS
        },
        "combined": {
            "amplitudes": [float(v) for v in comb_amps],
            "ul95_median": combined_entry["ul95_median"],
            "targeted": combined_entry["targeted"],
        },
    }

    # --- headline pooled g2(0) ---
    pooled = {}
    for pair in est.SAME_MESA_PAIRS:
        N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
        g, s, n_used = est.pooled_zero_delay(
            N1.astype(np.float64), N2.astype(np.float64), b1, b2, mask=mask
        )
        pooled[pair] = {
            "g2": g, "sem": s, "z_vs_1": (g - 1.0) / s, "bins": n_used,
        }
    out["pooled"] = pooled

    # --- surrogate diagnostic figure ---
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)
    axes[0].semilogy(SCAN_PERIODS_PS, per_pair["L"]["amps"], lw=1.0, color="tab:blue",
                     label="pair L |A|")
    axes[0].semilogy(SCAN_PERIODS_PS, per_pair["R"]["amps"], lw=1.0, color="tab:red",
                     label="pair R |A|")
    axes[0].semilogy(SCAN_PERIODS_PS, comb_amps, lw=1.2, color="tab:purple",
                     label="anti-phase L+R |A|")
    axes[0].semilogy(SCAN_PERIODS_PS, ul95_comb, lw=0.9, ls="--", color="gray",
                     label="combined UL95 = |A|+2σ")
    for period in PRIOR_PERIODS_PS:
        axes[0].axvline(period, color="black", ls=":", lw=0.7)
    axes[0].set_xscale("log")
    axes[0].set_xticks([3, 10, 30, 100, 200])
    axes[0].set_xticklabels(["3", "10", "30", "100", "200"])
    axes[0].set_xlabel("period [ps]")
    axes[0].set_ylabel("g2-1 sinusoid amplitude")
    axes[0].set_title(f"Oscillation scan, {width} px bins")
    axes[0].legend(fontsize=8)
    axes[1].hist(surr["amp_comb"], bins=30, color="gray", alpha=0.7,
                 label="surrogate null")
    axes[1].axvline(obs["amp_comb"], color="tab:red", lw=1.4,
                    label=f"observed {obs['amp_comb']:.2f}")
    axes[1].set_xlabel("max over periods |A|/σ (combined)")
    axes[1].set_title(f"anti-phase oscillation null, p = {verdicts['amp_comb']['p_value']:.3f}")
    axes[1].legend(fontsize=8)
    figure.savefig(out_dir / f"fig2c_science_tests_bin{width}px.png", dpi=150)
    plt.close(figure)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--widths", type=int, nargs="+", default=[4, 2])
    parser.add_argument("--surrogates", type=int, default=N_SURROGATES)
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args()

    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "window_ps": list(WINDOW_PS),
        "early_ps": list(EARLY_PS),
        "late_ps": list(LATE_PS),
        "n_surrogates": int(args.surrogates),
        "criteria_doc": "docs/26-08-04-fig2c.md (pre-registered 2026-08-06)",
        "per_width": {},
    }
    for width in args.widths:
        print(f"evaluating {width} px bins with {args.surrogates} surrogates...",
              flush=True)
        report["per_width"][str(width)] = evaluate_width(
            args.out_dir, width, args.surrogates, args.seed + width
        )
        w = report["per_width"][str(width)]
        for pair in est.SAME_MESA_PAIRS:
            p = w["pooled"][pair]
            print(f"  {width} px pair {pair}: pooled g2(0) = {p['g2']:.4f} +/- "
                  f"{p['sem']:.4f} (z = {p['z_vs_1']:+.2f})")
        for key in ("amp_comb", "r_antiphase"):
            v = w["tests"][key]
            print(f"  {width} px {key}: observed {v['observed']:.3f}, "
                  f"null q95 {v['null_q95']:.3f}, p = {v['p_value']:.3f}")

    with open(args.out_dir / "fig2c_science_report.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"wrote {args.out_dir / 'fig2c_science_report.json'}")


if __name__ == "__main__":
    main()
