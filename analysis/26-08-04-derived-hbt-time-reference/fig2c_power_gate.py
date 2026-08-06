"""G3 statistical-power gate for the 26-08-04 Fig. 2(c) g2 replica.

Pre-registered in docs/adiyatullin2017-fig2c.md (G3): before any full
computation, quantify what this dataset can resolve. All inputs come from
the Fig. 2(a) accumulated traces (drafts/26-08-06_fig2a_cfd50/hbt_traces.npz):

  * in-band photon-event traces give the per-bin TOTAL event rates
    m_c(tau) (signal + residual detector-noise events) per band c;
  * out-of-band sub-window photon traces give the residual noise density
    rho(tau) per x-pixel, hence the background expectation beta_c(tau);
  * signal rates s_c = m_c - beta_c.

Products
  1. analytic per-bin SEM(g2)(t) per same-mesa pair (L = B1 x B3,
     R = B2 x B4) at 1/2/4 px binning, with a Poisson plug-in variance
     (the paper-style sample-moment SEM degenerates to 0 in
     coincidence-free bins at this occupancy; its ensemble rms is checked
     against Monte Carlo instead);
  2. Monte Carlo validation: Poisson cubes at the measured rates through
     fig2c_g2_estimator; MC spread vs analytic SEM (criterion: <10%);
  3. window-pooled g2(0) precision (jackknife, MC-checked);
  4. minimum resolvable g2-1 sinusoid amplitude (95%, phase-quadrature
     weighted fit on the g2 curve with quadratic detrend, period scan
     3-200 ps, bin-averaging attenuation included) vs the paper's 0.2-1.2;
  5. frame counts required to cover the paper magnitudes.

Gate verdict (docs/adiyatullin2017-fig2c.md section 5, G3): clear pass if
the minimum resolvable amplitude covers the paper's lower edge (0.2);
borderline up to 1.2; beyond that the target degrades to an upper-limit
report and the pipeline stops for user confirmation.
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

DEFAULT_NPZ = Path(
    "E:/Documents/Codex/joris-homodyne/drafts/26-08-06_fig2a_cfd50/hbt_traces.npz"
)
DEFAULT_OUT_DIR = Path("drafts/26-08-06_fig2c")

BIN_WIDTHS = (1, 2, 4)
WINDOW_PS = (0.0, 820.0)  # Fig. 2(a) convention
PAPER_AMPLITUDE_RANGE = (0.2, 1.2)  # g2-1: early peaks to late peaks
# Period scan grid identical to the Fig. 2(a) z-trace scan.
PERIODS_PS = np.concatenate([np.arange(3.0, 20.0, 0.5), np.arange(20.0, 200.0, 2.0)])
# Expected-null 95% upper limit factor on a two-quadrature amplitude:
# E|A| + 2 sigma with E|A| = 1.253 sigma (Rayleigh mean).
UL95_FACTOR = 1.253 + 2.0


def load_rates(npz_path):
    """Per-band total rates, residual-noise beta and signal rates per 1-px bin.

    Returns dict with tau_ps (1-px grid), total (4, n), beta (4, n),
    signal (4, n) in events/frame/bin, the noise density rho per x-px, its
    flatness chi2/ndf inside the analysis window, and the kept frame count.
    """
    data = np.load(npz_path, allow_pickle=True)
    tau_ps = np.asarray(data["time_ps"], dtype=np.float64)  # 1-px grid
    photon_sub = np.asarray(data["traces_photon"], dtype=np.float64)
    edges = np.asarray(data["subwindow_edges"], dtype=np.int64)
    meta = json.loads(str(data["meta_json"]))
    frames = int(meta["n_frames_used"])

    band_traces = np.asarray(data["band_traces_photon"], dtype=np.float64)
    band_traces = np.nan_to_num(band_traces, nan=0.0)

    # Out-of-band sub-windows: fully outside every band window.
    outside = np.ones(edges.size - 1, dtype=bool)
    for lo, hi in est.BAND_WINDOWS:
        outside &= ~((edges[:-1] >= lo) & (edges[1:] <= hi))
    outside_width = float((edges[1:] - edges[:-1])[outside].sum())
    rho = np.nan_to_num(photon_sub[outside], nan=0.0).sum(axis=0) / outside_width

    widths = est.BAND_WINDOWS[:, 1] - est.BAND_WINDOWS[:, 0]
    beta = rho[None, :] * widths[:, None]
    signal = band_traces - beta

    mask = (tau_ps >= WINDOW_PS[0]) & (tau_ps <= WINDOW_PS[1])
    # Flatness of the noise along tau, tested on the raw out-of-band event
    # counts per 1-px bin (Poisson chi2; 1.0 = consistent with flat).
    counts_out = np.nan_to_num(photon_sub[outside], nan=0.0).sum(axis=0) * frames
    expected = counts_out[mask].mean()
    chi2_ndf = float(
        np.sum((counts_out[mask] - expected) ** 2)
        / max(expected, 1e-15)
        / (mask.sum() - 1)
    )
    return {
        "tau_ps": tau_ps,
        "total": band_traces,
        "beta": beta,
        "signal": signal,
        "rho": rho,
        "rho_chi2_ndf_window": chi2_ndf,
        "frames": frames,
        "mask": mask,
    }


def rebin(values, width):
    """Sum consecutive 1-px bins into width-px bins (trailing partial dropped)."""
    n = values.shape[-1] // width * width
    return values[..., :n].reshape(*values.shape[:-1], n // width, width).sum(axis=-1)


def sem_plugin(m1, m2, beta1, beta2, frames):
    """Analytic per-bin SEM of the background-corrected g2 under the null.

    Plug-in Poisson variance (independence): var(n1 n2) = m1 m2 (1+m1+m2)
    with the total means m; derivatives of the corrected ratio as in
    fig2c_g2_estimator.g2_zero_delay, evaluated at g2 = 1.
    """
    s1 = m1 - beta1
    s2 = m2 - beta2
    denom = s1 * s2
    with np.errstate(divide="ignore", invalid="ignore"):
        var_product = m1 * m2 * (1.0 + m1 + m2) / frames
        var_m1 = m1 / frames
        var_m2 = m2 / frames
        d_product = 1.0 / denom
        d_m1 = -beta2 / denom - 1.0 / s1
        d_m2 = -beta1 / denom - 1.0 / s2
        variance = (
            d_product**2 * var_product + d_m1**2 * var_m1 + d_m2**2 * var_m2
        )
        sem = np.sqrt(variance)
    sem[(s1 <= 0.0) | (s2 <= 0.0)] = np.nan
    return sem


def analytic_pooled(m1, m2, beta1, beta2, frames):
    """Analytic pooled g2(0) SEM: coincidence-count dominated ratio estimate."""
    s1 = m1 - beta1
    s2 = m2 - beta2
    var_a = (m1 * m2 * (1.0 + m1 + m2) / frames).sum()  # numerator variance
    var_d = ((s2**2 * m1 + s1**2 * m2) / frames).sum()  # denominator variance
    d = (s1 * s2).sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        sem = np.sqrt(var_a / d**2 + var_d / d**2)  # g2 = 1 under the null
    return float(sem)


def amplitude_sigma(t_ps, sem, periods=PERIODS_PS, bin_ps=None):
    """Two-quadrature sinusoid amplitude sigma vs period on the g2 curve.

    Weighted least squares with design [1, t, t2, cos, sin] (quadratic
    detrend as in the Fig. 2(a) scan); the reported sigma is the per-
    quadrature standard deviation of the amplitude components. Bin
    averaging attenuates a sinusoid by sinc(pi dt/T); the design matrix
    uses the attenuated basis so the sensitivity degrades accordingly.
    """
    finite = np.isfinite(sem) & (sem > 0.0)
    t = t_ps[finite]
    w = 1.0 / sem[finite] ** 2
    sigma = np.full(periods.size, np.nan)
    for k, period in enumerate(periods):
        attenuation = 1.0
        if bin_ps is not None:
            attenuation = np.sinc(bin_ps / period)  # sin(pi x)/(pi x), x = dt/T
        basis = attenuation * np.exp(2j * np.pi * t / period)
        design = np.column_stack(
            [np.ones_like(t), t, t**2, basis.real, basis.imag]
        )
        xt_w = design.T * w
        try:
            cov = np.linalg.inv(xt_w @ design)
        except np.linalg.LinAlgError:
            continue
        sigma[k] = np.sqrt(max(cov[3, 3], cov[4, 4]))
    return sigma


def monte_carlo(total, beta, frames, reps, tau_m, bin_ps, seed=7):
    """Poisson cubes at the measured rates through the real estimator.

    total/beta: (4, bins) per-frame rates for the masked bins. Per-bin g2 at
    this occupancy is a discrete, heavy-tailed statistic (mostly exactly 0,
    rare huge spikes), so no per-bin Gaussian SEM is meaningful; the gate
    therefore validates the COHERENT statistics directly:

      * pooled g2(0): MC spread of the pooled estimate vs mean jackknife SEM;
      * sinusoid amplitude: per-replica two-quadrature weighted fit of the
        per-bin g2 curve (plug-in analytic weights, quadratic detrend, same
        period scan), giving the null distribution of |A| per period; the
        mean of (|A| + 2 sigma_A) over replicas is the MC counterpart of the
        analytic UL95, and q95(|A|) is the detection threshold.

    Returns (pooled_stats, amp_stats) keyed by pair.
    """
    rng = np.random.default_rng(seed)
    bins = total.shape[1]
    pooled = {pair: [] for pair in est.SAME_MESA_PAIRS}
    pooled_jk = {pair: [] for pair in est.SAME_MESA_PAIRS}
    amps = {pair: [] for pair in est.SAME_MESA_PAIRS}
    analytic_sem = {}
    mats = {}
    for pair, (a, b) in est.SAME_MESA_PAIRS.items():
        analytic_sem[pair] = sem_plugin(
            total[a], total[b], beta[a], beta[b], frames
        )
        mats[pair] = fit_matrices(tau_m, analytic_sem[pair], bin_ps)
    for _ in range(reps):
        cube = rng.poisson(total[None, :, :], size=(frames, 4, bins)).astype(np.float64)
        for pair in est.SAME_MESA_PAIRS:
            N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
            g2, _ = est.g2_zero_delay(N1, N2, b1, b2)
            p, s, _ = est.pooled_zero_delay(N1, N2, b1, b2)
            pooled[pair].append(p)
            pooled_jk[pair].append(s)
            if np.all(np.isfinite(g2)):
                amps[pair].append(fit_amplitudes_from(g2, mats[pair]))
            else:
                # Rare at 2/4 px, common at 1 px: some bins have non-finite
                # g2 (corrected signal fluctuated to <= 0); fall back to the
                # per-replica mask, mirroring what the real-data fit will do.
                amps[pair].append(
                    fit_amplitudes(tau_m, g2, analytic_sem[pair], bin_ps)
                )
    pooled_stats = {
        pair: {
            "mc_spread": float(np.std(pooled[pair], ddof=1)),
            "mean_jackknife": float(np.mean(pooled_jk[pair])),
            "mean_g2": float(np.mean(pooled[pair])),
        }
        for pair in est.SAME_MESA_PAIRS
    }
    amp_stats = {}
    for pair in est.SAME_MESA_PAIRS:
        a = np.array(amps[pair])  # (reps, periods)
        sigma_a = amplitude_sigma(tau_m, analytic_sem[pair], bin_ps=bin_ps)
        ul95_replicas = a + 2.0 * sigma_a[None, :]
        amp_stats[pair] = {
            "ul95_mc_mean": np.nanmean(ul95_replicas, axis=0),
            "amp_q95": np.nanquantile(a, 0.95, axis=0),
            "amp_median": np.nanmedian(a, axis=0),
        }
    return pooled_stats, amp_stats


def fit_amplitudes(t_ps, g2, sem, bin_ps, periods=PERIODS_PS):
    """Two-quadrature sinusoid amplitude per period on one g2 curve.

    Weighted least squares with design [1, t, t2, attenuated cos/sin]; the
    per-bin weights are the plug-in analytic SEMs (fixed, not
    sample-dependent), which keeps the fit well-behaved where the per-bin
    g2 has heavy tails. Bins with non-finite g2 or SEM are dropped.
    """
    finite = np.isfinite(g2) & np.isfinite(sem) & (sem > 0.0)
    t = t_ps[finite]
    y = g2[finite]
    w = 1.0 / sem[finite] ** 2
    out = np.full(periods.size, np.nan)
    for k, period in enumerate(periods):
        attenuation = np.sinc(bin_ps / period)
        basis = attenuation * np.exp(2j * np.pi * t / period)
        design = np.column_stack(
            [np.ones_like(t), t, t**2, basis.real, basis.imag]
        )
        xt_w = design.T * w
        try:
            cov = np.linalg.inv(xt_w @ design)
        except np.linalg.LinAlgError:
            continue
        beta_hat = cov @ (xt_w @ y)
        out[k] = np.hypot(beta_hat[3], beta_hat[4])
    return out


def fit_matrices(t_ps, sem, bin_ps, periods=PERIODS_PS):
    """Precompute (cov, XtW, finite-mask) per period for fixed plug-in weights.

    The weights depend only on the analytic SEMs, so everything except the
    final projection onto the data vector is shared across Monte Carlo
    replicas.
    """
    finite = np.isfinite(sem) & (sem > 0.0)
    t = t_ps[finite]
    w = 1.0 / sem[finite] ** 2
    mats = []
    for period in periods:
        attenuation = np.sinc(bin_ps / period)
        if abs(attenuation) < 0.05:
            # Bin averaging kills the basis (period ~ bin width / k): the
            # weighted problem is numerically degenerate and carries no
            # sensitivity; skip (NaN) instead of inverting near-singular.
            mats.append(None)
            continue
        basis = attenuation * np.exp(2j * np.pi * t / period)
        design = np.column_stack(
            [np.ones_like(t), t, t**2, basis.real, basis.imag]
        )
        xt_w = design.T * w
        try:
            cov = np.linalg.inv(xt_w @ design)
        except np.linalg.LinAlgError:
            mats.append(None)
            continue
        mats.append((cov, xt_w, finite))
    return mats


def fit_amplitudes_from(g2, mats):
    """Per-period |A| of one g2 curve given precomputed fit matrices."""
    out = np.full(len(mats), np.nan)
    for k, entry in enumerate(mats):
        if entry is None:
            continue
        cov, xt_w, finite = entry
        y = g2[finite]
        if not np.all(np.isfinite(y)):
            continue
        beta_hat = cov @ (xt_w @ y)
        out[k] = np.hypot(beta_hat[3], beta_hat[4])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--reps", type=int, default=24)
    parser.add_argument("--figures-only", action="store_true",
                        help="load the existing report JSON and only re-render "
                             "figures (no Monte Carlo)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.figures_only:
        rates = load_rates(args.npz)
        with open(out_dir / "power_gate_report.json") as handle:
            report = json.load(handle)
        make_figures(rates, report, out_dir)
        return

    rates = load_rates(args.npz)
    tau_ps = rates["tau_ps"]
    frames = rates["frames"]
    mask = rates["mask"]
    print(f"frames kept: {frames}; analysis window {WINDOW_PS} ps")
    print(f"noise density rho: window mean {rates['rho'][mask].mean():.3e} "
          f"events/frame/x-px/1-px bin, flatness chi2/ndf = {rates['rho_chi2_ndf_window']:.2f}")

    report = {
        "npz": str(args.npz),
        "frames": frames,
        "window_ps": list(WINDOW_PS),
        "rho_window_mean": float(rates["rho"][mask].mean()),
        "rho_flatness_chi2_ndf": rates["rho_chi2_ndf_window"],
        "reps": int(args.reps),
        "periods_ps": [float(v) for v in PERIODS_PS],
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "per_width": {},
    }

    # Per-band rates inside the window (1-px bins, for the report).
    for k, name in enumerate(est.BAND_NAMES):
        report.setdefault("band_rates_window", {})[name] = {
            "total_per_frame": float(rates["total"][k][mask].sum()),
            "beta_per_frame": float(rates["beta"][k][mask].sum()),
            "signal_per_frame": float(rates["signal"][k][mask].sum()),
        }
        print(f"{name}: total {rates['total'][k][mask].sum():.3f}, noise "
              f"{rates['beta'][k][mask].sum():.3f}, signal "
              f"{rates['signal'][k][mask].sum():.3f} events/frame in window")

    for width in BIN_WIDTHS:
        # Rebin, then apply the window mask on the rebinned grid.
        total_r = rebin(rates["total"], width)
        beta_r = rebin(rates["beta"], width)
        tau_r = rebin(tau_ps, width) / width  # mean tau per rebinned bin
        mask_r = (tau_r >= WINDOW_PS[0]) & (tau_r <= WINDOW_PS[1])
        total_m = total_r[:, mask_r]
        beta_m = beta_r[:, mask_r]
        tau_m = tau_r[mask_r]
        bin_ps = width * est.PS_PER_PX
        i_paper = int(np.argmin(np.abs(PERIODS_PS - 9.84)))

        width_report = {
            "bin_width_px": width,
            "bin_width_ps": float(bin_ps),
            "n_bins_window": int(mask_r.sum()),
            "pairs": {},
        }
        sigmas = {}
        for pair, (a, b) in est.SAME_MESA_PAIRS.items():
            m1, m2 = total_m[a], total_m[b]
            b1, b2 = beta_m[a], beta_m[b]
            sem = sem_plugin(m1, m2, b1, b2, frames)
            pooled_sem = analytic_pooled(m1, m2, b1, b2, frames)
            sigma_a = amplitude_sigma(tau_m, sem, bin_ps=bin_ps)
            sigmas[pair] = sigma_a
            ul95 = UL95_FACTOR * sigma_a
            width_report["pairs"][pair] = {
                "sem_median": float(np.nanmedian(sem)),
                "sem_min": float(np.nanmin(sem)),
                "sem_max": float(np.nanmax(sem)),
                "frac_bins_sem_gt_0.5": float(np.mean(sem > est.MAX_DISPLAY_SEM)),
                "pooled_sem_analytic": pooled_sem,
                "pooled_ul95_analytic": float(2.0 * pooled_sem),
                "ul95_median": float(np.nanmedian(ul95)),
                "ul95_min": float(np.nanmin(ul95)),
                "ul95_max": float(np.nanmax(ul95)),
                "ul95_at_9p84ps": float(ul95[i_paper]),
                "ul95_vs_period": [float(v) for v in ul95],
                "frames_for_ul95_0.2": float(
                    frames * (np.nanmedian(ul95) / PAPER_AMPLITUDE_RANGE[0]) ** 2
                ),
                "frames_for_ul95_1.0": float(
                    frames * (np.nanmedian(ul95) / 1.0) ** 2
                ),
            }
            print(f"\n{width} px bins ({bin_ps:.2f} ps), pair {pair}: "
                  f"median per-bin SEM {np.nanmedian(sem):.2f}, "
                  f"pooled SEM {pooled_sem:.3f}, "
                  f"median 95% min-resolvable amplitude {np.nanmedian(ul95):.2f}, "
                  f"at 9.84 ps {ul95[i_paper]:.2f}")

        # Anti-phase combined (L and R independent, opposite sign):
        # inverse-variance combination of the two amplitude sigmas.
        sigma_comb = 1.0 / np.sqrt(
            1.0 / sigmas["L"] ** 2 + 1.0 / sigmas["R"] ** 2
        )
        ul95_comb = UL95_FACTOR * sigma_comb
        width_report["combined_antiphase"] = {
            "ul95_median": float(np.nanmedian(ul95_comb)),
            "ul95_min": float(np.nanmin(ul95_comb)),
            "ul95_max": float(np.nanmax(ul95_comb)),
            "ul95_at_9p84ps": float(ul95_comb[i_paper]),
            "ul95_vs_period": [float(v) for v in ul95_comb],
            "frames_for_ul95_0.2": float(
                frames * (np.nanmedian(ul95_comb) / PAPER_AMPLITUDE_RANGE[0]) ** 2
            ),
        }
        print(f"  combined anti-phase L+R: median 95% min-resolvable "
              f"{np.nanmedian(ul95_comb):.2f}, at 9.84 ps {ul95_comb[i_paper]:.2f}")

        # Monte Carlo validation at this width (masked bins only): the
        # coherent statistics (pooled g2, sinusoid amplitude) are checked
        # directly; per-bin g2 is heavy-tailed at this occupancy and gets
        # no per-bin Gaussian validation by design.
        pooled_stats, amp_stats = monte_carlo(
            total_m, beta_m, frames, args.reps, tau_m, bin_ps
        )
        for pair, (a, b) in est.SAME_MESA_PAIRS.items():
            analytic_sem = sem_plugin(
                total_m[a], total_m[b], beta_m[a], beta_m[b], frames
            )
            sigma_a = amplitude_sigma(tau_m, analytic_sem, bin_ps=bin_ps)
            analytic_ul95 = UL95_FACTOR * sigma_a
            mc_ul95 = amp_stats[pair]["ul95_mc_mean"]
            ratio = mc_ul95 / analytic_ul95
            # The doc criterion: MC reproduces the prediction within 10%.
            validated = bool(0.9 <= np.nanmedian(ratio) <= 1.1)
            pooled_ratio = (
                pooled_stats[pair]["mean_jackknife"] / pooled_stats[pair]["mc_spread"]
            )
            width_report["pairs"][pair]["mc"] = {
                "ul95_mc_over_analytic_median": float(np.nanmedian(ratio)),
                "ul95_mc_over_analytic_max": float(np.nanmax(ratio)),
                "ul95_mc_over_analytic_vs_period": [float(v) for v in ratio],
                "mc_validated": validated,
                "amp_q95_median": float(np.nanmedian(amp_stats[pair]["amp_q95"])),
                "pooled_mc_spread": pooled_stats[pair]["mc_spread"],
                "pooled_mean_jackknife": pooled_stats[pair]["mean_jackknife"],
                "pooled_jk_over_mc": float(pooled_ratio),
                "pooled_mean_g2": pooled_stats[pair]["mean_g2"],
            }
            print(f"  MC({args.reps} reps) {pair}: UL95 MC/analytic median "
                  f"{np.nanmedian(ratio):.3f} (max {np.nanmax(ratio):.3f}) "
                  f"-> {'validated' if validated else 'NOT validated'}; "
                  f"pooled jackknife/MC {pooled_ratio:.3f}; "
                  f"pooled mean g2 {pooled_stats[pair]['mean_g2']:.3f}")
        width_report["combined_antiphase"]["mc_validated"] = bool(
            width_report["pairs"]["L"]["mc"]["mc_validated"]
            and width_report["pairs"]["R"]["mc"]["mc_validated"]
        )
        report["per_width"][str(width)] = width_report

    # Gate verdict on the best (smallest) median UL95 across MC-validated
    # widths/pairs, including the anti-phase combined sensitivity.
    candidates = [
        (w["pairs"][p]["ul95_median"], int(width), p)
        for width, w in report["per_width"].items()
        for p in est.SAME_MESA_PAIRS
        if w["pairs"][p]["mc"]["mc_validated"]
    ]
    candidates += [
        (w["combined_antiphase"]["ul95_median"], int(width), "L+R")
        for width, w in report["per_width"].items()
        if w["combined_antiphase"]["mc_validated"]
    ]
    best = min(candidates, key=lambda item: item[0])
    verdict = (
        "pass" if best[0] <= PAPER_AMPLITUDE_RANGE[0]
        else "borderline" if best[0] <= PAPER_AMPLITUDE_RANGE[1]
        else "upper-limit-only"
    )
    report["gate"] = {
        "paper_amplitude_range": list(PAPER_AMPLITUDE_RANGE),
        "best_median_ul95": best[0],
        "best_width_px": best[1],
        "best_pair": best[2],
        "verdict": verdict,
        "rule": (
            "pass if min-resolvable 95% amplitude <= 0.2 (paper lower edge); "
            "borderline up to 1.2; else upper-limit-only (docs/adiyatullin2017-fig2c.md G3)"
        ),
    }
    print(f"\nGATE VERDICT: {verdict} (best median UL95 {best[0]:.2f} "
          f"at {best[1]} px, pair {best[2]})")

    with open(out_dir / "power_gate_report.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"wrote {out_dir / 'power_gate_report.json'}")

    make_figures(rates, report, out_dir)


def make_figures(rates, report, out_dir):
    """Diagnostic figures; every quantity labelled, per repo figure rules."""
    tau_ps = rates["tau_ps"]
    mask = rates["mask"]

    # 1. Rates overview: band signals with noise shading + noise flatness.
    figure, axes = plt.subplots(2, 1, figsize=(11.0, 7.5), constrained_layout=True)
    colors = ["tab:blue", "tab:red", "tab:orange", "tab:green"]
    for k, name in enumerate(est.BAND_NAMES):
        axes[0].plot(tau_ps[mask], rates["signal"][k][mask], lw=0.8,
                     color=colors[k], label=f"{name} signal")
        axes[0].plot(tau_ps[mask], rates["beta"][k][mask], lw=0.6, ls="--",
                     color=colors[k], alpha=0.6, label=f"{name} noise beta")
    axes[0].set_ylabel("events / frame / 1-px bin")
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=8, ncol=2)
    axes[0].set_title("In-band total vs residual-noise event rates (photon traces)")
    axes[1].plot(tau_ps[mask], rates["rho"][mask], lw=0.7, color="black")
    axes[1].set_ylabel("rho [events / frame / x-px / bin]")
    axes[1].set_xlabel("time relative to t0 [ps]")
    axes[1].set_title(
        f"Residual-noise density along tau; flatness chi2/ndf = "
        f"{rates['rho_chi2_ndf_window']:.2f} in window"
    )
    figure.savefig(out_dir / "power_gate_rates.png", dpi=150)
    plt.close(figure)

    # 2. Per-bin SEM(t) medians per width and pair (bar-style hlines).
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True,
                                constrained_layout=True)
    for axis, pair in zip(axes, est.SAME_MESA_PAIRS):
        for width, color in zip(BIN_WIDTHS, ["0.3", "tab:blue", "tab:red"]):
            wr = report["per_width"][str(width)]
            axis.axhline(wr["pairs"][pair]["sem_median"], color=color, ls=":",
                         lw=1.2,
                         label=f"{width} px median {wr['pairs'][pair]['sem_median']:.2f}")
        axis.axhline(est.MAX_DISPLAY_SEM, color="gray", ls="--", lw=1.0,
                     label=f"MAX_DISPLAY_SEM = {est.MAX_DISPLAY_SEM}")
        axis.set_yscale("log")
        axis.set_title(f"pair {pair} (arms {'B1xB3' if pair == 'L' else 'B2xB4'})")
        axis.legend(fontsize=8)
        axis.set_xticks([])
    axes[0].set_ylabel("per-bin SEM(g2) median [log]")
    figure.suptitle("Analytic per-bin SEM(g2) medians vs bin width (window 0-820 ps)")
    figure.savefig(out_dir / "power_gate_sem_medians.png", dpi=150)
    plt.close(figure)

    # 3. Min-resolvable amplitude summary vs width + frames scaling.
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    widths = [int(w) for w in BIN_WIDTHS]
    for pair, color in zip(est.SAME_MESA_PAIRS, ["tab:blue", "tab:red"]):
        ul = [report["per_width"][str(w)]["pairs"][pair]["ul95_median"] for w in widths]
        axes[0].plot(widths, ul, "o-", color=color, label=f"pair {pair}")
    combined = [
        report["per_width"][str(w)]["combined_antiphase"]["ul95_median"] for w in widths
    ]
    axes[0].plot(widths, combined, "s-", color="tab:purple",
                 label="anti-phase combined L+R")
    axes[0].axhspan(*PAPER_AMPLITUDE_RANGE, color="gray", alpha=0.2,
                    label="paper amplitude range 0.2-1.2")
    axes[0].axhline(PAPER_AMPLITUDE_RANGE[0], color="gray", ls="--", lw=0.8)
    axes[0].set_xticks(widths)
    axes[0].set_xlabel("bin width [px]")
    axes[0].set_ylabel("median min-resolvable |A| (95%)")
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=8)
    axes[0].set_title("Coherent 3-200 ps sinusoid sensitivity (median over periods)")
    for pair, color in zip(est.SAME_MESA_PAIRS, ["tab:blue", "tab:red"]):
        frames_needed = [
            report["per_width"][str(w)]["pairs"][pair]["frames_for_ul95_0.2"]
            for w in widths
        ]
        axes[1].plot(widths, np.array(frames_needed) / 1000.0, "o-", color=color,
                     label=f"pair {pair}")
    frames_comb = [
        report["per_width"][str(w)]["combined_antiphase"]["frames_for_ul95_0.2"]
        for w in widths
    ]
    axes[1].plot(widths, np.array(frames_comb) / 1000.0, "s-", color="tab:purple",
                 label="anti-phase combined L+R")
    axes[1].axhline(report["frames"] / 1000.0, color="black", ls="--",
                    label=f"this dataset ({report['frames']} frames)")
    axes[1].set_xticks(widths)
    axes[1].set_xlabel("bin width [px]")
    axes[1].set_ylabel("frames needed for 95% sensitivity to |A| = 0.2 [k]")
    axes[1].set_yscale("log")
    axes[1].legend(fontsize=8)
    axes[1].set_title("Frame-count scaling to cover the paper's lower edge")
    gate = report["gate"]
    figure.suptitle(f"G3 gate verdict: {gate['verdict']} "
                    f"(best median UL95 {gate['best_median_ul95']:.2f})")
    figure.savefig(out_dir / "power_gate_sensitivity.png", dpi=150)
    plt.close(figure)

    # 4. MC validation summary.
    figure, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), constrained_layout=True)
    for pair, color in zip(est.SAME_MESA_PAIRS, ["tab:blue", "tab:red"]):
        ratios = [
            report["per_width"][str(w)]["pairs"][pair]["mc"]["ul95_mc_over_analytic_median"]
            for w in widths
        ]
        axes[0].plot(widths, ratios, "o-", color=color, label=f"pair {pair}")
    axes[0].axhline(1.0, color="black", lw=0.8)
    axes[0].axhspan(0.9, 1.1, color="gray", alpha=0.2, label="+/-10% reference")
    axes[0].set_yscale("log")
    axes[0].set_ylim(0.8, 3.0)  # 1-px ratios (~1e12) clip away by design
    axes[0].set_xticks(widths)
    axes[0].set_xlabel("bin width [px]")
    axes[0].set_ylabel("MC UL95 / analytic UL95 (median over periods)")
    axes[0].legend(fontsize=8)
    axes[0].set_title("Sinusoid-scan sensitivity: Monte Carlo vs analytic\n"
                      "(log scale; 1 px explodes = tail-dominated, NOT validated)")
    for pair, color in zip(est.SAME_MESA_PAIRS, ["tab:blue", "tab:red"]):
        ratios = [
            report["per_width"][str(w)]["pairs"][pair]["mc"]["pooled_jk_over_mc"]
            for w in widths
        ]
        axes[1].plot(widths, ratios, "o-", color=color, label=f"pair {pair}")
    axes[1].axhline(1.0, color="black", lw=0.8)
    axes[1].axhspan(0.9, 1.1, color="gray", alpha=0.2)
    axes[1].set_xticks(widths)
    axes[1].set_xlabel("bin width [px]")
    axes[1].set_ylabel("pooled jackknife SEM / MC spread")
    axes[1].legend(fontsize=8)
    axes[1].set_title("Pooled g2(0): jackknife honesty")
    figure.savefig(out_dir / "power_gate_mc_validation.png", dpi=150)
    plt.close(figure)

    # 5. Min-resolvable amplitude vs period (fig2a-style UL scan curves).
    periods = np.asarray(report["periods_ps"], dtype=np.float64)
    figure, axes = plt.subplots(1, 3, figsize=(16.5, 4.8), sharey=True,
                                constrained_layout=True)
    for axis, width in zip(axes, widths):
        wr = report["per_width"][str(width)]
        for pair, color in zip(est.SAME_MESA_PAIRS, ["tab:blue", "tab:red"]):
            axis.plot(periods, wr["pairs"][pair]["ul95_vs_period"], lw=1.2,
                      color=color, label=f"pair {pair}")
        axis.plot(periods, wr["combined_antiphase"]["ul95_vs_period"], lw=1.4,
                  color="tab:purple", label="anti-phase L+R")
        axis.axhspan(*PAPER_AMPLITUDE_RANGE, color="gray", alpha=0.2)
        axis.axvline(9.84, color="black", ls=":", lw=0.9,
                     label="9.84 ps (paper Josephson period)")
        axis.axvline(37.2, color="gray", ls=":", lw=0.9,
                     label="37.2 ps (k-space fringe period)")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xticks([3, 10, 30, 100, 200])
        axis.set_xticklabels(["3", "10", "30", "100", "200"])
        validated = {
            p: wr["pairs"][p]["mc"]["mc_validated"] for p in est.SAME_MESA_PAIRS
        }
        axis.set_title(f"{width} px bins ({wr['bin_width_ps']:.2f} ps) — MC: "
                       + ", ".join(
                           f"{p} {'ok' if v else 'NOT ok'}"
                           for p, v in validated.items()
                       ))
        axis.set_xlabel("period [ps]")
    axes[0].set_ylabel("min-resolvable g2-1 amplitude (95%)")
    axes[0].legend(fontsize=8)
    figure.suptitle("Phase-quadrature sinusoid sensitivity vs period "
                    "(quadratic detrend, bin-averaging attenuation included; "
                    "gray band = paper amplitude 0.2-1.2)")
    figure.savefig(out_dir / "power_gate_ul95_vs_period.png", dpi=150)
    plt.close(figure)
    print(f"wrote figures to {out_dir}")


if __name__ == "__main__":
    main()
