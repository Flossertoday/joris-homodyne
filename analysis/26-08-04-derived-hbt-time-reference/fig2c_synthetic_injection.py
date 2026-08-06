"""G7 sensitivity validation for the 26-08-04 Fig. 2(c) g2 pipeline.

Injects a paper-style g2 modulation into synthetic count cubes and checks
recovery through the same estimator + oscillation fit used on the real
data. Injection model (per frame, per bin): independent Poisson counts per
band at the MEASURED total rates (signal + residual noise), plus correlated
same-bin pairs across each same-mesa arm pair at rate s_a(t) s_b(t) A(t),
giving E[g2(0)(t)] - 1 = A(t) exactly after the beta correction. The
continuous A(t) = env(t) * A_pk * (1 +/- sin(2 pi t / T)) / 2 (growing
envelope env = 0.75 -> 1.0, opposite signs for L and R = anti-phase) is
bin-averaged by the measurement; the injection applies the sinc(dt/T)
attenuation to the oscillating part, so the sinc-correcting fit recovers
the true continuous amplitude.

Amplitude estimators: the two-quadrature magnitude |z| is Rayleigh-inflated
(E|z| ~ sqrt(A^2 + 2 sigma^2)), so bias is evaluated on the SIGNED
quadrature projection Im(z) (the injection phase is +90 deg by
construction), a linear unbiased estimator; the ensemble mean over replicas
carries SEM sigma/sqrt(N_rep).

Two runs:
  * strong (A_pk ramping 6 -> 8, T = 9.84 ps): 16 replicas; criteria
    (docs/adiyatullin2017-fig2c.md G7): recovered period within 2% of the
    injected value (per replica), ensemble amplitude bias < 10%;
  * paper strength (A_pk 1.5 -> 2): 16 replicas; ensemble recovery vs
    target and per-replica significance |z|/sigma_A vs the two-quadrature
    null q95 = 2.45 — documents the response at the paper's amplitudes
    (the G3 gate says this sits at our sensitivity edge).

Scope: synthesis is at the count-cube level. The raw-frame chain
(background subtraction, CFD50 t0, event detection) was validated
end-to-end in the Fig. 2(a) work (synthetic_recovery_cfd50: t0 bias
-0.10 px, std 0.33 px; period recovery 9.87 vs 9.84 ps injected); the
residual t0 error is much smaller than the 4 px bins used here.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fig2c_g2_estimator as est  # noqa: E402
from fig2c_power_gate import (  # noqa: E402
    fit_matrices,
    sem_plugin,
)
from fig2c_science_tests import fit_complex  # noqa: E402

DEFAULT_OUT_DIR = Path("drafts/26-08-06_fig2c")
WINDOW_PS = (0.0, 820.0)
T_INJECT_PS = 9.84
BIN_WIDTH = 4
RAYLEIGH_Q95 = 2.45  # two-quadrature null q95 of |z|/sigma


def measured_rates(out_dir, width):
    """Total rates m_c(t) and beta_c(t) per band per bin from the real cube."""
    results = np.load(out_dir / f"g2_results_bin{width}px.npz")
    t_ps = np.asarray(results["t_ps"], dtype=np.float64)
    mask = (t_ps >= WINDOW_PS[0]) & (t_ps <= WINDOW_PS[1])
    total = np.asarray(results["mean_counts"], dtype=np.float64)[:, mask]
    beta = np.asarray(results["beta"], dtype=np.float64)[:, mask]
    return t_ps[mask], total, beta


def inject(total, beta, t_ps, a_pk0, a_pk1, period, frames, rng):
    """Synthetic cube (frames, 4, bins) with the paper-style modulation.

    Returns the cube and the per-pair BINNED true curves (what the
    measurement should see), which the noiseless target fit uses.
    """
    bins = total.shape[1]
    bin_ps = float(np.median(np.diff(t_ps)))
    attenuation = np.sinc(bin_ps / period)
    cube = rng.poisson(total[None, :, :], size=(frames, 4, bins)).astype(np.float64)
    ramp = (t_ps - t_ps[0]) / (t_ps[-1] - t_ps[0])
    envelope = 0.75 + 0.25 * ramp
    a_pk = a_pk0 + (a_pk1 - a_pk0) * ramp
    phase = np.sin(2.0 * np.pi * t_ps / period)
    signal = np.clip(total - beta, 0.0, None)
    binned_curves = {}
    for pair, sign in (("L", +1.0), ("R", -1.0)):
        a, b = est.SAME_MESA_PAIRS[pair]
        binned = envelope * a_pk * (1.0 + sign * attenuation * phase) / 2.0
        binned_curves[pair] = binned
        rate = signal[a] * signal[b] * binned  # pair rate per frame per bin
        pairs = rng.poisson(rate[None, :], size=(frames, bins))
        cube[:, a, :] += pairs
        cube[:, b, :] += pairs
    return cube, binned_curves


def fit_at_period(g2, sem, t_ps, period):
    """Complex amplitude and per-quadrature sigma at one exact period."""
    bin_ps = float(np.median(np.diff(t_ps)))
    mats = fit_matrices(t_ps, sem, bin_ps, np.array([period]))
    cov = mats[0][0]
    sigma_a = float(np.sqrt(max(cov[3, 3], cov[4, 4])))
    z = fit_complex(t_ps, g2, sem, bin_ps, mats, np.array([period]))[0]
    return z, sigma_a


def scan_period(g2, sem, t_ps, period_target=T_INJECT_PS):
    """Recovered period: argmax |A| on a 0.02 ps fine scan around the target."""
    bin_ps = float(np.median(np.diff(t_ps)))
    fine = np.arange(period_target / 1.2, period_target * 1.2, 0.02)
    mats = fit_matrices(t_ps, sem, bin_ps, fine)
    z = fit_complex(t_ps, g2, sem, bin_ps, mats, fine)
    return float(fine[int(np.nanargmax(np.abs(z)))])


def analyze(cube, beta, t_ps, frames, period_target=T_INJECT_PS):
    """Same chain as the real data: g2, plug-in SEM, exact-period fit, scan."""
    out = {}
    for pair in est.SAME_MESA_PAIRS:
        N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
        N1 = N1.astype(np.float64)
        N2 = N2.astype(np.float64)
        g2, _ = est.g2_zero_delay(N1, N2, b1, b2)
        sem = sem_plugin(N1.mean(axis=0), N2.mean(axis=0), b1, b2, frames)
        z, sigma_a = fit_at_period(g2, sem, t_ps, period_target)
        out[pair] = {
            "z": z,
            "sigma_a": sigma_a,
            "period": scan_period(g2, sem, t_ps, period_target),
        }
    return out


def targets_from(binned_curves, sems, t_ps, period_target=T_INJECT_PS):
    """Noiseless target: same exact-period fit on the true binned curves."""
    targets = {}
    for pair, binned in binned_curves.items():
        z, _ = fit_at_period(1.0 + binned, sems[pair], t_ps, period_target)
        targets[pair] = float(z.imag)  # injection phase is +90 deg
    return targets


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--replicas", type=int, default=16)
    parser.add_argument("--seed", type=int, default=23)
    args = parser.parse_args()

    t_ps, total, beta = measured_rates(args.out_dir, BIN_WIDTH)
    with open(args.out_dir / "g2_summary.json") as handle:
        frames = json.load(handle)["n_frames_kept"]
    print(f"rates from the real {BIN_WIDTH} px cube: {frames} frames, "
          f"{t_ps.size} bins in window")

    rng = np.random.default_rng(args.seed)
    report = {"bin_width_px": BIN_WIDTH, "t_inject_ps": T_INJECT_PS,
              "frames": frames, "created_utc":
              datetime.now(timezone.utc).isoformat(timespec="seconds")}

    def run_ensemble(a_pk0, a_pk1, label):
        recs = {pair: [] for pair in est.SAME_MESA_PAIRS}
        targets = {pair: [] for pair in est.SAME_MESA_PAIRS}
        for rep in range(args.replicas):
            cube, binned_curves = inject(total, beta, t_ps, a_pk0, a_pk1,
                                         T_INJECT_PS, frames, rng)
            sems = {}
            for pair in est.SAME_MESA_PAIRS:
                N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
                sems[pair] = sem_plugin(
                    N1.mean(axis=0), N2.mean(axis=0), b1, b2, frames
                )
            tgt = targets_from(binned_curves, sems, t_ps)
            analysis = analyze(cube, beta, t_ps, frames)
            for pair in est.SAME_MESA_PAIRS:
                entry = analysis[pair]
                recs[pair].append({
                    "z": complex(entry["z"]),
                    "sigma_a": entry["sigma_a"],
                    "period": entry["period"],
                })
                targets[pair].append(tgt[pair])
            print(f"  {label} replica {rep}: " + ", ".join(
                f"{pair} Im(z)={analysis[pair]['z'].imag:.2f} "
                f"(target {tgt[pair]:.2f}), |z|/σ="
                f"{abs(analysis[pair]['z']) / analysis[pair]['sigma_a']:.2f}, "
                f"period {analysis[pair]['period']:.2f} ps"
                for pair in est.SAME_MESA_PAIRS
            ), flush=True)
        summary = {}
        for pair in est.SAME_MESA_PAIRS:
            z = np.array([r["z"] for r in recs[pair]])
            sigma_a = float(np.mean([r["sigma_a"] for r in recs[pair]]))
            periods = np.array([r["period"] for r in recs[pair]])
            tgt = float(np.mean(targets[pair]))
            mean_signed = float(z.imag.mean())
            ensemble_sem = float(z.imag.std(ddof=1) / np.sqrt(z.size))
            significance = np.abs(z) / sigma_a
            summary[pair] = {
                "target_amplitude": tgt,
                "ensemble_mean_signed_amplitude": mean_signed,
                "ensemble_sem": ensemble_sem,
                "amplitude_bias": mean_signed / tgt - 1.0,
                "single_shot_sigma_a": sigma_a,
                "period_within_2pct": int(
                    np.sum(np.abs(periods - T_INJECT_PS) / T_INJECT_PS < 0.02)
                ),
                "period_range": [float(periods.min()), float(periods.max())],
                "frac_replicas_significant": float(
                    np.mean(significance > RAYLEIGH_Q95)
                ),
                "mean_abs_z_over_sigma": float(significance.mean()),
            }
        return summary

    print("strong injection (A_pk 6 -> 8):", flush=True)
    report["strong"] = {
        "a_pk_ramp": [6.0, 8.0],
        "replicas": int(args.replicas),
        "pairs": run_ensemble(6.0, 8.0, "strong"),
    }
    print("paper-strength injection (A_pk 1.5 -> 2):", flush=True)
    report["paper_strength"] = {
        "a_pk_ramp": [1.5, 2.0],
        "replicas": int(args.replicas),
        "pairs": run_ensemble(1.5, 2.0, "paper"),
    }

    with open(args.out_dir / "fig2c_injection_report.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"wrote {args.out_dir / 'fig2c_injection_report.json'}")


if __name__ == "__main__":
    main()
