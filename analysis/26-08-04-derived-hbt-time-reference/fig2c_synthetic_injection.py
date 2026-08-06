"""G7 sensitivity validation for the 26-08-04 Fig. 2(c) g2 pipeline.

Injects a paper-style g2 modulation into synthetic count cubes and checks
recovery through the same estimator + oscillation scan used on the real
data. Injection model (per frame, per 1-px bin): independent Poisson counts
per band at the MEASURED total rates (signal + residual noise), plus
correlated same-bin pairs across each same-mesa arm pair at rate
s_a(t) s_b(t) A(t), giving E[g2(0)(t)] - 1 = A(t) exactly after the beta
correction. A(t) = env(t) * A_pk * (1 +/- sin(2 pi t / T)) / 2 with a
growing envelope env(t) = 0.75 + 0.25 t/T_window (paper: peaks grow as the
population decays) and opposite signs for the L and R pairs (anti-phase).

Two runs:
  * paper strength (A_pk ramping 1.5 -> 2, T = 9.84 ps): 16 replicas; the
    ensemble mean of the recovered amplitude tests the <10% bias criterion
    (single-shot scatter at this strength is ~30%, set by the G3 gate);
  * strong (A_pk ramping 6 -> 8): 1 replica; single-shot period recovery
    <2% and amplitude bias <10%.

Scope: synthesis is at the count-cube level. The raw-frame chain (background
subtraction, CFD50 t0, event detection) was validated end-to-end in the
Fig. 2(a) work (synthetic_recovery_cfd50: t0 bias -0.10 px, std 0.33 px;
period recovery 9.87 vs 9.84 ps injected); the residual t0 error is much
smaller than the 4 px bins used here, so the cube-level injection isolates
the correlation pipeline that is new in this task.
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
    PERIODS_PS,
    fit_matrices,
    sem_plugin,
)
from fig2c_science_tests import fit_complex  # noqa: E402

DEFAULT_OUT_DIR = Path("drafts/26-08-06_fig2c")
WINDOW_PS = (0.0, 820.0)
T_INJECT_PS = 9.84
BIN_WIDTH = 4


def measured_rates(out_dir, width):
    """Total rates m_c(t) and beta_c(t) per band per bin from the real cube."""
    results = np.load(out_dir / f"g2_results_bin{width}px.npz")
    t_ps = np.asarray(results["t_ps"], dtype=np.float64)
    mask = (t_ps >= WINDOW_PS[0]) & (t_ps <= WINDOW_PS[1])
    total = np.asarray(results["mean_counts"], dtype=np.float64)[:, mask]
    beta = np.asarray(results["beta"], dtype=np.float64)[:, mask]
    return t_ps[mask], total, beta


def inject(total, beta, t_ps, a_pk0, a_pk1, period, frames, rng):
    """Synthetic cube (frames, 4, bins) with a paper-style g2 modulation.

    Also returns the true per-pair g2-1 curves A_L(t), A_R(t) (for the
    noiseless target amplitude) and the fundamental envelope.
    """
    bins = total.shape[1]
    cube = rng.poisson(total[None, :, :], size=(frames, 4, bins)).astype(np.float64)
    ramp = (t_ps - t_ps[0]) / (t_ps[-1] - t_ps[0])
    envelope = 0.75 + 0.25 * ramp
    a_pk = a_pk0 + (a_pk1 - a_pk0) * ramp
    phase = np.sin(2.0 * np.pi * t_ps / period)
    signal = np.clip(total - beta, 0.0, None)
    true_curves = {}
    for pair, sign in (("L", +1.0), ("R", -1.0)):
        a, b = est.SAME_MESA_PAIRS[pair]
        # Fundamental-component amplitude of A(t): a_pk * envelope / 2.
        amplitude = envelope * a_pk * (1.0 + sign * phase) / 2.0
        true_curves[pair] = amplitude
        rate = signal[a] * signal[b] * amplitude  # pair rate per frame per bin
        pairs = rng.poisson(rate[None, :], size=(frames, bins))
        cube[:, a, :] += pairs
        cube[:, b, :] += pairs
    fundamental = envelope * a_pk / 2.0
    return cube, fundamental, true_curves


def target_amplitudes(true_curves, sems, t_ps, period_target=T_INJECT_PS):
    """Noiseless target: fit the true g2-1 curve (+1) with the same machinery.

    Separates estimator bias from the slow-envelope model mismatch: the
    recovered amplitude of noisy data is compared against what the same fit
    returns on the noiseless injected curve.
    """
    bin_ps = float(np.median(np.diff(t_ps)))
    targets = {}
    for pair, amplitude in true_curves.items():
        g2_true = 1.0 + amplitude
        sem = sems[pair]
        mats = fit_matrices(t_ps, sem, bin_ps)
        z = fit_complex(t_ps, g2_true, sem, bin_ps, mats)
        near = (PERIODS_PS >= period_target / 1.3) & (
            PERIODS_PS <= period_target * 1.3
        )
        k = int(np.nanargmax(np.where(near, np.abs(z), np.nan)))
        targets[pair] = {
            "period": float(PERIODS_PS[k]),
            "amplitude": float(abs(z[k])),
        }
    return targets


def analyze(cube, beta, t_ps, frames):
    """Same chain as the real data: g2 curves, plug-in SEM, complex fit."""
    bin_ps = float(np.median(np.diff(t_ps)))
    out = {}
    for pair in est.SAME_MESA_PAIRS:
        N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
        g2, _ = est.g2_zero_delay(
            N1.astype(np.float64), N2.astype(np.float64), b1, b2
        )
        sem = sem_plugin(N1.mean(axis=0), N2.mean(axis=0), b1, b2, frames)
        mats = fit_matrices(t_ps, sem, bin_ps)
        z = fit_complex(t_ps, g2, sem, bin_ps, mats)
        out[pair] = {"g2": g2, "sem": sem, "z": z}
    return out


def recover(analysis, period_target=T_INJECT_PS):
    """Recovered period (peak |z| within +-30% of target) and amplitude."""
    rec = {}
    for pair, entry in analysis.items():
        amps = np.abs(entry["z"])
        near = (PERIODS_PS >= period_target / 1.3) & (PERIODS_PS <= period_target * 1.3)
        k = int(np.nanargmax(np.where(near, amps, np.nan)))
        rec[pair] = {"period": float(PERIODS_PS[k]), "amplitude": float(amps[k])}
    return rec


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

    def run_replica(a_pk0, a_pk1):
        cube, fundamental, true_curves = inject(
            total, beta, t_ps, a_pk0, a_pk1, T_INJECT_PS, frames, rng
        )
        analysis = analyze(cube, beta, t_ps, frames)
        sems = {pair: entry["sem"] for pair, entry in analysis.items()}
        targets = target_amplitudes(true_curves, sems, t_ps)
        rec = recover(analysis)
        return rec, targets

    # --- paper-strength ensemble ---
    recs = []
    targets_all = []
    for rep in range(args.replicas):
        rec, targets = run_replica(1.5, 2.0)
        recs.append(rec)
        targets_all.append(targets)
        print(f"  replica {rep}: L {rec['L']['amplitude']:.2f} @ "
              f"{rec['L']['period']:.1f} ps (target {targets['L']['amplitude']:.2f}), "
              f"R {rec['R']['amplitude']:.2f} @ {rec['R']['period']:.1f} ps "
              f"(target {targets['R']['amplitude']:.2f})", flush=True)
    amp_l = np.array([r["L"]["amplitude"] for r in recs])
    amp_r = np.array([r["R"]["amplitude"] for r in recs])
    tgt_l = np.array([t["L"]["amplitude"] for t in targets_all])
    tgt_r = np.array([t["R"]["amplitude"] for t in targets_all])
    per_l = np.array([r["L"]["period"] for r in recs])
    per_r = np.array([r["R"]["period"] for r in recs])
    report["paper_strength"] = {
        "a_pk_ramp": [1.5, 2.0],
        "replicas": int(args.replicas),
        "target_amplitude_L": float(tgt_l.mean()),
        "target_amplitude_R": float(tgt_r.mean()),
        "recovered_L_mean": float(amp_l.mean()),
        "recovered_L_std": float(amp_l.std(ddof=1)),
        "recovered_R_mean": float(amp_r.mean()),
        "recovered_R_std": float(amp_r.std(ddof=1)),
        "period_L_unique": sorted(set(per_l.tolist())),
        "period_R_unique": sorted(set(per_r.tolist())),
        "amplitude_bias_L": float(amp_l.mean() / tgt_l.mean() - 1.0),
        "amplitude_bias_R": float(amp_r.mean() / tgt_r.mean() - 1.0),
    }
    print(f"paper strength: target amplitudes L {tgt_l.mean():.2f} / "
          f"R {tgt_r.mean():.2f}; recovered L {amp_l.mean():.2f} +/- "
          f"{amp_l.std(ddof=1):.2f}, R {amp_r.mean():.2f} +/- {amp_r.std(ddof=1):.2f}")

    # --- strong single-shot ---
    rec, targets = run_replica(6.0, 8.0)
    report["strong"] = {
        "a_pk_ramp": [6.0, 8.0],
        "target_amplitude_L": targets["L"]["amplitude"],
        "target_amplitude_R": targets["R"]["amplitude"],
        "recovered": rec,
        "period_rel_error_L": abs(rec["L"]["period"] - T_INJECT_PS) / T_INJECT_PS,
        "period_rel_error_R": abs(rec["R"]["period"] - T_INJECT_PS) / T_INJECT_PS,
        "amplitude_bias_L": rec["L"]["amplitude"] / targets["L"]["amplitude"] - 1.0,
        "amplitude_bias_R": rec["R"]["amplitude"] / targets["R"]["amplitude"] - 1.0,
    }
    print(f"strong: targets L {targets['L']['amplitude']:.2f} / "
          f"R {targets['R']['amplitude']:.2f}; recovered "
          f"L {rec['L']['amplitude']:.2f} @ {rec['L']['period']:.1f} ps, "
          f"R {rec['R']['amplitude']:.2f} @ {rec['R']['period']:.1f} ps")

    with open(args.out_dir / "fig2c_injection_report.json", "w") as handle:
        json.dump(report, handle, indent=2)
    print(f"wrote {args.out_dir / 'fig2c_injection_report.json'}")


if __name__ == "__main__":
    main()
