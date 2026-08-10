"""Render Figure: injection recovery — real-shape synthetic validation + g2 summary.

Output: analysis/26-08-10-project-report/figures/fig_injection_recovery.png
"""

import json
import sys
import time
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt

WORKTREE_ROOT = Path(__file__).resolve().parent.parent.parent
ANALYSIS_DIR = WORKTREE_ROOT / "analysis" / "26-08-04-derived-hbt-time-reference"
sys.path.insert(0, str(ANALYSIS_DIR))
import extract_hbt_traces as eht  # noqa: E402

DATA_H5 = Path(
    "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/"
    "26-08-04_hbt_time_reference.h5"
)
INJECTION_JSON = Path(
    "E:/Documents/Codex/joris-homodyne/drafts/26-08-06_fig2c/fig2c_injection_report.json"
)
OUT_PATH = Path(__file__).resolve().parent / "figures" / "fig_injection_recovery.png"

N_TEMPLATE_FRAMES = 400
N_FRAMES = 1500
PERIOD_PX = 10.0
AMPL_INJECTED = 0.5
MEAN_PHOTONS_PER_BAND = 60
BACKGROUND = 105.0
SEED = 42

NY, NX = 1016, 1344
BANDS = [(380, 450), (570, 680), (810, 880), (1010, 1110)]
LR_OF_BAND = [0, 1, 0, 1]
ROI_R0, ROI_R1, ROI_C0, ROI_C1 = eht.T0_ROI
T0_NOMINAL = 160.0
PS_PER_PX = eht.PS_PER_PX

rng = np.random.default_rng(SEED)
rows_grid = np.arange(NY)


def load_reference_patches():
    with h5py.File(DATA_H5, "r") as handle:
        background = np.asarray(handle[eht.BACKGROUND_DATASET], dtype=np.float64)
        frames = handle[eht.FRAMES_DATASET]
        picks = np.linspace(0, frames.shape[0] - 1, N_TEMPLATE_FRAMES).astype(int)
        patches, cfd50 = [], []
        for index in picks:
            frame_bg = np.asarray(frames[int(index)], dtype=np.float64) - background
            t50 = eht.estimate_t0_cfd(frame_bg, roi=eht.T0_ROI)
            if not np.isfinite(t50):
                continue
            patch = frame_bg[ROI_R0:ROI_R1, ROI_C0:ROI_C1].copy()
            patches.append(patch)
            cfd50.append(t50)
    return patches, np.asarray(cfd50)


def deposit_patch(frame, patch, patch_cfd50, target_cfd50):
    shift = target_cfd50 - patch_cfd50
    n_rows = ROI_R1 - ROI_R0
    src = np.arange(n_rows, dtype=np.float64) - shift
    out = np.empty_like(patch)
    for col in range(patch.shape[1]):
        out[:, col] = np.interp(
            src, np.arange(n_rows, dtype=np.float64), patch[:, col], left=0.0, right=0.0
        )
    frame[ROI_R0:ROI_R1, ROI_C0:ROI_C1] += out
    return frame


def envelope(t0_true):
    env = np.zeros(NY)
    start = int(t0_true + 20)
    env[start:950] = 1.0
    return env


def make_frame(patches, patch_cfd50):
    jitter = rng.normal(0, 8) + (rng.choice([0.0, -25.0], p=[0.7, 0.3]))
    t0_true = T0_NOMINAL + jitter
    pick = int(rng.integers(0, len(patches)))
    lam = np.full((NY, NX), BACKGROUND)
    deposit_patch(lam, patches[pick] * rng.uniform(0.7, 1.3), patch_cfd50[pick], t0_true)
    np.clip(lam, 0.0, None, out=lam)
    frame = rng.poisson(lam).astype(np.float64)
    total_signal = 0
    for bi, (x0, x1) in enumerate(BANDS):
        tau = rows_grid - t0_true
        sign = 1.0 if LR_OF_BAND[bi] == 0 else -1.0
        profile = envelope(t0_true) * (
            1 + sign * AMPL_INJECTED * np.cos(2 * np.pi * tau / PERIOD_PX)
        )
        n_ph = rng.poisson(MEAN_PHOTONS_PER_BAND)
        total_signal += int(n_ph)
        if n_ph == 0:
            continue
        probs = profile / profile.sum()
        ph_rows = rng.choice(rows_grid, size=n_ph, p=probs)
        ph_cols = rng.integers(x0 + 3, x1 - 3, size=n_ph)
        for r, c in zip(ph_rows, ph_cols):
            counts = rng.uniform(150, 800)
            frame[r, c] += counts * 0.4
            frame[r, c + 1] += counts * 0.25
            frame[r + 1, c] += counts * 0.2
            frame[r + 1, c + 1] += counts * 0.15
    return frame, t0_true, total_signal


def run_synthetic_validation():
    patches, patch_cfd50 = load_reference_patches()
    tau_grid = eht.TAU_GRID
    edges = eht.SUBWINDOW_EDGES
    acc_photon = np.zeros((len(edges) - 1, len(tau_grid)))
    acc_count = np.zeros(len(tau_grid))
    err_cfd50, err_centroid = [], []
    background = np.full((NY, NX), BACKGROUND)
    total_signal = 0

    t_start = time.time()
    n_fail = 0
    for i in range(N_FRAMES):
        frame, t0_true, sig = make_frame(patches, patch_cfd50)
        total_signal += sig
        frame_bg = frame - background
        t0 = eht.estimate_t0_cfd(frame_bg, roi=eht.T0_ROI)
        t0c = eht.estimate_t0(frame_bg, roi=eht.T0_ROI)
        if np.isfinite(t0):
            err_cfd50.append(t0 - t0_true)
        if np.isfinite(t0c):
            err_centroid.append(t0c - t0_true)
        if not np.isfinite(t0):
            n_fail += 1
            continue
        rows_ev, cols_ev = eht.detect_events(
            frame_bg, threshold=25, min_size=2, col_range=(350, 1150)
        )
        if len(rows_ev):
            tau_ev = rows_ev - t0
            sw = np.clip(np.searchsorted(edges, cols_ev) - 1, 0, len(edges) - 2)
            bin_idx = np.floor(tau_ev - tau_grid[0] + 0.5).astype(np.int64)
            ok = (bin_idx >= 0) & (bin_idx < len(tau_grid))
            np.add.at(acc_photon, (sw[ok], bin_idx[ok]), 1)
        valid = (tau_grid + t0 >= 0) & (tau_grid + t0 <= NY - 1)
        acc_count[valid] += 1
        if (i + 1) % 300 == 0:
            print(f"synthetic {i + 1}/{N_FRAMES} frames, {time.time() - t_start:.0f}s", flush=True)

    err_cfd50 = np.asarray(err_cfd50)
    err_centroid = np.asarray(err_centroid)
    total_detected = int(acc_photon.sum())

    t_ps = tau_grid * PS_PER_PX

    def band_sw(x0, x1):
        return [i for i in range(len(edges) - 1) if edges[i] >= x0 and edges[i + 1] <= x1]

    L = sum(acc_photon[band_sw(*BANDS[i])].sum(axis=0) for i in [0, 2])
    R = sum(acc_photon[band_sw(*BANDS[i])].sum(axis=0) for i in [1, 3])
    mean_L = L / np.maximum(acc_count, 1)
    mean_R = R / np.maximum(acc_count, 1)
    m = (t_ps > 40) & (t_ps < 800) & (acc_count > 100)
    z = (mean_L - mean_R) / (mean_L + mean_R)

    # FFT peak period
    tt = t_ps[m]
    zz = z[m] - z[m].mean()
    freq = np.fft.rfftfreq(len(zz), d=PS_PER_PX)
    spec = np.abs(np.fft.rfft(zz * np.hanning(len(zz))))
    peak = np.argmax(spec[1:]) + 1
    fft_period = 1.0 / freq[peak]
    fft_amp = 2 * spec[peak] / len(zz)

    # Least-squares sine+cosine fit at the TRUE injected period.
    true_period = PERIOD_PX * PS_PER_PX
    omega = 2 * np.pi / true_period
    X = np.column_stack([np.sin(omega * tt), np.cos(omega * tt)])
    coeffs, *_ = np.linalg.lstsq(X, zz, rcond=None)
    a_sin, a_cos = coeffs
    fitted_amp = float(np.hypot(a_sin, a_cos))
    resid = zz - X @ coeffs
    n = zz.size
    sigma = float(np.sqrt(np.sum(resid**2) / max(n - 2, 1)))
    XtX_inv = np.linalg.inv(X.T @ X)
    var_sin, var_cos = sigma**2 * np.diag(XtX_inv)
    cov_sc = sigma**2 * XtX_inv[0, 1]
    amp_sigma = float(np.sqrt(
        (a_sin**2 * var_sin + a_cos**2 * var_cos + 2 * a_sin * a_cos * cov_sc)
        / max(fitted_amp**2, 1e-12)
    ))

    # Expected (diluted) amplitude from the simulation's own signal/noise
    # counts in the same fit window. We estimate the signal oscillation
    # amplitude in each arm by a sine+cosine fit at the true period, and the
    # noise floor as the mean total count minus the signal offset. The predicted
    # z amplitude is A_predicted = A_injected * S0 / (S0 + N0), equivalently
    # (amp_L + amp_R) / (mean_L_total + mean_R_total).
    n_bins = int(m.sum())
    tt = t_ps[m]
    omega = 2 * np.pi / true_period
    X = np.column_stack([np.sin(omega * tt), np.cos(omega * tt)])

    def detrend_and_amp(trace):
        t = np.arange(n_bins)
        poly = np.polyfit(t, trace, 1)
        det = trace - np.polyval(poly, t)
        c, *_ = np.linalg.lstsq(X, det, rcond=None)
        return float(np.hypot(c[0], c[1]))

    amp_L = detrend_and_amp(mean_L[m])
    amp_R = detrend_and_amp(mean_R[m])
    mean_total = float(mean_L[m].mean() + mean_R[m].mean())
    expected_amplitude = float(
        (amp_L + amp_R) / max(mean_total, 1e-12)
    )
    dilution = expected_amplitude / AMPL_INJECTED

    return {
        "err_cfd50": err_cfd50,
        "err_centroid": err_centroid,
        "n_fail": n_fail,
        "t_ps": t_ps,
        "z": z,
        "m": m,
        "fft_period": fft_period,
        "fft_amp": fft_amp,
        "expected_amplitude": expected_amplitude,
        "dilution": dilution,
        "total_signal": total_signal,
        "total_detected": total_detected,
        "fitted_amplitude": fitted_amp,
        "fitted_amplitude_sigma": amp_sigma,
        "fitted_sin_coeff": a_sin,
        "fitted_cos_coeff": a_cos,
        "true_period": true_period,
    }


def render_g2_summary(axis):
    with open(INJECTION_JSON, "r", encoding="utf-8") as handle:
        data = json.load(handle)

    groups = [
        ("strong\n(peak g2−1 = 6→8)", "strong", [6.0, 8.0]),
        ("paper-strength\n(1.5→2)", "paper_strength", [1.5, 2.0]),
    ]
    colors = {"L": "tab:blue", "R": "tab:red"}
    offsets = {"L": -0.2, "R": 0.2}

    xticks, xlabels = [], []
    x_base = 0.25
    for label, key, ramp in groups:
        g = data[key]
        pairs = g["pairs"]
        n_replicas = g["replicas"]
        frac = {
            side: f"{int(round(pairs[side]['frac_replicas_significant'] * n_replicas))}/{n_replicas} "
                  f"({100*pairs[side]['frac_replicas_significant']:.1f}%)"
            for side in ("L", "R")
        }

        for side in ("L", "R"):
            p = pairs[side]
            target = abs(p["target_amplitude"])
            rec = abs(p["ensemble_mean_signed_amplitude"])
            sem = p["ensemble_sem"]
            bias = p["amplitude_bias"] * 100
            x = x_base + offsets[side]
            axis.errorbar(
                x, rec, yerr=sem,
                fmt="o", color=colors[side], markersize=9,
                capsize=4, capthick=1.5, elinewidth=1.5, zorder=4,
                label=f"{side} recovered" if x_base == 0.25 else "",
            )
            axis.plot([x - 0.12, x + 0.12], [target, target],
                      color=colors[side], lw=1.2, linestyle="--", zorder=2)
            # Annotation: place below points for strong L to stay inside axes,
            # above for others.
            txt = f"{side}: bias {bias:+.1f}%, det {frac[side]}"
            if side == "L" and x_base == 0.25:
                y_txt = rec - sem - 0.55
                va = "top"
            else:
                y_txt = rec + sem + 0.35
                va = "bottom"
            axis.text(x, y_txt, txt, ha="center", va=va,
                      fontsize=8, color=colors[side])

        xticks.append(x_base)
        xlabels.append(label)
        x_base += 1.5

    axis.set_xticks(xticks)
    axis.set_xticklabels(xlabels)
    axis.set_xlim(-0.1, 2.1)
    axis.set_ylabel("|amplitude|")
    axis.set_title("g^(2) injection recovery (bin 4 px, anti-phase)")
    axis.set_ylim(0, 8)

    # legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="gray",
               markersize=8, label="ensemble recovered |A| ± SEM"),
        Line2D([0], [0], color="gray", lw=1.2, linestyle="--", label="dashed = target |A|"),
        Line2D([0], [0], color="tab:blue", lw=2, label="L"),
        Line2D([0], [0], color="tab:red", lw=2, label="R"),
    ]
    axis.legend(handles=legend_elements, loc="upper right", fontsize=8)


def render(synth, out_path=OUT_PATH):
    fig = plt.figure(figsize=(16, 9.33), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.25, 1.0])
    ax_z = fig.add_subplot(gs[0, :])
    ax_t0 = fig.add_subplot(gs[1, 0])
    ax_g2 = fig.add_subplot(gs[1, 1])

    t_ps = synth["t_ps"]
    z = synth["z"]
    m = synth["m"]
    expected_amp = synth["expected_amplitude"]
    fitted_amp = synth["fitted_amplitude"]
    fitted_sigma = synth["fitted_amplitude_sigma"]
    period = synth["fft_period"]
    true_period = synth["true_period"]

    # ---- Panel (i): recovered z(t) ----
    ax_z.plot(t_ps[m], z[m], ".", ms=2, color="tab:purple", alpha=0.6, rasterized=True,
              label="recovered z(t)")
    sine_overlay = expected_amp * np.sin(2 * np.pi * t_ps[m] / true_period)
    ax_z.plot(t_ps[m], sine_overlay, "-", color="black", lw=1.4,
              label=f"expected diluted sine (A={expected_amp:.3f}, true={AMPL_INJECTED:.2f})")
    a_sin = synth["fitted_sin_coeff"]
    a_cos = synth["fitted_cos_coeff"]
    fit_overlay = a_sin * np.sin(2 * np.pi * t_ps[m] / true_period) + a_cos * np.cos(2 * np.pi * t_ps[m] / true_period)
    ax_z.plot(t_ps[m], fit_overlay, "--", color="tab:cyan", lw=1.4,
              label=f"LSQ sine+cos fit at true period (A={fitted_amp:.3f})")
    ax_z.set_xlim(t_ps[m].min(), t_ps[m].max())
    ax_z.set_xlabel("time after reference pulse [ps]")
    ax_z.set_ylabel("z = (L − R) / (L + R)")
    ax_z.set_title(
        f"real-shape synthetic validation: recovered z(t); "
        f"fitted A={fitted_amp:.3f} ± {fitted_sigma:.3f}, expected {expected_amp:.3f} "
        f"(dilution S/(S+N)); FFT peak period = {period:.2f} ps (true {true_period:.2f} ps)"
    )
    ax_z.legend(loc="upper right", fontsize=8)

    # ---- Panel (ii): t0 estimator errors ----
    err_cfd = synth["err_cfd50"]
    err_cen = synth["err_centroid"]
    bins = np.linspace(-3, 3, 61)
    ax_t0.hist(
        err_cfd, bins=bins, histtype="step", lw=1.6, color="tab:green",
        label=f"CFD50 − true  bias={err_cfd.mean():+.2f}  std={err_cfd.std():.2f} px"
    )
    ax_t0.hist(
        err_cen, bins=bins, histtype="step", lw=1.6, color="tab:orange",
        label=f"centroid − true  bias={err_cen.mean():+.2f}  std={err_cen.std():.2f} px"
    )
    ax_t0.axvline(0, color="black", lw=0.8, linestyle="--")
    ax_t0.set_xlabel("t0 estimator error [px]")
    ax_t0.set_ylabel("frames")
    ax_t0.set_title("per-frame time-origin estimator errors")
    ax_t0.legend(loc="upper right", fontsize=8)

    # ---- Panel (iii): g2 injection summary ----
    render_g2_summary(ax_g2)

    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    synth = run_synthetic_validation()
    print(
        f"CFD50    bias={synth['err_cfd50'].mean():+.3f} px  "
        f"std={synth['err_cfd50'].std():.3f} px  n={synth['err_cfd50'].size}"
    )
    print(
        f"centroid bias={synth['err_centroid'].mean():+.3f} px  "
        f"std={synth['err_centroid'].std():.3f} px  n={synth['err_centroid'].size}"
    )
    print(f"FFT peak period {synth['fft_period']:.2f} ps, z amp {synth['fft_amp']:.3f}")
    print(
        f"dilution: signal={synth['total_signal']}, detected={synth['total_detected']}, "
        f"expected A={synth['expected_amplitude']:.3f}"
    )
    print(
        f"LSQ fit at true period: A={synth['fitted_amplitude']:.3f} ± "
        f"{synth['fitted_amplitude_sigma']:.3f}"
    )
    discrepancy = abs(synth['fitted_amplitude'] - synth['expected_amplitude'])
    if discrepancy > 2 * synth['fitted_amplitude_sigma']:
        raise RuntimeError(
            f"Fitted amplitude {synth['fitted_amplitude']:.3f} ± "
            f"{synth['fitted_amplitude_sigma']:.3f} disagrees with expected "
            f"{synth['expected_amplitude']:.3f} by {discrepancy:.3f} (>2σ). "
            "STOPPING as requested."
        )
    render(synth)


if __name__ == "__main__":
    main()
