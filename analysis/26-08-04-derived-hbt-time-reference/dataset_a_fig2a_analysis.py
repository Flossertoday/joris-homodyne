"""Fig. 2(a)-style z analyses on the dataset-A traces npz.

Ports the ad-hoc tmp/ scripts of the original fig2a work
(fig2a_z_crossarm_cfd50.py, fig2a_upper_limit_cfd50.py) into one
parameterized script, and adds the noise-dilution correction from the
dataset-A noise model: the measured z amplitude is the true population
imbalance amplitude times S/(S+N) per arm, so the floor-corrected upper
limit is UL_true = UL_measured / dilution.

Outputs a JSON summary + one figure (cross-arm panel and UL scan).
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from hbt_z_targeted_fits import PAPER_PERIOD_PS, targeted_fit  # noqa: E402
from kspace_fringe_priors import measured_periods_ps  # noqa: E402

DEFAULT_NPZ = Path("drafts/26-08-08_dataset_a/hbt_traces_dataset_a.npz")
DEFAULT_TABLE = (
    Path(__file__).resolve().parent
    / "detected-peaks/26-08-04_hbt_time_reference_dataset_a.h5"
)
WINDOW_PS = (20.0, 800.0)
SCAN_PERIODS_PS = np.concatenate(
    [np.arange(3.0, 20.0, 0.1), np.arange(20.0, 200.0, 0.5)]
)


def z_arm(tr_i, se_i, tr_j, se_j):
    z = (tr_i - tr_j) / (tr_i + tr_j)
    sz = (
        2 * np.sqrt((tr_j * se_i) ** 2 + (tr_i * se_j) ** 2)
        / (tr_i + tr_j) ** 2
    )
    return z, sz


def combined(t, bands, sems, mask):
    z1, sz1 = z_arm(bands[0], sems[0], bands[1], sems[1])
    z2, sz2 = z_arm(bands[2], sems[2], bands[3], sems[3])
    w = 1 / sz1[mask] ** 2 + 1 / sz2[mask] ** 2
    z = (z1[mask] / sz1[mask] ** 2 + z2[mask] / sz2[mask] ** 2) / w
    sz = np.sqrt(1 / w)
    zr = z - np.polyval(np.polyfit(t[mask], z, 2, w=1 / sz), t[mask])
    return t[mask], zr, sz, (z1, sz1), (z2, sz2)


def crossarm_r(tt, arm1, arm2, mask1, rng, n_null=200):
    """Zero-lag Pearson r of the two arms' detrended z, with a Gaussian
    null built on the per-bin SEMs."""
    (z1, sz1), (z2, sz2) = arm1, arm2
    d1 = z1[mask1] - np.polyval(np.polyfit(tt, z1[mask1], 2, w=1 / sz1[mask1]), tt)
    d2 = z2[mask1] - np.polyval(np.polyfit(tt, z2[mask1], 2, w=1 / sz2[mask1]), tt)
    r = float(np.corrcoef(d1, d2)[0, 1])
    null = []
    for _ in range(n_null):
        n1 = rng.normal(0, sz1[mask1])
        n2 = rng.normal(0, sz2[mask1])
        null.append(np.corrcoef(n1, n2)[0, 1])
    null = np.asarray(null)
    p = float((np.abs(null) >= abs(r)).mean())
    return r, p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE,
                        help="dataset-A h5 (for the noise floors)")
    parser.add_argument("--out-dir", type=Path,
                        default=Path("drafts/26-08-08_dataset_a"))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    data = np.load(args.npz)
    t = np.asarray(data["time_ps"], dtype=np.float64)
    bands = np.asarray(data["band_traces_photon"], dtype=np.float64)
    sems = np.asarray(data["band_sem_photon"], dtype=np.float64)
    mask = (t > WINDOW_PS[0]) & (t < WINDOW_PS[1])

    with h5py.File(args.table, "r") as handle:
        floor = np.asarray(handle["noise_model/floor_per_band"])
        n_kept = int(handle.attrs["n_frames_kept"])
        band_windows = np.asarray(handle.attrs["band_windows"], dtype=np.float64)

    # dilution S/(S+N) per arm from the noise model and the measured traces
    widths = band_windows[:, 1] - band_windows[:, 0] + 1.0
    noise_rate = floor * widths  # events/frame/bin(1px)
    signal_rate = np.array(
        [np.nanmean(bands[k][mask]) - noise_rate[k] for k in range(4)]
    )
    dilution_arm = {
        "arm1 (B1,B2)": float(
            (signal_rate[0] + signal_rate[1])
            / (signal_rate[0] + signal_rate[1]
               + noise_rate[0] + noise_rate[1])
        ),
        "arm2 (B3,B4)": float(
            (signal_rate[2] + signal_rate[3])
            / (signal_rate[2] + signal_rate[3]
               + noise_rate[2] + noise_rate[3])
        ),
    }
    dilution = float(np.mean(list(dilution_arm.values())))

    tt, zr, sz, arm1, arm2 = combined(t, bands, sems, mask)
    rng = np.random.default_rng(0)
    r_cross, p_cross = crossarm_r(tt, arm1, arm2, np.ones(tt.size, bool), rng)

    # sinusoid UL scan on the combined detrended z
    amps = np.zeros(SCAN_PERIODS_PS.size)
    uls = np.zeros(SCAN_PERIODS_PS.size)
    for k, period in enumerate(SCAN_PERIODS_PS):
        result = targeted_fit(tt, zr, sz, period)
        amps[k] = result["amplitude"]
        uls[k] = result["ul95"]
    priors = {}
    for period in (PAPER_PERIOD_PS, *measured_periods_ps()):
        priors[f"{period:.2f}"] = targeted_fit(tt, zr, sz, float(period))

    summary = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "npz": str(args.npz),
        "window_ps": list(WINDOW_PS),
        "n_frames_kept": n_kept,
        "noise_rate_per_bin_per_band": noise_rate.tolist(),
        "signal_rate_per_bin_per_band": signal_rate.tolist(),
        "dilution_per_arm": dilution_arm,
        "crossarm_r": r_cross,
        "crossarm_p": p_cross,
        "sinusoid_ul95_median": float(np.median(uls)),
        "sinusoid_ul95_worst": float(uls.max()),
        "sinusoid_ul95_median_dilution_corrected": float(np.median(uls) / dilution),
        "targeted_priors": priors,
        "z_sem_median": float(np.median(sz)),
    }

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    axes[0].semilogy(SCAN_PERIODS_PS, amps, ".", ms=3, label="fitted |A|")
    axes[0].semilogy(SCAN_PERIODS_PS, uls, "-", lw=1.2,
                     label="|A| + 2σ (≈95% UL)")
    axes[0].set_xlabel("period [ps]")
    axes[0].set_ylabel("z oscillation amplitude")
    axes[0].set_title("dataset A: z(t) sinusoid amplitude vs period")
    axes[0].legend(fontsize=8)
    axes[1].plot(tt, zr, ".", ms=2, color="magenta")
    axes[1].fill_between(tt, -sz, sz, alpha=0.3, color="gray")
    axes[1].set_xlabel("t [ps]")
    axes[1].set_title(
        f"combined detrended z ±1 SEM; cross-arm r={r_cross:.3f} (p={p_cross:.3f})"
    )
    fig.savefig(args.out_dir / "dataset_a_z_analysis.png", dpi=150)
    plt.close(fig)

    with open(args.out_dir / "dataset_a_fig2a_analysis.json", "w") as stream:
        json.dump(summary, stream, indent=2)
    print(json.dumps({k: v for k, v in summary.items() if k != "targeted_priors"},
                     indent=2))
    print("targeted priors:", json.dumps(priors, indent=2))


if __name__ == "__main__":
    main()
