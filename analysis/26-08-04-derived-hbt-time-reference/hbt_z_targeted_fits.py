"""Targeted sinusoid fits on the HBT z traces at the prior periods.

Produces the |A| / 95% UL numbers quoted at the prior periods in
docs/26-08-04-fig2a.md and docs/26-08-06-hbt-supplement.md: inverse-variance
combined z of the two mesas (photon or analog band traces), quadratic
detrend, two-quadrature weighted fit at each exact prior period.

Priors: the 9.84 ps paper Josephson period plus the measured k-space fringe
periods read from results/kspace_fringe_oscillation.json via
kspace_fringe_priors — do not hardcode them (the periods scale with the
synchroscan window; 37.2/41.8 -> 74.4/83.7 ps on 2026-08-07).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from kspace_fringe_priors import measured_periods_ps  # noqa: E402

PAPER_PERIOD_PS = 9.84  # Adiyatullin 2017 Josephson period (resonant driving)


def combined_z(t, bands, sems, mask):
    """Inverse-variance combination of the two mesas' z traces on `mask`.

    Mesa arms are bands (0, 1) and (2, 3) (L = B1+B3, R = B2+B4 per-mesa
    sums enter as the two z components), matching the fig2a-era estimator.
    """

    def z_arm(i, j):
        tr_i, se_i = bands[i], sems[i]
        tr_j, se_j = bands[j], sems[j]
        z = (tr_i - tr_j) / (tr_i + tr_j)
        sz = (
            2 * np.sqrt((tr_j * se_i) ** 2 + (tr_i * se_j) ** 2)
            / (tr_i + tr_j) ** 2
        )
        return z, sz

    z1, sz1 = z_arm(0, 1)
    z2, sz2 = z_arm(2, 3)
    w = 1 / sz1[mask] ** 2 + 1 / sz2[mask] ** 2
    z = (z1[mask] / sz1[mask] ** 2 + z2[mask] / sz2[mask] ** 2) / w
    return t[mask], z, np.sqrt(1 / w)


def detrend_quadratic(tt, z, sz):
    c = np.polyfit(tt, z, 2, w=1 / sz)
    return z - np.polyval(c, tt)


def targeted_fit(tt, zr, sz, period):
    """Two-quadrature weighted sinusoid amplitude at one exact period."""
    X = np.column_stack(
        [np.cos(2 * np.pi * tt / period), np.sin(2 * np.pi * tt / period)]
    )
    W = 1 / sz**2
    XtW = X.T * W
    cov = np.linalg.inv(XtW @ X)
    beta = cov @ (XtW @ zr)
    amp = float(np.hypot(*beta))
    sigma = float(np.sqrt(max(cov[0, 0], cov[1, 1])))
    return {"amplitude": amp, "sigma": sigma, "ul95": amp + 2.0 * sigma}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--npz",
        type=Path,
        required=True,
        help="hbt_traces npz with band_traces_/band_sem_ arrays",
    )
    parser.add_argument(
        "--source", choices=("photon", "analog", "both"), default="both"
    )
    parser.add_argument(
        "--window", type=float, nargs=2, default=(20.0, 800.0), metavar=("T0", "T1")
    )
    parser.add_argument(
        "--periods",
        type=float,
        nargs="+",
        default=None,
        help="explicit fit periods [ps]; default: 9.84 paper period plus the "
        "measured k-space fringe periods from the tracked results JSON",
    )
    args = parser.parse_args()

    periods = args.periods or [PAPER_PERIOD_PS, *measured_periods_ps()]
    d = np.load(args.npz)
    t = np.asarray(d["time_ps"], dtype=np.float64)
    mask = (t > args.window[0]) & (t < args.window[1])
    sources = ("photon", "analog") if args.source == "both" else (args.source,)
    for src in sources:
        bands = d[f"band_traces_{src}"]
        sems = d[f"band_sem_{src}"]
        tt, z, sz = combined_z(t, bands, sems, mask)
        zr = detrend_quadratic(tt, z, sz)
        print(
            f"[{src}] window {args.window[0]:.0f}-{args.window[1]:.0f} ps, "
            f"{int(mask.sum())} 1-px bins"
        )
        for period in periods:
            r = targeted_fit(tt, zr, sz, period)
            print(
                f"  P={period:7.2f} ps: |A|={r['amplitude']:.4f}, "
                f"sigma={r['sigma']:.4f}, 95% UL={r['ul95']:.4f}"
            )


if __name__ == "__main__":
    main()
