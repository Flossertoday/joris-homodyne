"""ROI-interval correctness: band emission morphology vs method windows.

The method comparison in compare_band_window_methods.py scores detection
count and S/sqrt(S+N); this script answers the more direct question: which
method places the window EDGES at the right columns, judged against the
measured emission morphology?

Inputs per batch:
- photon-event column histogram (bands_report.json; Poisson-clean but only
  inside detection columns 350-1150);
- full-column analog profile (compute_full_column_profiles.py; covers all
  1344 columns, lumpy fixed-pattern baseline, but no truncation).

Per band we define the morphological EMISSION SUPPORT: the contiguous region
around the band peak where the smoothed excess over the local floor exceeds
max(3 x Poisson sigma of the local floor, 1.5% of the band peak excess).
The analog profile gives the same walk independently (and beyond the
detection-range edges), cross-checking the support and measuring how much
emission the 350/1150 detection range clips.

Per method per batch we then report:
- detection: n bands (gate);
- coverage: fraction of the support's total excess inside the windows;
- floor fraction: local-floor counts inside the windows / total counts
  inside the windows;
- edge deviation: window edge minus support edge per band side (px;
  negative = edge inside support, tail deliberately excluded; positive =
  edge outside support, floor swallowed);
- S/sqrt(S+N) FOM (from compare_band_window_methods).

Figures (--out-dir):
- roi_fullrange_overview.png  full 1344-column analog profile with the
  detection range, reference-spot region, emission supports and windows
- roi_band_morphology.png     per-band zoom panels (4 bands x 2 batches):
  event histogram (sqrt) + scaled analog profile + local floor + support
  shading + one window strip per method

Run from the repository/worktree root:

    python analysis/26-08-08-band-window-methods/render_roi_morphology.py \
        --main-report drafts/26-08-06_fig2a_cfd50/bands_report.json \
        --supp-report drafts/26-08-06_fig2a_supplement/bands_report.json \
        --profiles-dir drafts/26-08-08_band_window_methods \
        --out-dir drafts/26-08-08_band_window_methods
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

import compare_band_window_methods as cw

REFERENCE_SPOT_COLS = (140, 270)
SUPPORT_MIN_FRAC = 0.015
SUPPORT_SMOOTH = 3.0


def emission_support(resid, floor, center, det_edges=None, edge_margin=3):
    """Contiguous emission support around a band center.

    resid: smoothed (hist - floor); floor: local floor per column.
    Walks outward from the peak while resid > max(3*sqrt(floor),
    SUPPORT_MIN_FRAC * peak). Returns (x0, x1, truncated_left,
    truncated_right, threshold) in full-frame column coordinates.

    det_edges caps the walk (for the event histogram, which only exists
    inside detection cols). Smoothing bleeds the zero outside-region into
    the first few in-range columns, so a support edge landing within
    edge_margin of the boundary is also flagged as truncated — the flag
    means "the walk may have been stopped by the detection boundary, check
    the full-range analog profile".
    """
    c = int(round(center))
    lo_search = max(0, c - 60)
    peak_i = lo_search + int(np.argmax(resid[lo_search : c + 61]))
    peak = resid[peak_i]
    thr = max(3.0 * np.sqrt(max(floor[peak_i], 1.0)), SUPPORT_MIN_FRAC * peak)
    lo, hi = peak_i, peak_i
    while lo - 1 >= 0 and resid[lo - 1] > thr:
        lo -= 1
    while hi + 1 < resid.size and resid[hi + 1] > thr:
        hi += 1
    trunc_l = trunc_r = False
    if det_edges is not None:
        if lo <= det_edges[0] + edge_margin:
            lo, trunc_l = det_edges[0], True
        if hi + 1 >= det_edges[1] - edge_margin:
            hi, trunc_r = det_edges[1] - 1, True
    return lo, hi + 1, trunc_l, trunc_r, thr


def analog_residual(profile):
    baseline = ndimage.median_filter(profile, size=201, mode="nearest")
    return profile - baseline, baseline


def method_metrics(hist, floor, windows, supports):
    """Coverage / floor-fraction / edge deviations of a window set."""
    if not windows:
        return {"n_bands": 0}
    res = np.clip(hist - floor, 0.0, None)
    tot_in = sum(hist[x0:x1].sum() for (x0, x1) in windows)
    n_in = sum(floor[x0:x1].sum() for (x0, x1) in windows)
    cov_num = cov_den = 0.0
    edge_dev = []
    for (x0, x1), (s0, s1, *_rest) in zip(windows, supports):
        inter0, inter1 = max(x0, s0), min(x1, s1)
        if inter1 > inter0:
            cov_num += res[inter0:inter1].sum()
        cov_den += res[s0:s1].sum()
        edge_dev.append([x0 - s0, x1 - s1])
    out = {
        "n_bands": len(windows),
        "coverage": float(cov_num / max(cov_den, 1e-9)),
        "floor_fraction": float(n_in / max(tot_in, 1e-9)),
    }
    if len(windows) == len(supports):
        out["edge_dev_px"] = edge_dev
    return out


def evaluate_batch(hist, floor, results, ref_windows, batch_name):
    refs = [0.5 * (a + b - 1) for (a, b) in ref_windows]
    res_smooth = ndimage.gaussian_filter1d(
        cw._mask_edge_hot(hist, floor) - floor, SUPPORT_SMOOTH, mode="nearest"
    )
    supports = []
    for c in refs:
        supports.append(
            emission_support(res_smooth, floor, c, det_edges=cw.DET_RANGE)
        )
    metrics = {}
    for key, r in results.items():
        metrics[key] = method_metrics(hist, floor, r["windows"], supports)
    metrics["reference"] = method_metrics(hist, floor, ref_windows, supports)
    return {
        "batch": batch_name,
        "supports": [
            {"x0": s[0], "x1": s[1], "trunc_left": s[2], "trunc_right": s[3],
             "threshold": s[4]}
            for s in supports
        ],
        "metrics": metrics,
    }, supports


# ----------------------------------------------------------------------
# figures
# ----------------------------------------------------------------------

METHOD_ROWS = [
    ("subwin5", "sub-window 5σ (incumbent)"),
    ("subwin3", "sub-window 3σ"),
    ("gauss4-mu2s", "4-Gaussian μ±2σ"),
    ("gauss4-rule", "4-Gaussian k≤3 rule"),
    ("otsu", "Otsu"),
    ("findpeaks", "find_peaks"),
    ("matched", "matched filter"),
    ("fdr-snr", "FDR + S/N (proposed)"),
    ("snr-aperture", "matched + S/N"),
]

ROW_COLORS = {
    "subwin5": "tab:red", "subwin3": "tab:red",
    "gauss4-mu2s": "tab:purple", "gauss4-rule": "tab:purple",
    "otsu": "tab:brown", "findpeaks": "tab:orange", "matched": "tab:orange",
    "fdr-snr": "tab:green", "snr-aperture": "tab:green",
    "reference": "black",
}


def fig_fullrange(batches_info, out_dir):
    fig, axes = plt.subplots(len(batches_info), 1, figsize=(13.5, 4.2 * len(batches_info)),
                             sharex=True, squeeze=False)
    for ax, (name, info) in zip(axes.ravel(), batches_info.items()):
        prof = info["analog"]
        cols = np.arange(prof.size)
        ax.plot(cols, np.sqrt(np.clip(prof, 0, None)), color="black", lw=0.7)
        ax.axvspan(*cw.DET_RANGE, color="tab:blue", alpha=0.06)
        ax.axvspan(*REFERENCE_SPOT_COLS, color="tab:red", alpha=0.10)
        for s in info["supports"]:
            ax.axvspan(s[0], s[1], color="gray", alpha=0.15)
        for (x0, x1) in info["ref_windows"]:
            ax.axvline(x0, color="tab:red", ls="--", lw=0.8)
            ax.axvline(x1, color="tab:red", ls="--", lw=0.8)
        # keep the band/baseline structure readable: the reference-spot peak
        # is allowed off-scale
        masked = np.sqrt(np.clip(prof, 0, None)).copy()
        masked[REFERENCE_SPOT_COLS[0] : REFERENCE_SPOT_COLS[1]] = np.nan
        top = np.nanmax(masked)
        ax.set_ylim(0, top * 1.08)
        spot_peak = np.sqrt(max(prof[REFERENCE_SPOT_COLS[0] : REFERENCE_SPOT_COLS[1]].max(), 0))
        ax.annotate(
            f"reference spot peak off-scale (sqrt ~ {spot_peak:.0f})",
            xy=(0.5 * sum(REFERENCE_SPOT_COLS), top * 0.9), fontsize=7.5,
            ha="center", color="tab:red",
        )
        ax.set_ylabel("sqrt(counts)")
        ax.set_title(
            f"{name} batch — full-frame analog profile (bg-subtracted, "
            f"rows 32-986); red = reference-spot region, blue = detection "
            f"columns {cw.DET_RANGE[0]}-{cw.DET_RANGE[1]}, gray = emission support",
            fontsize=9.5,
        )
        ax.set_xlim(0, 1344)
    axes.ravel()[-1].set_xlabel("detector column [pixel]")
    fig.suptitle(
        f"What the detection column range {cw.DET_RANGE[0]}-{cw.DET_RANGE[1]} "
        f"contains and excludes",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = out_dir / "roi_fullrange_overview.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def fig_band_morphology(batches_info, out_dir):
    fig, axes = plt.subplots(
        4, 2, figsize=(13.5, 13.2), squeeze=False,
        gridspec_kw={"hspace": 0.42, "wspace": 0.14},
    )
    band_names = ["B1", "B2", "B3", "B4"]
    rows = METHOD_ROWS + [("reference", "definitive (ref)")]
    for col, (name, info) in enumerate(batches_info.items()):
        hist, floor = info["hist"], info["floor"]
        resid_analog_scaled = info["analog_res_scaled"]
        for row, (bn, s) in enumerate(zip(band_names, info["supports"])):
            ax = axes[row][col]
            s0, s1 = s[0], s[1]
            z0, z1 = max(cw.DET_RANGE[0], s0 - 70), min(cw.DET_RANGE[1], s1 + 70)
            cols = np.arange(z0, z1)
            ymax = float(np.sqrt(np.clip(hist[z0:z1], 0, None)).max())
            ax.plot(cols, np.sqrt(np.clip(hist[z0:z1], 0, None)),
                    color="black", lw=1.0, label="event histogram")
            ax.plot(cols, np.sqrt(floor[z0:z1]), color="tab:blue", lw=1.0,
                    label="local floor")
            if resid_analog_scaled is not None:
                ana = np.clip(
                    ndimage.gaussian_filter1d(resid_analog_scaled, 2.0) + floor,
                    0, None)
                ax.plot(cols, np.sqrt(ana[z0:z1]), color="tab:cyan", lw=1.0,
                        alpha=0.85, label="analog profile (scaled, smoothed)")
            ax.axvspan(s0, s1, color="gray", alpha=0.18, label="emission support")
            ana_s = info.get("analog_supports")
            if ana_s is not None:
                a0, a1 = ana_s[row][0], ana_s[row][1]
                ax.axvline(a0, color="tab:cyan", ls="--", lw=1.1,
                           label="analog support edges" if row == 0 else None)
                ax.axvline(a1, color="tab:cyan", ls="--", lw=1.1)
            if s[2]:
                ax.axvline(cw.DET_RANGE[0], color="tab:red", lw=1.4)
            if s[3]:
                ax.axvline(cw.DET_RANGE[1], color="tab:red", lw=1.4)
            # method window strips in the reserved negative-y band of the panel
            ax.axhline(0, color="0.6", lw=0.5)
            band_center = 0.5 * (s0 + s1)
            for i, (key, label) in enumerate(rows):
                y = -(0.10 + 0.062 * i) * ymax
                ws = (info["ref_windows"] if key == "reference"
                      else info["results"][key]["windows"])
                color = ROW_COLORS[key]
                # match this panel's band by nearest window center (indexing
                # by row misaligns when a method returns partial windows)
                match = None
                for (x0, x1) in ws:
                    if abs(0.5 * (x0 + x1 - 1) - band_center) < 80:
                        match = (x0, x1)
                        break
                if match is not None:
                    x0, x1 = match
                    lo, hi = max(x0, z0), min(x1, z1)
                    if lo < hi:
                        ax.plot([lo, hi], [y, y], color=color, lw=3.2,
                                solid_capstyle="butt", clip_on=False)
                    if x0 < z0:  # window extends beyond the panel: arrow it
                        ax.plot([z0], [y], marker="<", color=color, ms=4,
                                clip_on=False)
                    if x1 > z1:
                        ax.plot([z1], [y], marker=">", color=color, ms=4,
                                clip_on=False)
                else:
                    ax.text(z0 + 3, y, "(no window)", color=color, fontsize=6,
                            va="center")
            ax.set_ylim(-(0.10 + 0.062 * len(rows)) * ymax, ymax * 1.05)
            ax.set_yticks([0, ymax])
            ax.set_yticklabels(["0", f"{ymax:.0f}"], fontsize=7)
            ax.set_xlim(z0, z1)
            trunc = []
            if s[2]:
                trunc.append("LEFT EDGE AT RANGE BOUNDARY")
            if s[3]:
                trunc.append("RIGHT EDGE AT RANGE BOUNDARY")
            title = f"{bn}  support ({s0},{s1})"
            if row == 0:
                title = f"{name} batch  " + title
            if trunc:
                title += "  [" + "; ".join(trunc) + "]"
            ax.set_title(title, fontsize=9)
            if row == 3:
                ax.set_xlabel("detector column [pixel]")
        axes[0][col].legend(fontsize=7, loc="upper right")
    handles = [
        plt.Line2D([0], [0], color=ROW_COLORS[key], lw=3.2, label=label)
        for key, label in rows
    ]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(
        "Per-band emission morphology and method ROI edges "
        "(gray = emission support; strips below zero = window spans, "
        "top-to-bottom in legend reading order)",
        fontsize=11, y=0.998,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.975))
    path = out_dir / "roi_band_morphology.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--main-report", default="drafts/26-08-06_fig2a_cfd50/bands_report.json")
    ap.add_argument("--supp-report", default="drafts/26-08-06_fig2a_supplement/bands_report.json")
    ap.add_argument("--profiles-dir", default="drafts/26-08-08_band_window_methods")
    ap.add_argument("--out-dir", default="drafts/26-08-08_band_window_methods")
    ap.add_argument("--det-range", type=int, nargs=2, default=None,
                    help="detection column range (default: as archived, 350 1150)")
    args = ap.parse_args()
    if args.det_range is not None:
        cw.configure(args.det_range)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdir = Path(args.profiles_dir)
    batches = {}
    for name, report, mu0, refs in (
        ("main", args.main_report, [415.0, 625.0, 845.0, 1060.0], cw.MAIN_REF_WINDOWS),
        ("supplement", args.supp_report, [395.0, 610.0, 815.0, 1035.0], cw.SUPP_MANUAL_WINDOWS),
    ):
        profile_path = pdir / f"full_column_profile_{name}.npy"
        if not profile_path.exists():
            print(f"skip {name}: {profile_path} not found")
            continue
        hist = cw.load_hist(report)
        floor = cw.local_floor(hist)
        results = cw.run_methods(hist, floor, mu0)
        analog = np.load(profile_path)
        batches[name] = {
            "hist": hist, "floor": floor, "results": results,
            "analog": analog, "ref_windows": refs,
        }
    if not batches:
        raise SystemExit("no batch profiles available")

    summary = {}
    for name, b in batches.items():
        ev, supports = evaluate_batch(
            b["hist"], b["floor"], b["results"], b["ref_windows"], name)
        summary[name] = ev
        b["supports"] = supports
        # analog overlay scaled to the event residual over all supports
        ana_res, ana_base = analog_residual(b["analog"])
        b["analog_res"] = ana_res
        ev_res = np.clip(b["hist"] - b["floor"], 0, None)
        mask = np.zeros(b["hist"].size, bool)
        for (s0, s1, *_r) in supports:
            mask[max(0, s0 - 10) : s1 + 10] = True
        mask &= (np.arange(b["hist"].size) >= cw.DET_RANGE[0]) & (
            np.arange(b["hist"].size) < cw.DET_RANGE[1])
        scale = (ev_res[mask] * ana_res[mask]).sum() / max(
            (ana_res[mask] ** 2).sum(), 1e-9)
        b["analog_res_scaled"] = ana_res * scale
        # analog-side support beyond the detection range (full-frame walk)
        ana_supports = []
        refs = [0.5 * (a + bb - 1) for (a, bb) in b["ref_windows"]]
        ana_smooth = ndimage.gaussian_filter1d(ana_res, SUPPORT_SMOOTH)
        for c in refs:
            ana_supports.append(
                emission_support(ana_smooth, ana_base, c, det_edges=None))
        summary[name]["analog_supports"] = [
            {"x0": s[0], "x1": s[1], "threshold": s[4]} for s in ana_supports
        ]
        b["analog_supports"] = ana_supports
        print(f"== {name}")
        print("  event supports :", [(s[0], s[1], s[2], s[3]) for s in supports])
        print("  analog supports:", [(s[0], s[1]) for s in ana_supports])
        for key, _ in METHOD_ROWS + [("reference", "ref")]:
            m = ev["metrics"][key]
            if m["n_bands"] == 4:
                print(
                    f"  {key:12s} cov={m['coverage']:.3f} floor%={m['floor_fraction']:.3f} "
                    f"edge_dev={['%+d/%+d' % tuple(d) for d in m.get('edge_dev_px', [])]}"
                )
            else:
                print(f"  {key:12s} n={m['n_bands']} (failed)")

    with open(out_dir / "roi_morphology_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    paths = [fig_fullrange(batches, out_dir), fig_band_morphology(batches, out_dir)]
    for p in paths:
        print("wrote", p)


if __name__ == "__main__":
    main()
