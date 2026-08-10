"""Figures and synthetic benchmark for compare_band_window_methods.py."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import compare_band_window_methods as cw


def _shade(ax, windows, color, alpha=0.22, lw=0):
    for (x0, x1) in windows:
        ax.axvspan(x0, x1, color=color, alpha=alpha, lw=lw)


def _sqrt_ax_setup(ax, hist, floor, xlim=None):
    if xlim is None:
        xlim = (cw.DET_RANGE[0] - 30, cw.DET_RANGE[1] + 30)
    cols = np.arange(hist.size)
    ax.plot(cols, np.sqrt(np.clip(hist, 0, None)), color="black", lw=0.7)
    ax.plot(cols, np.sqrt(floor), color="tab:blue", lw=1.0)
    ax.set_xlim(*xlim)
    # ignore the hot-column spikes when scaling the compressed axis
    inner = np.arange(cw.DET_RANGE[0] + 5, cw.DET_RANGE[1] - 5)
    inner = np.setdiff1d(inner, list(cw.EDGE_HOT))
    ax.set_ylim(0, np.sqrt(hist[inner].max()) * 1.14)


def fig_method_grid(batches, out_dir):
    rows = cw.METHOD_ORDER
    n_rows = len(rows)
    fig, axes = plt.subplots(
        n_rows, 2, figsize=(13.5, 2.05 * n_rows), sharex=True,
        squeeze=False,
    )
    batch_titles = {
        "main": "main batch (20 000 frames)",
        "supplement": "supplement batch (72 000 frames, structured glow)",
    }
    for col, name in enumerate(("main", "supplement")):
        b = batches[name]
        hist, floor = b["hist"], b["floor"]
        refs = b["ref_windows"]
        for row, (key, label) in enumerate(rows):
            ax = axes[row][col]
            _sqrt_ax_setup(ax, hist, floor)
            r = b["results"][key]
            if "threshold" in r and key in ("subwin5", "gauss4-mu2s",
                                            "gauss4-rule", "otsu"):
                # only thresholds in counts/col units are drawn (the matched
                # filter threshold lives in response units and would mislead)
                ax.axhline(np.sqrt(r["threshold"]), color="tab:red", ls="--",
                           lw=0.9)
            if key == "fdr-snr":
                _shade(ax, r.get("fdr_runs", []), color="tab:orange", alpha=0.18)
            n = len(r["windows"])
            ok = n == 4
            _shade(ax, r["windows"],
                   color="tab:green" if ok else "tab:red",
                   alpha=0.28 if ok else 0.18)
            for (x0, x1) in refs:
                ax.axvline(x0, color="gray", ls=":", lw=0.8)
                ax.axvline(x1, color="gray", ls=":", lw=0.8)
            if row == 0:
                ax.set_title(batch_titles[name], fontsize=11)
            if col == 0:
                ax.set_ylabel(f"{label}\nsqrt(events/col)", fontsize=8.5)
            extra = ""
            if "error" in r:
                extra = f"  [{r['error']}]"
            ax.text(
                0.01, 0.97,
                f"n={n}  " + ", ".join(f"({x0},{x1})" for x0, x1 in r["windows"]) + extra,
                transform=ax.transAxes, fontsize=7.0, va="top",
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=1),
            )
    for ax in axes[-1]:
        ax.set_xlabel("detector column [pixel]")
    handles = [
        Line2D([0], [0], color="black", lw=0.7, label="event histogram"),
        Line2D([0], [0], color="tab:blue", lw=1.0, label="common local floor"),
        Line2D([0], [0], color="tab:red", ls="--", lw=0.9, label="method threshold"),
        plt.Rectangle((0, 0), 1, 1, fc="tab:green", alpha=0.28, label="windows (n=4)"),
        plt.Rectangle((0, 0), 1, 1, fc="tab:red", alpha=0.18, label="windows (n≠4: failure)"),
        plt.Rectangle((0, 0), 1, 1, fc="tab:orange", alpha=0.18, label="FDR raw runs"),
        Line2D([0], [0], color="gray", ls=":", lw=0.8, label="reference edges"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8,
               frameon=False, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle(
        "Band-window localization methods on the photon-event column histogram"
        " (sqrt-compressed axis)",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0.02, 1, 0.995))
    path = out_dir / "band_windows_method_grid.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def fig_robustness(batches, out_dir):
    def _n(r):
        if "n_raw" in r:
            return r["n_raw"]
        if "fdr_runs" in r:
            return len(r["fdr_runs"])
        return len(r["windows"])

    sweeps = [
        ("sub-window kσ", np.arange(2.0, 7.01, 0.5),
         lambda h, f, p: cw.method_subwindow(h, k=p)),
        ("find_peaks prominence kσ", np.arange(1.0, 10.01, 0.5),
         lambda h, f, p: cw.method_findpeaks(h, f, k=p)),
        ("matched filter kσ", np.arange(1.0, 10.01, 0.5),
         lambda h, f, p: cw.method_matched(h, f, k=p)),
        ("FDR q", np.logspace(-6, np.log10(0.2), 13),
         lambda h, f, p: cw.method_fdr_snr(h, f, q=p)),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=True)
    colors = {"main": ("tab:blue", "-", "main (20k)"),
              "supplement": ("tab:red", "--", "supplement (72k)")}
    for ax, (title, params, fn) in zip(axes.ravel(), sweeps):
        for name in ("main", "supplement"):
            b = batches[name]
            ns = []
            for p in params:
                try:
                    r = fn(b["hist"], b["floor"], p)
                    ns.append(_n(r))
                except Exception:
                    ns.append(0)
            c, ls, lab = colors[name]
            ax.plot(params, ns, color=c, ls=ls, marker=".", ms=4, lw=1.2,
                    label=lab)
        ax.axhline(4, color="gray", lw=0.8, ls=":")
        ax.set_ylim(0, 10)
        ax.set_yticks(range(0, 11, 2))
        ax.set_title(title, fontsize=10)
        ax.set_ylabel("bands detected")
        ax.grid(alpha=0.3)
        if "q" in title:
            ax.set_xscale("log")
            ax.set_xlabel("FDR level q")
        else:
            ax.set_xlabel("threshold multiple k")
    axes[0][0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Detection robustness: parameter plateau yielding the four bands")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path = out_dir / "band_windows_robustness.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _growth_curves(hist, floor, centers):
    hs = np.arange(15, 101)
    curves = []
    for c in centers:
        fom = []
        for h in hs:
            x0, x1 = int(round(c - h)), int(round(c + h)) + 1
            tot = hist[x0:x1].sum()
            n = floor[x0:x1].sum()
            s = tot - n
            fom.append(s / np.sqrt(max(tot, 1.0)) if s > 0 else np.nan)
        curves.append(np.asarray(fom))
    return hs, curves


def fig_snr_growth(batches, summary, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), sharey=True)
    band_names = ["B1", "B2", "B3", "B4"]
    for ax, name in zip(axes, ("main", "supplement")):
        b = batches[name]
        centers = [0.5 * (x0 + x1 - 1) for (x0, x1) in b["ref_windows"]]
        hs, curves = _growth_curves(b["hist"], b["floor"], centers)
        arr = np.vstack(curves)
        best = np.nanmax(arr, axis=1, keepdims=True)
        for curve, bn, pk in zip(arr, band_names, best.ravel()):
            ax.plot(hs, curve / pk, lw=1.2, label=bn)
        for key, color, lab in (
            ("subwin5", "tab:red", "incumbent 5σ"),
            ("gauss4-mu2s", "tab:purple", "μ±2σ"),
            ("fdr-snr", "tab:green", "proposed FDR+S/N"),
        ):
            ws = b["results"][key]["windows"]
            if ws:
                hw = [0.5 * (x1 - x0) for (x0, x1) in ws]
                ax.plot(hw, [1.0] * len(hw), "v", color=color, ms=7,
                        label=f"{lab} widths", clip_on=False)
        ax.set_ylim(0.9, 1.02)
        ax.set_xlabel("window half-width [pixel]")
        ax.set_title(f"{name} batch", fontsize=10)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("S/√(S+N)  (normalized per band)")
    axes[0].legend(fontsize=7.5, loc="lower right")
    fig.suptitle("Aperture growth curves: where each method's window width lands")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = out_dir / "band_windows_snr_growth.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# --------------------------------------------------------------------------
# synthetic benchmark
# --------------------------------------------------------------------------

def _band_amplitudes(hist, floor, centers):
    ys = cw.smooth(hist, 2.0)
    amps = []
    for c in centers:
        lo, hi = int(c) - 10, int(c) + 11
        amps.append(float(np.clip(ys[lo:hi] - floor[lo:hi], 0, None).max()))
    return amps


def _truncate_det_range(profile):
    """Events are only recorded inside the detection column range, so the
    synthetic profile is zero outside it (as in the real histograms)."""
    p = np.array(profile, dtype=float)
    p[: cw.DET_RANGE[0]] = 0.0
    p[cw.DET_RANGE[1]:] = 0.0
    return p


def synthetic_regimes(batches):
    main, supp = batches["main"], batches["supplement"]
    mu_main = np.array([415.0, 625.0, 845.0, 1060.0])
    mu_supp = np.array([395.0, 610.0, 815.0, 1035.0])
    sigma = 36.0
    regimes = {}
    for tag, batch, mus in (("flat", main, mu_main),
                            ("struct-main", main, mu_main),
                            ("struct-supp", supp, mu_supp)):
        floor = batch["floor"]
        amps = _band_amplitudes(batch["hist"], floor, mus)
        profile = floor.copy()
        if tag == "flat":
            profile = np.full_like(floor, np.median(floor[cw.det_cols()]))
        for a, mu in zip(amps, mus):
            profile = profile + a * np.exp(-0.5 * ((np.arange(cw.N_COLS) - mu) / sigma) ** 2)
        regimes[tag] = {"profile": _truncate_det_range(profile),
                        "mus": mus, "amps": amps}
    return regimes


def _amplify_structure(regimes):
    f = regimes["struct-supp"]["profile"]
    mus, amps = regimes["struct-supp"]["mus"], regimes["struct-supp"]["amps"]
    base = f.copy()
    for a, mu in zip(amps, mus):
        base = base - a * np.exp(-0.5 * ((np.arange(cw.N_COLS) - mu) / 36.0) ** 2)
    floor_mean = np.median(base[cw.det_cols()])
    floor2 = floor_mean + 2.0 * (base - floor_mean)
    prof = np.clip(floor2, 1.0, None)
    for a, mu in zip(amps, mus):
        prof = prof + a * np.exp(-0.5 * ((np.arange(cw.N_COLS) - mu) / 36.0) ** 2)
    regimes["struct-supp-x2"] = {"profile": _truncate_det_range(prof),
                                 "mus": mus, "amps": amps}


SYN_METHODS = [
    ("subwin5", lambda h, f, mu: cw.method_subwindow(h, k=5.0)),
    ("gauss4-mu2s", lambda h, f, mu: cw.method_gauss4(h, list(mu), rule="mu2sigma")),
    ("otsu", lambda h, f, mu: cw.method_otsu(h)),
    ("findpeaks", lambda h, f, mu: cw.method_findpeaks(h, f, k=5.0)),
    ("matched", lambda h, f, mu: cw.method_matched(h, f, k=5.0)),
    ("fdr-snr", lambda h, f, mu: cw.method_fdr_snr(h, f, q=0.01)),
    ("snr-aperture", lambda h, f, mu: cw.method_snr_aperture(h, f)),
]


def run_synthetic(batches, reps, seed=20260808):
    regimes = synthetic_regimes(batches)
    _amplify_structure(regimes)
    rng = np.random.default_rng(seed)
    results = {}
    for tag, reg in regimes.items():
        profile = reg["profile"]
        mus = reg["mus"]
        # oracle windows on the noiseless profile
        floor_true = profile.copy()
        for a, mu in zip(reg["amps"], mus):
            floor_true = floor_true - a * np.exp(
                -0.5 * ((np.arange(cw.N_COLS) - mu) / 36.0) ** 2)
        oracle = []
        for mu in mus:
            w = cw._snr_optimal_window(profile, floor_true, mu)
            oracle.append(w)
        oracle_fom = cw.score_windows(profile, floor_true, oracle, mus)["fom"]
        stats_by_method = {k: {"success": 0, "n4": 0, "rmse": [], "fom_ratio": []}
                           for k, _ in SYN_METHODS}
        for _ in range(reps):
            hist = rng.poisson(np.clip(profile, 0, None)).astype(float)
            floor = cw.local_floor(hist)
            for key, fn in SYN_METHODS:
                try:
                    r = fn(hist, floor, mus)
                    ws = r["windows"]
                except Exception:
                    ws = []
                st = stats_by_method[key]
                if len(ws) == 4:
                    st["n4"] += 1
                    centers = np.array([0.5 * (x0 + x1 - 1) for (x0, x1) in ws])
                    if np.all(np.abs(centers - mus) < 30):
                        st["success"] += 1
                        st["rmse"].append(np.sqrt(np.mean((centers - mus) ** 2)))
                        sc = cw.score_windows(hist, floor, ws, mus)
                        st["fom_ratio"].append(sc["fom"] / oracle_fom)
        results[tag] = {
            "methods": {
                k: {
                    "n4_rate": st["n4"] / reps,
                    "success_rate": st["success"] / reps,
                    "center_rmse_median": float(np.median(st["rmse"])) if st["rmse"] else None,
                    "fom_ratio_median": float(np.median(st["fom_ratio"])) if st["fom_ratio"] else None,
                }
                for k, st in stats_by_method.items()
            },
            "mus": mus.tolist(),
        }
        print(f"  synthetic {tag}: " + "; ".join(
            f"{k} ok={stats_by_method[k]['success']}/{reps}" for k, _ in SYN_METHODS))
    return results


def fig_synthetic(syn_results, out_dir):
    tags = list(syn_results.keys())
    methods = [k for k, _ in SYN_METHODS]
    x = np.arange(len(methods))
    width = 0.8 / len(tags)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for i, tag in enumerate(tags):
        res = syn_results[tag]["methods"]
        succ = [res[m]["success_rate"] for m in methods]
        rmse = [res[m]["center_rmse_median"] or np.nan for m in methods]
        axes[0].bar(x + i * width - 0.4 + width / 2, succ, width, label=tag)
        axes[1].bar(x + i * width - 0.4 + width / 2, rmse, width, label=tag)
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("success rate (4 bands, centers within 30 px)")
    axes[1].set_ylabel("center RMSE on successes [pixel]")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(methods, rotation=25, ha="right", fontsize=8)
        ax.grid(alpha=0.3, axis="y")
    axes[1].set_ylim(0, max(8.0, np.nanmax([
        syn_results[t]["methods"][m]["center_rmse_median"] or 0
        for t in tags for m in methods]) * 1.15))
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(tags), fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.995))
    fig.suptitle("Synthetic recovery benchmark (Poisson realizations of "
                 "realistic structured baselines + 4 known bands)", y=0.93)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    path = out_dir / "band_windows_synthetic.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def make_all_figures(batches, summary, out_dir, synthetic_reps):
    paths = [fig_method_grid(batches, out_dir),
             fig_robustness(batches, out_dir),
             fig_snr_growth(batches, summary, out_dir)]
    syn = run_synthetic(batches, synthetic_reps)
    summary["synthetic"] = syn
    paths.append(fig_synthetic(syn, out_dir))
    import json

    class _NpEnc(json.JSONEncoder):
        def default(self, o):
            if isinstance(o, (np.floating, np.integer)):
                return o.item()
            return super().default(o)

    with open(out_dir / "band_windows_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, cls=_NpEnc)
    for p in paths:
        print("wrote", p)
