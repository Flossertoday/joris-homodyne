"""Band-window localization method comparison (26-08-04 HBT main + supplement).

Compares candidate algorithms that turn the photon-event column histogram into
band integration windows, on both 26-08-04 HBT batches:

- main batch (20000 frames, drafts/26-08-06_fig2a_cfd50/bands_report.json):
  the incumbent sub-window 5-sigma count-threshold method works;
- supplement batch (72000 frames, drafts/26-08-06_fig2a_supplement/
  bands_report.json): the incumbent fails (structured out-of-band glow) and
  windows had to be fixed manually on the aligned sub-window profile
  (docs/26-08-06-hbt-supplement.md section 4).

Candidate methods:

- subwin5 / subwin3   incumbent sub-window count threshold (BAND_SIGMA = 5,
                      and the 3-sigma variant);
- gauss4-mu2s         constant + 4 Gaussians fit, windows = mu +/- 2 sigma
                      (26-07-31 acceptance-window semantics);
- gauss4-rule         the actual 26-07-31 rule: largest k <= 3 (0.5 steps)
                      with disjoint mu +/- k sigma intervals;
- otsu                Otsu threshold on the smoothed histogram;
- findpeaks           scipy find_peaks on the baseline-subtracted histogram,
                      windows = half-prominence widths;
- matched             SExtractor-style matched filtering (Gaussian kernel on
                      the baseline-subtracted histogram), threshold on the
                      filtered response;
- fdr-snr             proposed: local sigma-clipped baseline, per-column
                      Poisson p-values with Benjamini-Hochberg FDR control,
                      then S/sqrt(S+N)-optimal width per band (aperture
                      photometry growth curve);
- snr-aperture        aperture-photometry sibling: centers from the matched
                      filter, same S/N-optimal width scan.

Scoring uses one common, method-independent local baseline (iterative
sigma-clipped rolling mean calibrated against out-of-band columns) and the
S/sqrt(S+N) figure of merit already used for the 2026-08-08 window-width scan
(docs/26-08-08-dataset-a.md section 1). Reference window sets: the incumbent
windows on the main batch and the manual windows on the supplement batch.

Outputs (--out-dir):
- band_windows_method_grid.png   per-method panels on sqrt-compressed axis
- band_windows_robustness.png    detection plateau vs threshold parameter
- band_windows_snr_growth.png    S/sqrt(S+N) growth curves + method widths
- band_windows_synthetic.png     synthetic-recovery benchmark
- band_windows_summary.json      machine-readable scores

Run (from the repository/worktree root; JSON paths are repo-relative):

    .\\.venv\\Scripts\\python.exe analysis\\26-08-04-derived-hbt-time-reference\\compare_band_window_methods.py \
        --main-report drafts/26-08-06_fig2a_cfd50/bands_report.json \
        --supp-report drafts/26-08-06_fig2a_supplement/bands_report.json \
        --out-dir drafts/26-08-08_band_window_methods
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage, optimize, signal, stats

DET_RANGE = (350, 1150)  # detect_events column range (both batches)
EDGE_HOT = (350, 351, 1147, 1148)  # known hot edge columns (dataset-A doc)
N_COLS = 1344

MAIN_REF_WINDOWS = [(380, 450), (570, 680), (810, 880), (1010, 1110)]
SUPP_MANUAL_WINDOWS = [(350, 440), (540, 680), (770, 860), (970, 1100)]

SUBWINDOW_EDGES = np.arange(DET_RANGE[0], DET_RANGE[1] + 1, 10)  # 80 x 10 px


# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------

def load_hist(report_path):
    det = json.load(open(report_path))["band_detection"]
    return np.asarray(det["event_col_hist"], dtype=np.float64)


def det_cols():
    return np.arange(DET_RANGE[0], DET_RANGE[1])


def local_floor(hist, win=201, n_iter=3):
    """Iterative sigma-clipped local baseline for the event column histogram.

    Rolling 20th-percentile start; columns more than 3 sigma (Poisson) above
    the current baseline are masked as band columns; the baseline is the
    rolling mean of the surviving floor columns (median-filter smoothed).
    Returns the per-column baseline over the full frame (outside the
    detection range it is a constant extrapolation).
    """
    x = np.arange(hist.size)
    det = (x >= DET_RANGE[0]) & (x < DET_RANGE[1])
    y = hist.astype(np.float64)
    base = ndimage.percentile_filter(y, percentile=20, size=win, mode="nearest")
    mask = det.copy()
    for _ in range(n_iter):
        smooth = ndimage.median_filter(base, size=51, mode="nearest")
        lam = np.maximum(smooth, 1.0)
        mask = det & (y <= lam + 3.0 * np.sqrt(lam))
        # rolling mean of surviving columns via masked cumulative sums
        vals = np.where(mask, y, np.nan)
        base = _rolling_nanmean(vals, win)
    # fill any remaining NaN (fully-masked stretches) by interpolation
    good = ~np.isnan(base)
    base = np.interp(x, x[good], base[good])
    base[~det] = np.interp(x[~det], x[det], base[det])
    return base


def _rolling_nanmean(vals, win):
    n = vals.size
    half = win // 2
    out = np.full(n, np.nan)
    cs_v = np.nancumsum(np.where(np.isnan(vals), 0.0, vals))
    cs_c = np.nancumsum(~np.isnan(vals))
    lo = np.clip(np.arange(n) - half, 0, n - 1)
    hi = np.clip(np.arange(n) + half + 1, 0, n)
    cnt = _window_sum(cs_c, lo, hi)
    tot = _window_sum(cs_v, lo, hi)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = tot / cnt
    return out


def _window_sum(cs, lo, hi):
    cs = np.concatenate(([0.0], cs))
    return cs[np.clip(hi, 0, cs.size - 1)] - cs[np.clip(lo, 0, cs.size - 1)]


def smooth(y, sigma=3.0):
    return ndimage.gaussian_filter1d(y, sigma, mode="nearest")


def runs_from_mask(mask, cols, bridge=10, min_width=20, drop_edges=True,
                   bridge_mask=None):
    """Contiguous True runs -> (x0, x1) windows; bridge short gaps.

    bridge: gaps of False columns of length <= bridge are filled when
    bridge_mask (default: all gap columns) holds for the whole gap.
    """
    m = mask.copy()
    if bridge > 0:
        i = 0
        n = m.size
        while i < n:
            if not m[i]:
                j = i
                while j < n and not m[j]:
                    j += 1
                gap = j - i
                if 0 < i and j < n and gap <= bridge:
                    ok = bridge_mask[i:j].all() if bridge_mask is not None else True
                    if ok:
                        m[i:j] = True
                i = j
            else:
                i += 1
    windows = []
    start = None
    for i in range(m.size + 1):
        if i < m.size and m[i] and start is None:
            start = i
        elif (i == m.size or not m[i]) and start is not None:
            x0, x1 = int(cols[start]), int(cols[i - 1]) + 1
            if x1 - x0 >= min_width:
                windows.append((x0, x1))
            start = None
    if drop_edges:
        windows = [
            (x0, x1) for (x0, x1) in windows
            if x0 > DET_RANGE[0] + 1 and x1 < DET_RANGE[1] - 1
        ]
    return windows


# --------------------------------------------------------------------------
# candidate methods: hist -> dict(windows=[(x0, x1)], plot extras...)
# --------------------------------------------------------------------------

def method_subwindow(hist, k=5.0, bridge_k=2.0):
    """Incumbent: per-sub-window counts > low-40% baseline + k*sigma."""
    edges = SUBWINDOW_EDGES
    n_cols = edges[1] - edges[0]
    sub = np.array([hist[e0:e1].sum() for e0, e1 in zip(edges[:-1], edges[1:])])
    ordered = np.sort(sub)
    low = ordered[: max(2, int(np.ceil(0.4 * ordered.size)))]
    median, sigma_robust = float(np.mean(low)), float(np.std(low, ddof=1))
    sigma = max(sigma_robust, np.sqrt(max(median * n_cols, 1.0)))
    thr = median + k * sigma
    member = sub > thr
    bridge = sub > median + bridge_k * sigma
    for w in range(1, member.size - 1):
        if not member[w] and member[w - 1] and member[w + 1] and bridge[w]:
            member[w] = True
    windows = []
    start = None
    for w in range(member.size + 1):
        if w < member.size and member[w] and start is None:
            start = w
        elif (w == member.size or not member[w]) and start is not None:
            windows.append((int(edges[start]), int(edges[w])))
            start = None
    windows = [(x0, x1) for (x0, x1) in windows
               if x0 > int(edges[0]) and x1 < int(edges[-1])]
    return {"windows": windows, "threshold": thr / n_cols,
            "baseline": median / n_cols}


def _four_gaussians(x, c, *params):
    out = np.full_like(x, c, dtype=float)
    for i in range(4):
        a, mu, s = params[3 * i: 3 * i + 3]
        out = out + a * np.exp(-0.5 * ((x - mu) / s) ** 2)
    return out


def _fit_four_gaussians(hist, mu0):
    cols = det_cols()
    y = hist[cols]
    p0 = [float(np.percentile(y, 10))]
    for mu in mu0:
        p0 += [float(y.max()) * 0.8, float(mu), 30.0]
    lo, hi = [0.0], [np.inf]
    for mu in mu0:
        lo += [0.0, mu - 60, 3.0]
        hi += [np.inf, mu + 60, 80.0]
    popt, _ = optimize.curve_fit(
        _four_gaussians, cols.astype(float), y, p0=p0, bounds=(lo, hi),
        sigma=np.sqrt(np.maximum(y, 1.0)), maxfev=120000,
    )
    return popt


def method_gauss4(hist, mu0, rule="mu2sigma"):
    """26-07-31-style constant + 4 Gaussians fit windows."""
    try:
        popt = _fit_four_gaussians(hist, mu0)
    except Exception as exc:  # fit failure is a method outcome
        return {"windows": [], "error": f"fit failed: {exc}"}
    mus = popt[2::3]
    sigs = popt[3::3]
    if rule == "mu2sigma":
        ks = [2.0] * 4
    else:  # largest k <= 3 (0.5 steps) with disjoint intervals
        for k in np.arange(3.0, 0.9, -0.5):
            iv = [(mu - k * s, mu + k * s) for mu, s in zip(mus, sigs)]
            if all(iv[i][1] < iv[i + 1][0] for i in range(3)):
                ks = [k] * 4
                break
        else:
            ks = [1.0] * 4
    windows = [
        (int(round(mu - k * s)), int(round(mu + k * s)))
        for mu, s, k in zip(mus, sigs, ks)
    ]
    return {"windows": windows, "mus": mus.tolist(), "sigmas": sigs.tolist(),
            "k": ks}


def _otsu_threshold(values, n_bins=256):
    hist, edges = np.histogram(values, bins=n_bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    w = hist.astype(np.float64)
    w0 = np.cumsum(w)
    w1 = w0[-1] - w0
    m0 = np.cumsum(w * centers) / np.maximum(w0, 1e-12)
    m1 = (np.cumsum(w * centers)[-1] - np.cumsum(w * centers)) / np.maximum(w1, 1e-12)
    between = w0 * w1 * (m0 - m1) ** 2
    return centers[int(np.argmax(between))]


def method_otsu(hist, smooth_sigma=3.0):
    y = smooth(hist, smooth_sigma)
    cols = det_cols()
    thr = _otsu_threshold(y[cols])
    mask = y[cols] > thr
    windows = runs_from_mask(mask, cols, bridge=10, min_width=20)
    return {"windows": windows, "threshold": thr}


def _mask_edge_hot(hist, floor):
    """Replace the known hot edge columns by the local floor value (defect
    masking, same practice as the dataset-A curation) so peak-based methods
    do not latch onto the 1-px spikes at 350/351/1147/1148.
    """
    y = hist.copy()
    y[list(EDGE_HOT)] = floor[list(EDGE_HOT)]
    return y


def method_findpeaks(hist, floor, k=5.0, smooth_sigma=3.0, n_bands=4):
    """find_peaks on the floor-subtracted histogram.

    The prominence floor is Poisson-aware (per-column floor noise scaled by
    the smoothing-kernel energy), not a MAD over the segment — the latter is
    contaminated by band flanks and the structured floor. Four bands are
    expected on physics grounds (two mesas x two HBT arms), so the n_bands
    most prominent peaks are kept; windows = half-prominence widths.
    """
    cols = det_cols()
    res = _mask_edge_hot(hist, floor) - floor
    ys = smooth(res, smooth_sigma)
    seg = ys[cols]
    var_eff = np.median(np.maximum(floor[cols], 1.0)) / (
        2 * np.sqrt(np.pi) * smooth_sigma)
    prom_min = k * np.sqrt(var_eff)
    peaks, props = signal.find_peaks(seg, prominence=prom_min, distance=150)
    n_raw = peaks.size
    if peaks.size > n_bands:
        keep = np.argsort(props["prominences"])[-n_bands:]
        peaks = np.sort(peaks[keep])
    if peaks.size == 0:
        return {"windows": [], "prominence": float(prom_min), "n_raw": 0}
    widths = signal.peak_widths(seg, peaks, rel_height=0.5)
    windows = []
    for p, w in zip(peaks, widths[0]):
        x0 = int(round(cols[p] - w / 2))
        x1 = int(round(cols[p] + w / 2))
        windows.append((max(x0, DET_RANGE[0] + 2), min(x1, DET_RANGE[1] - 4)))
    windows.sort()
    return {"windows": windows, "prominence": float(prom_min),
            "n_raw": int(n_raw)}


def method_matched(hist, floor, k=5.0, sigma_band=36.0, n_bands=4):
    """SExtractor-style matched filtering: Gaussian kernel matched to the band
    width on the floor-subtracted histogram; the n_bands strongest response
    peaks above a Poisson-aware height become bands; windows = half-height
    segments of each response peak.
    """
    cols = det_cols()
    res = smooth(_mask_edge_hot(hist, floor) - floor, 2.0)
    half = int(3 * sigma_band)
    t = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (t / sigma_band) ** 2)
    kernel /= kernel.sum()
    resp = signal.fftconvolve(res, kernel, mode="same")[cols]
    var_resp = np.median(np.maximum(floor[cols], 1.0)) * float((kernel ** 2).sum())
    med = np.median(resp)
    thr = med + k * np.sqrt(var_resp)
    peaks, _ = signal.find_peaks(resp, height=thr, distance=150)
    n_raw = peaks.size
    if peaks.size > n_bands:
        keep = np.argsort(resp[peaks])[-n_bands:]
        peaks = np.sort(peaks[keep])
    windows = []
    for p in peaks:
        level = med + 0.5 * (resp[p] - med)
        lo, hi = p, p
        while lo > 0 and resp[lo - 1] > level:
            lo -= 1
        while hi < resp.size - 1 and resp[hi + 1] > level:
            hi += 1
        x0, x1 = int(cols[lo]), int(cols[hi]) + 1
        windows.append((max(x0, DET_RANGE[0] + 2), min(x1, DET_RANGE[1] - 4)))
    windows.sort()
    return {"windows": windows, "threshold": float(thr), "n_raw": int(n_raw)}


def _poisson_fdr_columns(hist, floor, q):
    """Per-column Poisson p-values vs local floor, Benjamini-Hochberg mask."""
    cols = det_cols()
    lam = np.maximum(floor[cols], 1e-6)
    n = hist[cols]
    p = stats.poisson.sf(n - 1.0, lam)
    tested = ~np.isin(cols, EDGE_HOT)
    pt = p[tested]
    order = np.argsort(pt)
    m = pt.size
    k_star = 0
    for rank, idx in enumerate(order, start=1):
        if pt[idx] <= (rank / m) * q:
            k_star = rank
    sig = np.zeros(cols.size, bool)
    if k_star > 0:
        cutoff = pt[order[k_star - 1]]
        sig[tested] = p[tested] <= cutoff
    loose = p <= 0.05
    return cols, sig, loose


def _snr_optimal_window(hist, floor, center, h_min=15, h_max=100, cap=None):
    """Aperture-photometry growth curve: half-width maximizing S/sqrt(S+N).

    The scan never includes the detection-range edge hot columns
    (350/351/1147/1148, known defects), so windows stay inside
    [DET_RANGE[0]+2, DET_RANGE[1]-4].
    """
    if cap is not None:
        h_max = min(h_max, cap)
    best = None
    for h in range(h_min, h_max + 1):
        x0 = int(round(center - h))
        x1 = int(round(center + h)) + 1
        if x0 < DET_RANGE[0] + 2 or x1 > DET_RANGE[1] - 4:
            continue
        tot = hist[x0:x1].sum()
        n = floor[x0:x1].sum()
        s = tot - n
        if s <= 0:
            continue
        fom = s / np.sqrt(tot)
        if best is None or fom > best[0]:
            best = (fom, x0, x1, s, n)
    if best is None:
        return None
    return (best[1], best[2])


def method_fdr_snr(hist, floor, q=0.01):
    """Proposed: FDR column detection + S/N-optimal width refinement."""
    cols, sig, loose = _poisson_fdr_columns(hist, floor, q)
    runs = runs_from_mask(sig, cols, bridge=10, min_width=20,
                          bridge_mask=loose)
    windows = []
    for (x0, x1) in runs:
        seg = hist[x0:x1] - floor[x0:x1]
        seg = np.clip(seg, 0.0, None)
        center = x0 + float(np.average(np.arange(seg.size), weights=seg)) \
            if seg.sum() > 0 else 0.5 * (x0 + x1)
        w = _snr_optimal_window(hist, floor, center)
        windows.append(w if w is not None else (x0, x1))
    windows.sort()
    return {"windows": windows, "fdr_runs": runs}


def method_snr_aperture(hist, floor, sigma_band=36.0, n_bands=4):
    """Matched-filter centers + S/N-optimal widths (aperture photometry)."""
    cols = det_cols()
    res = smooth(_mask_edge_hot(hist, floor) - floor, 2.0)
    half = int(3 * sigma_band)
    t = np.arange(-half, half + 1)
    kernel = np.exp(-0.5 * (t / sigma_band) ** 2)
    kernel /= kernel.sum()
    resp = signal.fftconvolve(res, kernel, mode="same")[cols]
    peaks, _ = signal.find_peaks(resp, distance=150)
    if peaks.size == 0:
        return {"windows": []}
    top = peaks[np.argsort(resp[peaks])[-n_bands:]]
    centers = np.sort(cols[top]).astype(float)
    windows = []
    for i, c in enumerate(centers):
        cap = None
        if i + 1 < centers.size:
            cap = int((centers[i + 1] - c) / 2) - 5
        w = _snr_optimal_window(hist, floor, c, cap=cap)
        if w is not None:
            windows.append(w)
    windows.sort()
    return {"windows": windows, "centers": centers.tolist()}


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------

def score_windows(hist, floor, windows, ref_centers):
    """Common, method-independent scores for a window set."""
    out = {"n_bands": len(windows), "windows": [list(w) for w in windows]}
    if not windows:
        return out
    centers = np.array([0.5 * (x0 + x1 - 1) for (x0, x1) in windows])
    out["centers"] = centers.tolist()
    tot = sum(hist[x0:x1].sum() for (x0, x1) in windows)
    n = sum(floor[x0:x1].sum() for (x0, x1) in windows)
    s = tot - n
    out["S"] = float(s)
    out["N"] = float(n)
    out["fom"] = float(s / np.sqrt(max(tot, 1.0)))
    if len(windows) == len(ref_centers):
        ref = np.asarray(sorted(ref_centers), dtype=float)
        out["center_offsets"] = (centers - ref).tolist()
        out["center_rmse"] = float(np.sqrt(np.mean((centers - ref) ** 2)))
    return out


def ref_centers_for(batch):
    if batch == "main":
        return [0.5 * (a + b - 1) for (a, b) in MAIN_REF_WINDOWS]
    return [0.5 * (a + b - 1) for (a, b) in SUPP_MANUAL_WINDOWS]


# --------------------------------------------------------------------------
# main driver
# --------------------------------------------------------------------------

METHOD_ORDER = [
    ("subwin5", "sub-window 5σ (incumbent)"),
    ("gauss4-mu2s", "4-Gaussian fit μ±2σ (26-07-31)"),
    ("gauss4-rule", "4-Gaussian fit, largest k≤3 rule"),
    ("otsu", "Otsu threshold"),
    ("findpeaks", "find_peaks prominence 5σ"),
    ("matched", "matched filter 5σ"),
    ("fdr-snr", "Poisson FDR + S/N-optimal width (proposed)"),
    ("snr-aperture", "matched centers + S/N-optimal width"),
]


def run_methods(hist, floor, mu0):
    res = {}
    res["subwin5"] = method_subwindow(hist, k=5.0)
    res["subwin3"] = method_subwindow(hist, k=3.0)
    g = method_gauss4(hist, mu0, rule="mu2sigma")
    res["gauss4-mu2s"] = g
    res["gauss4-rule"] = method_gauss4(hist, mu0, rule="largest_k")
    res["otsu"] = method_otsu(hist)
    res["findpeaks"] = method_findpeaks(hist, floor, k=5.0)
    res["matched"] = method_matched(hist, floor, k=5.0)
    res["fdr-snr"] = method_fdr_snr(hist, floor, q=0.01)
    res["snr-aperture"] = method_snr_aperture(hist, floor)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--main-report", default="drafts/26-08-06_fig2a_cfd50/bands_report.json")
    ap.add_argument("--supp-report", default="drafts/26-08-06_fig2a_supplement/bands_report.json")
    ap.add_argument("--out-dir", default="drafts/26-08-08_band_window_methods")
    ap.add_argument("--synthetic-reps", type=int, default=50)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    batches = {
        "main": {
            "hist": load_hist(args.main_report),
            "mu0": [415.0, 625.0, 845.0, 1060.0],
            "ref_windows": MAIN_REF_WINDOWS,
            "ref_label": "incumbent windows",
        },
        "supplement": {
            "hist": load_hist(args.supp_report),
            "mu0": [395.0, 610.0, 815.0, 1035.0],
            "ref_windows": SUPP_MANUAL_WINDOWS,
            "ref_label": "manual windows",
        },
    }

    summary = {"batches": {}}
    for name, b in batches.items():
        hist = b["hist"]
        floor = local_floor(hist)
        b["floor"] = floor
        refs = ref_centers_for("main" if name == "main" else "supplement")
        results = run_methods(hist, floor, b["mu0"])
        b["results"] = results
        scored = {}
        for key, r in results.items():
            scored[key] = dict(r)
            scored[key].update(score_windows(hist, floor, r["windows"], refs))
        summary["batches"][name] = {
            "scores": scored,
            "ref_windows": b["ref_windows"],
        }
        print(f"== {name}: floor range {floor[det_cols()].min():.0f}.."
              f"{floor[det_cols()].max():.0f} events/col")
        for key, _ in METHOD_ORDER:
            s = scored[key]
            print(f"  {key:12s} n={s['n_bands']} fom={s.get('fom', float('nan')):.1f} "
                  f"rmse={s.get('center_rmse', float('nan'))} windows={s['windows']}")
        s3 = scored["subwin3"]
        print(f"  {'subwin3':12s} n={s3['n_bands']} fom={s3.get('fom', float('nan')):.1f} "
              f"windows={s3['windows']}")

    with open(out_dir / "band_windows_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    from compare_band_window_figures import make_all_figures
    make_all_figures(batches, summary, out_dir, args.synthetic_reps)


if __name__ == "__main__":
    main()
