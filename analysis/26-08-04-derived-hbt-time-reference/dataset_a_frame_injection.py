"""End-to-end injection validation THROUGH the dataset-A curator.

The strongest answer to "is the null a pipeline artifact": inject a known
paper-style signal into REAL frames (real noise, real t0 wander, real
background), run the exact production curation chain
(build_dataset_a.detect_events_features + curate_events), and recover the
signal from the curated events only.

Two modes:
  * z: antiphase L/R population oscillation (Fig. 2(a) target). Per-frame
    per-band injected rates s_k(tau) = s_k * (1 +/- A sin(2 pi t/T)), + for
    the L bands (B1, B3), - for R (B2, B4). Recovery: band traces ->
    combined z -> targeted two-quadrature fit at the injected period.
  * bunching: same-tau cross-arm event pairs (Fig. 2(c) target). Pair rate
    per tau bin s_a s_b A(t) with A(t) = env * A_pk (1 +/- sin(2 pi t/T))/2,
    envelope 0.75 -> 1.0, opposite signs for the L and R pairs — the
    cube-level injection model of fig2c_synthetic_injection.py lifted to
    raw frames. Recovery: 4-px count cube -> beta-corrected g2 -> targeted
    fit at the injected period.

Injected blobs are REAL event patches (15x15 px, bg-subtracted, clipped at
0) sampled around production-table event centroids, so the injected PSF is
the measured one by construction. Injection times anchor on the production
per-frame t0; the curator recomputes t0 on the modified frames and the max
|t0 shift| is reported (blobs never touch the reference-pulse ROI).

Also reported: injected-event retention through detection+curation (match
within 4 px), which directly answers whether the curator discards signal.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_dataset_a as bda  # noqa: E402
import fig2c_g2_estimator as est  # noqa: E402
from dataset_a_g2 import beta_from_noise_model, build_cube  # noqa: E402
from extract_hbt_traces import (  # noqa: E402
    DEFAULT_H5,
    EVENT_COL_RANGE,
    FRAMES_DATASET,
    N_ROWS,
    PS_PER_PX,
    estimate_t0_cfd,
    load_background,
)
from fig2c_power_gate import sem_plugin  # noqa: E402
from hbt_z_targeted_fits import targeted_fit  # noqa: E402

OLD_TABLE = (
    Path("E:/Documents/Codex/joris-homodyne")
    / "analysis/26-08-04-derived-hbt-time-reference/detected-peaks"
    / "26-08-04_hbt_time_reference_events.h5"
)
DATASET_A = bda.Path(
    Path(__file__).resolve().parent
    / "detected-peaks/26-08-04_hbt_time_reference_dataset_a.h5"
)
BANDS = bda.DEFAULT_BANDS
TAU_WIN_PX = (0.0, 833.0)  # 0-820 ps
BLOB_HALF = 7  # 15x15 patches
MATCH_PX = 4.0


def build_blob_library(h5_path, background, old_table_path, n_blobs=400, seed=3):
    """15x15 bg-subtracted patches around production-table event centroids."""
    rng = np.random.default_rng(seed)
    with h5py.File(old_table_path, "r") as handle:
        x = np.asarray(handle["events/x"])
        y = np.asarray(handle["events/y"])
        amp = np.asarray(handle["events/amplitude"])
        frame = np.asarray(handle["events/frame"])
    lo, hi = np.percentile(amp, [25, 95])
    sel = np.flatnonzero((amp >= lo) & (amp <= hi))
    pick = rng.choice(sel, size=min(n_blobs, sel.size), replace=False)
    blobs = []
    with h5py.File(h5_path, "r") as handle:
        frames = handle[FRAMES_DATASET]
        for index in pick:
            fr = np.asarray(frames[int(frame[index])], dtype=np.float64) - background
            r, c = int(round(y[index])), int(round(x[index]))
            if r - BLOB_HALF < 0 or r + BLOB_HALF + 1 > fr.shape[0]:
                continue
            if c - BLOB_HALF < 0 or c + BLOB_HALF + 1 > fr.shape[1]:
                continue
            patch = fr[r - BLOB_HALF:r + BLOB_HALF + 1,
                       c - BLOB_HALF:c + BLOB_HALF + 1]
            blobs.append(np.clip(patch, 0.0, None))
    return blobs


def draw_positions(rng, band_index, n, t0, period_ps, amplitude,
                   ramp_sign=0.0):
    """Injected (row, col) positions for n events in one band of one frame.

    tau uniform in the window, modulated acceptance for the oscillation;
    rows anchored on the frame's production t0.
    """
    lo, hi = BANDS[band_index]
    cols = rng.uniform(lo + BLOB_HALF + 1, hi - BLOB_HALF - 1, size=n)
    taus = np.empty(0)
    while taus.size < n:  # rejection sampling of 1 +/- A sin
        cand = rng.uniform(*TAU_WIN_PX, size=4 * n)
        t_ps = cand * PS_PER_PX
        weight = 1.0 + ramp_sign * amplitude * np.sin(2 * np.pi * t_ps / period_ps)
        accept = rng.uniform(0.0, 1.0 + abs(amplitude), size=cand.size) < weight
        taus = np.concatenate([taus, cand[accept]])
    taus = taus[:n]
    rows = t0 + taus
    ok = (rows > BLOB_HALF + 1) & (rows < 1016 - BLOB_HALF - 1)
    return rows[ok], cols[ok]


def inject_z(frame_bg, t0, signal_rates, amplitude, period_ps, blobs, rng):
    """Antiphase population oscillation; returns list of (row, col, band)."""
    injected = []
    for k in range(4):
        sign = +1.0 if k in (0, 2) else -1.0  # L bands +, R bands -
        mean_n = signal_rates[k]  # windowed per-frame per-band mean
        n = rng.poisson(mean_n)
        if n == 0:
            continue
        rows, cols = draw_positions(rng, k, n, t0, period_ps,
                                    amplitude, ramp_sign=sign)
        for r, c in zip(rows, cols):
            blob = blobs[rng.integers(len(blobs))]
            stamp = _stamp(frame_bg, blob, r, c)
            if stamp is not None:
                injected.append((stamp[0], stamp[1], k))
    return injected


def inject_bunching(frame_bg, t0, signal_rates, a_pk, period_ps, blobs, rng):
    """Same-tau cross-arm pairs with growing envelope; returns injected list."""
    injected = []
    ramp = (TAU_WIN_PX[1] - TAU_WIN_PX[0])
    for pair, (a, b), sign in (("L", (0, 2), +1.0), ("R", (1, 3), -1.0)):
        # per-frame expected pairs: sum_tau s_a s_b A(tau)
        grid = np.arange(*TAU_WIN_PX, 1.0)
        t_ps = grid * PS_PER_PX
        env = 0.75 + 0.25 * (grid - grid[0]) / ramp
        a_t = env * a_pk * (1.0 + sign * np.sin(2 * np.pi * t_ps / period_ps)) / 2
        sa = signal_rates[a]  # windowed per-frame signal rate
        sb = signal_rates[b]
        win = TAU_WIN_PX[1] - TAU_WIN_PX[0]
        rate = (sa / win) * (sb / win) * a_t  # pair rate per 1-px bin
        n = rng.poisson(rate.sum())
        if n == 0:
            continue
        probs = rate / rate.sum()
        taus = rng.choice(grid, size=n, p=probs)
        for tau in taus:
            row = t0 + tau
            if row <= BLOB_HALF + 1 or row >= 1016 - BLOB_HALF - 1:
                continue
            for band_index in (a, b):
                lo, hi = BANDS[band_index]
                col = rng.uniform(lo + BLOB_HALF + 1, hi - BLOB_HALF - 1)
                blob = blobs[rng.integers(len(blobs))]
                stamp = _stamp(frame_bg, blob, row, col)
                if stamp is not None:
                    injected.append((stamp[0], stamp[1], band_index))
    return injected


def _stamp(frame_bg, blob, row, col):
    """Add a blob patch centered at (row, col); returns the stamped center."""
    r, c = int(round(row)), int(round(col))
    r0, r1 = r - BLOB_HALF, r + BLOB_HALF + 1
    c0, c1 = c - BLOB_HALF, c + BLOB_HALF + 1
    if r0 < 0 or c0 < 0 or r1 > frame_bg.shape[0] or c1 > frame_bg.shape[1]:
        return None
    frame_bg[r0:r1, c0:c1] += blob
    return float(r), float(c)


def match_retention(injected_rows, injected_cols, det_rows, det_cols,
                    match_px=MATCH_PX):
    """Fraction of injected blobs with a detected event within match_px."""
    retained = 0
    det = np.stack([det_rows, det_cols], axis=1) if det_rows.size else np.empty((0, 2))
    for r, c in zip(injected_rows, injected_cols):
        if det.size == 0:
            continue
        if np.any((np.abs(det[:, 0] - r) <= match_px)
                  & (np.abs(det[:, 1] - c) <= match_px)):
            retained += 1
    return retained


def run(h5_path, dataset_a_path, old_table_path, out_dir, mode, amplitude,
        period_ps, sample_step, max_frames, seed=11):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    background = load_background(h5_path)

    with h5py.File(dataset_a_path, "r") as handle:
        floor = np.asarray(handle["noise_model/floor_per_band"])
        frame_t0_old = np.asarray(handle["frame_t0"])
        n_kept_total = int(handle.attrs["n_frames_kept"])
        band_totals = np.asarray(
            [int(np.sum(np.asarray(handle["events/band"]) == k)) for k in range(4)]
        )
        band_windows = np.asarray(handle.attrs["band_windows"], dtype=np.float64)
    widths = band_windows[:, 1] - band_windows[:, 0] + 1.0
    noise_rate = floor * widths  # events/frame over the full tau grid
    total_rate = band_totals / n_kept_total
    # per-frame per-band signal rate inside the 0-820 ps injection window:
    # both signal and noise are tau-flat, so scale the whole-frame excess by
    # the window's fraction of the reachable rows (N_ROWS)
    win_frac = (TAU_WIN_PX[1] - TAU_WIN_PX[0]) / N_ROWS
    signal_rates = np.clip(total_rate - noise_rate, 0.0, None) * win_frac
    print(f"signal rates per band per frame (window): {signal_rates.round(3)}",
          flush=True)

    blobs = build_blob_library(h5_path, background, old_table_path)
    print(f"blob library: {len(blobs)} patches", flush=True)

    rng = np.random.default_rng(seed)
    frame_indices_all = np.arange(0, frame_t0_old.size, sample_step)
    if max_frames is not None:
        frame_indices_all = frame_indices_all[:max_frames]
    keep_frames = frame_indices_all[np.isfinite(frame_t0_old[frame_indices_all])]

    store = {"frame": [], "x": [], "y": [], "tau": [], "band": []}
    n_injected = 0
    n_retained = 0
    t0_shifts = []
    noise_col_hist = np.zeros(EVENT_COL_RANGE[1] - EVENT_COL_RANGE[0], dtype=np.int64)
    import time
    started = time.time()
    with h5py.File(h5_path, "r") as handle:
        frames = handle[FRAMES_DATASET]
        for done, index in enumerate(keep_frames):
            frame_bg = np.asarray(frames[int(index)], dtype=np.float64) - background
            t0_old = frame_t0_old[int(index)]
            if mode == "z":
                injected = inject_z(frame_bg, t0_old, signal_rates, amplitude,
                                    period_ps, blobs, rng)
            else:
                injected = inject_bunching(frame_bg, t0_old, signal_rates,
                                           amplitude, period_ps, blobs, rng)
            t0_new = estimate_t0_cfd(frame_bg, roi=bda.T0_ROI)
            if not np.isfinite(t0_new):
                continue
            t0_shifts.append(abs(t0_new - t0_old))
            feats = bda.detect_events_features(frame_bg, bda.DETECTION_THRESHOLD)
            keep_mask, band, reason = bda.curate_events(feats)
            if injected:
                inj = np.asarray(injected)
                n_injected += inj.shape[0]
                n_retained += match_retention(
                    inj[:, 0], inj[:, 1], feats["rows"], feats["cols"])
            col_idx = np.clip(
                np.floor(feats["cols"]).astype(np.int64) - EVENT_COL_RANGE[0],
                0, noise_col_hist.size - 1)
            np.add.at(noise_col_hist, col_idx[band < 0], 1)
            if keep_mask.any():
                store["frame"].append(np.full(int(keep_mask.sum()), done,
                                              dtype=np.int32))
                store["x"].append(feats["cols"][keep_mask])
                store["y"].append(feats["rows"][keep_mask])
                store["tau"].append(feats["rows"][keep_mask] - t0_new)
                store["band"].append(band[keep_mask].astype(np.int8))
            if (done + 1) % 500 == 0:
                print(f"  {done + 1}/{keep_frames.size} frames, "
                      f"{time.time() - started:.0f}s", flush=True)

    events = {key: np.concatenate(value) if value else np.empty(0)
              for key, value in store.items()}
    n_frames = keep_frames.size
    frame_t0 = frame_t0_old[keep_frames]
    frame_indices = np.arange(n_frames)
    col_rate = noise_col_hist / n_frames
    member = np.zeros(col_rate.size, dtype=bool)
    for lo, hi in BANDS:
        member |= (np.arange(col_rate.size) + EVENT_COL_RANGE[0] >= lo) & \
                  (np.arange(col_rate.size) + EVENT_COL_RANGE[0] <= hi)
    col_rate[member] = 0.0
    floor_inj, _, _ = bda.fit_band_floors(col_rate, bands=BANDS)

    report = {
        "mode": mode,
        "amplitude_injected": amplitude,
        "period_ps": period_ps,
        "n_frames": n_frames,
        "n_injected_blobs": int(n_injected),
        "n_retained": int(n_retained),
        "retention": float(n_retained / max(n_injected, 1)),
        "t0_shift_max_px": float(np.max(t0_shifts)) if t0_shifts else None,
        "t0_shift_median_px": float(np.median(t0_shifts)) if t0_shifts else None,
        "events_kept": int(events["frame"].size),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }

    if mode == "z":
        report["recovery"] = analyze_z(events, frame_t0, amplitude, period_ps)
    else:
        report["recovery"] = analyze_bunching(
            events, frame_t0, frame_indices, floor_inj, band_windows,
            amplitude, period_ps, n_frames)

    with open(out_dir / f"injection_{mode}_report.json", "w") as stream:
        json.dump(report, stream, indent=2)
    print(json.dumps(report, indent=2), flush=True)
    return report


def analyze_z(events, frame_t0, amplitude, period_ps):
    """Band traces -> combined z -> targeted fit at the injected period."""
    from extract_hbt_traces import mean_sem

    n_tau = bda.TAU_GRID.size
    kept_t0 = frame_t0
    valid = (bda.TAU_GRID[None, :] + kept_t0[:, None] >= 0.0) & (
        bda.TAU_GRID[None, :] + kept_t0[:, None] <= 1015)
    count = np.broadcast_to(
        valid.astype(np.float64).sum(axis=0), (4, n_tau)).copy()
    bin_index = np.floor(events["tau"] - bda.TAU_GRID[0]).astype(np.int64)
    ok = (bin_index >= 0) & (bin_index < n_tau)
    cell = events["frame"][ok].astype(np.int64) * (4 * n_tau) + (
        events["band"][ok].astype(np.int64) * n_tau + bin_index[ok])
    uniq, c = np.unique(cell, return_counts=True)
    cell_index = uniq % (4 * n_tau)
    acc_sum = np.bincount(cell_index, weights=c, minlength=4 * n_tau)
    acc_sumsq = np.bincount(cell_index, weights=c**2, minlength=4 * n_tau)
    band_means, band_sems = mean_sem(
        acc_sum.reshape(4, n_tau), acc_sumsq.reshape(4, n_tau), count)
    t = bda.TAU_GRID.astype(np.float64) * PS_PER_PX
    mask = (t > 20.0) & (t < 800.0)

    def z_arm(i, j):
        tr_i, se_i, tr_j, se_j = band_means[i], band_sems[i], band_means[j], band_sems[j]
        with np.errstate(invalid="ignore", divide="ignore"):
            z = (tr_i - tr_j) / (tr_i + tr_j)
            sz = (2 * np.sqrt((tr_j * se_i) ** 2 + (tr_i * se_j) ** 2)
                  / (tr_i + tr_j) ** 2)
        return z, sz

    z1, sz1 = z_arm(0, 1)
    z2, sz2 = z_arm(2, 3)
    with np.errstate(invalid="ignore", divide="ignore"):
        w = 1 / sz1[mask] ** 2 + 1 / sz2[mask] ** 2
        z = (z1[mask] / sz1[mask] ** 2 + z2[mask] / sz2[mask] ** 2) / w
        sz = np.sqrt(1 / w)
    tt = t[mask]
    # sparse-frame runs can have NaN/inf bins (zero valid frames or zero
    # counts); drop them before the detrend/fit
    fin = np.isfinite(z) & np.isfinite(sz) & (sz > 0)
    tt, z, sz = tt[fin], z[fin], sz[fin]
    zr = z - np.polyval(np.polyfit(tt, z, 2, w=1 / sz), tt)
    fit = targeted_fit(tt, zr, sz, period_ps)
    return {
        "targeted_fit_at_injected_period": fit,
        "n_bins_used": int(fin.size),
        "note": "recovered amplitude is in measured z space (diluted by the "
                "residual noise floor, same as any real signal)",
    }


def analyze_bunching(events, frame_t0, frame_indices, floor_inj, band_windows,
                     a_pk, period_ps, n_frames):
    """4-px cube -> beta-corrected g2 -> targeted fit at injected period."""
    bin_width = 4
    cube, nbins, frames = build_cube(events, frame_t0, frame_indices, bin_width)
    beta = beta_from_noise_model(floor_inj, band_windows, nbins, bin_width)
    bin_ps = bin_width * est.PS_PER_PX
    t_ps = (est.TAU_RANGE_PX[0] + (np.arange(nbins) + 0.5) * bin_width) * est.PS_PER_PX
    window = (t_ps >= 0.0) & (t_ps <= 820.0)
    out = {}
    for pair in est.SAME_MESA_PAIRS:
        N1, N2, b1, b2 = est.pair_counts(cube, beta, pair)
        N1 = N1.astype(np.float64)
        N2 = N2.astype(np.float64)
        g2, _ = est.g2_zero_delay(N1, N2, b1, b2)
        m1 = N1.mean(axis=0)
        m2 = N2.mean(axis=0)
        sem_plug = sem_plugin(m1, m2, b1, b2, frames)
        pooled, pooled_sem, _ = est.pooled_zero_delay(N1, N2, b1, b2, mask=window)
        sel = window & np.isfinite(g2) & np.isfinite(sem_plug) & (sem_plug > 0)
        attenuation = np.sinc(bin_ps / period_ps)
        basis = attenuation * np.exp(2j * np.pi * t_ps[sel] / period_ps)
        design = np.column_stack(
            [np.ones(int(sel.sum())), t_ps[sel], t_ps[sel] ** 2,
             basis.real, basis.imag])
        w = 1.0 / sem_plug[sel] ** 2
        xt_w = design.T * w
        cov = np.linalg.inv(xt_w @ design)
        beta_hat = cov @ (xt_w @ g2[sel])
        sigma_a = float(np.sqrt(max(cov[3, 3], cov[4, 4])))
        out[pair] = {
            "pooled_g2": pooled,
            "pooled_sem": pooled_sem,
            "fit_amplitude_at_period": float(np.hypot(beta_hat[3], beta_hat[4])),
            "fit_signed_quadrature": float(beta_hat[4]),
            "sigma_a": sigma_a,
            "injected_A_pk": a_pk,
        }
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5", type=Path, default=DEFAULT_H5)
    parser.add_argument("--dataset-a", type=Path, default=DATASET_A)
    parser.add_argument("--old-table", type=Path, default=OLD_TABLE)
    parser.add_argument("--out-dir", type=Path,
                        default=Path("drafts/26-08-08_dataset_a/injection"))
    parser.add_argument("--mode", choices=("z", "bunching"), required=True)
    parser.add_argument("--amplitude", type=float, default=0.5,
                        help="z: z-space antiphase amplitude; bunching: A_pk")
    parser.add_argument("--period-ps", type=float, default=9.84)
    parser.add_argument("--sample-step", type=int, default=7)
    parser.add_argument("--max-frames", type=int, default=None)
    args = parser.parse_args()
    run(args.h5, args.dataset_a, args.old_table, args.out_dir, args.mode,
        args.amplitude, args.period_ps, args.sample_step, args.max_frames)


if __name__ == "__main__":
    main()
