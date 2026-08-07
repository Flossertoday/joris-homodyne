"""Shared g2(0)(t) estimator core for the 26-08-04 Fig. 2(c) replica.

Implements the same Supplementary Note 1 estimator of Adiyatullin et al.,
Nat. Commun. 8:1329 (2017) already validated for the 26-07-31 dataset in
``analysis/26-07-31-derived-single-sweep/compute_g2_matrix.py``, adapted to
the 26-08-04 HBT condition:

  * classes are the four emission bands B1..B4 (fixed windows from the Fig.
    2(a) pipeline); the same-mesa pairs are L = B1 x B3 and R = B2 x B4
    (the two HBT arms of each mesa);
  * events are binned in the jitter-corrected coordinate tau = y - t0
    (per-frame CFD50 reference-pulse time zero), not in absolute rows;
  * the residual background beta is estimated from events inside the
    detection column range but outside the four bands (frames were already
    measured-background-subtracted before event detection, so beta only
    carries detector-noise clusters that pass the event threshold).

Only the quantities Fig. 2(c) needs are implemented: the zero-delay diagonal
g2(0)(t) with the paper-style three-term propagated SEM, and the
window-pooled g2(0) with a leave-one-frame-out jackknife SEM. The full 2-D
map g2(t1, t2) (paper Fig. 3) and cross-mesa correlations are out of scope
(docs/adiyatullin2017-fig2c.md, G9). Diagonal-only keeps the per-frame
algebra O(frames x bins) with no FFT, which matters at 20000 frames.

Conventions match compute_g2_matrix.py: N1, N2 are (frames, bins) per-frame
count matrices; beta1, beta2 are per-frame background expectations per bin;
bins with non-positive background-corrected intensity are NaN.
"""
import os

import numpy as np

from streak_time_axis import ps_per_px

# Acquisition metadata (metadata-conventions.toml [streak_camera.time_range_0804]):
# the HBT condition is single-sweep, 1 ns window over 1016 rows.
PS_PER_PX = ps_per_px("01_HBT_condensate_Josephson_regime")


def _band_windows() -> np.ndarray:
    """Definitive band windows from the Fig. 2(a) pipeline (docs/26-08-04-fig2a.md).

    The 2026-08-06 supplement acquisition (New folder (10), 72000 frames) shows
    the same four-band geometry but with slightly shifted/wider bump edges; set
    FIG2C_BAND_WINDOWS to a comma list of 2N column edges (e.g.
    "370,450,540,680,800,870,980,1100") to override for that dataset.
    """
    override = os.environ.get("FIG2C_BAND_WINDOWS")
    if override:
        values = [float(item) for item in override.split(",") if item.strip()]
        if len(values) % 2 or not values:
            raise ValueError(f"FIG2C_BAND_WINDOWS needs an even number of edges: {override!r}")
        return np.asarray(values, dtype=np.float64).reshape(-1, 2)
    return np.array(
        [[380, 450], [570, 680], [810, 880], [1010, 1110]], dtype=np.float64
    )


BAND_WINDOWS = _band_windows()
BAND_NAMES = ("B1", "B2", "B3", "B4")
# Same-mesa pairs correlate the two independent HBT arms of one mesa:
# L = B1 x B3, R = B2 x B4.
SAME_MESA_PAIRS = {"L": (0, 2), "R": (1, 3)}

# Event-detection column range of the Fig. 2(a) pipeline.
EVENT_COL_RANGE = (350, 1150)

# Aligned tau window (px) covering the full sweep plus margin, as in
# extract_hbt_traces.TAU_GRID.
TAU_RANGE_PX = (-64, 1080)

# Bins whose zero-delay SEM exceeds this are marked statistics-limited in
# figures (26-07-31 precedent; analogue of the paper's SNR<20 gray band).
MAX_DISPLAY_SEM = 0.5


def assign_bands(x: np.ndarray, windows=BAND_WINDOWS) -> np.ndarray:
    """Band index per event (-1 outside all windows; windows disjoint)."""
    band = np.full(len(x), -1, dtype=np.int64)
    for index, (lo, hi) in enumerate(windows):
        band[(x >= lo) & (x <= hi)] = index
    return band


def aligned_bin_index(
    tau: np.ndarray, bin_width: int, tau_lo: float = TAU_RANGE_PX[0]
) -> tuple[np.ndarray, int]:
    """Integer bin index per aligned coordinate tau = y - t0 (-1 outside)."""
    nbins = int(np.ceil((TAU_RANGE_PX[1] - TAU_RANGE_PX[0]) / bin_width))
    index = np.floor((tau - tau_lo) / bin_width).astype(np.int64)
    index[(index < 0) | (index >= nbins)] = -1
    return index, nbins


def count_cube_aligned(
    frame: np.ndarray,
    band: np.ndarray,
    tau: np.ndarray,
    frame_count: int,
    bin_width: int,
) -> tuple[np.ndarray, int]:
    """(frame_count, 4, nbins) per-frame event counts per band and tau bin."""
    bin_index, nbins = aligned_bin_index(tau, bin_width)
    valid = (band >= 0) & (bin_index >= 0)
    flat = (frame[valid].astype(np.int64) * len(BAND_NAMES) + band[valid]) * nbins
    flat = flat + bin_index[valid]
    cube = np.bincount(flat, minlength=frame_count * len(BAND_NAMES) * nbins)
    return cube.reshape(frame_count, len(BAND_NAMES), nbins).astype(np.int32), nbins


def background_beta_aligned(
    tau: np.ndarray,
    x: np.ndarray,
    frame_count: int,
    bin_width: int,
    windows=BAND_WINDOWS,
    col_range=EVENT_COL_RANGE,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Expected residual-background events per frame, band and tau bin.

    The density rho(tau) is measured on the events inside the detection
    column range but outside every band window, resolved along the aligned
    coordinate; beta_c(tau) = rho(tau) * x-width of band c. Returns
    (beta[4, nbins], rho[nbins], flatness chi2/ndf over bins).
    """
    bin_index, nbins = aligned_bin_index(tau, bin_width)
    outside = np.ones(len(x), dtype=bool)
    for lo, hi in windows:
        outside &= ~((x >= lo) & (x <= hi))
    inside_range = (x >= col_range[0]) & (x <= col_range[1])
    keep = outside & inside_range & (bin_index >= 0)
    outside_width = float(col_range[1] - col_range[0])
    for lo, hi in windows:
        outside_width -= float(hi - lo)
    hist = np.bincount(bin_index[keep], minlength=nbins).astype(np.float64)
    rho = hist / (frame_count * outside_width)
    widths = windows[:, 1] - windows[:, 0]
    beta = rho[None, :] * widths[:, None]
    expected = hist.mean()
    chi2_ndf = float(
        np.sum((hist - expected) ** 2) / max(expected, 1e-12) / max(nbins - 1, 1)
    )
    return beta, rho, chi2_ndf


def g2_zero_delay(
    N1: np.ndarray, N2: np.ndarray, beta1: np.ndarray, beta2: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Background-corrected g2(0)(t) diagonal and its propagated SEM.

    Same estimator as compute_g2_matrix.g2_map restricted to the diagonal
    t1 = t2: g2(t) = <s1 s2>(t) / (<s1>(t) <s2>(t)) with the signal moments
    corrected by the residual background expectations beta. SEM propagates
    the standard errors of <n1>, <n2> and <n1*n2> (the paper's approach).
    """
    N1 = np.asarray(N1, dtype=np.float64)
    N2 = np.asarray(N2, dtype=np.float64)
    frames = N1.shape[0]
    m1 = N1.mean(axis=0)
    m2 = N2.mean(axis=0)
    s1 = m1 - beta1
    s2 = m2 - beta2
    product = (N1 * N2).mean(axis=0)
    corrected = product - m1 * beta2 - beta1 * m2 + beta1 * beta2
    denom = s1 * s2
    with np.errstate(divide="ignore", invalid="ignore"):
        g2 = corrected / denom

    ddof = frames / (frames - 1.0)
    var_product = np.clip((N1**2 * N2**2).mean(axis=0) - product**2, 0.0, None)
    var_product = var_product / frames * ddof
    var_m1 = np.clip((N1**2).mean(axis=0) - m1**2, 0.0, None) / frames * ddof
    var_m2 = np.clip((N2**2).mean(axis=0) - m2**2, 0.0, None) / frames * ddof

    with np.errstate(divide="ignore", invalid="ignore"):
        d_product = 1.0 / denom
        d_m1 = -beta2 / denom - g2 / s1
        d_m2 = -beta1 / denom - g2 / s2
        variance = (
            d_product**2 * var_product + d_m1**2 * var_m1 + d_m2**2 * var_m2
        )
        sem = np.sqrt(variance)
    invalid = (s1 <= 0.0) | (s2 <= 0.0)
    g2[invalid] = np.nan
    sem[invalid] = np.nan
    return g2, sem


def pooled_zero_delay(
    N1: np.ndarray,
    N2: np.ndarray,
    beta1: np.ndarray,
    beta2: np.ndarray,
    mask: np.ndarray | None = None,
) -> tuple[float, float, int]:
    """Window-pooled g2(0) = sum_t <s1 s2> / sum_t <s1><s2> with jackknife SEM.

    The pooled ratio-of-sums is the tau = 0 entry of the tau profile in
    compute_g2_matrix.g2_tau_profile. mask selects the bins entering the
    window sums (default: all bins with positive corrected intensity). SEM
    is a leave-one-frame-out jackknife, which carries all covariances
    between the aggregated numerator and denominator. Returns
    (g2_pooled, sem, n_bins_used).
    """
    N1 = np.asarray(N1, dtype=np.float64)
    N2 = np.asarray(N2, dtype=np.float64)
    frames, bins = N1.shape
    m1 = N1.mean(axis=0)
    m2 = N2.mean(axis=0)
    s1 = m1 - beta1
    s2 = m2 - beta2
    product = (N1 * N2).mean(axis=0)
    corrected = product - m1 * beta2 - beta1 * m2 + beta1 * beta2
    denom = s1 * s2

    valid = np.isfinite(denom) & (denom > 0.0)
    if mask is not None:
        valid &= mask
    n_used = int(valid.sum())
    if n_used == 0:
        return float("nan"), float("nan"), 0

    # Leave-one-frame-out moments via the sum-subtract trick.
    sum1 = N1.sum(axis=0)
    sum2 = N2.sum(axis=0)
    sum12 = (N1 * N2).sum(axis=0)
    m1_loo = (sum1[None, :] - N1) / (frames - 1.0)
    m2_loo = (sum2[None, :] - N2) / (frames - 1.0)
    product_loo = (sum12[None, :] - N1 * N2) / (frames - 1.0)
    corrected_loo = (
        product_loo
        - m1_loo * beta2[None, :]
        - beta1[None, :] * m2_loo
        + (beta1 * beta2)[None, :]
    )
    denom_loo = (m1_loo - beta1[None, :]) * (m2_loo - beta2[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        numerator = corrected_loo[:, valid].sum(axis=1)
        denominator = denom_loo[:, valid].sum(axis=1)
        g2_loo = numerator / denominator
    g2 = float(corrected[valid].sum() / denom[valid].sum())
    sem = float(
        np.sqrt((frames - 1.0) / frames * np.sum((g2_loo - g2_loo.mean()) ** 2))
    )
    return g2, sem, n_used


def pair_counts(
    cube: np.ndarray, beta: np.ndarray, pair: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(N1, N2, beta1, beta2) for a same-mesa pair ('L' or 'R')."""
    first, second = SAME_MESA_PAIRS[pair]
    return cube[:, first, :], cube[:, second, :], beta[first], beta[second]
