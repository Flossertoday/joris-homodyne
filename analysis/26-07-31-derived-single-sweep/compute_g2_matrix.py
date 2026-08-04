"""Compute the time-resolved g2 correlation matrix from the HBT peak table.

Replicates the photon-counting g2 procedure of Adiyatullin et al., Nat.
Commun. 8:1329 (2017), Supplementary Note 1, on the 26-07-31 single-sweep
0.18-condition peak tables, and extends it to the full mesa correlation
matrix. The four spots (left to right) are condensate A x HBT arm 1,
condensate B x arm 1, condensate A x arm 2, condensate B x arm 2.

Per laser pulse (frame) and time bin t, n_c(t) is the number of detected
photons in class c. The second-order correlation is

    g2(t1, t2) = <n1(t1) * n2(t2)> / (<n1(t1)> * <n2(t2)>)

with the average over frames. Same-mesa correlations (AA, BB) correlate the
two independent HBT arms of one condensate, which makes the estimator
normally ordered (no self-coincidence term). Cross-mesa correlations (AB)
correlate the arm-summed counts of the two condensates; all terms are
cross-detector there as well. g2_BA(t1, t2) = g2_AB(t2, t1) by construction
and is not stored separately for the map; both orientations are stored for
the tau profile.

The window x-ranges contain a uniform detection background (~14% of kept
counts). With n = s + b and the background independent of the signal and of
the other arm, the signal moments are

    <s_i>        = <n_i> - beta_i
    <s1 * s2>    = <n1 * n2> - <n1>*beta2 - beta1*<n2> + beta1*beta2

where beta_c(t) is the expected background count per frame in class c and
bin t, estimated from the full peak table outside the four x-windows.
``--no-background-correction`` disables this (beta = 0).

Errors on the 2-D map propagate the standard errors of <n1>, <n2> and
<n1*n2> (the paper's approach; at our per-bin occupancies the neglected
covariance terms contribute <1%). Pooled quantities -- the tau profile
g2(tau) aggregated over t and, at tau = 0, the whole-window pooled g2(0) --
use a leave-one-frame-out jackknife, which captures all covariances.

Example:

    python compute_g2_matrix.py peaks_xrange2p0sigma.h5 peaks.h5 \
        --output-h5 g2.h5 --figure-dir results --time-bins-px 1 3 6
"""
import argparse
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fit_spot_x_ranges import load_peaks  # noqa: E402

# metadata-conventions.toml [streak_camera.time_range]: 1016 px span 1 ns.
WINDOW_NS = 1.0
# Interval order of the four spots: condensate (A or B) x HBT arm (1 or 2).
CLASS_NAMES = ("A1", "B1", "A2", "B2")
# Same-mesa pairs correlate the two independent HBT arms of one condensate.
SAME_MESA_PAIRS = {"AA": (0, 2), "BB": (1, 3)}
# Cross-mesa pairs correlate arm-summed condensate totals.
CROSS_MESA_PAIRS = {"AB": (0, 2, 1, 3)}  # A total = A1+A2, B total = B1+B2

DEFAULT_TIME_BINS_PX = (1, 3, 6)
# Bins whose zero-delay SEM exceeds this are marked as statistics-limited in
# figures (analogue of the gray low-SNR region in the paper's Fig. 2a).
MAX_DISPLAY_SEM = 0.5


def assign_classes(x: np.ndarray, intervals: np.ndarray) -> np.ndarray:
    """Class index per peak (-1 outside all intervals; intervals disjoint)."""
    cls = np.full(len(x), -1, dtype=np.int64)
    for index, (lo, hi) in enumerate(intervals):
        cls[(x >= lo) & (x <= hi)] = index
    return cls


def time_bin_index(y: np.ndarray, height: int, bin_width: int) -> tuple[np.ndarray, int]:
    """Integer time-bin index per peak and the bin count covering [0, height)."""
    nbins = int(np.ceil(height / bin_width))
    index = np.clip((y // bin_width).astype(np.int64), 0, nbins - 1)
    return index, nbins


def count_cube(
    frame: np.ndarray,
    cls: np.ndarray,
    bin_index: np.ndarray,
    nbins: int,
    frame_count: int,
) -> np.ndarray:
    """(frame_count, 4, nbins) per-frame photon counts per class and bin."""
    valid = cls >= 0
    flat = (frame[valid].astype(np.int64) * len(CLASS_NAMES) + cls[valid]) * nbins
    flat = flat + bin_index[valid]
    cube = np.bincount(flat, minlength=frame_count * len(CLASS_NAMES) * nbins)
    return cube.reshape(frame_count, len(CLASS_NAMES), nbins).astype(np.float64)


def background_beta(
    full_x: np.ndarray,
    full_y: np.ndarray,
    intervals: np.ndarray,
    detector_width: int,
    frame_count: int,
    bin_width: int,
    nbins: int,
    height: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Expected background counts per frame, class and time bin.

    Background density is estimated from full-table detections whose x falls
    outside every acceptance interval, resolved along y. Returns
    (beta[4, nbins], rho_y[nbins], y_flatness_chi2_ndf). rho_y is the
    background density per frame per x-pixel per y-bin; beta_c(t) =
    rho_y(t) * window width of class c. The chi2/ndf measures how flat the
    background is along y (1.0 = consistent with flat).
    """
    outside = np.ones(len(full_x), dtype=bool)
    for lo, hi in intervals:
        outside &= ~((full_x >= lo) & (full_x <= hi))
    outside_width = float(detector_width - np.sum(intervals[:, 1] - intervals[:, 0]))
    ybin, _ = time_bin_index(full_y[outside], height, bin_width)
    hist = np.bincount(ybin, minlength=nbins).astype(np.float64)
    rho_y = hist / (frame_count * outside_width)
    widths = intervals[:, 1] - intervals[:, 0]
    beta = rho_y[None, :] * widths[:, None]
    expected = hist.mean()
    chi2_ndf = float(np.sum((hist - expected) ** 2) / max(expected, 1e-12) / max(nbins - 1, 1))
    return beta, rho_y, chi2_ndf


def g2_map(
    N1: np.ndarray,
    N2: np.ndarray,
    beta1: np.ndarray,
    beta2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Background-corrected g2(t1, t2) map and its propagated SEM.

    N1, N2 are (frames, bins) count matrices; beta1, beta2 are the per-frame
    background expectations per bin. SEM propagates the standard errors of
    the three ensemble means (paper-style; covariance terms are negligible
    at our per-bin occupancies). Bins with non-positive corrected intensity
    are NaN.
    """
    frames = N1.shape[0]
    m1 = N1.mean(axis=0)
    m2 = N2.mean(axis=0)
    s1 = m1 - beta1
    s2 = m2 - beta2
    product = (N1.T @ N2) / frames
    corrected = (
        product
        - np.outer(m1, beta2)
        - np.outer(beta1, m2)
        + np.outer(beta1, beta2)
    )
    denom = np.outer(s1, s2)
    with np.errstate(divide="ignore", invalid="ignore"):
        g2 = corrected / denom

    ddof = frames / (frames - 1.0)
    var_product = np.clip(((N1 ** 2).T @ (N2 ** 2)) / frames - product ** 2, 0.0, None)
    var_product = var_product / frames * ddof
    var_m1 = np.clip((N1 ** 2).mean(axis=0) - m1 ** 2, 0.0, None) / frames * ddof
    var_m2 = np.clip((N2 ** 2).mean(axis=0) - m2 ** 2, 0.0, None) / frames * ddof

    with np.errstate(divide="ignore", invalid="ignore"):
        d_product = 1.0 / denom
        d_m1 = -beta2[None, :] / denom - g2 / s1[:, None]
        d_m2 = -beta1[:, None] / denom - g2 / s2[None, :]
        variance = (
            d_product ** 2 * var_product
            + d_m1 ** 2 * var_m1[:, None]
            + d_m2 ** 2 * var_m2[None, :]
        )
        sem = np.sqrt(variance)
    invalid = (s1[:, None] <= 0.0) | (s2[None, :] <= 0.0)
    g2[invalid] = np.nan
    sem[invalid] = np.nan
    return g2, sem


def _cross_corr_all(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Batched linear cross-correlation: out[f, tau] = sum_n a[f, n+tau] b[f, n].

    tau runs -(B-1)..B-1 along axis 1. Zero-padded FFT, so the result equals
    the linear (non-circular) correlation.
    """
    bins = a.shape[1]
    length = 1 << (2 * bins - 1).bit_length()
    spectrum_a = np.fft.rfft(a, length, axis=1)
    spectrum_b = np.fft.rfft(b, length, axis=1)
    full = np.fft.irfft(spectrum_a * np.conj(spectrum_b), length, axis=1)
    out = np.empty((a.shape[0], 2 * bins - 1), dtype=np.float64)
    out[:, bins - 1:] = full[:, :bins]  # tau = 0..B-1
    out[:, :bins - 1] = full[:, length - bins + 1:]  # tau = -(B-1)..-1
    return out


def g2_tau_profile(
    N1: np.ndarray,
    N2: np.ndarray,
    beta1: np.ndarray,
    beta2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pooled g2(tau) = sum_t <s1(t) s2(t+tau)> / sum_t <s1(t)> <s2(t+tau)>.

    Returns (tau_index, g2_tau, sem_tau, contributing_bins). The tau = 0
    entry is the whole-window pooled zero-delay g2(0). SEM is a
    leave-one-frame-out jackknife over the 1000 realizations, which carries
    all covariances between the aggregated numerator and denominator.
    """
    frames, bins = N1.shape
    m1 = N1.mean(axis=0)
    m2 = N2.mean(axis=0)
    s1 = m1 - beta1
    s2 = m2 - beta2

    coincidence = _cross_corr_all(N2, N1)  # per frame: sum_t N1(t) N2(t+tau)
    c_m1_b2 = _cross_corr_all(beta2[None, :], m1[None, :])[0]
    c_b1_m2 = _cross_corr_all(m2[None, :], beta1[None, :])[0]
    c_b1_b2 = _cross_corr_all(beta2[None, :], beta1[None, :])[0]
    numerator = coincidence.mean(axis=0) - c_m1_b2 - c_b1_m2 + c_b1_b2
    denominator = _cross_corr_all(s2[None, :], s1[None, :])[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        g2 = numerator / denominator

    # Leave-one-frame-out jackknife.
    m1_loo = (N1.sum(axis=0)[None, :] - N1) / (frames - 1.0)
    m2_loo = (N2.sum(axis=0)[None, :] - N2) / (frames - 1.0)
    coincidence_loo = (coincidence.sum(axis=0)[None, :] - coincidence) / (frames - 1.0)
    c_m1_b2_loo = _cross_corr_all(
        np.broadcast_to(beta2, (frames, bins)), m1_loo
    )
    c_b1_m2_loo = _cross_corr_all(
        m2_loo, np.broadcast_to(beta1, (frames, bins))
    )
    numerator_loo = coincidence_loo - c_m1_b2_loo - c_b1_m2_loo + c_b1_b2[None, :]
    denominator_loo = _cross_corr_all(m2_loo - beta2[None, :], m1_loo - beta1[None, :])
    with np.errstate(divide="ignore", invalid="ignore"):
        g2_loo = numerator_loo / denominator_loo
        # Extreme tau can leave some leave-one-out denominators at ~0 (inf/nan
        # g2_loo); those tau are masked below via the full-sample denominator.
        sem = np.sqrt(
            (frames - 1.0) / frames
            * np.sum((g2_loo - g2_loo.mean(axis=0)) ** 2, axis=0)
        )

    tau = np.arange(-(bins - 1), bins, dtype=np.int64)
    contributing = bins - np.abs(tau)
    g2 = np.where(denominator > 0.0, g2, np.nan)
    sem = np.where(denominator > 0.0, sem, np.nan)
    return tau, g2, sem, contributing


def pair_counts(
    cube: np.ndarray, beta: np.ndarray, pair: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(N1, N2, beta1, beta2) for a named pair ('AA', 'BB' or 'AB')."""
    if pair in SAME_MESA_PAIRS:
        first, second = SAME_MESA_PAIRS[pair]
        return cube[:, first, :], cube[:, second, :], beta[first], beta[second]
    a1, a2, b1, b2 = CROSS_MESA_PAIRS[pair]
    return (
        cube[:, a1, :] + cube[:, a2, :],
        cube[:, b1, :] + cube[:, b2, :],
        beta[a1] + beta[a2],
        beta[b1] + beta[b2],
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_results_h5(output_path: Path, per_width: dict, attributes: dict) -> Path:
    """Write the g2 results tree, refusing to overwrite; atomic via .partial."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing g2 file: {output_path}")
    partial_path = output_path.with_name(f"{output_path.name}.partial")
    if partial_path.exists():
        raise FileExistsError(
            f"Interrupted computation artifact already exists: {partial_path}"
        )
    try:
        with h5py.File(partial_path, "x") as output:
            for key, value in attributes.items():
                output.attrs[key] = value
            for width, results in per_width.items():
                group = output.create_group(f"bin_px_{width}")
                group.attrs["bin_width_px"] = int(width)
                group.attrs["bin_width_ps"] = float(results["bin_width_ps"])
                group.attrs["bin_count"] = int(results["nbins"])
                group.create_dataset("t_ps", data=results["t_ps"])
                group.create_dataset("mean_counts", data=results["mean_counts"])
                group.create_dataset("beta", data=results["beta"])
                group.create_dataset("intensity_A", data=results["intensity_A"])
                group.create_dataset("intensity_B", data=results["intensity_B"])
                group.create_dataset("imbalance_z", data=results["imbalance_z"])
                for pair in ("AA", "BB", "AB"):
                    entry = results["pairs"][pair]
                    for name in ("map", "map_sem"):
                        group.create_dataset(
                            f"g2_{pair}_{name}",
                            data=entry[name].astype(np.float32),
                            compression="gzip",
                        )
                    group.create_dataset(
                        f"g2_{pair}_zero_delay", data=entry["zero_delay"]
                    )
                    group.create_dataset(
                        f"g2_{pair}_zero_delay_sem", data=entry["zero_delay_sem"]
                    )
                    group.create_dataset(f"g2_{pair}_tau", data=entry["tau_g2"])
                    group.create_dataset(f"g2_{pair}_tau_sem", data=entry["tau_sem"])
                    group.attrs[f"g2_{pair}_pooled_zero_delay"] = float(
                        entry["pooled_zero_delay"]
                    )
                    group.attrs[f"g2_{pair}_pooled_zero_delay_sem"] = float(
                        entry["pooled_zero_delay_sem"]
                    )
                group.create_dataset(
                    "g2_BA_tau", data=results["pairs"]["AB"]["tau_g2"][::-1]
                )
                group.create_dataset(
                    "g2_BA_tau_sem", data=results["pairs"]["AB"]["tau_sem"][::-1]
                )
                group.create_dataset("tau_index", data=results["tau_index"])
                group.create_dataset("tau_ps", data=results["tau_ps"])
                group.create_dataset("tau_contributing_bins", data=results["tau_contrib"])
        if output_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing g2 file: {output_path}"
            )
        partial_path.rename(output_path)
    except BaseException:
        if partial_path.exists():
            partial_path.unlink()
        raise
    return output_path


def _check_figure_path(output_path: Path) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing figure: {output_path}"
        )
    return output_path


def save_intensity_figure(results: dict, title: str, output_path: Path) -> Path:
    t_ps = results["t_ps"]
    figure, (axis, z_axis) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )
    axis.plot(t_ps, results["intensity_A"], color="tab:blue", lw=0.9,
              label="condensate A (both arms, background-corrected)")
    axis.plot(t_ps, results["intensity_B"], color="tab:red", lw=0.9,
              label="condensate B (both arms, background-corrected)")
    axis.set_ylabel("counts / frame / bin")
    axis.set_title(title)
    axis.legend(fontsize=8)
    z_axis.plot(t_ps, results["imbalance_z"], color="tab:purple", lw=0.9)
    z_axis.axhline(0.0, color="gray", lw=0.7)
    # |z| > 1 is unphysical and only appears where the corrected intensity
    # denominator fluctuates to ~0; clip the display.
    z_axis.set_ylim(-1.1, 1.1)
    z_axis.set_ylabel("z = (A-B)/(A+B)")
    z_axis.set_xlabel("time [ps]")
    figure.savefig(_check_figure_path(output_path), dpi=150)
    plt.close(figure)
    return output_path


def save_zero_delay_figure(results: dict, title: str, output_path: Path) -> Path:
    t_ps = results["t_ps"]
    figure, axes = plt.subplots(
        3, 1, figsize=(9, 9), sharex=True, constrained_layout=True
    )
    for axis, pair, color in zip(
        axes, ("AA", "BB", "AB"), ("tab:blue", "tab:red", "tab:green")
    ):
        entry = results["pairs"][pair]
        g2 = entry["zero_delay"]
        sem = entry["zero_delay_sem"]
        axis.plot(t_ps, g2, color=color, lw=0.7, label=f"$g^{{(2)}}_{{{pair}}}(0)(t)$")
        axis.fill_between(t_ps, g2 - sem, g2 + sem, color=color, alpha=0.25, lw=0)
        limited = np.isfinite(sem) & (sem > MAX_DISPLAY_SEM)
        if np.any(limited):
            axis.fill_between(
                t_ps, 0.0, 1.0, where=limited,
                transform=axis.get_xaxis_transform(),
                color="gray", alpha=0.15, lw=0,
                label=f"SEM > {MAX_DISPLAY_SEM:g} (statistics-limited)",
            )
        axis.axhline(1.0, color="black", lw=0.8, ls="--")
        # Clip to the physically meaningful range: per-bin values beyond it
        # are noise artefacts from near-zero corrected intensities.
        axis.set_ylim(0.0, 3.0)
        pooled = entry["pooled_zero_delay"]
        pooled_sem = entry["pooled_zero_delay_sem"]
        axis.set_ylabel(f"$g^{{(2)}}_{{{pair}}}(0)$")
        axis.legend(fontsize=8, loc="upper right")
        axis.set_title(
            f"pooled $g^{{(2)}}_{{{pair}}}(0)$ = {pooled:.3f} +/- {pooled_sem:.3f}",
            fontsize=9,
        )
    axes[0].set_title(
        axes[0].get_title() + "\n" + title + "\ndisplay clipped to [0, 3]",
        fontsize=10,
    )
    axes[-1].set_xlabel("time [ps]")
    figure.savefig(_check_figure_path(output_path), dpi=150)
    plt.close(figure)
    return output_path


def save_map_figure(results: dict, title: str, output_path: Path) -> Path:
    t_ps = results["t_ps"]
    figure, axes = plt.subplots(
        1, 3, figsize=(16, 5.2), constrained_layout=True
    )
    for axis, pair in zip(axes, ("AA", "BB", "AB")):
        g2 = results["pairs"][pair]["map"]
        finite = np.isfinite(g2)
        span = (
            float(np.nanpercentile(np.abs(g2[finite] - 1.0), 90))
            if np.any(finite)
            else 1.0
        )
        # Cap the display span: beyond |g2-1| = 3 the map only shows noise,
        # and an unclipped percentile span would flatten all structure.
        span = min(max(span, 0.05), 3.0)
        norm = TwoSlopeNorm(vmin=1.0 - span, vcenter=1.0, vmax=1.0 + span)
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("lightgray")
        image = axis.imshow(
            np.ma.masked_invalid(g2),
            origin="lower",
            extent=[t_ps[0], t_ps[-1], t_ps[0], t_ps[-1]],
            cmap=cmap,
            norm=norm,
            aspect="equal",
        )
        axis.set_xlabel("$t_1$ [ps]")
        axis.set_ylabel("$t_2$ [ps]")
        axis.set_title(f"$g^{{(2)}}_{{{pair}}}(t_1, t_2)$")
        colorbar = figure.colorbar(image, ax=axis, shrink=0.85)
        colorbar.set_label(
            f"RdBu_r, linear around 1, range [{1.0 - span:.2f}, {1.0 + span:.2f}]",
            fontsize=8,
        )
    figure.suptitle(title, fontsize=10)
    figure.savefig(_check_figure_path(output_path), dpi=150)
    plt.close(figure)
    return output_path


def save_tau_figure(results: dict, title: str, output_path: Path) -> Path:
    tau_ps = results["tau_ps"]
    # Display the central half of the tau range: towards the window edges the
    # number of contributing (t, t+tau) bin pairs collapses and the pooled
    # estimate is pure noise.
    half = 0.5 * np.max(np.abs(tau_ps))
    shown = np.abs(tau_ps) <= half
    figure, axes = plt.subplots(
        3, 1, figsize=(9, 9), sharex=True, constrained_layout=True
    )
    for axis, pair, color in zip(
        axes, ("AA", "BB", "AB"), ("tab:blue", "tab:red", "tab:green")
    ):
        entry = results["pairs"][pair]
        g2 = entry["tau_g2"]
        sem = entry["tau_sem"]
        axis.plot(tau_ps[shown], g2[shown], color=color, lw=0.8,
                  label=f"$g^{{(2)}}_{{{pair}}}(\\tau)$ (pooled over t)")
        axis.fill_between(
            tau_ps[shown], (g2 - sem)[shown], (g2 + sem)[shown],
            color=color, alpha=0.25, lw=0,
        )
        axis.axhline(1.0, color="black", lw=0.8, ls="--")
        axis.axvline(0.0, color="gray", lw=0.7)
        center_sem = entry["pooled_zero_delay_sem"]
        axis.set_ylim(1.0 - 6 * center_sem - 0.05, 1.0 + 6 * center_sem + 0.05)
        axis.set_ylabel(f"$g^{{(2)}}_{{{pair}}}(\\tau)$")
        axis.legend(fontsize=8, loc="upper right")
    axes[0].set_title(
        title + "\n$g^{(2)}_{BA}(\\tau) = g^{(2)}_{AB}(-\\tau)$; "
        "central half of the tau range, y = 1 +/- 6 x pooled SEM",
        fontsize=10,
    )
    axes[-1].set_xlabel("$\\tau = t_1 - t_2$ [ps]")
    figure.savefig(_check_figure_path(output_path), dpi=150)
    plt.close(figure)
    return output_path


def compute_width(
    frame: np.ndarray,
    cls: np.ndarray,
    y: np.ndarray,
    beta: np.ndarray,
    frame_count: int,
    height: int,
    bin_width: int,
    correct_background: bool,
) -> dict:
    """All g2 products for one time-bin width."""
    bin_index, nbins = time_bin_index(y, height, bin_width)
    cube = count_cube(frame, cls, bin_index, nbins, frame_count)
    used_beta = beta if correct_background else np.zeros_like(beta)
    ps_per_bin = WINDOW_NS * 1000.0 / height * bin_width
    t_ps = (np.arange(nbins) + 0.5) * ps_per_bin

    pairs = {}
    for pair in ("AA", "BB", "AB"):
        N1, N2, beta1, beta2 = pair_counts(cube, used_beta, pair)
        g2, sem = g2_map(N1, N2, beta1, beta2)
        tau, tau_g2, tau_sem, tau_contrib = g2_tau_profile(N1, N2, beta1, beta2)
        zero = nbins - 1  # index of tau = 0
        pairs[pair] = {
            "map": g2,
            "map_sem": sem,
            "zero_delay": np.diag(g2),
            "zero_delay_sem": np.diag(sem),
            "tau_g2": tau_g2,
            "tau_sem": tau_sem,
            "pooled_zero_delay": float(tau_g2[zero]),
            "pooled_zero_delay_sem": float(tau_sem[zero]),
        }

    mean_counts = cube.mean(axis=0)
    intensity_A = mean_counts[0] + mean_counts[2] - used_beta[0] - used_beta[2]
    intensity_B = mean_counts[1] + mean_counts[3] - used_beta[1] - used_beta[3]
    with np.errstate(divide="ignore", invalid="ignore"):
        imbalance = (intensity_A - intensity_B) / (intensity_A + intensity_B)
    return {
        "nbins": nbins,
        "bin_width_ps": ps_per_bin,
        "t_ps": t_ps,
        "mean_counts": mean_counts,
        "beta": used_beta,
        "intensity_A": intensity_A,
        "intensity_B": intensity_B,
        "imbalance_z": imbalance,
        "pairs": pairs,
        "tau_index": tau,
        "tau_ps": tau * ps_per_bin,
        "tau_contrib": tau_contrib,
    }


def print_width_summary(width: int, results: dict, beta_fraction: np.ndarray) -> None:
    print(
        f"\nbin width {width} px ({results['bin_width_ps']:.2f} ps), "
        f"{results['nbins']} bins; background fraction per class: "
        + ", ".join(f"{name} {frac:.1%}" for name, frac in zip(CLASS_NAMES, beta_fraction))
    )
    for pair in ("AA", "BB", "AB"):
        entry = results["pairs"][pair]
        print(
            f"  g2_{pair}: pooled g2(0) = {entry['pooled_zero_delay']:.4f} "
            f"+/- {entry['pooled_zero_delay_sem']:.4f} (jackknife SEM)"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute the time-resolved g2 correlation matrix "
        "(AA, BB cross-arm; AB cross-mesa) from the x-filtered HBT peak "
        "table, replicating Adiyatullin et al. 2017 Supplementary Note 1."
    )
    parser.add_argument("filtered_h5", type=Path,
                        help="x-range-filtered peak table (carries x_range_intervals)")
    parser.add_argument("full_h5", type=Path,
                        help="full peak table (frame_shape attribute + background)")
    parser.add_argument("--output-h5", type=Path, default=None)
    parser.add_argument("--figure-dir", type=Path, default=None)
    parser.add_argument("--time-bins-px", type=int, nargs="+",
                        default=list(DEFAULT_TIME_BINS_PX),
                        help="time bin widths in pixels (default: %(default)s)")
    parser.add_argument("--no-background-correction", action="store_true",
                        help="set beta = 0 (skip the uniform-background "
                        "correction estimated from the full peak table)")
    arguments = parser.parse_args()

    filtered_attrs, frame, y, x, _ = load_peaks(arguments.filtered_h5)
    full_attrs, full_frame, full_y, full_x, _ = load_peaks(arguments.full_h5)
    frame_shape = tuple(int(dim) for dim in full_attrs["frame_shape"])
    height, detector_width = frame_shape
    frame_count = int(full_attrs["frame_count"])
    intervals = np.asarray(filtered_attrs["x_range_intervals"], dtype=np.float64)
    cls = assign_classes(x, intervals)
    unassigned = int((cls < 0).sum())
    condition_name = filtered_attrs.get("condition_name", "unknown condition")
    print(
        f"{condition_name}: {len(x)} peaks, {frame_count} frames, "
        f"frame {frame_shape[0]}x{frame_shape[1]} px, "
        f"{unassigned} peaks outside the four x-ranges"
    )

    correct_background = not arguments.no_background_correction
    per_width = {}
    for width in arguments.time_bins_px:
        beta, rho_y, chi2_ndf = background_beta(
            full_x, full_y, intervals, detector_width, frame_count, width,
            int(np.ceil(height / width)), height,
        )
        results = compute_width(
            frame, cls, y, beta, frame_count, height, width, correct_background
        )
        per_width[width] = results
        beta_fraction = results["beta"].sum(axis=1) / np.clip(
            results["mean_counts"].sum(axis=1), 1e-12, None
        )
        print(
            f"background y-flatness chi2/ndf = {chi2_ndf:.2f} "
            f"({results['nbins']} bins)"
        )
        print_width_summary(width, results, beta_fraction)

    if arguments.output_h5 is not None:
        attributes = {
            "condition_name": condition_name,
            "central_dot_strength": float(filtered_attrs["central_dot_strength"]),
            "frame_count": frame_count,
            "frame_shape": np.asarray(frame_shape),
            "window_ns": WINDOW_NS,
            "time_bins_px": np.asarray(arguments.time_bins_px),
            "background_correction": correct_background,
            "background_model": (
                "uniform detection background estimated per y-bin from "
                "full-table detections outside the four x-ranges"
            ),
            "class_names": ",".join(CLASS_NAMES),
            "pair_definitions": (
                "AA = A arm1 x A arm2; BB = B arm1 x B arm2; "
                "AB = A total x B total; g2_BA(t1,t2) = g2_AB(t2,t1)"
            ),
            "x_range_intervals": intervals,
            "source_filtered_h5": Path(arguments.filtered_h5).name,
            "source_filtered_h5_sha256": _sha256(arguments.filtered_h5),
            "source_full_h5": Path(arguments.full_h5).name,
            "source_full_h5_sha256": _sha256(arguments.full_h5),
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "script": Path(__file__).name,
        }
        print("saved", write_results_h5(arguments.output_h5, per_width, attributes))

    if arguments.figure_dir is not None:
        for width, results in per_width.items():
            stem = f"26-07-31_{condition_name}_g2"
            title = (
                f"{condition_name} — bin width {width} px "
                f"({results['bin_width_ps']:.2f} ps), {frame_count} frames"
            )
            figures = [
                (save_intensity_figure, f"{stem}_intensities_bin{width}px.png"),
                (save_zero_delay_figure, f"{stem}_zero_delay_bin{width}px.png"),
                (save_map_figure, f"{stem}_map_bin{width}px.png"),
                (save_tau_figure, f"{stem}_tau_bin{width}px.png"),
            ]
            for saver, name in figures:
                print("saved", saver(results, title, arguments.figure_dir / name))
