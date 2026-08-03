"""Fit the summed peak x-profile with four Gaussians and derive x-acceptance ranges.

Loads the per-frame peak table produced by ``extract_condition_peaks.py``,
sums all detected peaks onto detector axis 2 (1-pixel histogram, i.e. the
vertical axis is compressed away), and fits the profile with a constant
background plus four Gaussians ordered left to right:

    spot 1: condensate A, arm A
    spot 2: condensate B, arm A
    spot 3: condensate A, arm B
    spot 4: condensate B, arm B

Only peaks whose x falls inside one of the four fitted spot ranges can be
real spot data. The criterion is mu_i +/- k*sigma_i, starting from the
conservative k = 3 (about 99.7% containment per Gaussian). If the resulting
intervals of adjacent spots overlap, k is automatically reduced in steps
until the four intervals are disjoint (see disjoint_intervals); the
acceptance region is the union of the four disjoint intervals. The script
writes a diagnostic figure, a fit CSV, and the x-filtered peak table for
downstream analysis.

Example:

    python fit_spot_x_ranges.py peaks.h5 --figure fit.png --fit-csv fit.csv \
        --filtered-h5 peaks_xfiltered.h5 --sigma-multiple 3.0
"""
import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

matplotlib.use("Agg")
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract_condition_peaks import write_peaks_h5  # noqa: E402

N_SPOTS = 4
# Left-to-right beam-splitter assignment of the four spots.
SPOT_LABELS = (
    "condensate A arm A",
    "condensate B arm A",
    "condensate A arm B",
    "condensate B arm B",
)

DEFAULT_SIGMA_MULTIPLE = 3.0
DEFAULT_SMOOTH_SIGMA = 3.0
DEFAULT_MIN_SEPARATION = 30
SIGMA_MULTIPLE_STEP = 0.5
SIGMA_MULTIPLE_MIN = 1.0


def load_peaks(peaks_path: Path) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (attributes, frame, y, x, mass) from an extraction HDF5 file."""
    with h5py.File(peaks_path, "r") as handle:
        attributes = dict(handle.attrs)
        frame = np.asarray(handle["frame"])
        y = np.asarray(handle["y"])
        x = np.asarray(handle["x"])
        mass = np.asarray(handle["mass"])
    return attributes, frame, y, x, mass


def x_histogram(x: np.ndarray, width: int) -> tuple[np.ndarray, np.ndarray]:
    """1-pixel-bin histogram of subpixel x on the detector grid [0, width)."""
    centers = np.arange(width, dtype=np.float64)
    counts, _ = np.histogram(x, bins=np.arange(width + 1, dtype=np.float64) - 0.5)
    return centers, counts.astype(np.float64)


def initial_centers(
    smoothed: np.ndarray,
    n_spots: int = N_SPOTS,
    min_separation: int = DEFAULT_MIN_SEPARATION,
) -> list[int]:
    """Bin indices of the n strongest local maxima, greedily separated."""
    available = np.ones(len(smoothed), dtype=bool)
    # Local maxima only: strictly above both neighbours (edges count too).
    is_maximum = np.zeros(len(smoothed), dtype=bool)
    is_maximum[1:-1] = (smoothed[1:-1] > smoothed[:-2]) & (
        smoothed[1:-1] >= smoothed[2:]
    )
    is_maximum[0] = smoothed[0] >= smoothed[1] if len(smoothed) > 1 else True
    is_maximum[-1] = smoothed[-1] > smoothed[-2] if len(smoothed) > 1 else True

    centers = []
    for _ in range(n_spots):
        candidate = np.where(available & is_maximum, smoothed, -np.inf)
        index = int(np.argmax(candidate))
        if not np.isfinite(candidate[index]):
            raise ValueError(
                f"Could not find {n_spots} separated local maxima in the "
                "x-profile; check min_separation and the input data"
            )
        centers.append(index)
        lo = max(0, index - min_separation)
        hi = min(len(smoothed), index + min_separation + 1)
        available[lo:hi] = False
    return sorted(centers)


def four_gaussians(
    x: np.ndarray,
    background: float,
    a1: float, mu1: float, sigma1: float,
    a2: float, mu2: float, sigma2: float,
    a3: float, mu3: float, sigma3: float,
    a4: float, mu4: float, sigma4: float,
) -> np.ndarray:
    """Constant background plus four Gaussians (peak-height amplitudes)."""
    result = np.full(np.shape(x), background, dtype=np.float64)
    for amplitude, mu, sigma in (
        (a1, mu1, sigma1),
        (a2, mu2, sigma2),
        (a3, mu3, sigma3),
        (a4, mu4, sigma4),
    ):
        result = result + amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    return result


def fit_profile(
    centers: np.ndarray,
    counts: np.ndarray,
    min_separation: int = DEFAULT_MIN_SEPARATION,
    smooth_sigma: float = DEFAULT_SMOOTH_SIGMA,
) -> tuple[np.ndarray, np.ndarray, list[int]]:
    """Fit background + four Gaussians; return (popt, pcov, initial centers).

    The four mu bounds are the midpoints between neighbouring initial
    centers, which enforces the left-to-right spot order.
    """
    smoothed = gaussian_filter1d(counts, smooth_sigma)
    mu0 = initial_centers(smoothed, N_SPOTS, min_separation)
    background0 = float(np.median(counts))
    width = len(counts)

    p0 = [background0]
    lower = [0.0]
    upper = [np.inf]
    for i, center in enumerate(mu0):
        amplitude0 = max(smoothed[center] - background0, 1.0)
        mu_lo = 0.0 if i == 0 else 0.5 * (mu0[i - 1] + mu0[i])
        mu_hi = float(width - 1) if i == N_SPOTS - 1 else 0.5 * (mu0[i] + mu0[i + 1])
        p0.extend([amplitude0, float(center), 12.0])
        lower.extend([0.0, mu_lo, 2.0])
        upper.extend([np.inf, mu_hi, 80.0])

    weights = np.sqrt(np.clip(counts, 1.0, None))
    popt, pcov = curve_fit(
        four_gaussians,
        centers,
        counts,
        p0=p0,
        bounds=(lower, upper),
        sigma=weights,
        absolute_sigma=False,
        maxfev=20000,
    )
    return popt, pcov, mu0


def unpack_components(popt: np.ndarray, pcov: np.ndarray) -> tuple[dict, list[dict]]:
    """Split the fit vector into a background dict and one dict per spot."""
    stderr = np.sqrt(np.diag(pcov))
    background = {"value": float(popt[0]), "std": float(stderr[0])}
    spots = []
    for i in range(N_SPOTS):
        base = 1 + 3 * i
        spots.append(
            {
                "label": SPOT_LABELS[i],
                "amplitude": float(popt[base]),
                "amplitude_std": float(stderr[base]),
                "mu": float(popt[base + 1]),
                "mu_std": float(stderr[base + 1]),
                "sigma": float(popt[base + 2]),
                "sigma_std": float(stderr[base + 2]),
            }
        )
    return background, spots


def spot_intervals(popt: np.ndarray, sigma_multiple: float) -> np.ndarray:
    """(N_SPOTS, 2) acceptance intervals [mu - k*sigma, mu + k*sigma]."""
    intervals = []
    for i in range(N_SPOTS):
        base = 1 + 3 * i
        mu, sigma = float(popt[base + 1]), float(popt[base + 2])
        intervals.append([mu - sigma_multiple * sigma, mu + sigma_multiple * sigma])
    return np.array(intervals, dtype=np.float64)


def intervals_are_disjoint(intervals: np.ndarray) -> bool:
    """True when no two intervals share any x (touching edges allowed).

    The spots are ordered left to right, so only adjacent pairs can overlap.
    """
    return bool(np.all(intervals[:-1, 1] <= intervals[1:, 0]))


def disjoint_intervals(
    popt: np.ndarray,
    sigma_multiple: float,
    step: float = SIGMA_MULTIPLE_STEP,
    minimum: float = SIGMA_MULTIPLE_MIN,
) -> tuple[float, np.ndarray]:
    """Largest k <= sigma_multiple (in ``step`` decrements) with disjoint intervals.

    Open problem (not yet solved): when the mu +/- k*sigma intervals of two
    adjacent spots overlap, the shared x-region contains light from both
    spots and we currently have no rule for assigning a peak there to one
    spot or the other. Until such a separation rule exists, we avoid the
    problem by shrinking k below the overlap threshold and accept losing the
    outer Gaussian tails. If even ``minimum`` still overlaps (pathological
    fit), that multiple is returned unchanged and the caller must warn.
    """
    k = float(sigma_multiple)
    while k > minimum:
        intervals = spot_intervals(popt, k)
        if intervals_are_disjoint(intervals):
            return k, intervals
        k = round(k - step, 10)
    return minimum, spot_intervals(popt, minimum)


def acceptance_mask(x: np.ndarray, intervals: np.ndarray) -> np.ndarray:
    """Boolean mask: x inside any of the acceptance intervals (union)."""
    mask = np.zeros(len(x), dtype=bool)
    for lo, hi in intervals:
        mask |= (x >= lo) & (x <= hi)
    return mask


def write_fit_csv(
    output_path: Path,
    background: dict,
    spots: list[dict],
    intervals: np.ndarray,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing fit CSV: {output_path}"
        )
    header = (
        "component,label,amplitude,amplitude_std,mu,mu_std,"
        "sigma,sigma_std,range_lo,range_hi"
    )
    lines = [header]
    lines.append(
        f"background,,{background['value']:.4f},{background['std']:.4f},,,,,,"
    )
    for i, spot in enumerate(spots):
        lines.append(
            f"spot{i + 1},{spot['label']},"
            f"{spot['amplitude']:.4f},{spot['amplitude_std']:.4f},"
            f"{spot['mu']:.4f},{spot['mu_std']:.4f},"
            f"{spot['sigma']:.4f},{spot['sigma_std']:.4f},"
            f"{intervals[i, 0]:.4f},{intervals[i, 1]:.4f}"
        )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def save_fit_figure(
    centers: np.ndarray,
    counts: np.ndarray,
    popt: np.ndarray,
    pcov: np.ndarray,
    intervals: np.ndarray,
    sigma_multiple: float,
    title: str,
    output_path: Path,
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing figure: {output_path}"
        )
    background, spots = unpack_components(popt, pcov)
    model = four_gaussians(centers, *popt)

    figure, (axis, residual_axis) = plt.subplots(
        2, 1, figsize=(11, 7.5), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )
    axis.step(centers, counts, where="mid", color="black", lw=0.8,
              label="peak x histogram (1 px bins)")
    axis.plot(centers, model, color="red", lw=1.6, label="background + 4 Gaussians")
    axis.axhline(
        background["value"], color="gray", ls="--", lw=0.9,
        label=f"background = {background['value']:.1f}",
    )
    colors = ("tab:blue", "tab:orange", "tab:green", "tab:purple")
    for i, spot in enumerate(spots):
        component = spot["amplitude"] * np.exp(
            -0.5 * ((centers - spot["mu"]) / spot["sigma"]) ** 2
        )
        axis.plot(centers, background["value"] + component, color=colors[i],
                  lw=1.0, alpha=0.85,
                  label=(
                      f"spot{i + 1} ({spot['label']}): "
                      f"$\\mu$={spot['mu']:.1f}, $\\sigma$={spot['sigma']:.1f}"
                  ))
        axis.axvspan(
            intervals[i, 0], intervals[i, 1], color=colors[i], alpha=0.10,
        )
    axis.set_ylabel("Detected peaks per 1-px bin")
    axis.set_title(title)
    axis.legend(loc="upper right", fontsize=8)

    residuals = counts - model
    residual_axis.step(centers, residuals, where="mid", color="black", lw=0.7)
    residual_axis.axhline(0.0, color="red", lw=0.9)
    for lo, hi in intervals:
        residual_axis.axvspan(lo, hi, color="gray", alpha=0.08)
    residual_axis.set_xlabel("Detector axis 2 [pixel]")
    residual_axis.set_ylabel("residual")

    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


def print_summary(
    background: dict,
    spots: list[dict],
    intervals: np.ndarray,
    x: np.ndarray,
    sigma_multiple: float,
) -> None:
    print(f"\nbackground: {background['value']:.2f} +/- {background['std']:.2f}")
    for i, spot in enumerate(spots):
        inside = int(((x >= intervals[i, 0]) & (x <= intervals[i, 1])).sum())
        print(
            f"spot{i + 1} ({spot['label']}): "
            f"amplitude {spot['amplitude']:.1f} +/- {spot['amplitude_std']:.1f}, "
            f"mu {spot['mu']:.2f} +/- {spot['mu_std']:.2f}, "
            f"sigma {spot['sigma']:.2f} +/- {spot['sigma_std']:.2f}, "
            f"range [{intervals[i, 0]:.1f}, {intervals[i, 1]:.1f}] "
            f"({inside} peaks)"
        )
    accepted = int(acceptance_mask(x, intervals).sum())
    print(
        f"acceptance: union of {N_SPOTS} mu +/- {sigma_multiple:g}*sigma "
        f"intervals keeps {accepted} of {len(x)} peaks "
        f"({accepted / len(x) * 100:.1f}%)"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Fit the summed peak x-profile with background + four "
        "Gaussians and derive the x-acceptance ranges (union of "
        "mu +/- k*sigma intervals)."
    )
    parser.add_argument("peaks_h5", type=Path)
    parser.add_argument("--figure", type=Path, default=None)
    parser.add_argument("--fit-csv", type=Path, default=None)
    parser.add_argument("--filtered-h5", type=Path, default=None)
    parser.add_argument(
        "--sigma-multiple", type=float, default=DEFAULT_SIGMA_MULTIPLE,
        help="half-width of each spot's acceptance interval in units of its "
        "fitted sigma (default: %(default)s, conservative). If adjacent "
        "intervals overlap, the multiple is reduced in steps of "
        f"{SIGMA_MULTIPLE_STEP:g} until they are disjoint; see "
        "--allow-overlap.",
    )
    parser.add_argument(
        "--allow-overlap", action="store_true",
        help="keep the requested --sigma-multiple even when adjacent "
        "intervals overlap (overlap assignment is an unsolved problem; the "
        "default is to shrink the multiple until intervals are disjoint)",
    )
    parser.add_argument("--smooth-sigma", type=float, default=DEFAULT_SMOOTH_SIGMA)
    parser.add_argument("--min-separation", type=int, default=DEFAULT_MIN_SEPARATION)
    arguments = parser.parse_args()

    attributes, frame, y, x, mass = load_peaks(arguments.peaks_h5)
    width = int(np.asarray(attributes["frame_shape"])[1])
    centers, counts = x_histogram(x, width)
    popt, pcov, mu0 = fit_profile(
        centers,
        counts,
        min_separation=arguments.min_separation,
        smooth_sigma=arguments.smooth_sigma,
    )
    background, spots = unpack_components(popt, pcov)
    requested_multiple = arguments.sigma_multiple
    if arguments.allow_overlap:
        used_multiple = requested_multiple
        intervals = spot_intervals(popt, used_multiple)
    else:
        used_multiple, intervals = disjoint_intervals(popt, requested_multiple)

    condition_name = attributes.get("condition_name", "unknown condition")
    print(f"{condition_name}: {len(x)} peaks, initial centers {mu0}")
    if used_multiple != requested_multiple:
        print(
            f"requested mu +/- {requested_multiple:g}*sigma intervals "
            f"overlap; reduced to {used_multiple:g}*sigma so the four "
            "acceptance intervals are disjoint"
        )
    print_summary(background, spots, intervals, x, used_multiple)

    if arguments.fit_csv is not None:
        print("saved", write_fit_csv(arguments.fit_csv, background, spots, intervals))
    if arguments.figure is not None:
        title = (
            f"{condition_name} — summed peak x-profile, four-Gaussian fit\n"
            f"acceptance = mu +/- {used_multiple:g}*sigma"
        )
        if used_multiple != requested_multiple:
            title += (
                f" (reduced from {requested_multiple:g}*sigma to avoid overlap)"
            )
        print(
            "saved",
            save_fit_figure(
                centers, counts, popt, pcov, intervals,
                used_multiple, title, arguments.figure,
            ),
        )
    if arguments.filtered_h5 is not None:
        mask = acceptance_mask(x, intervals)
        rows = np.column_stack([frame[mask], y[mask], x[mask], mass[mask]])
        filtered_attributes = {
            "condition_name": condition_name,
            "central_dot_strength": float(attributes["central_dot_strength"]),
            "source_peaks_h5": Path(arguments.peaks_h5).name,
            "x_range_filter": (
                f"union of {N_SPOTS} mu +/- {used_multiple:g}*sigma "
                "intervals from the four-Gaussian x-profile fit"
            ),
            "sigma_multiple_requested": requested_multiple,
            "sigma_multiple_used": used_multiple,
            "x_range_intervals": np.asarray(intervals, dtype=np.float64),
            "spot_labels": ",".join(SPOT_LABELS),
            "peak_count": int(mask.sum()),
            "source_peak_count": len(x),
            "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "script": Path(__file__).name,
        }
        print(
            "saved",
            write_peaks_h5(arguments.filtered_h5, rows, filtered_attributes),
        )
