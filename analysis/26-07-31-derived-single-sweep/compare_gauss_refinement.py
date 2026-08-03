"""Compare Crocker-Grier centroid positions with 2D Gaussian least-squares fits.

Audit script for the peak-position estimator choice. For a reproducible random
sample of frames from one condition, it locates peaks with the same
``tp.locate`` call used by ``find_frame_peaks.py`` (Crocker-Grier bandpass +
iterative centroid) and re-fits every peak with ``tp.refine_leastsq``
(circular 2D Gaussian + constant background, fitted on the raw frame). It
writes one CSV row per peak, prints the summary statistics quoted in
``docs/peak-localization.md``, and optionally renders a summary figure and a
zoom grid of the worst outliers.

Reproduce the documented 50-frame, 1518-peak comparison with:

    python compare_gauss_refinement.py <h5> <csv> \\
        --central-dot-strength 0.18 --n-frames 50 --seed 20260803 \\
        --figure <summary.png> --outliers-figure <outliers.png>
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
import trackpy as tp

# Same detector-calibrated defaults as the production peak finder, imported so
# this audit cannot silently drift away from the parameters it validates.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_frame_peaks import (  # noqa: E402
    DEFAULT_DIAMETER,
    DEFAULT_MINMASS,
    DEFAULT_PERCENTILE,
    load_condition_frame,
)


def condition_frame_count(input_path: Path, central_dot_strength: float) -> int:
    """Frame count of the condition selected exactly as load_condition_frame."""
    with h5py.File(input_path, "r") as handle:
        conditions = handle["conditions"]
        matches = [
            name
            for name, condition in conditions.items()
            if np.isclose(
                float(condition.attrs["central_dot_strength"]),
                central_dot_strength,
                rtol=0.0,
                atol=1e-12,
            )
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one condition with central_dot_strength="
                f"{central_dot_strength}, found {len(matches)}"
            )
        return int(conditions[matches[0]]["frames"].shape[0])


def select_frame_indices(n_total: int, n_frames: int, seed: int) -> list[int]:
    """Sorted reproducible uniform sample without replacement."""
    rng = np.random.default_rng(seed)
    return sorted(rng.choice(n_total, n_frames, replace=False).tolist())


def compare_frame(frame: np.ndarray, diameter: int, minmass: float,
                  percentile: float) -> tuple[np.ndarray, float]:
    """Return per-peak rows for one frame and the number of failed fits.

    Row columns: y_centroid, x_centroid, y_gauss, x_gauss, dy, dx, mass, ecc.
    """
    frame = np.asarray(frame, dtype=np.float64)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="No maxima survived.*", category=UserWarning
        )
        centroid = tp.locate(
            frame, diameter=diameter, minmass=minmass, percentile=percentile
        )
    gauss = tp.refine_leastsq(centroid, frame, diameter=diameter)
    if len(gauss) != len(centroid):
        raise RuntimeError("refine_leastsq changed the number of features")
    rows = np.column_stack(
        [
            centroid["y"],
            centroid["x"],
            gauss["y"],
            gauss["x"],
            gauss["y"] - centroid["y"],
            gauss["x"] - centroid["x"],
            centroid["mass"],
            centroid["ecc"],
        ]
    )
    # Failed fits keep the initial coordinates and get cost = NaN.
    failed = float(np.isnan(gauss["cost"].to_numpy()).sum())
    return rows, failed


def nearest_neighbor_distances(frames: np.ndarray, y: np.ndarray,
                               x: np.ndarray) -> np.ndarray:
    """Distance to the nearest other peak within the same frame."""
    from scipy.spatial import cKDTree

    nn = np.full(len(frames), np.inf)
    for f in np.unique(frames):
        sel = np.where(frames == f)[0]
        if len(sel) > 1:
            pts = np.column_stack([y[sel], x[sel]])
            dist, _ = cKDTree(pts).query(pts, k=2)
            nn[sel] = dist[:, 1]
    return nn


def print_stats(frames: np.ndarray, rows: np.ndarray, failed: float,
                radius: int) -> None:
    """Print the audit summary. ``radius`` is the refine_leastsq position
    clamp (diameter//2; trackpy/refine/least_squares.py validate_bounds limits
    position shifts to the mask radius)."""
    dy, dx, mass, ecc = rows[:, 4], rows[:, 5], rows[:, 6], rows[:, 7]
    dist = np.hypot(dy, dx)
    at_bound = (np.abs(np.abs(dy) - radius) < 1e-6) & (
        np.abs(np.abs(dx) - radius) < 1e-6
    )
    near_corner = (np.abs(dy) > radius - 0.5) & (np.abs(dx) > radius - 0.5)
    bulk = dist <= 0.2
    out1 = dist > 1.0
    nn = nearest_neighbor_distances(frames, rows[:, 0], rows[:, 1])

    print(f"\npeaks compared: {len(rows)} (failed fits: {int(failed)})")
    print(
        f"distance [px]: median {np.median(dist):.4f}  mean {dist.mean():.4f}  "
        f"p90 {np.percentile(dist, 90):.4f}  p99 {np.percentile(dist, 99):.4f}"
        f"  max {dist.max():.4f}"
    )
    print(f"systematic: mean dy {dy.mean():+.4f} px, mean dx {dx.mean():+.4f} px")
    print(
        f"share > 0.1 px: {np.mean(dist > 0.1) * 100:.2f}%   "
        f"> 0.2 px: {np.mean(dist > 0.2) * 100:.2f}%   "
        f"> 1 px: {np.mean(out1) * 100:.2f}%"
    )
    print(
        f"bulk (dist<=0.2 px, {bulk.sum()} peaks): median "
        f"{np.median(dist[bulk]):.4f}  p90 {np.percentile(dist[bulk], 90):.4f}"
        f"  max {dist[bulk].max():.4f}  mean dy {dy[bulk].mean():+.4f}  "
        f"mean dx {dx[bulk].mean():+.4f}"
    )
    print(
        f"at clamp boundary (|dy|=|dx|={radius}): {at_bound.sum()} peaks; "
        f"near-corner (|dy| and |dx| > {radius - 0.5}) among dist>1 px: "
        f"{(out1 & near_corner).sum()} of {out1.sum()}"
    )
    print(
        f"ecc median: bulk {np.median(ecc[bulk]):.3f}  outliers(>0.2 px) "
        f"{np.median(ecc[~bulk]):.3f}"
    )
    print(
        f"mass median: bulk {np.median(mass[bulk]):.0f}  outliers(>0.2 px) "
        f"{np.median(mass[~bulk]):.0f}"
    )
    print(
        f"nearest-neighbour median: bulk {np.median(nn[bulk]):.1f} px  "
        f"outliers(>1 px) {np.median(nn[out1]):.1f} px;  with neighbour < "
        f"{2 * radius + 2} px: bulk {(nn[bulk] < 2 * radius + 2).sum()}, "
        f"outliers {(nn[out1] < 2 * radius + 2).sum()}"
    )


def save_summary_figure(frames: np.ndarray, rows: np.ndarray,
                        input_path: Path, central_dot_strength: float,
                        condition_name: str, n_frames: int,
                        output_path: Path) -> None:
    dy, dx, mass = rows[:, 4], rows[:, 5], rows[:, 6]
    dist = np.hypot(dy, dx)
    i_worst = int(np.argmax(dist))
    w = rows[i_worst]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)

    ax = axes[0, 0]
    ax.scatter(dx, dy, s=4, alpha=0.4)
    ax.add_patch(
        plt.Circle((0, 0), 0.1, fill=False, color="red", lw=1, label="0.1 px")
    )
    ax.set_xlabel("dx = x_gauss - x_centroid [pixel]")
    ax.set_ylabel("dy = y_gauss - y_centroid [pixel]")
    ax.set_title("2D offset (all peaks)")
    ax.legend(loc="upper right")
    ax.set_aspect("equal")
    lim = max(0.25, np.abs(np.c_[dx, dy]).max() * 1.1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    ax = axes[0, 1]
    ax.hist(dist, bins=np.linspace(0, max(0.5, dist.max() * 1.05), 60))
    ax.set_yscale("log")
    ax.set_xlabel("|gauss - centroid| [pixel]")
    ax.set_ylabel("peaks")
    ax.set_title(f"offset distance (median {np.median(dist):.3f} px)")

    ax = axes[1, 0]
    ax.scatter(mass, dist, s=4, alpha=0.4)
    ax.set_xscale("log")
    ax.set_xlabel("mass (bandpassed integrated brightness)")
    ax.set_ylabel("|gauss - centroid| [pixel]")
    ax.set_title("offset vs brightness")

    ax = axes[1, 1]
    _, frame = load_condition_frame(
        input_path, central_dot_strength, int(frames[i_worst])
    )
    yc, xc = w[0], w[1]
    r = 8
    y0, y1 = int(yc) - r, int(yc) + r + 1
    x0, x1 = int(xc) - r, int(xc) + r + 1
    ax.imshow(
        frame[y0:y1, x0:x1],
        origin="upper",
        cmap="magma",
        interpolation="nearest",
        extent=[x0 - 0.5, x1 - 0.5, y1 - 0.5, y0 - 0.5],
    )
    ax.plot(w[1], w[0], "c+", ms=11, mew=1.4,
            label=f"centroid ({w[0]:.2f}, {w[1]:.2f})")
    ax.plot(w[3], w[2], "rx", ms=9, mew=1.4,
            label=f"gauss ({w[2]:.2f}, {w[3]:.2f})")
    ax.set_xticks(np.arange(x0, x1, 2))
    ax.set_yticks(np.arange(y0, y1, 2))
    ax.grid(color="gray", lw=0.3, alpha=0.5)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"worst peak: frame {int(frames[i_worst])}, "
                 f"offset {dist[i_worst]:.2f} px")

    fig.suptitle(
        f"{condition_name}, {n_frames} sampled frames: 2D Gaussian fit vs "
        "Crocker-Grier centroid"
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print("saved", output_path)


def save_outliers_figure(frames: np.ndarray, rows: np.ndarray,
                         input_path: Path, central_dot_strength: float,
                         condition_name: str, output_path: Path,
                         n_show: int = 6) -> None:
    dist = np.hypot(rows[:, 4], rows[:, 5])
    worst = np.argsort(dist)[-n_show:][::-1]
    ncols = 3
    nrows = int(np.ceil(n_show / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.5 * nrows),
                             constrained_layout=True)
    cache: dict[int, np.ndarray] = {}
    for ax, i in zip(np.atleast_1d(axes.flat), worst):
        f = int(frames[i])
        if f not in cache:
            cache[f] = load_condition_frame(
                input_path, central_dot_strength, f
            )[1]
        frame = cache[f]
        yc, xc, yg, xg = rows[i, 0], rows[i, 1], rows[i, 2], rows[i, 3]
        r = 8
        y0, y1 = int(yc) - r, int(yc) + r + 1
        x0, x1 = int(xc) - r, int(xc) + r + 1
        ax.imshow(
            frame[y0:y1, x0:x1],
            origin="upper",
            cmap="magma",
            interpolation="nearest",
            extent=[x0 - 0.5, x1 - 0.5, y1 - 0.5, y0 - 0.5],
        )
        ax.plot(xc, yc, "c+", ms=11, mew=1.4, label="centroid")
        ax.plot(xg, yg, "rx", ms=9, mew=1.4, label="gauss fit")
        ax.set_xticks(np.arange(x0, x1, 4))
        ax.set_yticks(np.arange(y0, y1, 4))
        ax.grid(color="gray", lw=0.3, alpha=0.5)
        ax.set_title(
            f"frame {f}  |Δ|={dist[i]:.2f} px  mass={rows[i, 6]:.0f}",
            fontsize=9,
        )
        ax.legend(loc="upper right", fontsize=7)
    fig.suptitle(
        f"{condition_name} — the {n_show} largest centroid-vs-Gaussian "
        "offsets (raw frames)"
    )
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print("saved", output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare Crocker-Grier centroid positions with 2D "
        "Gaussian least-squares fits on a reproducible frame sample."
    )
    parser.add_argument("input_h5", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--central-dot-strength", type=float, required=True)
    parser.add_argument("--n-frames", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260803,
                        help="RNG seed of the documented 50-frame sample")
    parser.add_argument("--diameter", type=int, default=DEFAULT_DIAMETER)
    parser.add_argument("--minmass", type=float, default=DEFAULT_MINMASS)
    parser.add_argument("--percentile", type=float, default=DEFAULT_PERCENTILE)
    parser.add_argument("--figure", type=Path, default=None,
                        help="optional 4-panel summary figure")
    parser.add_argument("--outliers-figure", type=Path, default=None,
                        help="optional zoom grid of the largest offsets")
    arguments = parser.parse_args()

    n_total = condition_frame_count(
        arguments.input_h5, arguments.central_dot_strength
    )
    indices = select_frame_indices(n_total, arguments.n_frames, arguments.seed)
    print(f"selected {len(indices)} of {n_total} frames (seed "
          f"{arguments.seed}): {' '.join(map(str, indices))}", flush=True)

    all_frames = []
    all_rows = []
    failed = 0.0
    t0 = time.time()
    for n, idx in enumerate(indices):
        condition_name, frame = load_condition_frame(
            arguments.input_h5, arguments.central_dot_strength, idx
        )
        rows, failed_frame = compare_frame(
            frame, arguments.diameter, arguments.minmass, arguments.percentile
        )
        all_frames.append(np.full(len(rows), idx))
        all_rows.append(rows)
        failed += failed_frame
        if (n + 1) % 10 == 0:
            print(f"{n + 1}/{len(indices)} frames, {time.time() - t0:.0f}s",
                  flush=True)

    frame_col = np.concatenate(all_frames)
    rows = np.concatenate(all_rows)
    radius = arguments.diameter // 2

    arguments.output_csv.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        arguments.output_csv,
        np.column_stack([frame_col, rows]),
        delimiter=",",
        header="frame,y_centroid,x_centroid,y_gauss,x_gauss,dy,dx,mass,ecc",
        comments="",
        fmt="%.4f",
    )
    print("saved", arguments.output_csv)
    print_stats(frame_col, rows, failed, radius)

    if arguments.figure is not None:
        save_summary_figure(
            frame_col,
            rows,
            arguments.input_h5,
            arguments.central_dot_strength,
            condition_name,
            len(indices),
            arguments.figure,
        )
    if arguments.outliers_figure is not None:
        save_outliers_figure(
            frame_col,
            rows,
            arguments.input_h5,
            arguments.central_dot_strength,
            condition_name,
            arguments.outliers_figure,
        )
