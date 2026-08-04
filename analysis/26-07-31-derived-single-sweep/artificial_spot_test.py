"""Artificial-spot (ADDSTAR-style) localization test on real frames.

Injects synthetic spots with known subpixel positions into real frames of one
condition, re-runs the production Crocker-Grier peak finder, and compares
recovered positions against ground truth. This is the astronomy community's
standard artificial-star test (Stetson's ADDSTAR for DAOPHOT) adapted to
single-sweep streak frames: injection happens into REAL frames so the
background, read noise, impulse noise and crowding are the measured ones, and
injected spot counts are kept small relative to the real peak population
(~30/frame) so the frames' character does not change.

Spot model, matched to the measured 0.18-condition population:

- Elliptical 2D Gaussian with sigma ~ U(2.4, 3.0) px (measured FWHM 6-7 px),
  eccentricity ~ U(0, 0.6) (measured median ecc ~0.26 with an elongated tail)
  and uniform random orientation.
- Target bandpassed mass log-uniform in --mass-min..--mass-max, spanning the
  detection threshold (minmass 300) up to the brightest observed peaks
  (~22000). The raw amplitude that yields the target bandpassed mass is solved
  exactly per spot by computing the bandpassed mass of a unit-integral profile
  (no noise), so the target is accurate by construction.
- Detector noise: pixel values are not photon counts (background 105 counts
  with read noise sigma=3, i.e. variance/mean = 0.086 instead of Poisson 1).
  Injected spots are therefore sampled as Poisson(counts * GAIN)/GAIN with
  GAIN = 105/9 = 11.7, which gives the injected photons the same
  variance-per-count as the measured background instead of overstating their
  shot noise as raw Poisson sampling would.

Outputs one CSV row per injection (truth + matched detection + Savin-Doyle ep
of the detection) and a four-panel summary figure: (a) position-error scatter,
(b) RMS radial error vs detected mass with the median ep curve overlaid (the
cross-check of the theoretical error model against Monte-Carlo truth),
(c) detection completeness vs target mass, (d) radial-error distribution.

Example:

    python artificial_spot_test.py frames.h5 spots.csv --central-dot-strength 0.18
"""
import argparse
import sys
import time
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_frame_peaks import (  # noqa: E402
    DEFAULT_DIAMETER,
    DEFAULT_MINMASS,
    DEFAULT_PERCENTILE,
    find_peaks,
)
from localization_uncertainty import (  # noqa: E402
    bandpass_image,
    frame_background_noise,
    static_ep,
)
from extract_condition_peaks import match_condition_name  # noqa: E402

# Detector calibration, same measured values as find_frame_peaks.py:
# background 105 counts, read noise sigma=3, so variance-per-count is 9/105.
BACKGROUND = 105.0
READ_NOISE = 3.0
VARIANCE_GAIN = BACKGROUND / READ_NOISE**2

PATCH_HALF = 10  # injection patch is (2*PATCH_HALF+1)^2
MARGIN = 12  # injected centers stay this far from frame borders
MIN_SEPARATION = 30.0  # between injected spots of the same frame

CSV_HEADER = (
    "frame,y_true,x_true,mass_target,sigma,ecc,angle_deg,total_counts,"
    "matched,y_det,x_det,dy,dx,mass_det,ep_det,frame_noise"
)


def sample_frame_indices(rng, n_select, frame_count):
    """Sorted reproducible subset of frame indices, without replacement."""
    return np.sort(rng.choice(frame_count, size=n_select, replace=False))


def gaussian_profile(grid_y, grid_x, center, sigma_minor, sigma_major, angle):
    """Unit-integral elliptical Gaussian on the given pixel grid."""
    cy, cx = center
    dy = grid_y - cy
    dx = grid_x - cx
    cos_a, sin_a = np.cos(angle), np.sin(angle)
    along = dx * cos_a + dy * sin_a
    across = -dx * sin_a + dy * cos_a
    profile = np.exp(
        -0.5 * ((along / sigma_major) ** 2 + (across / sigma_minor) ** 2)
    )
    return profile / profile.sum()


def unit_bandpassed_mass(sigma_minor, sigma_major, diameter=DEFAULT_DIAMETER):
    """Bandpassed mass of a centered unit-integral spot (noiseless)."""
    side = 2 * PATCH_HALF + 1
    grid_y, grid_x = np.mgrid[0:side, 0:side]
    profile = gaussian_profile(
        grid_y, grid_x, (PATCH_HALF, PATCH_HALF), sigma_minor, sigma_major, 0.0
    )
    bandpassed = bandpass_image(profile, diameter=diameter)
    radius = diameter // 2
    mask = (grid_y - PATCH_HALF) ** 2 + (grid_x - PATCH_HALF) ** 2 <= radius**2
    return float(bandpassed[mask].sum())


def draw_injections(
    rng,
    n_spots,
    frame_shape,
    mass_range,
    sigma_range=(2.4, 3.0),
    ecc_max=0.6,
):
    """Draw positions and spot parameters; solve amplitude from target mass."""
    injections = []
    attempts = 0
    while len(injections) < n_spots and attempts < 100 * n_spots:
        attempts += 1
        y = rng.uniform(MARGIN, frame_shape[0] - MARGIN)
        x = rng.uniform(MARGIN, frame_shape[1] - MARGIN)
        if any(
            (y - other["y"]) ** 2 + (x - other["x"]) ** 2 < MIN_SEPARATION**2
            for other in injections
        ):
            continue
        mass_target = float(
            np.exp(rng.uniform(np.log(mass_range[0]), np.log(mass_range[1])))
        )
        sigma = rng.uniform(*sigma_range)
        ecc = rng.uniform(0.0, ecc_max)
        sigma_minor = sigma * np.sqrt(1.0 - ecc)
        sigma_major = sigma / np.sqrt(1.0 - ecc)
        angle = rng.uniform(0.0, np.pi)
        total_counts = mass_target / unit_bandpassed_mass(sigma_minor, sigma_major)
        injections.append(
            {
                "y": y,
                "x": x,
                "mass_target": mass_target,
                "sigma": sigma,
                "ecc": ecc,
                "angle": angle,
                "total_counts": total_counts,
            }
        )
    return injections


def inject_spots(frame, injections, rng, variance_gain=VARIANCE_GAIN):
    """Return a float copy of the frame with noisy synthetic spots added."""
    injected = np.asarray(frame, dtype=np.float64).copy()
    for spot in injections:
        iy, ix = int(round(spot["y"])), int(round(spot["x"]))
        y0, y1 = iy - PATCH_HALF, iy + PATCH_HALF + 1
        x0, x1 = ix - PATCH_HALF, ix + PATCH_HALF + 1
        grid_y, grid_x = np.mgrid[y0:y1, x0:x1]
        sigma_minor = spot["sigma"] * np.sqrt(1.0 - spot["ecc"])
        sigma_major = spot["sigma"] / np.sqrt(1.0 - spot["ecc"])
        profile = gaussian_profile(
            grid_y,
            grid_x,
            (spot["y"], spot["x"]),
            sigma_minor,
            sigma_major,
            spot["angle"],
        )
        expected = spot["total_counts"] * profile
        injected[y0:y1, x0:x1] += rng.poisson(expected * variance_gain) / variance_gain
    return injected


def match_to_truth(injections, peaks, radius=2.0):
    """Match each injected spot to the nearest detected peak within radius."""
    matched = []
    for spot in injections:
        row = dict(spot)
        if len(peaks):
            distances = np.hypot(peaks[:, 0] - spot["y"], peaks[:, 1] - spot["x"])
            nearest = int(distances.argmin())
            if distances[nearest] <= radius:
                row.update(
                    matched=True,
                    y_det=peaks[nearest, 0],
                    x_det=peaks[nearest, 1],
                    dy=peaks[nearest, 0] - spot["y"],
                    dx=peaks[nearest, 1] - spot["x"],
                    mass_det=peaks[nearest, 2],
                )
                matched.append(row)
                continue
        row.update(
            matched=False,
            y_det=np.nan,
            x_det=np.nan,
            dy=np.nan,
            dx=np.nan,
            mass_det=np.nan,
        )
        matched.append(row)
    return matched


def run_injection_test(
    frames_path,
    condition_name,
    n_frames=200,
    seed=20260804,
    spots_per_frame=4,
    mass_range=(250.0, 30000.0),
    match_radius=2.0,
    diameter=DEFAULT_DIAMETER,
    minmass=DEFAULT_MINMASS,
    percentile=DEFAULT_PERCENTILE,
    progress_every=25,
):
    """Inject spots into sampled real frames and re-run the production finder.

    Returns (frame_indices, rows) where rows is a list of per-injection dicts
    with truth, match status, detection coordinates/mass, ep and frame noise.
    """
    rng = np.random.default_rng(seed)
    rows = []
    with h5py.File(frames_path, "r") as handle:
        frames = handle["conditions"][condition_name]["frames"]
        frame_indices = sample_frame_indices(rng, n_frames, frames.shape[0])
        started = time.time()
        for count, frame_index in enumerate(frame_indices):
            frame = np.asarray(frames[frame_index])
            injections = draw_injections(rng, spots_per_frame, frame.shape, mass_range)
            injected = inject_spots(frame, injections, rng)
            peaks = find_peaks(
                injected, diameter=diameter, minmass=minmass, percentile=percentile
            )
            noise = frame_background_noise(injected, diameter=diameter)
            matched = match_to_truth(injections, peaks, radius=match_radius)
            for row in matched:
                row["frame"] = int(frame_index)
                row["frame_noise"] = noise
                row["ep_det"] = (
                    static_ep([row["mass_det"]], noise, diameter=diameter)[0]
                    if row["matched"]
                    else np.nan
                )
                rows.append(row)
            if progress_every and (count + 1) % progress_every == 0:
                print(
                    f"{count + 1}/{len(frame_indices)} frames, "
                    f"{time.time() - started:.0f}s",
                    flush=True,
                )
    return frame_indices, rows


def rows_to_array(rows):
    """Flat (N, len(CSV_HEADER)) float array in CSV column order."""
    columns = []
    for row in rows:
        columns.append(
            [
                row["frame"],
                row["y"],
                row["x"],
                row["mass_target"],
                row["sigma"],
                row["ecc"],
                np.degrees(row["angle"]),
                row["total_counts"],
                float(row["matched"]),
                row["y_det"],
                row["x_det"],
                row["dy"],
                row["dx"],
                row["mass_det"],
                row["ep_det"],
                row["frame_noise"],
            ]
        )
    return np.asarray(columns, dtype=np.float64)


def write_csv(output_path, rows):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing CSV: {output_path}"
        )
    np.savetxt(
        output_path,
        rows_to_array(rows),
        delimiter=",",
        header=CSV_HEADER,
        comments="",
        fmt="%.4f",
    )
    return output_path


def save_summary_figure(rows, output_path, minmass=DEFAULT_MINMASS):
    """Four panels: error scatter, RMS-vs-mass vs ep, completeness, error dist."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing figure: {output_path}"
        )
    table = rows_to_array(rows)
    col = {name: index for index, name in enumerate(CSV_HEADER.split(","))}
    matched = table[table[:, col["matched"]] == 1.0]

    figure, axes = plt.subplots(2, 2, figsize=(11.5, 9.0), constrained_layout=True)

    axis = axes[0, 0]
    limit = max(1.0, np.nanpercentile(np.abs(matched[:, [col["dy"], col["dx"]]]), 99))
    axis.scatter(matched[:, col["dx"]], matched[:, col["dy"]], s=4, alpha=0.35,
                 color="tab:blue", rasterized=True)
    axis.add_patch(plt.Circle((0, 0), 0.2, fill=False, color="tab:orange", lw=1))
    axis.set_xlim(-limit, limit)
    axis.set_ylim(-limit, limit)
    axis.set_aspect("equal")
    axis.set_xlabel("dx = x_det - x_true [px]")
    axis.set_ylabel("dy = y_det - y_true [px]")
    axis.set_title(f"position errors of {len(matched)} matched spots "
                   "(orange circle: 0.2 px)")

    axis = axes[0, 1]
    mass_det = matched[:, col["mass_det"]]
    radial = np.hypot(matched[:, col["dy"]], matched[:, col["dx"]])
    min_count = 5 if len(matched) >= 40 else 1
    edges = np.logspace(np.log10(mass_det.min()), np.log10(mass_det.max()), 11)
    centers, rms_radial, median_ep = [], [], []
    for lower, upper in zip(edges[:-1], edges[1:]):
        selected = (mass_det >= lower) & (mass_det < upper)
        if selected.sum() < min_count:
            continue
        centers.append(np.sqrt(lower * upper))
        rms_radial.append(np.sqrt(np.mean(radial[selected] ** 2)))
        median_ep.append(np.median(matched[selected, col["ep_det"]]))
    axis.plot(centers, rms_radial, "o-", color="tab:blue",
              label="RMS radial error (Monte Carlo)")
    axis.plot(centers, median_ep, "s--", color="tab:orange",
              label="median ep (Savin–Doyle)")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("detected bandpassed mass [counts]")
    axis.set_ylabel("localization error [px]")
    axis.set_title("measured error vs theoretical ep")
    axis.legend(loc="lower left")

    axis = axes[1, 0]
    mass_target = table[:, col["mass_target"]]
    min_count = 5 if len(table) >= 40 else 1
    edges = np.logspace(np.log10(mass_target.min()), np.log10(mass_target.max()), 11)
    centers, fractions = [], []
    for lower, upper in zip(edges[:-1], edges[1:]):
        selected = (mass_target >= lower) & (mass_target < upper)
        if selected.sum() < min_count:
            continue
        centers.append(np.sqrt(lower * upper))
        fractions.append(table[selected, col["matched"]].mean())
    axis.plot(centers, fractions, "o-", color="tab:blue")
    axis.axvline(minmass, color="tab:red", lw=1.2)
    axis.text(minmass, 0.05, f" minmass={minmass:.0f}", rotation=90,
              va="bottom", ha="right", color="tab:red", fontsize=9)
    axis.set_xscale("log")
    axis.set_ylim(-0.02, 1.05)
    axis.set_xlabel("target bandpassed mass [counts]")
    axis.set_ylabel("recovered fraction")
    axis.set_title(f"detection completeness ({len(table)} injections)")

    axis = axes[1, 1]
    upper = np.percentile(radial, 99)
    axis.hist(radial[radial <= upper], bins=50, color="tab:blue", alpha=0.8)
    for value, label, color in [
        (np.median(radial), f"median {np.median(radial):.3f} px", "tab:orange"),
        (np.percentile(radial, 90), f"p90 {np.percentile(radial, 90):.3f} px", "tab:red"),
    ]:
        axis.axvline(value, color=color, lw=1.5)
        axis.text(value, axis.get_ylim()[1], f" {label}", rotation=90,
                  va="top", ha="right", color=color, fontsize=9)
    axis.set_xlabel("radial error |Δ| [px]")
    axis.set_ylabel("spots")
    axis.set_title(f"radial-error distribution (shown up to p99 = {upper:.3f} px)")

    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


def summarize(rows, minmass=DEFAULT_MINMASS):
    """One-block text summary of accuracy and completeness."""
    table = rows_to_array(rows)
    col = {name: index for index, name in enumerate(CSV_HEADER.split(","))}
    matched = table[table[:, col["matched"]] == 1.0]
    above = table[table[:, col["mass_target"]] >= minmass]
    dy, dx = matched[:, col["dy"]], matched[:, col["dx"]]
    radial = np.hypot(dy, dx)
    ratio = radial / matched[:, col["ep_det"]]
    return (
        f"{len(table)} injections, {len(matched)} matched "
        f"({len(matched) / len(table):.1%}); completeness above "
        f"minmass={minmass:.0f}: {(above[:, col['matched']] == 1.0).mean():.1%} "
        f"({int((above[:, col['matched']] == 1.0).sum())}/{len(above)})\n"
        f"bias: mean dy {dy.mean():+.4f} px, mean dx {dx.mean():+.4f} px; "
        f"RMS per axis {np.sqrt(np.mean(dy**2)):.4f}/{np.sqrt(np.mean(dx**2)):.4f} px\n"
        f"radial error: median {np.median(radial):.4f} px, "
        f"p90 {np.percentile(radial, 90):.4f} px, max {radial.max():.4f} px; "
        f"radial/ep ratio median {np.median(ratio):.2f} "
        f"(p16–p84 {np.percentile(ratio, 16):.2f}–{np.percentile(ratio, 84):.2f})"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ADDSTAR-style artificial-spot localization test: inject "
        "synthetic spots of known truth into real frames, re-detect, and "
        "measure localization error and completeness vs mass."
    )
    parser.add_argument("frames_h5", type=Path, help="normalized frames HDF5")
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--central-dot-strength", type=float, required=True)
    parser.add_argument("--n-frames", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260804)
    parser.add_argument("--spots-per-frame", type=int, default=4)
    parser.add_argument("--mass-min", type=float, default=250.0)
    parser.add_argument("--mass-max", type=float, default=30000.0)
    parser.add_argument("--match-radius", type=float, default=2.0)
    parser.add_argument("--figure", type=Path, default=None)
    arguments = parser.parse_args()

    with h5py.File(arguments.frames_h5, "r") as handle:
        condition_name = match_condition_name(handle, arguments.central_dot_strength)

    frame_indices, rows = run_injection_test(
        arguments.frames_h5,
        condition_name,
        n_frames=arguments.n_frames,
        seed=arguments.seed,
        spots_per_frame=arguments.spots_per_frame,
        mass_range=(arguments.mass_min, arguments.mass_max),
        match_radius=arguments.match_radius,
    )
    csv_path = write_csv(arguments.output_csv, rows)
    print(
        f"frames ({len(frame_indices)}, seed {arguments.seed}): "
        f"{','.join(str(int(i)) for i in frame_indices)}",
        flush=True,
    )
    print(summarize(rows), flush=True)
    print(f"per-injection rows -> {csv_path}", flush=True)
    if arguments.figure is not None:
        figure_path = save_summary_figure(rows, arguments.figure)
        print(f"summary figure -> {figure_path}", flush=True)
