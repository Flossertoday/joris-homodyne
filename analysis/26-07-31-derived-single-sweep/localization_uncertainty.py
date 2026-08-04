"""Per-peak localization uncertainty (Savin-Doyle static error) for one condition.

Adds an ``ep`` column to a previously extracted peak table. ``ep`` is the
Savin-Doyle static localization error (Savin & Doyle, Biophys. J. 88, 623
(2005)) as implemented in ``trackpy.uncertainty.static_error``:

    ep = (noise / mass) * noise_size * sqrt(sum of x^2 over the refine mask)

with ``mass`` the bandpassed integrated brightness stored by the production
peak finder and ``noise`` the per-pixel standard deviation of the bandpassed
image background, measured per frame.

Two implementation choices, both validated against the artificial-spot
Monte-Carlo in ``artificial_spot_test.py`` (see docs/peak-localization-error.md):

- trackpy's ``measure_noise`` is NOT used: it identifies background as the
  complement of the dilated nonzero bandpassed image, which degenerates on
  these sparse frames (about half the pixels stay nonzero after clipping), and
  it measures the raw-image standard deviation while ``mass`` is in bandpass
  units. Instead the noise is the robust (MAD) standard deviation of the
  UNCLIPPED bandpass image ``gaussian(sigma=1) - 7x7 boxcar``, i.e. exactly
  the image on which centroid and mass are computed, before trackpy clips it
  at the threshold. Spots and impulse noise cover ~0.1% of pixels, negligible
  for a median-based estimator.
- ``noise_size=1`` matches the bandpass lshort used by ``tp.locate``.

The output HDF5 carries the full provenance of both inputs.

Example:

    python localization_uncertainty.py peaks.h5 frames.h5 out_ep.h5 --figure ep.png
"""
import argparse
import hashlib
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd
import trackpy

matplotlib.use("Agg")
from matplotlib import pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_frame_peaks import DEFAULT_DIAMETER  # noqa: E402

from scipy.ndimage import gaussian_filter, uniform_filter  # noqa: E402
from trackpy.uncertainty import static_error  # noqa: E402
import trackpy  # noqa: E402

NOISE_SIZE = 1  # bandpass lshort used by tp.locate defaults


def bandpass_image(frame, diameter=DEFAULT_DIAMETER, noise_size=NOISE_SIZE):
    """Unclipped trackpy bandpass: Gaussian lowpass minus boxcar background.

    Identical to ``trackpy.preprocessing.bandpass(frame, noise_size, diameter)``
    except that the result is not clipped at the threshold, so the background
    keeps its symmetric noise distribution.
    """
    lowpass = gaussian_filter(np.asarray(frame, dtype=np.float64), noise_size)
    background = uniform_filter(lowpass, size=diameter)
    return lowpass - background


def robust_background_std(bandpassed):
    """Median-based standard deviation of a mostly-background image."""
    median = np.median(bandpassed)
    mad = np.median(np.abs(bandpassed - median))
    return 1.4826 * mad


def frame_background_noise(frame, diameter=DEFAULT_DIAMETER, noise_size=NOISE_SIZE):
    """Per-pixel background noise of the bandpassed image for one frame."""
    return robust_background_std(bandpass_image(frame, diameter, noise_size))


def static_ep(masses, noise, diameter=DEFAULT_DIAMETER, noise_size=NOISE_SIZE):
    """Savin-Doyle static error (px) for peaks of given bandpassed masses."""
    features = pd.DataFrame({"mass": np.asarray(masses, dtype=np.float64)})
    ep = static_error(features, float(noise), diameter, noise_size=noise_size)
    return np.asarray(ep, dtype=np.float64)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def compute_frame_noises(frames_path, condition_name, diameter, progress_every=100):
    """Bandpassed background noise for every frame of the condition."""
    with h5py.File(frames_path, "r") as handle:
        frames = handle["conditions"][condition_name]["frames"]
        noises = np.empty(frames.shape[0], dtype=np.float64)
        started = time.time()
        for index in range(frames.shape[0]):
            noises[index] = frame_background_noise(
                np.asarray(frames[index]), diameter=diameter
            )
            if progress_every and (index + 1) % progress_every == 0:
                print(
                    f"{index + 1}/{frames.shape[0]} frames, "
                    f"{time.time() - started:.0f}s",
                    flush=True,
                )
    return noises


def compute_peak_ep(peaks_path, frame_noises, diameter, noise_size=NOISE_SIZE):
    """One ep per peak, using the noise of the frame the peak was found in."""
    with h5py.File(peaks_path, "r") as handle:
        features = pd.DataFrame(
            {
                "frame": np.asarray(handle["frame"], dtype=np.int64),
                "mass": np.asarray(handle["mass"], dtype=np.float64),
            }
        )
    noise = pd.DataFrame(
        {
            "frame": np.arange(len(frame_noises), dtype=np.int64),
            "noise": np.asarray(frame_noises, dtype=np.float64),
        }
    ).set_index("frame")
    ep = static_error(features, noise, diameter, noise_size=noise_size)
    return np.asarray(ep, dtype=np.float64)


def write_ep_h5(output_path, peaks_path, ep, frame_noises, attributes):
    """Write peaks + ep + per-frame noise, refusing to overwrite; atomic."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing ep file: {output_path}"
        )
    partial_path = output_path.with_name(f"{output_path.name}.partial")
    if partial_path.exists():
        raise FileExistsError(
            f"Interrupted ep artifact already exists: {partial_path}"
        )
    try:
        with h5py.File(peaks_path, "r") as source, h5py.File(
            partial_path, "x"
        ) as output:
            for name in ("frame", "y", "x", "mass"):
                output.create_dataset(name, data=np.asarray(source[name]))
            output.create_dataset("ep", data=ep)
            output.create_dataset("frame_noise", data=frame_noises)
            for key, value in attributes.items():
                output.attrs[key] = value
        if output_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing ep file: {output_path}"
            )
        partial_path.rename(output_path)
    except BaseException:
        if partial_path.exists():
            partial_path.unlink()
        raise
    return output_path


def save_ep_figure(mass, ep, output_path, draw_seed=20260804, n_draw=5000):
    """Two-panel summary: ep vs mass (log-log) and the ep distribution."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing figure: {output_path}"
        )

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True)

    axis = axes[0]
    rng = np.random.default_rng(draw_seed)
    if len(mass) > n_draw:
        drawn = rng.choice(len(mass), n_draw, replace=False)
    else:
        drawn = np.arange(len(mass))
    axis.scatter(
        mass[drawn], ep[drawn], s=3, alpha=0.3, color="tab:blue", rasterized=True
    )
    edges = np.logspace(np.log10(mass.min()), np.log10(mass.max()), 13)
    centers, medians, lows, highs = [], [], [], []
    for lower, upper in zip(edges[:-1], edges[1:]):
        selected = ep[(mass >= lower) & (mass < upper)]
        if len(selected) < 5:
            continue
        centers.append(np.sqrt(lower * upper))
        medians.append(np.median(selected))
        lows.append(np.percentile(selected, 16))
        highs.append(np.percentile(selected, 84))
    axis.fill_between(centers, lows, highs, color="tab:orange", alpha=0.25,
                      label="p16–p84")
    axis.plot(centers, medians, color="tab:orange", lw=1.5, label="median")
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("bandpassed mass [counts]")
    axis.set_ylabel("ep (Savin–Doyle static error) [px]")
    axis.set_title(f"ep vs mass — {len(mass)} peaks "
                   f"(scatter shows {len(drawn)} drawn)")
    axis.legend(loc="lower left")

    axis = axes[1]
    upper = np.percentile(ep, 99)
    axis.hist(ep[ep <= upper], bins=60, color="tab:blue", alpha=0.8)
    for value, label, color in [
        (np.median(ep), f"median {np.median(ep):.3f} px", "tab:orange"),
        (np.percentile(ep, 90), f"p90 {np.percentile(ep, 90):.3f} px", "tab:red"),
    ]:
        axis.axvline(value, color=color, lw=1.5)
        axis.text(value, axis.get_ylim()[1], f" {label}", rotation=90,
                  va="top", ha="right", color=color, fontsize=9)
    axis.set_xlabel("ep [px]")
    axis.set_ylabel("peaks")
    axis.set_title(f"ep distribution (shown up to p99 = {upper:.3f} px)")

    figure.savefig(output_path, dpi=150)
    plt.close(figure)
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add per-peak Savin-Doyle static localization error (ep) "
        "to an extracted peak table, using per-frame bandpass background noise."
    )
    parser.add_argument("peaks_h5", type=Path, help="extracted peak table")
    parser.add_argument("frames_h5", type=Path, help="normalized frames HDF5")
    parser.add_argument("output_h5", type=Path)
    parser.add_argument("--diameter", type=int, default=DEFAULT_DIAMETER)
    parser.add_argument("--figure", type=Path, default=None)
    arguments = parser.parse_args()

    with h5py.File(arguments.peaks_h5, "r") as handle:
        condition_name = str(handle.attrs["condition_name"])
        peak_count = int(handle.attrs["peak_count"])
        source_h5 = str(handle.attrs["source_h5"])
        source_h5_sha256 = str(handle.attrs["source_h5_sha256"])

    frame_noises = compute_frame_noises(
        arguments.frames_h5, condition_name, arguments.diameter
    )
    ep = compute_peak_ep(arguments.peaks_h5, frame_noises, arguments.diameter)

    attributes = {
        "condition_name": condition_name,
        "diameter": arguments.diameter,
        "noise_size": NOISE_SIZE,
        "noise_method": "MAD standard deviation of the unclipped bandpass "
        "image (gaussian sigma=1 minus boxcar diameter) per frame",
        "trackpy_version": trackpy.__version__,
        "frame_count": len(frame_noises),
        "peak_count": peak_count,
        "peaks_h5": arguments.peaks_h5.name,
        "peaks_h5_sha256": _sha256(arguments.peaks_h5),
        "source_h5": source_h5,
        "source_h5_sha256": source_h5_sha256,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
    }
    created_path = write_ep_h5(
        arguments.output_h5, arguments.peaks_h5, ep, frame_noises, attributes
    )

    with h5py.File(arguments.peaks_h5, "r") as handle:
        mass = np.asarray(handle["mass"])
    print(
        f"{len(ep)} peaks from {condition_name}: "
        f"frame noise median {np.median(frame_noises):.3f} "
        f"(p5–p95 {np.percentile(frame_noises, 5):.3f}–"
        f"{np.percentile(frame_noises, 95):.3f}) bandpass counts; "
        f"ep median {np.median(ep):.4f} px "
        f"(p5–p95 {np.percentile(ep, 5):.4f}–{np.percentile(ep, 95):.4f}) "
        f"-> {created_path}",
        flush=True,
    )
    if arguments.figure is not None:
        figure_path = save_ep_figure(mass, ep, arguments.figure)
        print(f"ep summary figure -> {figure_path}", flush=True)
