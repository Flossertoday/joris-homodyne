"""Run sparse peak detection over every frame of one central-dot condition.

Applies the production Crocker-Grier peak finder from ``find_frame_peaks.py``
(same defaults, imported so this batch job cannot drift away from the audited
parameters) to all frames of the condition selected by its
``central_dot_strength`` attribute, and stores every detected peak as one
``(frame, y, x, mass)`` row in a single HDF5 file. That file is the raw input
for downstream analyses (x-profile fitting, correlations), so it records the
full detection provenance as root attributes.

Example:

    python extract_condition_peaks.py <in.h5> <out.h5> --central-dot-strength 0.18
"""
import argparse
import hashlib
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

# Same detector-calibrated defaults as the production peak finder, imported so
# this batch extraction cannot silently drift away from the parameters that
# were validated in docs/peak-localization.md.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from find_frame_peaks import (  # noqa: E402
    DEFAULT_DIAMETER,
    DEFAULT_MINMASS,
    DEFAULT_PERCENTILE,
    find_peaks,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def match_condition_name(handle: h5py.File, central_dot_strength: float) -> str:
    """Name of the single condition carrying the given central-dot strength."""
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
    return matches[0]


def extract_condition_peaks(
    input_path: Path,
    central_dot_strength: float,
    diameter: int = DEFAULT_DIAMETER,
    minmass: float = DEFAULT_MINMASS,
    percentile: float = DEFAULT_PERCENTILE,
    progress_every: int = 100,
) -> tuple[str, int, tuple[int, int], np.ndarray]:
    """Detect peaks in every frame.

    Returns (condition_name, frame_count, frame_shape, rows). Rows are an
    (N, 4) float64 array with columns frame, y, x, mass, ordered by frame.
    Frames legitimately containing no peaks contribute zero rows.
    """
    with h5py.File(input_path, "r") as handle:
        condition_name = match_condition_name(handle, central_dot_strength)
        frames = handle["conditions"][condition_name]["frames"]
        if frames.ndim != 3:
            raise ValueError(f"Expected a 3-D frames dataset, got shape {frames.shape}")

        frame_count = frames.shape[0]
        frame_shape = tuple(int(dim) for dim in frames.shape[1:])
        chunks = []
        started = time.time()
        for index in range(frame_count):
            peaks = find_peaks(
                np.asarray(frames[index]),
                diameter=diameter,
                minmass=minmass,
                percentile=percentile,
            )
            if len(peaks):
                chunks.append(
                    np.column_stack([np.full(len(peaks), index), peaks])
                )
            if progress_every and (index + 1) % progress_every == 0:
                print(
                    f"{index + 1}/{frame_count} frames, "
                    f"{time.time() - started:.0f}s",
                    flush=True,
                )

    if chunks:
        rows = np.concatenate(chunks).astype(np.float64)
    else:
        rows = np.empty((0, 4), dtype=np.float64)
    return condition_name, frame_count, frame_shape, rows


def write_peaks_h5(
    output_path: Path,
    rows: np.ndarray,
    attributes: dict,
) -> Path:
    """Write the peak table, refusing to overwrite; atomic via .partial."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing peaks file: {output_path}"
        )
    partial_path = output_path.with_name(f"{output_path.name}.partial")
    if partial_path.exists():
        raise FileExistsError(
            f"Interrupted extraction artifact already exists: {partial_path}"
        )
    try:
        with h5py.File(partial_path, "x") as output:
            output.create_dataset("frame", data=rows[:, 0].astype(np.int32))
            output.create_dataset("y", data=rows[:, 1])
            output.create_dataset("x", data=rows[:, 2])
            output.create_dataset("mass", data=rows[:, 3])
            for key, value in attributes.items():
                output.attrs[key] = value
        if output_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing peaks file: {output_path}"
            )
        partial_path.rename(output_path)
    except BaseException:
        if partial_path.exists():
            partial_path.unlink()
        raise
    return output_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Detect sparse peaks in every frame of one central-dot "
        "condition and store them as a (frame, y, x, mass) HDF5 table."
    )
    parser.add_argument("input_h5", type=Path)
    parser.add_argument("output_h5", type=Path)
    parser.add_argument("--central-dot-strength", type=float, required=True)
    parser.add_argument("--diameter", type=int, default=DEFAULT_DIAMETER)
    parser.add_argument("--minmass", type=float, default=DEFAULT_MINMASS)
    parser.add_argument("--percentile", type=float, default=DEFAULT_PERCENTILE)
    arguments = parser.parse_args()

    condition_name, frame_count, frame_shape, rows = extract_condition_peaks(
        arguments.input_h5,
        central_dot_strength=arguments.central_dot_strength,
        diameter=arguments.diameter,
        minmass=arguments.minmass,
        percentile=arguments.percentile,
    )

    attributes = {
        "condition_name": condition_name,
        "central_dot_strength": arguments.central_dot_strength,
        "diameter": arguments.diameter,
        "minmass": arguments.minmass,
        "percentile": arguments.percentile,
        "frame_count": frame_count,
        "peak_count": len(rows),
        "frame_shape": np.asarray(frame_shape, dtype=np.int64),
        "source_h5": arguments.input_h5.name,
        "source_h5_sha256": _sha256(arguments.input_h5),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": Path(__file__).name,
    }
    created_path = write_peaks_h5(arguments.output_h5, rows, attributes)
    print(
        f"Extracted {len(rows)} peaks from {condition_name} "
        f"({frame_count} frames) -> {created_path}",
        flush=True,
    )
