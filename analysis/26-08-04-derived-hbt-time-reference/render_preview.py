"""Build a human-browsable PNG preview tree of the 26-08-04 raw folders.

Creates hbt-time-reference-0804-preview/ at the project root, one subfolder
per raw folder under DATA/26-08-04-hbt-time-reference (1:1 mapping). Frames come from the verified
analysis copy DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5
(pixel-exact against the source MAT files); 07_Spectrum has no MAT payload and
is rendered from its raw .tif instead, with the original .spe/.tif copied
alongside.

Every frame is written as a 16-bit grayscale PNG through a strictly increasing
sqrt LUT defined on the folder's global intensity range, so dim background and
sparse bright events are both visible and intensity stays one-to-one with
detector counts (LUT approach ported from the deleted branch's
render_condition_frames.py). Each folder also gets _overview_sum.png, a
matplotlib figure with an explicit colorbar (magma, linear, p1-p99.9 clip) and
the vertical streak axis scaled to ps (1 ns window over 1016 px,
metadata-conventions.toml; row 0 rendered as the window start, in-data
direction calibration pending). The tree's README.md records the exact mapping.

Frame selection: 01_HBT_condensate_Josephson_regime contributes its first
frame per source MAT file (2000 PNGs, user decision 2026-08-05); every other
condition contributes all frames.
"""

import shutil
from pathlib import Path

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from PIL import Image


H5_PATH = Path("DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5")
RAW_ROOT = Path("DATA/26-08-04-hbt-time-reference")
PREVIEW_DIR = Path("hbt-time-reference-0804-preview")
SPECTRUM_RAW_DIR = RAW_ROOT / "07_Spectrum"
WINDOW_PS = 1000.0
SUM_BATCH = 50
HBT_FRAME_LIMIT = 100  # conditions with more frames get first-frame-per-file


def build_sqrt_lut(minimum: int, maximum: int) -> np.ndarray:
    """Strictly increasing LUT mapping [minimum, maximum] counts to 16-bit."""
    if not minimum < maximum:
        raise ValueError(f"Expected minimum < maximum, got {minimum} >= {maximum}")
    values = np.arange(minimum, maximum + 1, dtype=np.float64)
    lut = np.round(np.sqrt((values - minimum) / (maximum - minimum)) * 65535).astype(np.uint16)
    if np.any(np.diff(lut.astype(np.int64)) <= 0):
        raise ValueError("sqrt LUT is not strictly increasing")
    return lut


def sum_and_range(frames: h5py.Dataset) -> tuple[np.ndarray, int, int]:
    """Return (uint64 sum, global min, global max) over every frame."""
    summed = np.zeros(frames.shape[1:], dtype=np.uint64)
    minimum = maximum = None
    for start in range(0, frames.shape[0], SUM_BATCH):
        block = frames[start : start + SUM_BATCH]
        summed += block.sum(axis=0, dtype=np.uint64)
        low, high = int(block.min()), int(block.max())
        minimum = low if minimum is None else min(minimum, low)
        maximum = high if maximum is None else max(maximum, high)
    return summed, minimum, maximum


def save_overview_sum(summed: np.ndarray, raw_dir: str, frame_count: int, output_path: Path) -> None:
    lower, upper = np.percentile(summed, [1.0, 99.9])
    figure, axis = plt.subplots(figsize=(9.6, 7.2), constrained_layout=True)
    image = axis.imshow(
        summed,
        origin="lower",
        aspect="auto",
        cmap="magma",
        interpolation="nearest",
        extent=(0.0, float(summed.shape[1]), 0.0, WINDOW_PS),
        vmin=lower,
        vmax=upper,
    )
    axis.set_title(f"26-08-04 {raw_dir} — sum of {frame_count:,} streak frames")
    axis.set_xlabel("Detector axis 2 [pixel]")
    axis.set_ylabel("Time [ps]")
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label(
        f"Summed camera counts [a.u.]\n(linear, magma, clipped to p1–p99.9 = {lower:,.0f}–{upper:,.0f})"
    )
    figure.savefig(output_path, dpi=200)
    plt.close(figure)


def save_gray16(array: np.ndarray, output_path: Path) -> None:
    Image.fromarray(array).save(output_path)


def render_condition(condition: h5py.Group, target_dir: Path) -> dict:
    """Render one condition's folder; return per-folder stats for the README."""
    raw_dir = condition.attrs["source_block_directory"]
    frames = condition["frames"]
    frame_count = int(frames.shape[0])
    target_dir.mkdir(parents=True, exist_ok=True)

    summed, minimum, maximum = sum_and_range(frames)
    save_overview_sum(summed, raw_dir, frame_count, target_dir / "_overview_sum.png")
    lut = build_sqrt_lut(minimum, maximum)

    file_numbers = condition["source_files/file_number"][:]
    first_frames = condition["source_files/first_output_frame"][:]
    frame_counts = condition["source_files/frame_count"][:]
    written = 0
    for file_number, first_frame, count in zip(file_numbers, first_frames, frame_counts):
        if frame_count > HBT_FRAME_LIMIT:
            indices = [int(first_frame)]  # first frame per source file
        else:
            indices = [int(first_frame) + i for i in range(int(count))]
        for index in indices:
            in_file = index - int(first_frame)
            name = f"file_{int(file_number):05d}.png" if frame_count > HBT_FRAME_LIMIT else f"file_{int(file_number):05d}_frame_{in_file:02d}.png"
            save_gray16(lut[frames[index] - minimum], target_dir / name)
            written += 1

    stats = {
        "raw_dir": raw_dir,
        "frames_total": frame_count,
        "frames_rendered": written,
        "minimum": minimum,
        "maximum": maximum,
        "background": False,
    }
    if "background" in condition:
        background = condition["background"][:]
        low, high = float(background.min()), float(background.max())
        stretched = np.round((background - low) / (high - low) * 65535).astype(np.uint16)
        save_gray16(stretched, target_dir / "_background.png")
        stats["background"] = True
        stats["background_range"] = (low, high)
    return stats


def render_spectrum(target_dir: Path) -> dict:
    """Render 07_Spectrum from its raw TIFF and copy the original SPE/TIFF."""
    target_dir.mkdir(parents=True, exist_ok=True)
    tif_path = next(SPECTRUM_RAW_DIR.glob("*.tif"))
    spe_path = next(SPECTRUM_RAW_DIR.glob("*.spe"))
    array = np.array(Image.open(tif_path))
    minimum, maximum = int(array.min()), int(array.max())
    save_gray16(build_sqrt_lut(minimum, maximum)[array - minimum], target_dir / "spectrum.png")
    shutil.copy2(tif_path, target_dir / tif_path.name)
    shutil.copy2(spe_path, target_dir / spe_path.name)
    return {"raw_dir": "07_Spectrum", "shape": array.shape, "minimum": minimum, "maximum": maximum}


def write_readme(stats: list[dict], spectrum_stats: dict) -> None:
    lines = [
        "# 26-08-04 preview — human-browsable PNG mirror of the raw folders",
        "",
        "One subfolder per raw folder under `DATA/26-08-04-hbt-time-reference/` (1:1).",
        "Frames are read from the verified analysis copy",
        "`DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5`",
        "(pixel-exact against the source MAT files, verified 2026-08-05).",
        "Regenerate with `analysis/26-08-04-derived-hbt-time-reference/render_preview.py`.",
        "",
        "## Intensity mapping (frame PNGs)",
        "",
        "16-bit grayscale PNG, one shared strictly increasing mapping per folder:",
        "`PNG = round(65535 * sqrt((counts - min) / (max - min)))` with the folder's",
        "global [min, max] listed below. PNG row 0 is detector row 0 (top); the vertical",
        "axis is the streak time axis (1 ns window over 1016 px, row 0 rendered as the",
        "window start; in-data direction calibration pending). `_overview_sum.png` carries",
        "the explicit colorbar figure (magma, linear, p1–p99.9 clip, vertical axis in ps).",
        "",
        "| folder | frames in data | frame PNGs | count range [min, max] |",
        "| --- | --- | --- | --- |",
    ]
    for entry in stats:
        lines.append(
            f"| {entry['raw_dir']} | {entry['frames_total']:,} | {entry['frames_rendered']:,} "
            f"| {entry['minimum']}–{entry['maximum']} |"
        )
    lines += [
        "",
        "Frame selection: `01_HBT_condensate_Josephson_regime` contributes the first frame",
        "of each source MAT file (`file_NNNNN.png` = source file number); all other",
        "conditions contribute every frame (`file_NNNNN_frame_MM.png`).",
        "",
    ]
    for entry in stats:
        if entry.get("background"):
            low, high = entry["background_range"]
            lines.append(
                f"`{entry['raw_dir']}/_background.png`: measured background "
                f"(`structVar/StreakCamBackG`), linear stretch of [{low:.1f}, {high:.1f}] to 16-bit."
            )
    lines += [
        "",
        "`01_HBT_condensate_Josephson_regime/_overview_sum_bg_subtracted.png`: sum of all 20,000",
        "HBT frames minus 20,000×`StreakCamBackG`, two panels (full columns; reference-spot columns",
        "150–270 masked), magma, linear, explicit colorbars, vertical axis in ps. Unlike",
        "`_overview_sum.png` (raw counts, dominated by the reference spot), this paneling makes the",
        "four data bands (B1 x=380–450, B2 570–680, B3 810–880, B4 1010–1110, cyan overlays)",
        "directly visible. Generated by `analysis/26-08-04-derived-hbt-time-reference/accumulate_block_traces.py`.",
        "",
        "## 07_Spectrum",
        "",
        f"No MAT payload; `spectrum.png` renders the raw TIFF "
        f"({spectrum_stats['shape'][0]}×{spectrum_stats['shape'][1]} uint16, "
        f"sqrt LUT over {spectrum_stats['minimum']}–{spectrum_stats['maximum']}). "
        "The original `.spe`/`.tif` are copied alongside.",
        "",
        "## Frame-overlay samples",
        "",
        "`frame_overlays/` and `frame_overlays_new_folder_10/` are annotated samples of",
        "the HBT condition (raw frame, detected events kept/rejected by the Fig. 2(a)",
        "band x-ROIs, band x-ROIs overlaid) — not raw-folder mirrors. Regenerate with",
        "`analysis/26-08-04-derived-hbt-time-reference/render_frame_overlays.py` (main",
        "run) and `render_frame_overlays_new_folder_10.py` (supplement run); see the",
        "README inside each folder.",
        "",
        "This tree is bulk preview material for manual browsing. The numbered",
        "raw-folder mirrors can be deleted and regenerated at any time with",
        "`render_preview.py`.",
    ]
    (PREVIEW_DIR / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    PREVIEW_DIR.mkdir(parents=True, exist_ok=True)
    stats = []
    with h5py.File(H5_PATH, "r") as handle:
        conditions = handle["conditions"]
        ordered = sorted(conditions.values(), key=lambda group: group.attrs["source_block_directory"])
        for condition in ordered:
            entry = render_condition(condition, PREVIEW_DIR / condition.attrs["source_block_directory"])
            stats.append(entry)
            print(
                f"{entry['raw_dir']}: {entry['frames_rendered']} frame PNGs "
                f"(of {entry['frames_total']} frames, range {entry['minimum']}–{entry['maximum']})",
                flush=True,
            )
    spectrum_stats = render_spectrum(PREVIEW_DIR / "07_Spectrum")
    print(f"07_Spectrum: spectrum.png + copied spe/tif {spectrum_stats['shape']}", flush=True)
    write_readme(stats, spectrum_stats)
    print(f"README written to {PREVIEW_DIR / 'README.md'}", flush=True)


if __name__ == "__main__":
    main()
