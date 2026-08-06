"""Composite definitive figure: paper Fig. 2(a) original vs this work's replica.

Top-left: cropped original Fig. 2(a) of Adiyatullin et al. 2017 (tracked PNG
`fig2a_paper_original.png` next to this script, cropped at 300 dpi from the
tracked paper PDF). Top-right: the replica re-rendered from the extracted
trace archive via render_fig2a primitives (same band selection and rebin as
the standalone definitive render). The caption block below the panels states
the dataset, the band-to-mesa summation, the per-frame CFD50 time-zero
calibration, and the outcome.

Default output overwrites the definitive `fig2a_replica.png`/`.pdf` next to
this script. The caption is plain English (academic presentation use).
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import numpy as np
from matplotlib import pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import render_fig2a as rf  # noqa: E402

DEFAULT_NPZ = Path("drafts/26-08-06_fig2a_cfd50/hbt_traces.npz")
DEFAULT_PAPER_IMAGE = SCRIPT_DIR / "fig2a_paper_original.png"
DEFAULT_OUT = SCRIPT_DIR / "fig2a_replica.png"

CAPTION = (
    "Left: original Fig. 2(a) of Adiyatullin et al. 2017 (Nat. Commun. 8:1329) — accumulated emission intensity of the left (blue) and right\n"
    "(red) mesa on a logarithmic scale, showing anti-phase population oscillations; the gray band marks the low-confidence SNR<20 region.\n"
    "Right: this work — replication on the 26-08-04 HBT batch (20000 single-sweep frames; photon-event source, rebin 8 ≈ 7.9 ps/bin, peak-\n"
    "normalized). Each mesa is recorded simultaneously by both HBT arms; the plotted intensity sums its two arm bands (Left = B1+B3,\n"
    "Right = B2+B4). The per-frame time zero is calibrated by CFD50, leading-edge constant-fraction discrimination on the reference\n"
    "pulse (50%-of-peak crossing, sub-row interpolated), placing the origin ~2.7 ps earlier than the former half-max centroid definition;\n"
    "the measured StreakCamBackG background is subtracted. Outcome: no anti-phase oscillation — cross-arm z fluctuations are\n"
    "uncorrelated (r=0.03–0.10, p>0.14); the 95% upper limit on a phase-locked sinusoidal z amplitude is 0.008–0.019 over 3–200 ps\n"
    "periods (median 0.012; ~1–2% intensity modulation); the t≤150 ps left-side dip is an overdamped monotonic recovery, not an oscillation."
)


def render_comparison(npz_path: Path, paper_image: Path, output_path: Path,
                      bands_l, bands_r, rebin, t_min, t_max, snr_threshold):
    data = rf.load_band_traces(npz_path, "photon")
    time_ps, trace_l, trace_r, sem_l, sem_r, _match_factor = rf.build_lr_traces(
        data, bands_l, bands_r, t_min, t_max, rebin, "peak"
    )
    paper = mpimg.imread(paper_image)

    figure = plt.figure(figsize=(8.8, 5.4))
    grid = figure.add_gridspec(
        2, 2, height_ratios=[1.0, 0.78], width_ratios=[1.25, 1.0],
        hspace=0.18, wspace=0.24,
        left=0.03, right=0.99, top=0.90, bottom=0.05,
    )
    ax_paper = figure.add_subplot(grid[0, 0])
    ax_ours = figure.add_subplot(grid[0, 1])
    ax_caption = figure.add_subplot(grid[1, :])
    ax_caption.axis("off")

    ax_paper.imshow(paper)
    ax_paper.axis("off")
    ax_paper.set_anchor("N")
    ax_paper.set_title(
        "Adiyatullin et al. 2017 (Nat. Commun. 8:1329), Fig. 2(a)", fontsize=8.5
    )
    rf.draw_fig2a_panel(
        ax_ours, time_ps, trace_l, trace_r, sem_l, sem_r, snr_threshold,
        label_l=rf.band_label("Left", bands_l),
        label_r=rf.band_label("Right", bands_r),
    )
    ax_ours.set_title(
        "This work: 26-08-04 HBT, 20000 frames, CFD50-timed", fontsize=8.5, loc="left"
    )

    ax_caption.text(
        0.0, 0.98, CAPTION, fontsize=7.8, va="top", ha="left", linespacing=1.55,
    )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pdf_path = output_path.with_suffix(".pdf")
    figure.savefig(output_path, dpi=300)
    figure.savefig(pdf_path)
    plt.close(figure)
    return output_path, pdf_path


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Composite: paper Fig. 2(a) vs this work's CFD50-timed replica, "
                    "with an English caption recording the calibration and outcome."
    )
    parser.add_argument("--npz", type=Path, default=DEFAULT_NPZ)
    parser.add_argument("--paper-image", type=Path, default=DEFAULT_PAPER_IMAGE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bands-l", default="0,2")
    parser.add_argument("--bands-r", default="1,3")
    parser.add_argument("--rebin", type=int, default=8)
    parser.add_argument("--t-min", type=float, default=0.0)
    parser.add_argument("--t-max", type=float, default=820.0)
    parser.add_argument("--snr-threshold", type=float, default=0.0)
    args = parser.parse_args(argv)

    png_path, pdf_path = render_comparison(
        args.npz, args.paper_image, args.out,
        rf.parse_band_selection(args.bands_l), rf.parse_band_selection(args.bands_r),
        args.rebin, args.t_min, args.t_max, args.snr_threshold,
    )
    print(f"Wrote {png_path} and {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
