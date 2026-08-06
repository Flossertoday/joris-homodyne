"""Composite definitive figure: paper Fig. 2(c) original vs this work's replica.

Top-left: cropped original Fig. 2(c) of Adiyatullin et al. 2017 (tracked PNG
`fig2c_paper_original.png` next to this script, cropped at 300 dpi from the
tracked paper PDF with the shared Time (ps) axis). Top-right: this work's
g2(0)(t) panel from the g2 results archive via render_fig2c.draw_fig2c_panel.
The caption block below states the dataset, the estimator, the G3 gate
outcome (statistics-limited per-bin regime), and the final verdict.

Default output: `fig2c_replica.png`/`.pdf` next to this script. The caption
is plain English (academic presentation use).
"""

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import numpy as np
from matplotlib import pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import render_fig2c as rf  # noqa: E402

DEFAULT_OUT_DIR = Path("drafts/26-08-06_fig2c")
DEFAULT_PAPER_IMAGE = SCRIPT_DIR / "fig2c_paper_original.png"
DEFAULT_OUT = SCRIPT_DIR / "fig2c_replica.png"


def build_caption(width, pooled_l, pooled_r, ul95_text):
    return (
        "Left: original Fig. 2(c) of Adiyatullin et al. 2017 (Nat. Commun. 8:1329) — time-resolved zero-delay second-order correlation\n"
        "g(2)(0)(t) of the left (blue) and right (red) mesa; shaded areas are ±1 s.e.m.; the light statistics switches between Poissonian\n"
        "and bunched in phase with the Josephson oscillations, the two mesas in counterphase. Right: this work — same estimator\n"
        "(Suppl. Note 1) on the 26-08-04 HBT batch (19992 single-sweep frames, CFD50 time-zero, measured-background subtraction, "
        f"{width} px bins).\n"
        "The per-bin regime is statistics-limited: ~50% of detected events are residual detector-noise clusters and the per-bin\n"
        "occupancy is ~1e-3 events/frame, so per-bin g2 is a heavy-tailed statistic (shown for the record with plug-in ±1 s.e.m. bands\n"
        "validated by Monte Carlo; gray strip: s.e.m. > 0.5). Coherent endpoints: window-pooled "
        f"g2(0) = {pooled_l[0]:.2f} ± {pooled_l[1]:.2f} (L), {pooled_r[0]:.2f} ± {pooled_r[1]:.2f} (R) (jackknife s.e.m.);\n"
        f"95% upper limit on a reference-locked sinusoidal g2-1 oscillation (3-200 ps): {ul95_text}. No bunching oscillation at the\n"
        "paper's amplitudes (0.2-1.2) is resolvable below the upper edge of that range with 20k frames; ~3e5 frames would be needed\n"
        "for 95% sensitivity to |A| = 0.2 (G3 power gate, docs/26-08-04-fig2c.md)."
    )


def render_comparison(out_dir: Path, paper_image: Path, output_path: Path, width: int):
    data = np.load(out_dir / f"g2_results_bin{width}px.npz")
    with open(out_dir / "g2_summary.json") as handle:
        summary = json.load(handle)
    science_path = out_dir / "fig2c_science_report.json"
    science = json.load(open(science_path)) if science_path.exists() else None

    t_ps = np.asarray(data["t_ps"], dtype=np.float64)
    g2_l = np.asarray(data["g2_L"], dtype=np.float64)
    g2_r = np.asarray(data["g2_R"], dtype=np.float64)
    sem_l = np.asarray(data["sem_plugin_L"], dtype=np.float64)
    sem_r = np.asarray(data["sem_plugin_R"], dtype=np.float64)
    pairs = summary["per_width"][str(width)]["pairs"]
    pooled = {
        "L": (pairs["L"]["pooled_g2"], pairs["L"]["pooled_sem_jackknife"]),
        "R": (pairs["R"]["pooled_g2"], pairs["R"]["pooled_sem_jackknife"]),
    }
    if science is not None:
        combined = science["per_width"][str(width)]["oscillation"]["combined"]
        ul95_text = (
            f"median {combined['ul95_median']:.2f}, "
            + ", ".join(
                f"{v['ul95']:.2f} at {period}"
                for period, v in combined["targeted"].items()
            )
        )
    else:
        ul95_text = "see power-gate report (science tests pending)"

    paper = mpimg.imread(paper_image)
    figure = plt.figure(figsize=(9.6, 5.6))
    grid = figure.add_gridspec(
        2, 2, height_ratios=[1.0, 0.86], width_ratios=[1.0, 1.28],
        hspace=0.18, wspace=0.22,
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
        "Adiyatullin et al. 2017 (Nat. Commun. 8:1329), Fig. 2(c)", fontsize=8.5
    )
    rf.draw_fig2c_panel(ax_ours, t_ps, g2_l, sem_l, g2_r, sem_r, pooled=pooled)
    ax_ours.set_title(
        f"This work: 26-08-04 HBT, 19992 frames, {width} px bins", fontsize=8.5,
        loc="left",
    )

    ax_caption.text(
        0.0, 0.98, build_caption(width, pooled["L"], pooled["R"], ul95_text),
        fontsize=7.8, va="top", ha="left", linespacing=1.55,
    )

    output_path = Path(output_path)
    pdf_path = output_path.with_suffix(".pdf")
    figure.savefig(output_path, dpi=300)
    figure.savefig(pdf_path)
    plt.close(figure)
    return output_path, pdf_path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--paper-image", type=Path, default=DEFAULT_PAPER_IMAGE)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--width", type=int, default=4)
    args = parser.parse_args(argv)

    png_path, pdf_path = render_comparison(
        args.out_dir, args.paper_image, args.out, args.width
    )
    print(f"Wrote {png_path} and {pdf_path}", flush=True)


if __name__ == "__main__":
    main()
