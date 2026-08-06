"""Composite definitive figure: paper Fig. 2(a) original vs this work's replica.

Top-left: cropped original Fig. 2(a) of Adiyatullin et al. 2017 (tracked PNG
`fig2a_paper_original.png` next to this script, cropped at 300 dpi from the
tracked paper PDF). Top-right: the replica re-rendered from the extracted
trace archive via render_fig2a primitives (same band selection and rebin as
the standalone definitive render). The caption block below the panels states
the dataset, the per-frame CFD50 time-zero calibration, and the outcome.

Default output overwrites the definitive `fig2a_replica.png`/`.pdf` next to
this script. Requires a CJK-capable font (Microsoft YaHei on Windows) for the
Chinese caption; falls back to the default family if absent.
"""

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import numpy as np
from matplotlib import font_manager, pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import render_fig2a as rf  # noqa: E402

DEFAULT_NPZ = Path("drafts/26-08-06_fig2a_cfd50/hbt_traces.npz")
DEFAULT_PAPER_IMAGE = SCRIPT_DIR / "fig2a_paper_original.png"
DEFAULT_OUT = SCRIPT_DIR / "fig2a_replica.png"

CAPTION = (
    "左: 论文原图 Adiyatullin et al. 2017 (Nat. Commun. 8:1329) Fig. 2(a) —— 左 (蓝) / 右 (红) mesa 的累积发光强度, 对数轴, 反相交替振荡; 底部灰区为 SNR<20 的\n"
    "低置信区。右: 本工作对 26-08-04 HBT 批次 20000 帧 single-sweep 的复刻 (photon 事件源, rebin 8 ≈ 7.9 ps/bin, 峰值归一)。逐帧时间零点经 CFD50 校准: 参考\n"
    "脉冲前沿恒比定时 (leading-edge constant-fraction discrimination, 50% 峰高穿越, 亚行插值), 时间原点较旧 half-max 质心定义整体提前 ~2.7 ps; 背景为实测\n"
    "StreakCamBackG 扣除。结果: 未见反相交替振荡 —— 两臂 z 涨落相关不显著 (r=0.03–0.10, p>0.14), 3–200 ps 周期锁相正弦振幅 95% 上限 0.012–0.019\n"
    "(约 1–2% 强度调制); t≤150 ps 的 L 侧缺口为过阻尼式单调恢复。"
)


def cjk_font(name="Microsoft YaHei") -> font_manager.FontProperties | None:
    """Return the CJK font when available so the Chinese caption renders."""
    available = {f.name for f in font_manager.fontManager.ttflist}
    return font_manager.FontProperties(family=name) if name in available else None


def render_comparison(npz_path: Path, paper_image: Path, output_path: Path,
                      bands_l, bands_r, rebin, t_min, t_max, snr_threshold):
    data = rf.load_band_traces(npz_path, "photon")
    time_ps, trace_l, trace_r, sem_l, sem_r = rf.build_lr_traces(
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
        ax_ours, time_ps, trace_l, trace_r, sem_l, sem_r, snr_threshold
    )
    ax_ours.set_title(
        "This work: 26-08-04 HBT, 20000 frames, CFD50-timed", fontsize=8.5, loc="left"
    )

    font = cjk_font()
    ax_caption.text(
        0.0, 0.98, CAPTION, fontsize=7.8, va="top", ha="left", linespacing=1.55,
        fontproperties=font,
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
                    "with a Chinese caption recording the calibration and outcome."
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
