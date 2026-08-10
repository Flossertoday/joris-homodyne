"""Render Figure: results dashboard — upper limits and pooled g^(2)(0) summary.

Output: analysis/26-08-10-project-report/figures/fig_results_dashboard.png
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
from matplotlib import pyplot as plt

SCIENCE_JSON = Path(
    "E:/Documents/Codex/joris-homodyne/drafts/26-08-06_fig2c/fig2c_science_report.json"
)
SUPPLEMENT_SCIENCE_JSON = Path(
    "E:/Documents/Codex/joris-homodyne/drafts/26-08-06_fig2c_supplement/fig2c_science_report.json"
)
OUT_PATH = Path(__file__).resolve().parent / "figures" / "fig_results_dashboard.png"


def _read_g2_ul95_median(json_path):
    """Read per_width.4.oscillation.combined.ul95_median from the science report."""
    with open(json_path, "r", encoding="utf-8") as handle:
        data = json.load(handle)
    return float(data["per_width"]["4"]["oscillation"]["combined"]["ul95_median"])


def build_dashboard_data():
    """Return the dashboard data as a plain dict with provenance comments."""
    data = {
        "panel_a": {
            "title": "population oscillation upper limits (Fig. 2(a)-type)",
            "ylabel": "z-space sinusoid amplitude UL95 (median over 3–200 ps scan)",
            "y_categories": [
                "main batch, 19992 frames",
                "supplement batch, 71971 frames",
                "dataset A, main",
                "dataset A, supplement",
            ],
            "values": np.array([
                0.0123,   # DOCS/26-08-04-fig2a.md 2026-08-06, CFD50, main batch
                0.0084,   # DOCS/26-08-06-hbt-supplement.md §5, supplement batch
                0.0131,   # DOCS/26-08-08-dataset-a.md, dataset A main (raw)
                0.0083,   # DOCS/26-08-08-dataset-a.md, dataset A supplement (raw)
            ]),
            "corrected_values": np.array([
                np.nan,
                np.nan,
                0.0226,   # DOCS/26-08-08-dataset-a.md: dilution-corrected (S/(S+N))
                0.0161,   # DOCS/26-08-08-dataset-a.md: dilution-corrected
            ]),
            "band_min": 0.4,
            "band_max": 0.7,
            "band_label": "paper Fig. 2(a) anti-phase signal (measured z ≈ ±0.4–0.5)",
        },
        "panel_b": {
            "title": "g^(2) oscillation upper limits (Fig. 2(c)-type)",
            "ylabel": "|A| + 2σ UL95 (median over 3–200 ps, anti-phase combined, 4 px bins)",
            "y_categories": [
                "main batch",
                "supplement batch",
                "dataset A, main",
                "dataset A, supplement",
            ],
            "values": np.array([
                _read_g2_ul95_median(SCIENCE_JSON),  # drafts/26-08-06_fig2c/fig2c_science_report.json
                _read_g2_ul95_median(SUPPLEMENT_SCIENCE_JSON),  # drafts/26-08-06_fig2c_supplement/fig2c_science_report.json; sanity: 0.496
                0.58,    # DOCS/26-08-08-dataset-a.md
                0.488,   # DOCS/26-08-08-dataset-a.md
            ]),
            "band_min": 0.2,
            "band_max": 1.2,
            "band_label": "paper Fig. 2(c): g2−1 = 0.2–1.2",
        },
        "panel_c": {
            "title": "window-pooled g^(2)(0): the full two-condensate matrix (4 px bins, ±1 jackknife SEM)",
            "ylabel": "pooled g^(2)(0) in 0–820 ps window",
            "y_categories": [
                "main batch",
                "supplement batch",
                "dataset A, main",
                "dataset A, supplement",
            ],
            "L_values": np.array([
                0.310,   # drafts/26-08-06_fig2c/g2_summary.json main 4 px
                1.304,   # drafts/26-08-06_fig2c_supplement/g2_summary.json supplement 4 px
                0.531,   # DOCS/26-08-08-dataset-a.md local-floor beta
                1.256,   # DOCS/26-08-08-dataset-a.md supplement local-floor beta
            ]),
            "L_sem": np.array([
                0.380,
                0.320,
                0.258,
                0.291,
            ]),
            "R_values": np.array([
                1.114,
                0.990,
                1.065,
                0.952,
            ]),
            "R_sem": np.array([
                0.188,
                0.117,
                0.135,
                0.108,
            ]),
            "LR_values": np.array([
                0.659,   # drafts/26-08-11_g2_cross_mesa/main/g2_cross_summary.json 4 px
                1.110,   # drafts/26-08-11_g2_cross_mesa/supplement/g2_cross_summary.json 4 px
                0.770,   # drafts/26-08-11_g2_cross_mesa/dataset_a_main (local-floor beta)
                1.075,   # drafts/26-08-11_g2_cross_mesa/dataset_a_supplement (local-floor beta)
            ]),
            "LR_sem": np.array([
                0.265,
                0.190,
                0.179,
                0.174,
            ]),
            "RL_values": np.array([
                0.875,
                1.065,
                0.936,
                1.076,
            ]),
            "RL_sem": np.array([
                0.281,
                0.200,
                0.203,
                0.183,
            ]),
        },
    }

    # Sanity-check: supplement JSON value should be ~0.496.
    # (read at render time from drafts/26-08-06_fig2c_supplement/fig2c_science_report.json)

    return data


def render(data, out_path=OUT_PATH):
    fig, axes = plt.subplots(3, 1, figsize=(12, 14.67), constrained_layout=True)

    # ---- Panel A: population oscillation upper limits ----
    ax = axes[0]
    cats = data["panel_a"]["y_categories"]
    vals = data["panel_a"]["values"]
    corr = data["panel_a"]["corrected_values"]
    y = np.arange(len(cats))
    ax.scatter(vals, y, color="tab:blue", s=120, zorder=4, label="raw UL95")
    for i, (v, c) in enumerate(zip(vals, corr)):
        if np.isfinite(c):
            ax.scatter(c, i, color="tab:blue", s=120, facecolors="none",
                       edgecolors="tab:blue", linewidths=2, zorder=4, label="dilution-corrected" if i == 2 else "")
            ax.plot([v, c], [i, i], color="tab:blue", lw=1.2, linestyle="--", zorder=2)
    ax.axvspan(data["panel_a"]["band_min"], data["panel_a"]["band_max"],
               color="tab:orange", alpha=0.15, zorder=1)
    ax.text(
        0.95,
        0.35,
        data["panel_a"]["band_label"],
        ha="right", va="top", fontsize=9, color="darkorange", clip_on=False,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(cats)
    ax.set_xscale("log")
    ax.set_xlim(0.005, 1.0)
    ax.set_xlabel(data["panel_a"]["ylabel"])
    ax.set_title(data["panel_a"]["title"])
    ax.invert_yaxis()
    ax.grid(True, which="both", axis="x", ls=":", lw=0.5)
    ax.legend(loc="lower right", fontsize=8)

    # ---- Panel B: g2 oscillation upper limits ----
    ax = axes[1]
    cats = data["panel_b"]["y_categories"]
    vals = data["panel_b"]["values"]
    y = np.arange(len(cats))
    ax.scatter(vals, y, color="#2f4f4f", s=120, zorder=4)
    ax.axvspan(data["panel_b"]["band_min"], data["panel_b"]["band_max"],
               color="tab:orange", alpha=0.15, zorder=1)
    ax.text(
        0.5 * (data["panel_b"]["band_min"] + data["panel_b"]["band_max"]),
        0.35,
        data["panel_b"]["band_label"],
        ha="center", va="top", fontsize=9, color="darkorange", clip_on=False,
    )
    ax.text(
        0.98, 0.02,
        "limits cover only the upper half of the paper range",
        transform=ax.transAxes,
        ha="right", va="bottom", fontsize=8, color="#2f4f4f",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none"),
    )
    ax.set_yticks(y)
    ax.set_yticklabels(cats)
    ax.set_xscale("log")
    ax.set_xlim(0.1, 3.0)
    ax.set_xlabel(data["panel_b"]["ylabel"])
    ax.set_title(data["panel_b"]["title"])
    ax.invert_yaxis()
    ax.grid(True, which="both", axis="x", ls=":", lw=0.5)

    # ---- Panel C: window-pooled g2(0), full two-condensate matrix ----
    ax = axes[2]
    cats = data["panel_c"]["y_categories"]
    y = np.arange(len(cats))
    series = (
        ("L", "L = B1\u00d7B3", "o", "tab:blue", -0.27),
        ("R", "R = B2\u00d7B4", "s", "tab:red", -0.09),
        ("LR", "LR = B1\u00d7B4", "D", "tab:purple", 0.09),
        ("RL", "RL = B2\u00d7B3", "^", "tab:green", 0.27),
    )
    for key, label, marker, color, offset in series:
        ax.errorbar(y + offset, data["panel_c"][f"{key}_values"],
                    yerr=data["panel_c"][f"{key}_sem"],
                    fmt=marker, color=color, markersize=8, capsize=4,
                    label=label)
    ax.axhline(1.0, color="black", linestyle="--", lw=1.2, label="Poisson")
    ax.set_xticks(y)
    ax.set_xticklabels(cats)
    ax.set_ylabel(data["panel_c"]["ylabel"])
    ax.set_title(data["panel_c"]["title"])
    ax.set_ylim(-0.2, 2.5)
    ax.set_yticks([-0.2, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    ax.grid(True, axis="y", ls=":", lw=0.5)

    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"wrote {out_path}")


def main():
    data = build_dashboard_data()
    render(data)


if __name__ == "__main__":
    main()
