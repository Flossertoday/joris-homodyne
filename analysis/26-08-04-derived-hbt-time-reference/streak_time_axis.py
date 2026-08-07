"""Streak time-axis windows for the 26-08-04 batch, from one source of truth.

All 26-08-04 analysis scripts take the row -> time conversion from
metadata-conventions.toml [streak_camera.time_range_0804] instead of
hardcoding ps/px values: the synchroscan folders 02-06 span 2 ns (user
verification 2026-08-07), the single-sweep folders 01/08/09 span 1 ns
(26-07-31 notes), see docs/streak-time-axis.md.
"""

from pathlib import Path

import tomli

CONVENTIONS_PATH = Path(__file__).resolve().parents[2] / "metadata-conventions.toml"
FRAME_HEIGHT_PX = 1016


def _table() -> dict:
    with CONVENTIONS_PATH.open("rb") as handle:
        return tomli.load(handle)["streak_camera"]["time_range_0804"]


def window_ns(raw_folder: str) -> float:
    """Time window [ns] spanned by the full vertical frame of a raw folder."""
    table = _table()
    for group in ("synchroscan", "single_sweep"):
        if raw_folder in table[f"{group}_folders"]:
            return float(table[f"{group}_window_ns"])
    raise KeyError(f"no streak time window recorded for folder {raw_folder!r}")


def window_ps(raw_folder: str) -> float:
    return window_ns(raw_folder) * 1000.0


def ps_per_px(raw_folder: str) -> float:
    return window_ps(raw_folder) / FRAME_HEIGHT_PX


def folder_windows_ps() -> dict:
    """{raw_folder: window_ps} for every folder listed in the conventions."""
    table = _table()
    return {
        folder: float(table[f"{group}_window_ns"]) * 1000.0
        for group in ("synchroscan", "single_sweep")
        for folder in table[f"{group}_folders"]
    }


def default_window_ps() -> float:
    """Batch-default window [ps]: the single-sweep value (26-07-31 convention)."""
    return float(_table()["single_sweep_window_ns"]) * 1000.0
