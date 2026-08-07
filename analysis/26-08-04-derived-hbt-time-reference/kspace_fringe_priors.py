"""Single source for the measured k-space fringe oscillation periods.

`measure_kspace_fringe_oscillation.py` writes the tracked JSON
`results/kspace_fringe_oscillation.json`; analyses that need the measured
oscillation period as a prior must read it from there instead of hardcoding
(the 2026-08-07 window correction 37.2/41.8 -> 74.4/83.7 ps is the
cautionary tale: the periods scale with the synchroscan time window).
"""
import json
from pathlib import Path

RESULTS_JSON = Path(__file__).resolve().parent / "results" / "kspace_fringe_oscillation.json"


def measured_periods_ps(path=None):
    """(primary, secondary) fringe periods [ps] from the tracked results JSON."""
    results_path = Path(path) if path is not None else RESULTS_JSON
    with open(results_path, encoding="utf-8") as handle:
        data = json.load(handle)
    return tuple(float(v) for v in data["prior_periods_ps"])
