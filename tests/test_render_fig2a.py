import sys
import unittest
from pathlib import Path

import numpy as np

SCRIPT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))
import render_fig2a as rf  # noqa: E402


def make_data():
    """Four synthetic band traces with a built-in 2:1 L/R brightness ratio."""
    time_ps = np.arange(100, dtype=np.float64) * 0.984
    base = np.linspace(100.0, 10.0, 100)  # shared decaying shape, all positive
    traces = np.vstack([base, 0.5 * base, 0.8 * base, 0.4 * base])
    sems = 0.05 * traces
    return {"time_ps": time_ps, "traces": traces, "sems": sems}


class BuildLrTracesNormalizeTest(unittest.TestCase):
    def test_peak_scales_larger_peak_to_one_and_preserves_ratio(self):
        data = make_data()
        _, trace_l, trace_r, _, _, factor = rf.build_lr_traces(
            data, [0, 2], [1, 3], 0.0, None, 1, "peak"
        )
        self.assertIsNone(factor)
        self.assertAlmostEqual(max(trace_l.max(), trace_r.max()), 1.0)
        # peak normalization is a common factor: the L/R sum ratio is untouched
        self.assertAlmostEqual(trace_l.sum() / trace_r.sum(), 2.0)

    def test_none_returns_raw_traces(self):
        data = make_data()
        _, trace_l, trace_r, _, _, factor = rf.build_lr_traces(
            data, [0, 2], [1, 3], 0.0, None, 1, "none"
        )
        self.assertIsNone(factor)
        raw_l = data["traces"][[0, 2]].sum(axis=0)
        np.testing.assert_allclose(trace_l, raw_l)

    def test_matched_removes_brightness_offset_then_peak_normalizes(self):
        data = make_data()
        _, trace_l, trace_r, sem_l, sem_r, factor = rf.build_lr_traces(
            data, [0, 2], [1, 3], 0.0, None, 1, "matched"
        )
        self.assertAlmostEqual(factor, 2.0)  # sum(L)/sum(R) of the raw traces
        # after matching, both arms share the same window sum ...
        self.assertAlmostEqual(trace_l.sum() / trace_r.sum(), 1.0)
        # ... and the joint peak scaling still puts the larger peak at 1
        self.assertAlmostEqual(max(trace_l.max(), trace_r.max()), 1.0)
        # SEMs get exactly the same scaling as their traces
        raw_l, raw_sem_l = rf.combine_bands(data["traces"], data["sems"], [0, 2])
        raw_r, raw_sem_r = rf.combine_bands(data["traces"], data["sems"], [1, 3])
        peak = max(raw_l.max(), 2.0 * raw_r.max())  # matched, then peak scaling
        np.testing.assert_allclose(trace_l, raw_l / peak)
        np.testing.assert_allclose(sem_l, raw_sem_l / peak)
        np.testing.assert_allclose(trace_r, 2.0 * raw_r / peak)
        np.testing.assert_allclose(sem_r, 2.0 * raw_sem_r / peak)

    def test_matched_rejects_traces_without_common_positive_bins(self):
        data = make_data()
        data["traces"][1] = -1.0  # right arm fully non-positive
        data["traces"][3] = -1.0
        with self.assertRaises(ValueError):
            rf.build_lr_traces(data, [0, 2], [1, 3], 0.0, None, 1, "matched")


if __name__ == "__main__":
    unittest.main()
