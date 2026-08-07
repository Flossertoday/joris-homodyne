"""Tests for kspace_fringe_priors and hbt_z_targeted_fits."""
import sys
import unittest
from pathlib import Path

import numpy as np

SCRIPT_DIR = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))

import hbt_z_targeted_fits as hz  # noqa: E402
from kspace_fringe_priors import RESULTS_JSON, measured_periods_ps  # noqa: E402


class TestKspaceFringePriors(unittest.TestCase):
    def test_results_json_exists_and_carries_2ns_priors(self):
        self.assertTrue(RESULTS_JSON.exists(), f"missing {RESULTS_JSON}")
        self.assertEqual((74.4, 83.7), measured_periods_ps())

    def test_periods_are_positive(self):
        self.assertTrue(all(p > 0 for p in measured_periods_ps()))


class TestTargetedFit(unittest.TestCase):
    def test_recovers_injected_amplitude(self):
        rng = np.random.default_rng(5)
        tt = np.arange(20.0, 800.0, 0.984)
        amp, phase = 0.05, 0.7
        sz = np.full_like(tt, 0.01)
        z = amp * np.cos(2 * np.pi * tt / 74.4 + phase) + rng.normal(0, 0.01, tt.size)
        r = hz.targeted_fit(tt, z, sz, 74.4)
        self.assertAlmostEqual(amp, r["amplitude"], delta=0.01)
        self.assertGreater(r["sigma"], 0.0)
        self.assertAlmostEqual(r["amplitude"] + 2.0 * r["sigma"], r["ul95"])

    def test_null_amplitude_within_noise(self):
        rng = np.random.default_rng(7)
        tt = np.arange(20.0, 800.0, 0.984)
        sz = np.full_like(tt, 0.01)
        z = rng.normal(0, 0.01, tt.size)
        r = hz.targeted_fit(tt, z, sz, 74.4)
        self.assertLess(r["amplitude"], 3.0 * r["sigma"])


if __name__ == "__main__":
    unittest.main()
