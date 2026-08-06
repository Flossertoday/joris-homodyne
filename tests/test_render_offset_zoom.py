import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

SCRIPT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))
import render_offset_zoom as roz  # noqa: E402


class PermutationPvalueTest(unittest.TestCase):
    def test_perfect_positive_correlation(self):
        x = np.arange(40.0)
        r, p = roz.permutation_pvalue(x, x.copy(), n_perm=500, seed=3)
        self.assertAlmostEqual(r, 1.0)
        self.assertAlmostEqual(p, 1.0 / 501)  # no permutation can beat r=1

    def test_perfect_anticorrelation_gives_negative_r(self):
        x = np.arange(40.0)
        r, p = roz.permutation_pvalue(x, -x, n_perm=500, seed=3)
        self.assertAlmostEqual(r, -1.0)
        self.assertAlmostEqual(p, 1.0 / 501)

    def test_deterministic_for_fixed_seed(self):
        rng = np.random.default_rng(11)
        x, y = rng.normal(size=60), rng.normal(size=60)
        self.assertEqual(
            roz.permutation_pvalue(x, y, n_perm=200, seed=5),
            roz.permutation_pvalue(x, y, n_perm=200, seed=5),
        )


class MainSmokeTest(unittest.TestCase):
    """Run main() on a tiny synthetic npz; check outputs and the printed offset."""

    def test_smoke(self):
        time_ps = np.arange(200, dtype=np.float64) * 0.984
        rng = np.random.default_rng(2)
        base = np.full(200, 100.0)
        traces = np.vstack([base + rng.normal(0, 1, 200),
                            2 * base + rng.normal(0, 1, 200),
                            0.8 * base + rng.normal(0, 1, 200),
                            1.6 * base + rng.normal(0, 1, 200)])
        sems = np.ones_like(traces)
        with tempfile.TemporaryDirectory() as tmp:
            npz_path = Path(tmp) / "traces.npz"
            out_path = Path(tmp) / "zoom.png"
            np.savez(npz_path, time_ps=time_ps,
                     band_traces_photon=traces, band_sem_photon=sems)
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                roz.main([
                    "--npz", str(npz_path), "--out", str(out_path),
                    "--rebin", "1", "--zoom-min", "100",
                ])
            self.assertTrue(out_path.is_file())
            self.assertTrue(out_path.with_suffix(".pdf").is_file())
            expected_offset = (traces[0] + traces[2]).mean() - (traces[1] + traces[3]).mean()
            self.assertIn(f"offset c={expected_offset:+.5g}", stdout.getvalue())
            self.assertIn("demeaned fluctuation correlation", stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
