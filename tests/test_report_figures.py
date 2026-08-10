"""Smoke tests for the 26-08-10 project-report summary figure."""

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

REPORT_DIR = (
    Path(__file__).resolve().parent.parent / "analysis" / "26-08-10-project-report"
)
sys.path.insert(0, str(REPORT_DIR))
import render_report_summary as rrs  # noqa: E402


class TestRenderReportSummary(unittest.TestCase):
    def test_build_dashboard_data_structure(self):
        data = rrs.build_dashboard_data()
        self.assertIn("panel_a", data)
        self.assertIn("panel_b", data)
        self.assertIn("panel_c", data)

        for key in ("panel_a", "panel_b"):
            panel = data[key]
            self.assertIn("y_categories", panel)
            self.assertEqual(len(panel["y_categories"]), 4)
            self.assertIn("values", panel)
            self.assertEqual(len(panel["values"]), 4)
            for v in panel["values"]:
                self.assertTrue(np.isfinite(v))
                self.assertGreater(v, 0.0)

        panel_c = data["panel_c"]
        self.assertEqual(len(panel_c["y_categories"]), 4)
        for side in ("L_values", "R_values", "L_sem", "R_sem"):
            self.assertEqual(len(panel_c[side]), 4)
            for v in panel_c[side]:
                self.assertTrue(np.isfinite(v))
                self.assertGreaterEqual(v, 0.0)

        expected_labels = [
            "main batch",
            "supplement batch",
            "dataset A, main",
            "dataset A, supplement",
        ]
        self.assertEqual(data["panel_b"]["y_categories"], expected_labels)
        self.assertEqual(data["panel_c"]["y_categories"], expected_labels)

    def test_render_writes_png(self):
        data = rrs.build_dashboard_data()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "dashboard_smoke.png"
            rrs.render(data, out_path=out)
            self.assertTrue(out.exists())
            self.assertGreater(out.stat().st_size, 20 * 1024)


if __name__ == "__main__":
    unittest.main()
