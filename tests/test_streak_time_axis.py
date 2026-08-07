import sys
import unittest
from pathlib import Path

SCRIPT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))
import streak_time_axis as sta  # noqa: E402


class TestStreakTimeAxis(unittest.TestCase):
    def test_window_per_group(self):
        self.assertEqual(sta.window_ns("02_Synchroscan_Kspace_image"), 2.0)
        self.assertEqual(sta.window_ns("06_Background_image_synchroscan_mode"), 2.0)
        self.assertEqual(sta.window_ns("01_HBT_condensate_Josephson_regime"), 1.0)
        self.assertEqual(sta.window_ns("09_Background_and_streak_pulses_delay_calibration"), 1.0)

    def test_ps_per_px(self):
        self.assertAlmostEqual(sta.ps_per_px("02_Synchroscan_Kspace_image"), 2000.0 / 1016.0)
        self.assertAlmostEqual(sta.ps_per_px("01_HBT_condensate_Josephson_regime"), 1000.0 / 1016.0)

    def test_folder_windows_cover_the_eight_mat_folders(self):
        self.assertEqual(
            set(sta.folder_windows_ps()),
            {
                "01_HBT_condensate_Josephson_regime",
                "02_Synchroscan_Kspace_image",
                "03_Synchroscan_real_space_image",
                "04_Real_space_no_pulse_intensity_calibration",
                "05_Kspace_no_pulse_intensity_calibration",
                "06_Background_image_synchroscan_mode",
                "08_Focus_cross_sections",
                "09_Background_and_streak_pulses_delay_calibration",
            },
        )

    def test_default_is_single_sweep(self):
        self.assertEqual(sta.default_window_ps(), 1000.0)

    def test_unknown_folder_raises(self):
        with self.assertRaises(KeyError):
            sta.window_ns("07_Spectrum")


if __name__ == "__main__":
    unittest.main()
