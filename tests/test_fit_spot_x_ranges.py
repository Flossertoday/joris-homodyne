import importlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np


FIT_MODULE = "analysis.26-07-31-derived-single-sweep.fit_spot_x_ranges"
FIT_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analysis"
    / "26-07-31-derived-single-sweep"
    / "fit_spot_x_ranges.py"
)

WIDTH = 400
TRUE_BACKGROUND = 25.0
# (amplitude, mu, sigma) per spot, left to right.
TRUE_SPOTS = (
    (300.0, 60.0, 8.0),
    (500.0, 140.0, 10.0),
    (450.0, 260.0, 9.0),
    (350.0, 330.0, 12.0),
)


def load_fit_module():
    try:
        return importlib.import_module(FIT_MODULE)
    except ModuleNotFoundError as exc:
        raise AssertionError(f"{FIT_MODULE} is not implemented") from exc


def true_profile(centers):
    profile = np.full(len(centers), TRUE_BACKGROUND)
    for amplitude, mu, sigma in TRUE_SPOTS:
        profile = profile + amplitude * np.exp(
            -0.5 * ((centers - mu) / sigma) ** 2
        )
    return profile


def make_noisy_counts(seed=31):
    rng = np.random.default_rng(seed)
    centers = np.arange(WIDTH, dtype=np.float64)
    return centers, rng.poisson(true_profile(centers)).astype(np.float64)


class FitSpotXRangesTests(unittest.TestCase):
    def test_x_histogram_bins_subpixel_positions_on_pixel_grid(self):
        fit = load_fit_module()
        x = np.array([0.2, 0.49, 1.6, 5.49, 5.5, 7.49])
        centers, counts = fit.x_histogram(x, width=8)

        np.testing.assert_array_equal(centers, np.arange(8, dtype=np.float64))
        np.testing.assert_array_equal(counts, [2, 0, 1, 0, 0, 1, 1, 1])

    def test_initial_centers_find_four_ordered_maxima(self):
        fit = load_fit_module()
        centers, counts = make_noisy_counts()
        smoothed = fit.gaussian_filter1d(counts, 3.0)

        found = fit.initial_centers(smoothed)

        self.assertEqual(sorted(found), found)
        for index, (_, mu, _) in enumerate(TRUE_SPOTS):
            self.assertLess(abs(found[index] - mu), 3.0)

    def test_fit_profile_recovers_synthetic_four_gaussians(self):
        fit = load_fit_module()
        centers, counts = make_noisy_counts()

        popt, pcov, mu0 = fit.fit_profile(centers, counts)

        self.assertEqual(13, len(popt))
        self.assertLess(abs(popt[0] - TRUE_BACKGROUND), 5.0)
        for index, (_, mu, sigma) in enumerate(TRUE_SPOTS):
            base = 1 + 3 * index
            self.assertLess(abs(popt[base + 1] - mu), 1.5, f"mu {index}")
            self.assertLess(abs(popt[base + 2] - sigma), 1.5, f"sigma {index}")
        # The fitted order must be the physical left-to-right order.
        mus = [popt[1 + 3 * i + 1] for i in range(4)]
        self.assertEqual(sorted(mus), mus)
        self.assertTrue(np.all(np.diag(pcov) >= 0.0))
        self.assertEqual(4, len(mu0))

    def test_spot_intervals_are_mu_plus_minus_k_sigma(self):
        fit = load_fit_module()
        popt = np.array(
            [10.0, 100.0, 50.0, 5.0, 200.0, 150.0, 10.0,
             300.0, 250.0, 8.0, 400.0, 350.0, 6.0]
        )

        intervals = fit.spot_intervals(popt, 3.0)

        np.testing.assert_allclose(
            intervals,
            [[35.0, 65.0], [120.0, 180.0], [226.0, 274.0], [332.0, 368.0]],
        )

    def test_acceptance_mask_is_union_of_intervals(self):
        fit = load_fit_module()
        intervals = np.array([[10.0, 20.0], [30.0, 40.0]])
        x = np.array([5.0, 10.0, 15.0, 25.0, 35.0, 45.0])

        mask = fit.acceptance_mask(x, intervals)

        np.testing.assert_array_equal(
            mask, [False, True, True, False, True, False]
        )

    def test_cli_writes_figure_csv_and_filtered_peaks(self):
        fit = load_fit_module()
        rng = np.random.default_rng(37)
        # Draw peak x positions from the true mixture plus a flat background,
        # so the CLI exercises the same code path as the real dataset.
        n_signal = 6000
        n_background = 1500
        component = rng.choice(4, size=n_signal, p=[0.2, 0.3, 0.28, 0.22])
        x_signal = np.array(
            [
                rng.normal(TRUE_SPOTS[c][1], TRUE_SPOTS[c][2])
                for c in component
            ]
        )
        x = np.concatenate(
            [x_signal, rng.uniform(0, WIDTH, n_background)]
        )
        x = x[(x >= 0) & (x < WIDTH)]
        n = len(x)

        with tempfile.TemporaryDirectory() as tmp:
            peaks_path = Path(tmp) / "peaks.h5"
            figure_path = Path(tmp) / "fit.png"
            csv_path = Path(tmp) / "fit.csv"
            filtered_path = Path(tmp) / "filtered.h5"
            with h5py.File(peaks_path, "w") as handle:
                handle.create_dataset("frame", data=rng.integers(0, 10, n))
                handle.create_dataset("y", data=rng.uniform(0, 50, n))
                handle.create_dataset("x", data=x)
                handle.create_dataset("mass", data=rng.uniform(300, 3000, n))
                handle.attrs["condition_name"] = "central_dot_0p18"
                handle.attrs["central_dot_strength"] = 0.18
                handle.attrs["frame_shape"] = np.array([50, WIDTH])

            result = subprocess.run(
                [
                    sys.executable,
                    str(FIT_SCRIPT),
                    str(peaks_path),
                    "--figure",
                    str(figure_path),
                    "--fit-csv",
                    str(csv_path),
                    "--filtered-h5",
                    str(filtered_path),
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(0, result.returncode, result.stderr)
            self.assertEqual(b"\x89PNG\r\n\x1a\n", figure_path.read_bytes()[:8])

            csv_lines = csv_path.read_text().strip().splitlines()
            self.assertEqual(
                "component,label,amplitude,amplitude_std,mu,mu_std,"
                "sigma,sigma_std,range_lo,range_hi",
                csv_lines[0],
            )
            self.assertEqual(6, len(csv_lines))  # header + background + 4 spots

            with h5py.File(filtered_path, "r") as handle:
                filtered_x = np.asarray(handle["x"])
                intervals = np.asarray(handle.attrs["x_range_intervals"])
                self.assertEqual(
                    len(filtered_x), int(handle.attrs["peak_count"])
                )
                self.assertEqual(n, int(handle.attrs["source_peak_count"]))
            self.assertEqual((4, 2), intervals.shape)
            self.assertLess(len(filtered_x), n)
            mask = fit.acceptance_mask(filtered_x, intervals)
            self.assertTrue(bool(mask.all()))
            # The conservative 3-sigma union should keep most signal peaks.
            self.assertGreater(len(filtered_x) / n, 0.6)


if __name__ == "__main__":
    unittest.main()
