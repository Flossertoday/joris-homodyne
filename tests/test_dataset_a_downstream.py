import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

SCRIPT_DIR = (
    Path(__file__).resolve().parent.parent
    / "analysis"
    / "26-08-04-derived-hbt-time-reference"
)
sys.path.insert(0, str(SCRIPT_DIR))
import dataset_a_g2 as dag2  # noqa: E402
import dataset_a_traces as datr  # noqa: E402
import fig2c_g2_estimator as est  # noqa: E402


def synthetic_dataset_a(path, n_frames=8, seed=6):
    """Minimal dataset-A-schema h5: Poisson events in 4 bands, one t0 failure."""
    rng = np.random.default_rng(seed)
    bands = np.array([[380, 450], [570, 680], [810, 880], [1010, 1110]])
    frame_t0 = np.full(n_frames, 155.0)
    frame_t0[3] = np.nan  # t0 failure
    records = {k: [] for k in ("frame", "x", "tau", "band")}
    rates = [0.5, 0.8, 0.4, 0.9]  # events/frame/band over the tau grid
    for frame in range(n_frames):
        if frame == 3:
            continue
        for band in range(4):
            n = rng.poisson(rates[band])
            xs = rng.uniform(bands[band, 0], bands[band, 1], size=n)
            taus = rng.uniform(-60, 900, size=n)
            records["frame"].append(np.full(n, frame, dtype=np.int32))
            records["x"].append(xs)
            records["tau"].append(taus)
            records["band"].append(np.full(n, band, dtype=np.int8))
    events = {k: np.concatenate(v) for k, v in records.items()}
    with h5py.File(path, "w") as handle:
        group = handle.create_group("events")
        for name, data in events.items():
            group.create_dataset(name, data=data)
        group.create_dataset("y", data=events["tau"] + 155.0)
        handle.create_dataset("frame_t0", data=frame_t0)
        handle.create_dataset("frame_indices", data=np.arange(n_frames))
        handle.create_dataset("failed_frames", data=np.array([3], dtype=np.int32))
        nm = handle.create_group("noise_model")
        nm.create_dataset("floor_per_band", data=np.full(4, 1e-4))
        nm.create_dataset("rho_tau_global",
                          data=np.full(1144, 1e-4 / 450.0))
        nm.create_dataset("tau_grid", data=np.arange(-64, 1080, dtype=np.float64))
        handle.attrs["band_windows"] = bands
        handle.attrs["n_frames_kept"] = n_frames - 1
    return path, events, rates


class G2HelpersTest(unittest.TestCase):
    def test_beta_from_noise_model(self):
        floor = np.array([1e-3, 2e-3, 1.5e-3, 2.5e-3])
        windows = np.array([[380, 450], [570, 680], [810, 880], [1010, 1110]],
                           dtype=np.float64)
        beta = dag2.beta_from_noise_model(floor, windows, nbins=10, bin_width=4)
        self.assertEqual(beta.shape, (4, 10))
        # per-row density x inclusive member columns x bin width
        expected_b1 = 1e-3 / 1016 * 71 * 4
        self.assertAlmostEqual(beta[0, 0], expected_b1)
        self.assertTrue(np.allclose(beta[:, 0], beta[:, 5]))

    def test_build_cube_excludes_failed_frames(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, events, rates = synthetic_dataset_a(Path(tmp) / "a.h5")
            data = dag2.load_dataset_a(path)
            cube, nbins, frames = dag2.build_cube(
                data["events"], data["frame_t0"], data["frame_indices"], 4)
            self.assertEqual(frames, 7)  # one t0 failure excluded
            self.assertEqual(cube.shape[0], 7)
            self.assertEqual(int(cube.sum()), events["frame"].size)
            # per-band totals match the planted band index
            for band in range(4):
                self.assertEqual(
                    int(cube[:, band, :].sum()),
                    int((events["band"] == band).sum()),
                )


class TracesTest(unittest.TestCase):
    def test_traces_match_planted_rates(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, events, rates = synthetic_dataset_a(
                Path(tmp) / "a.h5", n_frames=2000, seed=9)
            result = datr.build_traces(path)
            t = result["time_ps"]
            window = (t >= 0) & (t <= 820)
            for band in range(4):
                measured = np.nanmean(result["band_traces_photon"][band][window])
                # planted rate uniform over tau in [-60, 900] px
                expected = rates[band] / 960.0
                self.assertAlmostEqual(measured, expected, delta=0.25 * expected)
            # band windows propagated
            self.assertEqual(result["band_windows"].shape, (4, 2))
            # SEM non-negative, with positive median in the window (sparse
            # synthetic bins can be exactly zero)
            sem = result["band_sem_photon"][:, window]
            self.assertTrue(np.all(sem >= 0))
            self.assertGreater(np.median(sem), 0)


if __name__ == "__main__":
    unittest.main()
