# Project memory

## Purpose and current direction

- This project studies nonresonantly pumped, optical-Stark-triggered polariton Josephson dynamics and the measurements that could support later quantum-statistical claims.
- The immediate gate is whether the previously accumulated left/right intensity oscillation survives single-sweep acquisition as a statistically resolvable per-kick trajectory.
- After that gate, prioritize time-resolved relative phase and the full \(g^{(2)}_{LL},g^{(2)}_{RR},g^{(2)}_{LR},g^{(2)}_{RL}\) matrix.
- CSI is optional post-processing of a calibrated full correlation dataset, not the project objective. Direct single-mode quadrature squeezing would require a phase-related LO and balanced homodyne chain that the current apparatus does not have.

## Stable interpretation rules

- Keep shot-resolved Josephson dynamics, inter-mode intensity correlation, CSI, single-mode quadrature squeezing, and two-mode squeezing/entanglement as separate claims.
- Photon bunching alone is not evidence of quadrature squeezing. A CSI violation rules out a positive-\(P\) classical-field model under the measurement assumptions, but does not by itself prove squeezing or entanglement; non-violation does not refute squeezing.
- One single-sweep frame is one realization, not a \(g^{(2)}\) matrix. Normally ordered correlations require many sparse photon-counting trials with same-event identity and independent HBT arms.
- Treat Stark-triggered oscillations conservatively as population dynamics consistent with Josephson exchange until counterphase response, weak common-mode response, parameter controls and relative-phase measurements exclude reservoir/detection transients.
- For the 2017 Adiyatullin paper, use the Zotero attachment named `Full Text PDF`. The experiment measured time-resolved intensity correlations; its squeezing curve was model-derived rather than a direct quadrature measurement.

## Data rules and current dataset

- Raw experimental files under `DATA/` are immutable and excluded from Git. Name raw directories `YY-MM-DD-<experiment>` and their distinct sibling analysis directories `YY-MM-DD-derived-<experiment>` rather than using a shared generic output folder.
- Scientific Python work must use a clone of the local `sci` environment in a non-OneDrive path. The checked-in Codex Local environment invokes `scripts/init_sandbox.ps1` when active, while `AGENTS.md` makes the same script the agent-owned fallback; live paths are recorded under ignored `tmp/`. Do not record machine-specific interpreter paths or transient package versions here.
- Current raw data: `DATA/26-07-31-single-sweep/`; verified lossless analysis copy: `DATA/26-07-31-derived-single-sweep/26-07-31_single_sweep.h5`.
- Treat an HDF5 `source_root` attribute as conversion-time provenance, not a live path after a directory rename; rebind embedded `source_files/relative_path` values to the current raw directory.
- This acquisition is streak-only. Only `StreakCamImage` is analysis pixel payload; Main/Second camera arrays and Cam2 state are excluded.
- `ZeroOrderPowerRatio` is the central-dot-strength field: the two conditions are 0.05 and 0.18, each with 1000 frames. The vertical frame axis is the streak time axis (26-07-31 window: 1 ns over 1016 px; `StreakRange` device setting 5 per the acquisition README); the horizontal axis stays uncalibrated detector pixels. Agreed metadata mappings live in the root `metadata-conventions.toml`.
- The single derived HDF5 is a multi-condition container: raw acquisition blocks 2 and 3 map to `central_dot_0p05` and `central_dot_0p18`, respectively; it is not one HDF5 per measurement condition.
- The 0.18 condition's element-wise sum over 1000 streak frames shows four bright bands. Its script, display rule and reviewed PNG live under the matching `analysis/26-07-31-derived-single-sweep/` directory.

## Progress and detail index

- The repository uses `main` as its stable integration branch. `DATA/`, the pre-Git cleanup archive, caches, and local environments are excluded from version control.
- Project entry point and layout: `README.md`.
- Project working and Git rules: `AGENTS.md`.
- Staged experiment and decision gates: `docs/experiment-plan.md`.
- Claim/evidence boundaries and equations: `docs/measurement-claims.md`.
- Literature conclusions and primary sources: `docs/literature.md`.
- Data format, environment and reproducibility rules: `docs/data.md`.
- Streak vertical-axis time-window investigation and axis conventions: `docs/streak-time-axis.md`, agreed quantities `metadata-conventions.toml`.
- Managed-sandbox initialization observations and requirements: `docs/sandbox-initialization.md`.
- Automatic sandbox environment, initializer and dependency contract: `.codex/environments/environment.toml`, `scripts/init_sandbox.ps1`, `requirements.txt`, `tests/test_init_sandbox.py`.
- Converter and tests: `scripts/convert_matlab_v73_to_analysis_hdf5.py`, `tests/test_convert_matlab_v73_to_analysis_hdf5.py`.
- Central-dot frame-sum analysis: `analysis/26-07-31-derived-single-sweep/README.md`, `analysis/26-07-31-derived-single-sweep/sum_streak_condition.py`, `analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_sum.png`.
- Single-frame sparse peak detection (trackpy Crocker-Grier): `analysis/26-07-31-derived-single-sweep/find_frame_peaks.py`, `tests/test_find_frame_peaks.py`, reviewed overlay `analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_frame0230_peaks.png`.
- Peak-localization methodology audit (pipeline visualization + 2D Gaussian fit cross-check, 50 frames seed 20260803): `analysis/26-07-31-derived-single-sweep/visualize_peak_pipeline.py`, `analysis/26-07-31-derived-single-sweep/compare_gauss_refinement.py`, write-up `docs/peak-localization.md`.
- Full 0.18-condition peak extraction (30542 peaks, 1000 frames; downstream raw data): `analysis/26-07-31-derived-single-sweep/extract_condition_peaks.py`, dataset + provenance `analysis/26-07-31-derived-single-sweep/detected-peaks/26-07-31_central_dot_0p18_peaks.h5`, `tests/test_extract_condition_peaks.py`.
- Four-Gaussian x-profile fit and x-acceptance ranges (auto-shrunk from 3σ to disjoint 2σ because same-arm spot intervals overlap and overlap assignment is unsolved; 15958/30542 peaks kept; spots left-to-right = A·armA, B·armA, A·armB, B·armB): `analysis/26-07-31-derived-single-sweep/fit_spot_x_ranges.py`, filtered dataset `analysis/26-07-31-derived-single-sweep/detected-peaks/26-07-31_central_dot_0p18_peaks_xrange2p0sigma.h5`, `tests/test_fit_spot_x_ranges.py`.
- Detector-grid density overview plus per-frame scatter panels of the full vs x-range-filtered peak tables (both PNGs live inside `detected-peaks/`; the `frame` column covers all 1000 frames in both tables; kept/rejected split imports `acceptance_mask` from the production fit): `analysis/26-07-31-derived-single-sweep/visualize_detected_peaks.py`, figures `analysis/26-07-31-derived-single-sweep/detected-peaks/26-07-31_central_dot_0p18_detected_peaks_density.png` and `analysis/26-07-31-derived-single-sweep/detected-peaks/26-07-31_central_dot_0p18_detected_peaks_frames.png`, `tests/test_visualize_detected_peaks.py`.
- g2 correlation matrix replicating Adiyatullin 2017 Suppl. Note 1 (AA/BB cross-arm, AB cross-mesa, background-corrected; bin widths 1/3/6 px; pooled g2(0) consistent with Poisson within ~1σ at 1000 frames — null result, statistics-limited): `analysis/26-07-31-derived-single-sweep/compute_g2_matrix.py`, results HDF5 `DATA/26-07-31-derived-single-sweep/g2/26-07-31_central_dot_0p18_g2.h5`, figures `analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_g2_*_bin{1,3,6}px.png`, `tests/test_compute_g2_matrix.py`.
- Exact derived-file integrity evidence: `DATA/26-07-31-derived-single-sweep/README.md`.
