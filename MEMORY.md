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
- Scientific Python work must use a clone of the local `sci` environment in a non-OneDrive path. Do not record machine-specific interpreter paths or transient package versions here.
- Current raw data: `DATA/26-07-31-single-sweep/`; verified lossless analysis copy: `DATA/26-07-31-derived-single-sweep/26-07-31_single_sweep.h5`.
- Treat an HDF5 `source_root` attribute as conversion-time provenance, not a live path after a directory rename; rebind embedded `source_files/relative_path` values to the current raw directory.
- This acquisition is streak-only. Only `StreakCamImage` is analysis pixel payload; Main/Second camera arrays and Cam2 state are excluded.
- `ZeroOrderPowerRatio` is the central-dot-strength field: the two conditions are 0.05 and 0.18, each with 1000 frames. Detector axes remain uncalibrated.

## Progress and detail index

- The workspace is initialized as a Git repository on `main`; no initial commit has been created yet. `DATA/`, the pre-Git cleanup archive, caches, and local environments are excluded from version control.
- Project entry point and layout: `README.md`.
- Project working and Git rules: `AGENTS.md`.
- Staged experiment and decision gates: `docs/experiment-plan.md`.
- Claim/evidence boundaries and equations: `docs/measurement-claims.md`.
- Literature conclusions and primary sources: `docs/literature.md`.
- Data format, environment and reproducibility rules: `docs/data.md`.
- Converter and tests: `scripts/convert_matlab_v73_to_analysis_hdf5.py`, `tests/test_convert_matlab_v73_to_analysis_hdf5.py`.
- Exact derived-file integrity evidence: `DATA/26-07-31-derived-single-sweep/README.md`.
