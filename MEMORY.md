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

## Data rules

- Raw experimental files under `DATA/` are immutable and excluded from Git. Name raw directories `YY-MM-DD-<experiment>` and their distinct sibling analysis directories `YY-MM-DD-derived-<experiment>`.
- Treat an HDF5 `source_root` attribute as conversion-time provenance, not a live path after a directory rename; rebind embedded `source_files/relative_path` values to the current raw directory.
- Both acquisitions are streak-only: only `StreakCamImage` is analysis pixel payload; Main/Second camera arrays and Cam2 state are excluded.
- The vertical frame axis is the streak time axis (1 ns over 1016 px, `StreakRange` 5); the horizontal axis stays uncalibrated detector pixels. Agreed metadata mappings live in `metadata-conventions.toml`; strict in-data time calibration is deferred to the 26-08-04 time-reference follow-up.
- Exact integrity evidence (sizes, hashes, pixel-exact verification) lives in each derived directory's `README.md`, not here.
- Scientific Python work uses the isolated `.venv` per `AGENTS.md` and `docs/data.md`; live paths stay in ignored `tmp/`.

## Datasets

- `DATA/26-07-31-single-sweep/` → `DATA/26-07-31-derived-single-sweep/26-07-31_single_sweep.h5`: two central-dot conditions by `ZeroOrderPowerRatio` (0.05, 0.18), 1000 frames each; multi-condition container with raw blocks 2 and 3 mapped to `central_dot_0p05`/`central_dot_0p18`. Analysis lives in `analysis/26-07-31-derived-single-sweep/`.
- `DATA/26-08-04-hbt-time-reference/` (renamed from export `26-08-04_dataSaved`; pristine zip keeps the original name) → `DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5`: schema_version 2, `--condition-mode directory` (all files share `ZeroOrderPowerRatio` 0.26). Nine sub-experiment conditions; main one is `hbt_condensate_josephson_regime` (20000 single-sweep HBT frames with the laser time reference on the side), the rest are 10-frame calibration sets (synchroscan k/real space, no-pulse intensity calibrations, backgrounds, focus cross-sections, delay calibration with preserved measured background). `07_Spectrum` has no MAT payload; acquisition notes stored as `metadata/readme`. Pixel-exact verification passed 2026-08-05.

## Progress and detail index

- Entry point and layout: `README.md`; working and Git rules: `AGENTS.md`; experiment gates: `docs/experiment-plan.md`; claim/evidence boundaries: `docs/measurement-claims.md`; literature: `docs/literature.md`; data/environment rules: `docs/data.md`; streak axis: `docs/streak-time-axis.md` + `metadata-conventions.toml`; sandbox: `docs/sandbox-initialization.md`.
- Converter: `scripts/convert_matlab_v73_to_analysis_hdf5.py` (+ tests). Sandbox contract: `scripts/init_sandbox.ps1`, `requirements.txt`.
- 26-07-31 analysis chain (all under `analysis/26-07-31-derived-single-sweep/`, figures in `results/` and `detected-peaks/`): frame-sum (0.18 condition shows four bright bands) → trackpy peak detection → pipeline audit + Gaussian cross-check (`docs/peak-localization.md`) → full extraction (30542 peaks) → four-Gaussian x-fit with disjoint 2σ acceptance (15958 kept; spots left-to-right A·armA, B·armA, A·armB, B·armB) → localization error (0.16–0.17 px RMS/axis, dominated by deterministic subpixel window-truncation bias; `docs/peak-localization-error.md`) → g2 matrix replicating Adiyatullin 2017 Suppl. Note 1: all pooled g2(0) consistent with Poisson within ~1.5σ at 1000 frames — null result, statistics-limited; same-mesa values hover +0.05–0.08, not a detection. **26-07-31 g2 analysis is closed; time-resolved correlation work moves to the pulse-calibrated 26-08-04 dataset.**
- 26-08-04 progress: verified derived HDF5 and human-browsable preview tree (`analysis/26-08-04-derived-hbt-time-reference/preview/`, regenerable via `render_preview.py`). Fig. 2(a) replica pipeline on the 20000-frame HBT condition is complete (CFD50 t0, band pairing L={B1,B3}/R={B2,B4}, accumulated traces, paper-style render, sensitivity validation; `docs/26-08-04-fig2a.md`): **reproduction closed 2026-08-06, unsuccessful** — no paper-style antiphase population alternation in the averaged single-sweep data, 95% upper limit ~1–2% intensity modulation, not a pipeline-sensitivity artifact. A later retry on the 26-07-31 dataset is planned but not scheduled. Fig. 2(c) time-resolved g2(0)(t) completed 2026-08-06 as an **upper-limit report** (`docs/26-08-04-fig2c.md`): the G3 power gate showed the per-bin regime is statistics-limited at any binning (≈50% of detected events are residual detector-noise clusters; per-bin occupancy ~1e-3/frame), so the deliverable is the pooled measurement (g2(0) = 0.31±0.38 L / 1.11±0.19 R, both Poisson-consistent), reference-locked oscillation 95% upper limits (median 0.81; 1.09 at 9.84 ps), and the G7 injection validation (amplitude bias <10% at strong injection; paper-strength signals only partially resolvable — matching the gate). All four pre-registered science tests null. Event table (275218 events) and count cubes archived in `drafts/26-08-06_fig2c/` for the later 2-D g2/cross-mesa entries; ~3e5 frames would be needed for 95% sensitivity to the paper's lower amplitude edge (0.2).
