# KM curve reconstruction study

Reproduce the manuscript tables and figures from the study's saved IPD with `make reproduce`. Model reconstruction is not rerun. An optional exact-coordinate Guyot replay has a documented numerical difference; it is separate from the validated frozen-IPD analysis.

This repository reproduces the study analyses from the original simulated data and the IPD actually produced during the study. The default pipeline requires no Claude/Anthropic calls, API keys or paid services.

## Data and verification

- Copied original scientific payloads without changing their contents: 68 reference datasets, 360 actual assigned figures and 543 completed reconstructed IPD files. `data/source_checksums.csv` records original repository-relative source paths, destination paths or ZIP members, sizes and SHA-256 values. Path-redacted administrative artifacts are identified separately as described below.
- `data/run_manifest.csv` retains all 550 in-scope registrations, including seven without completed IPD. It links batch, scenario, dataset, sample size, input format, condition, run ID, public manifest, code seed, original/recovery relationship, saved output, archived program, replay receipt and analysis code.
- Generated all three main tables, eleven supplementary tables and six figures. The default calculations use frozen IPD. The verifier checks all 528 scored study reconstructions and 986 numerical cells at the manuscript's displayed precision. The fifteen completed preliminary-screen IPDs are retained separately; the preliminary screen reports completion, not a new pooled accuracy analysis.
- Recomputed 48 direct first–repeat comparisons across 24 datasets. Repeats are not independent primary datasets.
- Verified saved source anchors, extraction classifications and 400 historical byte-identical LLM replay receipts. These receipt checks are not fresh executions of the model-generated programs.

The default analysis has ten maintained scripts. Two additional optional scripts separately replay automated PNG digitization and Guyot reconstruction; these responsibilities justify twelve files rather than combining different methods into one large script. Every script states its purpose, inputs, outputs and command.

## Scientific discrepancy: optional exact-coordinate replay

The optional replay uses the unchanged archived Guyot numerical algorithm. All 64 PNG-derived Guyot outputs were reproduced byte for byte. For the exact-coordinate experiment, 58 of 64 CSV files differed in final decimal digits; the largest observed-time difference was approximately `1.03e-13` months. Although small as time differences, these changes affect the strict maximum difference between right-continuous KM curves:

| Dataset | Metric | Frozen IPD result (probability) | Replay result (probability) |
|---|---|---:|---:|
| PP32_S1_D02 | maximum_A | 0.0139185276089289 | 0.0140203547584752 |
| PP32_S2_D03 | maximum_A | 0.0127620891704286 | 0.0126858485937197 |
| PP32_S3_D06 | maximum_A | 0.00154394433040783 | 0.00154298217589044 |

Across all 64 exact-coordinate datasets, mean maximum_A changes from **0.7355173683240882 pp** to **0.7355558439777132 pp** (displayed to four decimals: **0.7355** versus **0.7356**). Three of 832 rescored metrics differ by more than `1e-10` on their native scales. The cause has not been established. The historical exact-coordinate experiment used a different R environment; that is a candidate explanation, not a demonstrated cause.

Running `make guyot-replay` writes `build/guyot_replay/scoring_comparison.csv`, `validation.csv` and `summary.json`. It compares all 13 scores for every serialized reconstruction, in addition to checking records and event times, and **exits nonzero if a score differs by more than `1e-10`**. This deliberately detects the known difference instead of certifying equivalence from a small time-distance tolerance. Replayed files never replace frozen IPD.

The default frozen-IPD analysis reproduces the current, CSV-rescored S8 values. No data, event times, tolerances or reported results were changed to make the optional replay match. The manuscript reproduction remains based on the actual study outputs; bitwise regeneration of all exact-coordinate Guyot outputs is not claimed.

## Files

```text
README.md
Makefile                           # One command: make reproduce
.gitignore
environment.yml                    # Original pinned author environment
environment/conda-osx-arm64.lock   # Exact macOS ARM package lock
scripts/                           # 10 default scripts + 2 optional replay scripts
data/
  datasets.csv
  run_manifest.csv
  source_checksums.csv
  P4_prompt.txt
  figure_characteristic_analysis_plan.json
  reference_ipd/                   # Original CSV bytes, including preliminary screen
  input_figures/                   # Actual model inputs, not replacement renderings
  guyot_digitization/              # Actual coordinates/risk tables given to Guyot
  input_generation.zip             # Archived original input-generation sources
  guyot_workflow.zip               # Original digitizer, algorithm, driver and protocol
frozen_results/
  llm_ipd/
  guyot_ipd/
  reconstruction_programs.zip      # Archived programs, prompts, checks and audit evidence
  validation/                     # Saved unrounded authorities; verifier only
build/                             # Locally regenerated files; excluded from Git
```

ZIP files are research artifacts, not maintained helper modules. The model-generated programs, P4 and run-specific prompts are retained verbatim. Unreported alternative sensitivity CSVs, redundant diagnostic overlays and plain execution logs are omitted (490 generated artifacts); source programs, audit evidence and the actual model inputs remain available. `source_root=old_repository` denotes the historical research repository; `numeric_validation_reference` denotes numeric table tokens and their display precision, used only for validation. All source locations in the checksum ledger are relative to these named roots.

For public release, personal absolute directory prefixes in 207 archived JSON records were replaced with `old_repository` or `<LOCAL_HOME>`. One archived figure-generation helper now defaults to `Rscript` instead of a personal executable path; its calculations are unchanged. These 208 members are explicitly labeled public derivatives, not byte-identical originals. No IPD, assigned figure, model-generated program, prompt, numerical field, or historical hash was changed by this cleanup. The unredacted administrative records are excluded from the release.

In `data/source_checksums.csv`, `sha256` and `bytes` describe the distributed file or ZIP member. For a public derivative, `original_sha256` and `original_bytes` preserve the original fingerprint, and `transformation` records `private-path-redaction-v1` or `portable-rscript-path-v1`. Those three fields are blank for unchanged originals. Historical hashes inside archived admissions and `LOCKED.json` records still describe their original, unredacted targets; compare them with `original_sha256` when a target was redacted. They have not been rewritten to imply that the sanitized record existed during the study. The default verifier checks distributed payloads and restricts these transformations to administrative JSON and the single generator default; it never requires access to the old repository.

## Populations and run manifest

`primary64` selects one completed PDF and PNG300 run for each current primary dataset. `complete_case62` selects the earlier complete pairs after the startup recovery R003R1, excluding the two later PNG recoveries R014R1 and R064R1. `original61` also excludes the dataset whose PDF required R003R1. These are **pair inclusion flags**: they must not be used to erase the observed partner of a missing run when computing missing-output bounds.

`is_first_run`, `attempt`, `repeatability` and `pairing_key` keep the first and repeat reconstructions together. `assignment_status`, `original_run_id`, `replaces_run_id` and `replaced_by_run_id` retain recovery history. An empty `reconstruction_output` means no completed output exists, not zero error. The eight uniformly superseded early attempts are documented in the archived restart evidence and are outside the 550 current registrations.

The manifest is one row per registration, not one row per independent dataset:

| Field group | Meaning |
|---|---|
| `batch`, `scenario`, `dataset_id`, `n_per_arm` | Study cohort and original simulated dataset; scenario is S1–S4. |
| `figure_format`, `condition`, `experiment`, `method` | File format, treatment of the input figure, analysis family and reconstruction method. `exact_coordinates` is a separate information condition. |
| `run_id`, `effort`, `code_seed`, `attempt` | Recorded run identity, model effort, supplied program seed and repeated-run assignment. Blank historical values are not inferred. |
| `is_first_run`, `repeatability`, `pairing_key` | Links the first and repeated reconstructions within dataset and condition. |
| `complete`, `assignment_status`, `original_run_id`, `replaces_run_id`, `replaced_by_run_id` | Whether a saved IPD exists and whether this is an original or recovered assignment. |
| `primary64`, `complete_case62`, `original61` | Explicit pair-level population membership. |
| `reference_ipd`, `input_figure`, `reconstruction_output` | Repository-relative original reference data, actual supplied figure and saved reconstructed IPD. For exact-coordinate Guyot, `input_figure` identifies the corresponding study figure; its actual inputs are the step coordinates/risk counts specified below. |
| `program_archive`, `program_directory`, `prompt_member`, `public_manifest_member`, `replay_receipt_member` | Archive filename, source directory and ZIP-member paths for the actual program, prompt, supplied public metadata and saved replay evidence. |
| `analysis_code`, `source_assignment` | Maintained scorer and historical assignment path. |

There are 181 primary/repeat registrations (176 complete), 160 modified-figure registrations, 32 image-only PDF registrations, 33 plotting-style registrations (32 complete), 128 Guyot reconstructions and 16 preliminary effort-screen registrations (15 complete). The preliminary screen has four reference datasets; its absent generation seeds remain blank in `datasets.csv`. Completion findings remain reproducible from the original assignments and saved outputs.

The Guyot PNG run uses the actual saved digitization from its PNG. Exact-coordinate Guyot receives exact event-step coordinates and the displayed risk counts, but no censoring marks or total event count. Reference IPD is never passed to the PNG digitizer or PNG Guyot reconstruction.

## Statistical settings

TAAE integrates right-continuous KM step differences through 24 months and averages the two arm errors. Events and censors tied at a time use the same risk set, with the events applied before censor removal. No nonzero survival tail is extrapolated. Maximum error includes the endpoint at 24 months. The original duration threshold is `0.02 + 32 * .Machine$double.eps`.

Cox models use full follow-up, Efron ties and `timefix=FALSE`, with the original optimization settings. Raw curve errors are probabilities; tables and figures explicitly multiply by 100 when reporting percentage points. Relative standard-error errors remain fractions unless a display explicitly uses percent.

All bootstrap comparisons use 10,000 paired whole-dataset resamples and type-1 percentile endpoints. Current64 uses four scenarios with 16 datasets each. Historical missing-output analyses retain eight batch-by-scenario strata and eight assigned positions per stratum. Follow-up experiments use four scenarios with eight datasets; survminer uses four with four. Revision-limit sensitivity retains its original eight draws per available batch/scenario stratum, including strata with fewer retained datasets. Original R and NumPy RNG families, seeds, row ordering and contrast-specific seed offsets remain in the scripts. The primary descriptive IQR uses the original exclusive/type-6 convention; figure-condition summaries use their original NumPy/type-7 convention.

## Requirements and running the analysis

Validated on macOS Apple Silicon with Python 3.12.14 and R 4.5.2 / survival 3.8-3. The supplied lock and `environment.yml` are the original specifications. On 2026-10-16, all 212 locked package SHA-256 hashes were checked, a new environment was installed offline from those exact package archives, and `make reproduce` passed in an independent repository copy starting with no `build/` directory. This validates the environment and local pipeline; fresh downloads from the package server and other OS/CPU combinations have not been tested. No API key is needed.

To install the recorded environment into a separate environment name:

```sh
conda create --name km-curve-study --file environment/conda-osx-arm64.lock
conda activate km-curve-study
```

From this repository root, run the full analysis and validation:

```sh
make reproduce
```

`make verify` rechecks existing outputs and input hashes. `make test` runs the analytic KM/tie/endpoint checks. You can select explicit executables with `make reproduce PYTHON=/path/to/python RSCRIPT=/path/to/Rscript`. The default target runs the following steps sequentially:

```sh
export PYTHONDONTWRITEBYTECODE=1
python scripts/verify_reproduction.py --inputs
Rscript --vanilla scripts/score_survival_reconstructions.R --self-test
Rscript --vanilla scripts/score_survival_reconstructions.R
python scripts/calculate_uncertainty_and_events.py
Rscript --vanilla scripts/compare_formats_and_guyot.R
python scripts/analyze_exploratory_comparisons.py
python scripts/audit_reconstruction_methods.py
python scripts/write_manuscript_tables.py
python scripts/plot_study_design.py
python scripts/plot_reconstruction_example.py
python scripts/plot_accuracy_comparisons.py
python scripts/verify_reproduction.py --results
```

All six final PDFs were rendered and visually inspected. Read `build/tables/tables.html` for all tables, `build/figures/` for the six PDFs, `build/reported_statistics.json` for primary mean/median/IQR, direct-repeat means, image-only/plotting-style summaries, event counts, condition summaries and sensitivity results, and `build/verification.json` for the completed checks. TeX table exports are standalone data tables with full CSV companions; they are not a typeset replica of the manuscript tables. A full LaTeX manuscript build is outside this analysis pipeline. Private manuscript sources are not required by the analysis or verification commands.

## Result-to-code mapping

| Manuscript item | Calculation and rendering |
|---|---|
| Main Table 1 | `write_manuscript_tables.py`, using the run manifest |
| Main Tables 2–3; S1, S3, S4, S6 | `score_survival_reconstructions.R` → `compare_formats_and_guyot.R` → `write_manuscript_tables.py`; S6 also uses `analyze_exploratory_comparisons.py` |
| S2 | Direct repeated-IPD scoring in `score_survival_reconstructions.R` → table writer |
| S5, S9 | `audit_reconstruction_methods.py` → table writer |
| S7, S8, S10, S11 | R scoring + `calculate_uncertainty_and_events.py` → `analyze_exploratory_comparisons.py` → table writer |
| Main Figure 1 | `plot_study_design.py`, using actual saved input PDFs |
| Main Figures 2, 4; S1, S2 Figures | Corresponding calculated comparisons → `plot_accuracy_comparisons.py` |
| Main Figure 3 | `plot_reconstruction_example.py`, using original and reconstructed IPD plus the actual PNG |

## Maintained scripts

| Script | Coherent responsibility |
|---|---|
| `score_survival_reconstructions.R` | Right-continuous KM/TAAE, RMST and Cox scoring of frozen IPD; direct repeat comparisons and reference estimands. |
| `calculate_uncertainty_and_events.py` | Original Efron-Cox/Greenwood uncertainty and event-count calculations. |
| `compare_formats_and_guyot.R` | Original R paired-bootstrap comparisons: current64, historical62/61, image-only PDF and plotting style. |
| `analyze_exploratory_comparisons.py` | Original NumPy bootstrap for precision, exact-coordinate decomposition and figure conditions; missing-output and revision sensitivity. |
| `audit_reconstruction_methods.py` | Verify saved extraction classifications/source anchors and replay receipts; calculate plotting-template diagnostics. |
| `write_manuscript_tables.py` | Export all 14 tables and supporting prose statistics from calculated results. |
| `plot_study_design.py` | Figure 1 using the actual vector figures and assigned dataset subsets. |
| `plot_reconstruction_example.py` | Figure 3 using the original selection rule, actual PNG and reconstructed error steps. |
| `plot_accuracy_comparisons.py` | Figures 2, 4, S1 and S2 from calculated comparisons. S2 tick spacing is adjusted for legibility. |
| `verify_reproduction.py` | Check distributed payload hashes, documented public derivatives, manifest populations, unrounded authorities and published numerical cells. |
| `replay_guyot_reconstructions.R` (optional) | Execute the unchanged archived Guyot algorithm with each experiment's original information and compare rescored outputs. |
| `digitize_guyot_pngs.py` (optional) | Execute the unchanged PNG digitizer using actual PNGs, public sample sizes and the original pinned OCR runtime. |

Every script includes its inputs, outputs and standalone command. Archived programs are not imported into the default pipeline. Frozen full-precision authorities and table exports in `frozen_results/validation/` are read only by the verifier; they are never used to calculate replacement observations or manuscript estimates. Extraction-method labels are historical reviewed observations, checked against their saved code evidence rather than reclassified by a new model.

## Optional method replay

```sh
make guyot-replay
# Expected nonzero status for the exact-coordinate score difference described above.

python scripts/digitize_guyot_pngs.py --runtime /path/to/frozen/ocr/runtime
# For one dataset only, append: --dataset PP32_S1_D01
```

PNG digitization requires the original Node/Tesseract runtime layout: `node/bin/node`, `node/node_modules/` and `language/eng.traineddata`. The archived `environment/comparator-ocr-lock.json` records file hashes and Node version. The wrapper requires these exact files, performs no downloads and starts a fresh OCR worker for each image. The runtime is not bundled in the public repository. The original input payload for one dataset was successfully replayed; a fresh run of all 64 digitizations has not been validated. All 64 actual historical digitization payloads and Guyot IPDs are included, so neither OCR nor this optional replay is required for manuscript reproduction.

The model programs and historical replay receipts are available for extraction-method audit. A modern rerun of every archived model program is not part of `make reproduce`; no new model output is substituted for the study observations.

## Scope and publication notes

Unpublished manuscripts and their bibliography must remain local and must not be
committed, pushed, or included in repository archives without explicit author
approval. This repository is private; a future change to public visibility also
requires explicit author approval. Numerical validation references contain table
values only, not manuscript prose.

The local checkout uses the tracked pre-push guard. After a fresh clone, enable it
with `git config core.hooksPath .githooks`. It rejects the protected manuscript paths and filenames anywhere
in the history being pushed. Keep private source copies outside release packages;
`.gitignore` prevents routine staging but does not sanitize manually created ZIPs.

The repository includes the data underlying the current manuscript and necessary audit evidence. It excludes unrelated historical experiments, superseded result collections, full execution logs and generated `build/` files from version control. It does not execute model programs, reconstruct missing outputs, retune algorithms, alter manuscript statistics or publish to GitHub.

Original manuscript wording ambiguities previously identified in the audit have not been edited: the current64 intervals use four scenario strata, whereas an S10 caption refers to eight batch/scenario means; Figure 3 uses the recorded lower-median first-batch PNG example, although another sentence refers to largest errors. The actual analysis choices and populations are explicit here and in the scripts.

No repository-wide license or author/citation metadata has been invented. Authors should supply those before public release. The archived Guyot source retains its attribution to the 2012 supplementary algorithm and its stated CC BY 2.0 provenance.

During initial implementation, old-repository scientific-file checks and Git status showed no changes. CodeGraph automatically changed `.codegraph/daemon.log` and removed `.codegraph/daemon.pid` and `.codegraph/daemon.sock`. Those runtime side effects were reported; no further CodeGraph calls are needed to use this independent repository. The implementation integrity record is local in `build/old_repository_integrity.json`.
