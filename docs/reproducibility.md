# Reproducibility details

[Back to the README](../README.md)

## Data and analysis groups

The repository contains 68 simulated reference datasets, 360 assigned input figures, and 543 completed reconstructed IPD files. The main comparison uses 64 datasets; four additional datasets belong to the preliminary model-setting screen.

[`data/run_manifest.csv`](../data/run_manifest.csv) records 550 study runs, including seven without completed output. Each row represents a run. The dataset ID links runs using the same underlying data.

| Manifest field | Meaning |
|---|---|
| `primary64` | The main analysis: one completed PDF and PNG300 reconstruction per dataset, including documented recoveries. |
| `complete_case62` | The earlier 62 complete pairs, excluding the two later PNG recoveries, R014R1 and R064R1. |
| `original61` | Also excludes the dataset whose PDF required the startup recovery R003R1. |
| `is_first_run`, `attempt`, `repeatability`, `pairing_key` | Identify first and repeated reconstructions and their pairing. |
| `complete`, `assignment_status`, `original_run_id`, `replaces_run_id`, `replaced_by_run_id` | Record completion and recovery history. |
| `reference_ipd`, `input_figure`, `reconstruction_output` | Locate the reference records, assigned figure, and saved reconstruction. |
| `program_archive`, `program_directory`, `prompt_member`, `public_manifest_member`, `replay_receipt_member` | Locate the archived program, prompt, input information, and verification record. |

An empty `reconstruction_output` means that no completed output exists. Missing errors are not assigned zero. The population flags identify complete pairs; missing-output analyses also retain the observed partner of an incomplete pair.

The 550 runs comprise 181 main/repeat runs, 160 modified-figure runs, 32 image-only PDF runs, 33 plotting-style runs, 128 Guyot reconstructions, and 16 preliminary-screen runs. There are 48 first–repeat comparisons across 24 datasets. Repeated runs do not add independent datasets to the main analysis.

The eight attempts superseded by the initial budget restart are documented in the archived restart record and are outside these 550 runs. Generation seeds unavailable for the preliminary screen remain blank in [`data/datasets.csv`](../data/datasets.csv).

## Analysis settings

- **Curve error:** time-averaged absolute error (TAAE) measures the area between reconstructed and reference KM curves through 24 months, divided by 24, then averaged across the two arms. Curves are evaluated as right-continuous steps. Events tied with censoring use the same risk set, with events applied first. Maximum error includes the 24-month endpoint. A nonzero survival tail is not extrapolated.
- **Treatment effects:** Cox models use full follow-up, Efron ties, and `timefix=FALSE`. The scripts retain the recorded optimization settings.
- **Units:** stored curve errors are probabilities; displayed percentage-point errors are multiplied by 100. Relative standard-error errors are fractions unless explicitly displayed as percentages. The duration-above-threshold calculation uses `0.02 + 32 * .Machine$double.eps`.
- **Confidence intervals:** comparisons use 10,000 paired bootstrap samples of whole datasets and type-1 percentile endpoints. The main 64-dataset analysis samples within four scenario groups of 16 datasets. Historical missing-output analyses retain eight batch-by-scenario groups with eight assigned positions each. Follow-up experiments use four groups of eight; the alternative plotting-template experiment uses four groups of four. Revision-limit sensitivity retains eight draws per available batch/scenario group, including groups with fewer retained datasets.
- **Other settings:** the scripts retain the original random-number generators, seeds, row order, and comparison-specific seed offsets. Main descriptive interquartile ranges use type-6 quantiles; figure-condition summaries use type-7 quantiles.

For interpreting figure selection, the reconstruction example uses the recorded lower-median first-batch PNG case. The analysis scripts and manifest specify the dataset selections and populations used to generate each output.

## Verification and file provenance

`make reproduce` checks inputs before analysis and checks both inputs and results after generating the outputs. `make verify` repeats those checks on existing outputs. Validation covers 528 scored reconstructions and 986 numerical values at their reported precision. The 15 completed preliminary-screen outputs are retained for completion reporting.

The reference values in [`frozen_results/validation/`](../frozen_results/validation/) are used to check results. The analysis calculates its estimates from the saved participant records. The extraction-method audit also checks 400 historical records of identical program reruns; it does not execute those programs again.

[`data/source_checksums.csv`](../data/source_checksums.csv) records file sizes and SHA-256 hashes. For archived records whose personal paths were removed, `original_sha256`, `original_bytes`, and `transformation` document the original file and the edit. These changes cover 207 administrative JSON records and one figure-generation helper's Rscript path. Historical hashes in archived records refer to their original targets.

The software environment was validated on macOS Apple Silicon with Python 3.12.14, R 4.5.2, and survival 3.8-3. An independent installation from the recorded local package archives reproduced the analysis. Fresh package downloads and other platforms have not been tested. To select existing executables explicitly:

```sh
make reproduce PYTHON=/path/to/python RSCRIPT=/path/to/Rscript
```

## Optional reruns

These commands rerun parts of the reconstruction process. They are separate from the standard analysis of saved records.

### Guyot reconstruction

```sh
make guyot-replay
```

The PNG comparison uses the saved digitized curve coordinates and risk counts. The exact-coordinate experiment uses exact curve coordinates and displayed risk counts. Neither method receives censoring marks or the total event count; reference participant records are used to derive inputs only for the exact-coordinate experiment.

In the validation run, all 64 PNG-derived Guyot outputs matched the saved files byte for byte. In the exact-coordinate experiment, 58 of 64 CSV files differed in their final decimal digits, with a maximum time difference of approximately `1.03e-13` months. These differences changed three maximum-error measurements:

| Dataset | Arm A maximum error in saved records | Arm A maximum error after rerun |
|---|---:|---:|
| PP32_S1_D02 | 0.0139185276089289 | 0.0140203547584752 |
| PP32_S2_D03 | 0.0127620891704286 | 0.0126858485937197 |
| PP32_S3_D06 | 0.00154394433040783 | 0.00154298217589044 |

Values in this table are probabilities. Across the 64 datasets, mean arm A maximum error changes from **0.7355173683240882 to 0.7355558439777132 percentage points**: **0.7355 versus 0.7356** when rounded to four decimals. The cause remains unresolved; the historical experiment used a different R environment, which is one possible explanation.

The command checks 13 scores per reconstruction. Three of the 832 exact-coordinate scores exceed its `1e-10` tolerance, so the command currently exits with a failure status. Reports are written to `build/guyot_replay/scoring_comparison.csv`, `validation.csv`, and `summary.json`. New files are saved separately from the original reconstructions. The standard analysis uses the original saved CSV records, including the CSV-rescored S8 results.

### PNG digitization

```sh
python scripts/digitize_guyot_pngs.py --runtime /path/to/frozen/ocr/runtime
```

This optional command needs the original Node/Tesseract OCR software and language files, which are not bundled here. Their versions and hashes are recorded in `environment/comparator-ocr-lock.json` inside [`data/guyot_workflow.zip`](../data/guyot_workflow.zip). The required layout is `node/bin/node`, `node/node_modules/`, and `language/eng.traineddata`.

Add `--dataset PP32_S1_D01` to check one dataset. That dataset's original digitization was reproduced during validation; a fresh run of all 64 digitizations has not been validated. The standard analysis uses the saved digitization and reconstruction files.

## Notes for contributors

Keep manuscript text, bibliography files, and unpublished document exports outside version control. Numerical results, table fragments, and underlying research data are included in this repository and become visible if it is made public. Publishing or changing visibility requires the author's approval.

Enable the manuscript-file check after cloning:

```sh
git config core.hooksPath .githooks
```

The hook rejects protected manuscript filenames in the history being pushed. Check manually created archives separately. A repository-wide license has not yet been specified; the archived Guyot source retains its own attribution and license information.
