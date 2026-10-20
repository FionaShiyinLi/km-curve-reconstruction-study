# KM curve reconstruction study

This project studies how figure format affects the reconstruction of individual participant data (IPD) from Kaplan–Meier survival curves. The main comparison uses **64 simulated datasets**, each supplied as a native vector PDF and a 300-dpi PNG image. Additional experiments examine image resolution, compression, censoring marks, plotting style, and the conventional Guyot reconstruction method.

Use this repository to reproduce the study's tables and figures from the saved reconstructions. **No AI model access or API key is needed.**

## What is included?

| Folder | Contents |
|---|---|
| `data/` | Simulated participant records, the figures used in the experiments, and a list of study runs. |
| `frozen_results/` | Saved reconstructed records, model-generated programs, and reference results for checking the analysis. |
| `scripts/` | Code to calculate errors, compare methods, and create tables and figures. |
| `environment/` | The recorded software versions needed to reproduce the analysis. |

The repository contains the underlying research data and numerical results. Manuscript text and the bibliography are maintained separately.

## Run the analysis

The recorded environment has been tested on **macOS with Apple Silicon**, using Python 3.12.14 and R 4.5.2. Other platforms have not yet been validated.

With Conda installed, open a terminal in the repository folder and run:

```sh
conda create --name km-curve-study --file environment/conda-osx-arm64.lock
conda activate km-curve-study
make reproduce
```

This checks the input files, calculates the results, creates **14 tables and 6 figures**, and checks the results against the saved reference values.

## Find the results

All generated files are written to `build/`:

| Location | What to open |
|---|---|
| `build/tables/tables.html` | All tables in a browser. |
| `build/figures/` | Figures as PDF files. |
| `build/reported_statistics.json` | Numerical summaries used in the manuscript. |
| `build/verification.json` | Results of the reproducibility checks. |

To check existing outputs again, run `make verify`. To run the survival-calculation checks alone, use `make test`.

## Reproducibility notes

The main analysis uses the saved study reconstructions. An **optional rerun of the Guyot method with exact curve coordinates** produces numerical differences that change some maximum-error values. Its verification currently reports a failure. This limitation is separate from reproducing the tables and figures from the saved records.

See [Reproducibility details](docs/reproducibility.md) for dataset definitions, statistical settings, optional reruns, and the documented differences.
