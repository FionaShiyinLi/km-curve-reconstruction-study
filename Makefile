# Reproduce the published analyses from saved study data. No model/API calls.
# Run from the repository root after activating the recorded environment.
PYTHON ?= python
RSCRIPT ?= Rscript
export PYTHONDONTWRITEBYTECODE := 1
export MPLBACKEND := Agg
export MPLCONFIGDIR := $(CURDIR)/build/.matplotlib
.DEFAULT_GOAL := reproduce
.NOTPARALLEL:
.PHONY: reproduce verify check-inputs test guyot-replay

reproduce:
	$(PYTHON) scripts/verify_reproduction.py --inputs
	$(RSCRIPT) --vanilla scripts/score_survival_reconstructions.R --self-test
	$(RSCRIPT) --vanilla scripts/score_survival_reconstructions.R
	$(PYTHON) scripts/calculate_uncertainty_and_events.py
	$(RSCRIPT) --vanilla scripts/compare_formats_and_guyot.R
	$(PYTHON) scripts/analyze_exploratory_comparisons.py
	$(PYTHON) scripts/audit_reconstruction_methods.py
	$(PYTHON) scripts/write_manuscript_tables.py
	$(PYTHON) scripts/plot_study_design.py
	$(PYTHON) scripts/plot_reconstruction_example.py
	$(PYTHON) scripts/plot_accuracy_comparisons.py
	$(PYTHON) scripts/verify_reproduction.py --inputs --results

verify:
	$(PYTHON) scripts/verify_reproduction.py --inputs --results

check-inputs:
	$(PYTHON) scripts/verify_reproduction.py --inputs

test:
	$(RSCRIPT) --vanilla scripts/score_survival_reconstructions.R --self-test

# Optional, separate from manuscript reproduction. Writes a discrepancy report
# and exits nonzero if freshly reconstructed IPD changes a scientific score.
guyot-replay:
	$(RSCRIPT) --vanilla scripts/replay_guyot_reconstructions.R
