#!/usr/bin/env Rscript
# Post hoc comparison with Guyot (2012), Additional file 1, on prepared inputs.
# Run: Rscript --vanilla scripts/validate_guyot_implementation.R (repo root).
# Uses only existing data/guyot_workflow.zip, digitizations, reference IPD and
# run_manifest.csv. Writes 256 arm comparisons under ignored build/guyot_validation/.
# This tests the numerical core, not digitization, input preparation or CSV replay.

root <- normalizePath(".")
archive <- file.path(root, "data/guyot_workflow.zip")
out <- file.path(root, "build/guyot_validation")
source_root <- file.path(out, "source")
dir.create(source_root, recursive=TRUE, showWarnings=FALSE)
unzip(archive, files=c("R/reference.R", "R/pdf_png300_guyot.R", "R/guyot_2012_core.R"),
      exdir=source_root)

# Read the published script already archived with its attribution. Supply its
# inputs below and omit its intermediate file output; retain every calculation.
published <- readLines(unz(archive,
  "docs/proposals/pdf_png300_guyot_baseline_v1/guyot_2012_supplement.txt"), warn=FALSE)
start <- which(trimws(published) == "#Initialise vectors")[1]
end <- which(trimws(published) == "#Output IPD") - 1L
body <- published[start:end]
body <- body[!grepl("^write.table\\(", trimws(body))]
reference_code <- parse(text=body, keep.source=FALSE)
custom_code <- parse(file.path(source_root, "R/guyot_2012_core.R"), keep.source=FALSE)
stopifnot(identical(reference_code, custom_code))
setwd(source_root) # The frozen adapter sources its numerical core by relative path.
source("R/reference.R")
source("R/pdf_png300_guyot.R")

reference_fit <- function(coordinates, risk) {
  e <- new.env(parent=globalenv())
  e$t.S <- coordinates$time; e$S <- coordinates$survival
  e$t.risk <- risk$time; e$n.risk <- risk$displayed
  e$lower <- match(e$t.risk, e$t.S)
  e$upper <- c(e$lower[-1] - 1L, length(e$t.S))
  e$n.int <- length(e$n.risk); e$n.t <- length(e$t.S)
  e$tot.events <- "NA"; e$arm.id <- 1L
  eval(reference_code, e)
  list(records=data.frame(time=e$t.IPD, event=e$event.IPD),
       counts=data.frame(time=e$t.S, n_risk=e$n.hat[seq_len(e$n.t)], events=e$d, censors=e$cen),
       adjusted_risk=e$n.risk)
}

runs <- read.csv(file.path(root, "data/run_manifest.csv"), stringsAsFactors=FALSE)
runs <- runs[runs$method %in% c("Guyot", "GuyotExact"), ]
stopifnot(nrow(runs) == 128L, !anyDuplicated(runs$run_id),
          length(unique(runs$dataset_id)) == 64L)
comparisons <- list()
for (i in seq_len(nrow(runs))) {
  run <- runs[i, ]
  if (run$method == "Guyot") {
    extraction <- jsonlite::read_json(file.path(root, "data/guyot_digitization",
                                              paste0(run$dataset_id, ".json")))
    stopifnot(isTRUE(extraction$usable))
  } else {
    reference <- read.csv(file.path(root, run$reference_ipd))
    reference <- reference[c("group", "observed_time", "event")]
  }
  for (g in c("A", "B")) {
    if (run$method == "Guyot") {
      arm <- extraction$arms[[g]]
      coordinates <- as.data.frame(do.call(rbind, lapply(arm$coordinates, unlist)))
      names(coordinates) <- c("time", "survival")
      risk <- data.frame(time=unlist(arm$risk$time), n_risk=unlist(arm$risk$n_risk))
      n <- arm$n
    } else {
      observed <- reference[reference$group == g, ]
      coordinates <- step_coordinates(km_from_observed(observed$observed_time, observed$event))
      risk <- data.frame(time=seq(0, 36, 6), n_risk=risk_counts(observed$observed_time, observed$event))
      n <- nrow(observed)
    }
    fit <- guyot_2012(coordinates, risk, n)
    check <- reference_fit(fit$coordinates, fit$risk)
    records_match <- isTRUE(all.equal(fit$records, check$records, tolerance=0))
    counts_match <- isTRUE(all.equal(fit$counts, check$counts, tolerance=0))
    stopifnot(records_match, counts_match, identical(fit$risk$algorithm, check$adjusted_risk))
    comparisons[[length(comparisons) + 1L]] <- data.frame(
      dataset=run$dataset_id, method=run$method, arm=g,
      records_identical=records_match, allocations_identical=counts_match,
      last_used_anchor=max(fit$risk$time),
      endpoint_risk_rows_not_used=fit$endpoint_risk_rows_not_used)
  }
}
results <- do.call(rbind, comparisons)
stopifnot(nrow(results) == 256L)
anchors <- results[results$method == "Guyot" & results$last_used_anchor == 36, ]
write.csv(results, file.path(out, "arm_comparisons.csv"), row.names=FALSE)
summary <- list(arms_compared=nrow(results), numerical_body_identical=TRUE,
  reconstructed_records_and_allocations_identical=all(results$records_identical & results$allocations_identical),
  png_arms_retaining_month36=nrow(anchors), png_datasets_retaining_month36=length(unique(anchors$dataset)),
  exact_arms_retaining_month36=sum(results$method == "GuyotExact" & results$last_used_anchor == 36))
jsonlite::write_json(summary, file.path(out, "summary.json"), pretty=TRUE, auto_unbox=TRUE)
cat("PASS:", nrow(results), "arm inputs matched published calculations;",
    nrow(anchors), "PNG arms retained the month-36 anchor. Results: build/guyot_validation/\n")
