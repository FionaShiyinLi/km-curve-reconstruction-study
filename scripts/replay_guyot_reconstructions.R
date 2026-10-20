#!/usr/bin/env Rscript
# Purpose: optional deterministic replay of both original Guyot experiments.
# Inputs: data/guyot_workflow.zip, data/guyot_digitization/*.json (PNG workflow),
#         data/reference_ipd/*.csv (exact-coordinate experiment ONLY).
# Outputs: build/guyot_replay/{Guyot,GuyotExact}_<dataset>.csv, validation.csv,
#          scoring_comparison.csv and summary.json.
# Run: Rscript --vanilla scripts/replay_guyot_reconstructions.R (repository root).
# Does not feed original records, exact coordinates, censor marks or total events
# to the PNG comparator. The exact experiment receives only step vertices, risk
# counts at 0,6,...,36 and sample size. The original Guyot numerical source is
# unpacked unchanged. New CSVs are compared with, and never replace, frozen IPDs.
# The known exact-coordinate replay difference is documented in the README.
# A time-distance check alone is insufficient for strict KM peak errors. Every
# serialized IPD is rescored; a changed scientific score produces a nonzero exit.

root <- normalizePath(".")
stopifnot(file.exists(file.path(root,"data/run_manifest.csv")))
scorer <- new.env()
sys.source(file.path(root,"scripts/score_survival_reconstructions.R"),envir=scorer)
source_root <- file.path(root,"build/guyot_source")
dir.create(source_root,recursive=TRUE,showWarnings=FALSE)
unzip(file.path(root,"data/guyot_workflow.zip"),
      files=c("R/reference.R","R/pdf_png300_guyot.R","R/guyot_2012_core.R"),exdir=source_root)
setwd(source_root)
source("R/reference.R")
source("R/pdf_png300_guyot.R")
runs <- read.csv(file.path(root,"data/run_manifest.csv"),stringsAsFactors=FALSE)
runs <- runs[runs$method %in% c("Guyot","GuyotExact"), ]
target <- file.path(root,"build/guyot_replay")
dir.create(target,recursive=TRUE,showWarnings=FALSE)
validation <- list()
comparisons <- list()
for(i in seq_len(nrow(runs))) {
  run <- runs[i, ]
  if(run$method=="Guyot") {
    extraction <- jsonlite::read_json(file.path(root,"data/guyot_digitization",paste0(run$dataset_id,".json")),simplifyVector=FALSE)
    stopifnot(isTRUE(extraction$usable))
  } else {
    reference <- read.csv(file.path(root,run$reference_ipd),check.names=FALSE)
    reference <- reference[c("group","observed_time","event")]
  }
  reconstructed <- lapply(c("A","B"),function(g) {
    if(run$method=="Guyot") {
      arm <- extraction$arms[[g]]
      coordinates <- as.data.frame(do.call(rbind,lapply(arm$coordinates,unlist)))
      names(coordinates) <- c("time","survival")
      risk <- data.frame(time=unlist(arm$risk$time),n_risk=unlist(arm$risk$n_risk))
      n <- arm$n
    } else {
      observed <- reference[reference$group==g, ]
      coordinates <- step_coordinates(km_from_observed(observed$observed_time,observed$event))
      risk <- data.frame(time=seq(0,36,6),n_risk=risk_counts(observed$observed_time,observed$event))
      n <- nrow(observed)
    }
    fit <- guyot_2012(coordinates,risk,n)
    data.frame(group=g,observed_time=fit$records$time,event=fit$records$event)
  })
  output <- file.path(target,paste0(run$run_id,".csv"))
  write.csv(do.call(rbind,reconstructed),output,row.names=FALSE)
  # Compare serialized files, not the transient in-memory object (S8 precision).
  new <- read.csv(output);old <- read.csv(file.path(root,run$reconstruction_output))
  stopifnot(identical(new$group,old$group),identical(new$event,old$event),nrow(new)==nrow(old))
  error <- max(abs(new$observed_time-old$observed_time))
  if(error>1e-12) stop(paste("Guyot replay discrepancy:",run$run_id,error))
  validation[[i]] <- data.frame(run_id=run$run_id,records=nrow(new),max_time_difference=error)
  # Reference records enter only this AFTER-reconstruction scoring step for
  # PNG-derived Guyot. They never enter that method's reconstruction inputs.
  reference <- read.csv(file.path(root,run$reference_ipd),check.names=FALSE)
  old_scores <- scorer$pp_session(reference,old,list(A=TRUE,B=TRUE))$values
  new_scores <- scorer$pp_session(reference,new,list(A=TRUE,B=TRUE))$values
  difference <- new_scores-old_scores
  comparisons[[i]] <- data.frame(run_id=run$run_id,method=run$method,
    metric=names(old_scores),frozen=unname(old_scores),replay=unname(new_scores),
    difference=unname(difference),within_tolerance=is.finite(difference)&abs(difference)<=1e-10)
}
write.csv(do.call(rbind,validation),file.path(target,"validation.csv"),row.names=FALSE)
comparison <- do.call(rbind,comparisons)
write.csv(comparison,file.path(target,"scoring_comparison.csv"),row.names=FALSE)
discrepancies <- comparison[!comparison$within_tolerance, ]
jsonlite::write_json(list(reconstructions=nrow(runs),scored_values=nrow(comparison),
  changed_values=nrow(discrepancies),score_tolerance=1e-10,
  status=if(nrow(discrepancies)) "scientific_score_difference" else "equivalent_within_tolerance",
  R_version=R.version.string,survival_version=as.character(packageVersion("survival"))),
  file.path(target,"summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=NA)
cat("Replayed",nrow(runs),"Guyot reconstructions; frozen outputs unchanged.\n")
if(nrow(discrepancies)) {
  print(discrepancies,row.names=FALSE)
  stop("Replay changes scientific scores; see build/guyot_replay/scoring_comparison.csv and README. Frozen-IPD manuscript reproduction is separate.",call.=FALSE)
}
cat("PASS: all rescored values agree within the stated tolerance.\n")
