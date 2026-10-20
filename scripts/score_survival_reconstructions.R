#!/usr/bin/env Rscript
# Purpose: score frozen reconstructed IPD against observed reference IPD.
# Inputs: data/run_manifest.csv, data/reference_ipd/, frozen_results/*_ipd/.
# Outputs: build/scores.csv, repeat_direct.csv, reference_estimates.csv.
# Run from the repository root: Rscript --vanilla scripts/score_survival_reconstructions.R
# Validation only: append --self-test. No model execution or reconstruction.
# Numerical kernels retained from R/reference.R and R/pdf_png300_scoring.R.
# CSV parsing is deliberately R read.csv: exact event-time ties must not be
# snapped or repaired. Curve errors are probabilities, NOT percentage points.

validate_observed <- function(time, event) {
  if (!is.numeric(time) || !length(time) || length(time) != length(event) ||
      anyNA(time) || any(!is.finite(time)) || any(time < 0) ||
      !(is.numeric(event) || is.logical(event)) || anyNA(event) ||
      any(!event %in% c(0, 1))) {
    stop("Observed data require nonempty finite nonnegative times and binary event indicators.")
  }
}

km_from_observed <- function(time, event) {
  validate_observed(time, event)
  ord <- order(time)
  y <- time[ord]
  delta <- as.integer(event[ord])
  blocks <- rle(y)
  count <- blocks$lengths
  end <- cumsum(count)
  start <- c(1L, head(end, -1L) + 1L)
  cumulative_events <- c(0L, cumsum(delta))
  deaths <- cumulative_events[end + 1L] - cumulative_events[start]
  risk <- length(y) - start + 1L
  survival <- cumprod(1 - deaths / risk)
  data.frame(time = blocks$values, n_risk = risk, n_event = deaths,
             n_censor = count - deaths,
             survival_before = c(1, head(survival, -1L)), survival = survival)
}

rmst_from_observed <- function(time, event, tau = 24) {
  if (length(tau) != 1L || !is.finite(tau) || tau <= 0) stop("tau must be positive.")
  km <- km_from_observed(time, event)
  endpoint <- max(time)
  zero_tail <- any(km$time <= tau & km$survival == 0)
  supported <- endpoint >= tau || zero_tail
  if (!supported) {
    return(list(rmst = NA_real_, available = FALSE,
                reason = "unsupported_nonzero_tail", endpoint = endpoint, km = km))
  }
  # Each interval uses the post-drop value at its LEFT endpoint. Vertical
  # drops (including an event at tau) contribute zero area.
  inside <- km$time < tau
  boundaries <- c(0, km$time[inside], tau)
  heights <- c(1, km$survival[inside])
  list(rmst = sum(diff(boundaries) * heights), available = TRUE,
       reason = if (endpoint < tau) "supported_zero_tail" else "supported_followup",
       endpoint = endpoint, km = km)
}
# Proposal A extensions. The legacy score_dataset entry point remains unchanged.
# Zero-survival Greenwood terms use their continuous boundary value of zero.
arm_inference <- function(time, event, tau = 24) {
  a <- rmst_from_observed(time, event, tau)
  k <- a$km
  at <- which(k$time <= tau)
  s <- if (length(at)) k$survival[tail(at, 1)] else 1
  use <- which(k$n_event > 0 & k$time <= tau & k$n_risk > k$n_event)
  gw <- k$n_event[use] / (k$n_risk[use] * (k$n_risk[use] - k$n_event[use]))
  # Backward integration from each observation through tau, including censor rows.
  inside <- which(k$time <= tau)
  areas <- if (length(inside)) rev(cumsum(rev(
    diff(c(k$time[inside], tau)) * k$survival[inside]))) else numeric()
  variance <- if (a$available) sum(areas[match(use, inside)]^2 * gw) else NA_real_
  crossed <- which(k$survival <= .5)
  c(rmst = a$rmst, se_rmst = sqrt(variance),
    survival = if (a$available) s else NA_real_,
    se_survival = if (!a$available) NA_real_ else if (s == 0) 0 else s * sqrt(sum(gw)),
    median = if (length(crossed)) k$time[crossed[1]] else NA_real_,
    events = sum(event))
}

cox_reference <- function(records) {
  z <- data.frame(time = records$observed_time, event = records$event,
                  treatment = as.integer(records$group == "B"))
  warnings <- character()
  fail <- function(reason) list(log_hr = NA_real_, se_log_hr = NA_real_,
                                status = reason, warnings = warnings)
  if (!sum(z$event)) return(fail("no_events"))
  if (length(unique(z$treatment[z$event == 1])) < 2L)
    return(fail("events_in_one_group"))
  fit <- tryCatch(withCallingHandlers(
    survival::coxph(survival::Surv(time, event) ~ treatment, data = z, ties = "efron",
      control = survival::coxph.control(timefix = FALSE, eps = 1e-11, iter.max = 100)),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) e)
  if (inherits(fit, "error")) return(fail(paste0("cox_error: ", conditionMessage(fit))))
  if (any(grepl("infinite|converg", warnings, ignore.case = TRUE)))
    return(fail("nonfinite_or_nonconverged_cox"))
  beta <- unname(fit$coefficients[1]); se <- unname(sqrt(fit$var[1, 1]))
  if (!is.finite(beta) || !is.finite(se) || se <= 0) return(fail("unidentified_cox"))
  list(log_hr = beta, se_log_hr = se, status = "available", warnings = warnings)
}

pp_curve <- function(reference, reconstructed, tau=24) {
  r <- rmst_from_observed(reference$observed_time, reference$event, tau)
  empty <- c(taae=NA_real_, maximum=NA_real_, duration_above_002=NA_real_,
             signed_rmst_error=NA_real_, absolute_rmst_error=NA_real_)
  if (is.null(reconstructed)) return(list(values=empty, reason="missing_or_invalid_arm"))
  x <- rmst_from_observed(reconstructed$observed_time, reconstructed$event, tau)
  if (!r$available || !x$available) return(list(values=empty, reason="unsupported_nonzero_tail"))
  knots <- sort(unique(c(0,tau,r$km$time[r$km$time<=tau],x$km$time[x$km$time<=tau])))
  height <- function(k) c(1,k$survival)[findInterval(knots,k$time)+1L]
  delta <- height(x$km)-height(r$km); left <- head(delta,-1L); widths <- diff(knots)
  signed <- sum(widths*left)
  list(values=c(taae=sum(widths*abs(left))/tau, maximum=max(abs(delta)),
    duration_above_002=sum(widths[abs(left)>.02+32*.Machine$double.eps]),
    signed_rmst_error=signed, absolute_rmst_error=abs(signed)),reason="available")
}

pp_session <- function(reference, reconstructed, arms) {
  required <- c("group","observed_time","event")
  reference <- reference[required] # Never use latent T/E/D or generating hazards.
  # One malformed arm can make read.csv infer a text column for BOTH arms.
  # Convert values explicitly and validate each arm independently before fitting.
  if(!is.null(reconstructed)) {
    reconstructed <- reconstructed[required]
    for(field in c("observed_time","event"))
      reconstructed[[field]] <- suppressWarnings(as.numeric(reconstructed[[field]]))
    for(g in c("A","B")) if(isTRUE(arms[[g]])) {
      d <- reconstructed[reconstructed$group==g,]
      arms[[g]] <- isTRUE(tryCatch({validate_observed(d$observed_time,d$event);TRUE},error=function(e) FALSE))
    }
  }
  curves <- setNames(lapply(c("A","B"), function(g) pp_curve(reference[reference$group==g,],
    if (isTRUE(arms[[g]])) reconstructed[reconstructed$group==g,required] else NULL)),c("A","B"))
  values <- unlist(lapply(c("A","B"),function(g) setNames(curves[[g]]$values,paste0(names(curves[[g]]$values),"_",g))))
  values <- c(values,taae_mean=mean(values[c("taae_A","taae_B")]),
    absolute_rmst_difference_error=abs(values[["signed_rmst_error_B"]]-values[["signed_rmst_error_A"]]))
  cx <- cox_reference(reference)
  rec <- if (all(unlist(arms))) cox_reference(reconstructed[required]) else list(log_hr=NA_real_,status="missing_or_invalid_arm",warnings=character())
  values <- c(values,absolute_log_hr_error=abs(rec$log_hr-cx$log_hr))
  list(values=values,curve_status=vapply(curves,`[[`,character(1),"reason"),
       cox_status=c(reference=cx$status,reconstructed=rec$status),cox_warnings=c(cx$warnings,rec$warnings))
}

score_frozen_runs <- function() {
  if (!file.exists("data/run_manifest.csv")) stop("Run this script from the repository root.")
  runs <- read.csv("data/run_manifest.csv", stringsAsFactors=FALSE, check.names=FALSE)
  for (field in c("complete", "is_first_run", "primary64", "complete_case62", "original61", "repeatability"))
    runs[[field]] <- runs[[field]] == "True"
  runs <- runs[runs$complete & runs$experiment != "effort_screen", ]
  stopifnot(nrow(runs)==528L, !anyDuplicated(runs$run_id))
  references <- lapply(setNames(unique(runs$reference_ipd), unique(runs$dataset_id)), read.csv)
  records <- setNames(lapply(runs$reconstruction_output, read.csv), runs$run_id)
  scores <- lapply(seq_len(nrow(runs)), function(i) {
    a <- runs[i, ]
    score <- pp_session(references[[a$dataset_id]], records[[a$run_id]], list(A=TRUE, B=TRUE))
    if (any(!is.finite(score$values))) stop(paste("Unavailable score:", a$run_id))
    cbind(a, as.data.frame(as.list(score$values)), status="complete_scorable")
  })
  scores <- do.call(rbind, scores)
  repeat_runs <- runs[runs$repeatability, ]
  repeated <- lapply(seq_len(nrow(repeat_runs)), function(i) {
    a <- repeat_runs[i, ]
    first <- runs[runs$experiment=="primary" & runs$is_first_run &
                    runs$dataset_id==a$dataset_id & runs$condition==a$condition, ]
    stopifnot(nrow(first)==1L)
    values <- pp_session(records[[first$run_id]], records[[a$run_id]], list(A=TRUE,B=TRUE))$values
    cbind(data.frame(dataset=a$dataset_id, scenario=a$scenario, representation=a$condition,
                     batch=a$batch, first_session=first$run_id, repeat_session=a$run_id),
          as.data.frame(as.list(values)))
  })
  stopifnot(length(repeated)==48L, length(unique(repeat_runs$dataset_id))==24L)
  estimates <- lapply(names(references), function(id) {
    d <- references[[id]]; cx <- cox_reference(d)
    values <- c(log_hr=cx$log_hr, se_log_hr=cx$se_log_hr)
    for (g in c("A","B")) {
      arm <- d[d$group==g, ]; v <- arm_inference(arm$observed_time, arm$event)
      values <- c(values, setNames(v,paste0(names(v),"_",g)))
    }
    cbind(data.frame(dataset_id=id), as.data.frame(as.list(values)))
  })
  dir.create("build",showWarnings=FALSE)
  write.csv(scores,"build/scores.csv",row.names=FALSE,na="")
  write.csv(do.call(rbind,repeated),"build/repeat_direct.csv",row.names=FALSE,na="")
  write.csv(do.call(rbind,estimates),"build/reference_estimates.csv",row.names=FALSE,na="")
  cat("Scored",nrow(scores),"saved IPDs and",length(repeated),"direct first-repeat comparisons.\n")
}

self_test <- function() {
  near <- function(actual, expected) {
    stopifnot(length(actual) == length(expected),
              all(is.finite(actual)), all(abs(actual - expected) < 1e-13))
  }
  one_jump <- function(deaths, when) {
    data.frame(observed_time = c(rep(when, deaths), rep(30, 100 - deaths)),
               event = c(rep(1, deaths), rep(0, 100 - deaths)))
  }
  
  for (deaths in c(1, 10, 25, 50)) {
    height <- deaths / 100
    reference <- one_jump(deaths, 12)
    for (shift in c(-0.2, -1e-6, 1e-6, 0.2)) {
      # Both jumps lie strictly within the evaluation window, with supported tails.
      reconstructed <- one_jump(deaths, 12 + shift)
      result <- pp_curve(reference, reconstructed, tau = 24)$values
      near(result['taae'], height * abs(shift) / 24)
      near(result['maximum'], height)
      near(result['signed_rmst_error'], height * shift)
      near(result['duration_above_002'], if (height > 0.02) abs(shift) else 0)
    }
    near(pp_curve(reference, reference)$values, rep(0, 5))
  }
  
  # Endpoint jumps have zero area but can change the maximum on the closed window.
  endpoint <- pp_curve(one_jump(25, 24), one_jump(25, 25))$values
  near(endpoint[c('taae', 'maximum', 'duration_above_002')], c(0, 0.25, 0))
  # A shift crossing the window is clipped by the window, not by the input times.
  crossing <- pp_curve(one_jump(25, 23.9), one_jump(25, 24.1))$values
  near(crossing['taae'], 0.25 * 0.1 / 24)
  
  # The primary session outcome averages two arm errors; one unchanged arm halves it.
  a <- one_jump(25, 12); a$group <- 'A'
  b <- one_jump(10, 8); b$group <- 'B'
  reference <- rbind(a, b)
  reconstructed <- reference
  reconstructed$observed_time[reconstructed$group == 'A' & reconstructed$event == 1] <- 12.2
  result <- pp_session(reference, reconstructed, list(A = TRUE, B = TRUE))$values
  near(result['taae_A'], 0.25 * 0.2 / 24)
  near(result['taae_B'], 0)
  near(result['taae_mean'], 0.25 * 0.2 / 48)
  cat('PASS: 16 interior shifts, four identities, endpoint/window checks, and two-arm averaging\n')
  cat('Example: delta S = 0.01, delta t = 1e-6 month: TAAE = 4.166667e-10; peak = 0.01\n')
  tied <- km_from_observed(c(0,6,6,24),c(1,1,0,0))
  near(tied$survival,c(.75,.5,.5))
  stopifnot(!rmst_from_observed(c(1,2),c(0,0))$available,
            rmst_from_observed(c(1,2),c(1,1))$available)
  ref <- read.csv("data/reference_ipd/PP32_S1_D01.csv")
  near(pp_session(ref,ref,list(A=TRUE,B=TRUE))$values,rep(0,13))
  set.seed(1701)
  for (i in 1:20) {
    z <- data.frame(time=sample(c(0,1.25,3,6,12,24,30),100,TRUE),event=rbinom(100,1,.6))
    near(km_from_observed(z$time,z$event)$survival,
         survival::survfit(survival::Surv(time,event)~1,data=z)$surv)
  }
  cat("PASS: ties, right continuity, support, two-arm averaging, analytic shift, survival cross-check.\n")
}
if (sys.nframe()==0L) {
  if ("--self-test" %in% commandArgs(trailingOnly=TRUE)) self_test() else score_frozen_runs()
}
