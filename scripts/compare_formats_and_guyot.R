#!/usr/bin/env Rscript
# Purpose: paired format/Guyot comparisons, historical missing-output sensitivities,
# and image-only/alternative-template experiments, using the original R algorithms.
# Inputs: build/scores.csv (newly scored frozen IPDs).
# Outputs: build/{primary_pairs,contrasts_current64,historical62,historical61,
#          image_only_pairs,image_only_bootstrap,survminer_pairs,survminer_bootstrap}.csv.
# Run: Rscript --vanilla scripts/compare_formats_and_guyot.R (repository root).
# Primary64: four scenario strata of 16; historical62/61: eight batch x scenario
# strata, eight assigned positions including missing outputs. All B=10000,
# Mersenne-Twister/Inversion/Rejection, type-1 percentile endpoints. Missing
# records are NEVER replaced by zeros. The zero in absolute-error contrasts is
# the comparison to perfect reconstruction, not an imputation.

pp_bootstrap <- function(pairs, metrics, B=10000L, seed=2026092306L) {
  scenarios <- paste0("S",1:4)
  stopifnot(nrow(pairs)==32L,!anyDuplicated(pairs$dataset),all(table(factor(pairs$scenario,levels=scenarios))==8L))
  RNGkind("Mersenne-Twister","Inversion","Rejection");set.seed(seed)
  # All eight assigned datasets are resampled in each scenario, including missing
  # pairs. The same whole-dataset indices are shared by every metric and scope.
  picks <- lapply(scenarios,function(s) matrix(sample.int(8L,8L*B,replace=TRUE),nrow=8L))
  scopes <- list(all=scenarios,n100=c("S1","S2"),n1000=c("S3","S4"),
                 noncrossing=c("S1","S3"),crossing=c("S2","S4"))
  result <- list()
  for (metric in metrics) {
    means <- rep(NA_real_,4); complete <- integer(4); draws <- matrix(NA_real_,B,4)
    for (i in 1:4) {
      d <- pairs[pairs$scenario==scenarios[i],paste0(metric,"_PDF")]-pairs[pairs$scenario==scenarios[i],paste0(metric,"_PNG300")]
      complete[i] <- sum(is.finite(d))
      means[i] <- if(any(is.finite(d))) mean(d[is.finite(d)]) else NA_real_
      z <- matrix(d[picks[[i]]],nrow=8); counts <- colSums(is.finite(z)); z[!is.finite(z)] <- 0
      draws[,i] <- ifelse(counts>0,colSums(z)/counts,NA_real_)
    }
    for (scope in names(scopes)) {
      use <- match(scopes[[scope]],scenarios); estimate <- mean(means[use]); x <- rowMeans(draws[,use,drop=FALSE])
      undefined <- sum(!is.finite(x)); degenerate <- !undefined && length(unique(x))==1L
      status <- if(!is.finite(estimate)) "undefined_estimate" else if(undefined) "undefined_resamples" else if(degenerate) "degenerate_resamples" else "available"
      available <- status=="available"
      result[[length(result)+1L]] <- data.frame(metric=metric,scope=scope,estimate=estimate,
        complete_pairs=sum(complete[use]),
        pairs_S1=if(1 %in% use) complete[1] else NA_integer_,pairs_S2=if(2 %in% use) complete[2] else NA_integer_,
        pairs_S3=if(3 %in% use) complete[3] else NA_integer_,pairs_S4=if(4 %in% use) complete[4] else NA_integer_,
        lower=if(available) unname(quantile(x,.025,type=1)) else NA_real_,
        upper=if(available) unname(quantile(x,.975,type=1)) else NA_real_,
        standard_error=if(available) sd(x) else NA_real_,undefined_resamples=undefined,status=status,B=B,seed=seed)
    }
  }
  do.call(rbind,result)
}

pp_pairs <- function(scores) {
  primary <- scores[scores$primary_attempt,]
  if(nrow(primary)!=64L || anyDuplicated(primary[c("dataset","representation")])) stop("Invalid primary assignment table")
  metrics <- names(primary)[vapply(primary,is.numeric,logical(1))]
  metrics <- intersect(metrics,c("taae_mean","taae_A","taae_B","maximum_A","maximum_B",
    "duration_above_002_A","duration_above_002_B","absolute_rmst_error_A","absolute_rmst_error_B",
    "absolute_rmst_difference_error","absolute_log_hr_error"))
  cols <- c("dataset","scenario","status",metrics)
  p <- merge(primary[primary$representation=="PDF",cols],primary[primary$representation=="PNG300",cols],
    by=c("dataset","scenario"),suffixes=c("_PDF","_PNG300"),sort=TRUE)
  pdf <- is.finite(p$taae_mean_PDF); png <- is.finite(p$taae_mean_PNG300)
  p$pair_status <- ifelse(pdf & png,"both",ifelse(pdf,"PDF_only",ifelse(png,"PNG300_only","neither")))
  p$taae_difference <- p$taae_mean_PDF-p$taae_mean_PNG300
  p$lower_bound <- ifelse(pdf,p$taae_mean_PDF,0)-ifelse(png,p$taae_mean_PNG300,1)
  p$upper_bound <- ifelse(pdf,p$taae_mean_PDF,1)-ifelse(png,p$taae_mean_PNG300,0)
  list(pairs=p,metrics=metrics)
}
pooled_bootstrap <- function(original, replication, metrics, B=10000L, seed=2026092707L) {
  scenarios <- paste0('S', 1:4)
  for(p in list(original, replication)) stopifnot(nrow(p)==32L, !anyDuplicated(p$dataset),
    all(table(factor(p$scenario, levels=scenarios))==8L))
  stopifnot(!length(intersect(original$dataset, replication$dataset)))
  original$cohort <- 'original'; replication$cohort <- 'replication'
  p <- rbind(original, replication)
  strata <- split(seq_len(nrow(p)), interaction(p$cohort, p$scenario, drop=TRUE))
  RNGkind('Mersenne-Twister','Inversion','Rejection'); set.seed(seed)
  picks <- lapply(strata, function(ids) matrix(sample.int(8L,8L*B,replace=TRUE),nrow=8L))
  results <- list()
  for(metric in metrics) {
    delta <- p[[paste0(metric,'_PDF')]]-p[[paste0(metric,'_PNG300')]]
    means <- vapply(strata, function(ids) if(any(is.finite(delta[ids]))) mean(delta[ids][is.finite(delta[ids])]) else NA_real_, numeric(1))
    draws <- sapply(seq_along(strata), function(i) {
      d <- delta[strata[[i]]]; z <- matrix(d[picks[[i]]], nrow=8L)
      count <- colSums(is.finite(z)); z[!is.finite(z)] <- 0
      ifelse(count>0, colSums(z)/count, NA_real_)
    })
    x <- rowMeans(draws); estimate <- mean(means); undefined <- sum(!is.finite(x))
    status <- if(!is.finite(estimate)) 'undefined_estimate' else if(undefined) 'undefined_resamples' else if(length(unique(x))==1L) 'degenerate_resamples' else 'available'
    available <- status=='available'
    results[[length(results)+1L]] <- data.frame(metric=metric, scope='pooled64_secondary',
      estimate=estimate, complete_pairs=sum(is.finite(delta)),
      lower=if(available) unname(quantile(x,.025,type=1)) else NA_real_,
      upper=if(available) unname(quantile(x,.975,type=1)) else NA_real_,
      standard_error=if(available) sd(x) else NA_real_, undefined_resamples=undefined, status=status, B=B, seed=seed)
  }
  do.call(rbind, results)
}

# Post-result three-format aggregation. Source pdf_png300_scoring.R first.
pp_three_format_pairs <- function(scores) {
  formats <- c("PDF", "RASTER_PDF", "PNG300")
  metrics <- c("taae_mean", "taae_A", "taae_B", "maximum_A", "maximum_B",
    "duration_above_002_A", "duration_above_002_B", "absolute_rmst_error_A",
    "absolute_rmst_error_B", "absolute_rmst_difference_error", "absolute_log_hr_error")
  primary <- scores[scores$primary_attempt, ]
  if (nrow(primary) != 96L || anyDuplicated(primary[c("dataset", "representation")]) ||
      !all(formats %in% primary$representation) || !all(metrics %in% names(primary)))
    stop("Expected exactly 32 unique primary scores for each of three formats")
  if (!all(table(factor(primary$representation, levels=formats)) == 32L))
    stop("Unequal representation counts")
  datasets <- unique(primary[c("dataset", "scenario")])
  if (nrow(datasets) != 32L || anyDuplicated(datasets$dataset) ||
      !all(table(factor(datasets$scenario, levels=paste0("S", 1:4))) == 8L))
    stop("Expected eight independent datasets per scenario")
  pairs <- datasets
  for (format in formats) {
    part <- primary[primary$representation == format, c("dataset", "scenario", "status", metrics)]
    names(part)[-(1:2)] <- paste0(names(part)[-(1:2)], "_", format)
    pairs <- merge(pairs, part, by=c("dataset", "scenario"), sort=TRUE)
  }
  if (nrow(pairs) != 32L) stop("A three-format pair is missing")
  for (contrast in c("PDF_minus_RASTER_PDF", "RASTER_PDF_minus_PNG300")) {
    formats_here <- strsplit(contrast, "_minus_", fixed=TRUE)[[1]]
    first <- pairs[[paste0("taae_mean_", formats_here[1])]]
    second <- pairs[[paste0("taae_mean_", formats_here[2])]]
    pairs[[paste0("complete_", contrast)]] <- is.finite(first) & is.finite(second)
    pairs[[paste0("taae_difference_", contrast)]] <- first - second
    pairs[[paste0("lower_bound_", contrast)]] <- ifelse(is.finite(first), first, 0) -
      ifelse(is.finite(second), second, 1)
    pairs[[paste0("upper_bound_", contrast)]] <- ifelse(is.finite(first), first, 1) -
      ifelse(is.finite(second), second, 0)
  }
  list(pairs=pairs, metrics=metrics)
}

pp_three_format_bootstrap <- function(pairs, metrics, B=10000L, seed=2026092309L) {
  scenarios <- paste0("S", 1:4)
  contrasts <- c("PDF_minus_RASTER_PDF", "RASTER_PDF_minus_PNG300")
  scopes <- list(all=scenarios, n100=c("S1", "S2"), n1000=c("S3", "S4"),
    noncrossing=c("S1", "S3"), crossing=c("S2", "S4"),
    S1="S1", S2="S2", S3="S3", S4="S4")
  if (nrow(pairs) != 32L || anyDuplicated(pairs$dataset) ||
      !all(table(factor(pairs$scenario, levels=scenarios)) == 8L))
    stop("Expected eight assigned datasets in each scenario")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  picks <- lapply(scenarios, function(s) matrix(sample.int(8L, 8L*B, replace=TRUE), nrow=8L))
  result <- list()
  for (contrast in contrasts) {
    formats_here <- strsplit(contrast, "_minus_", fixed=TRUE)[[1]]
    for (metric in metrics) {
      means <- rep(NA_real_, 4L)
      complete <- integer(4L)
      draws <- matrix(NA_real_, B, 4L)
      bounds <- matrix(NA_real_, 4L, 2L)
      for (i in 1:4) {
        block <- pairs[pairs$scenario == scenarios[i], ]
        first <- block[[paste0(metric, "_", formats_here[1])]]
        second <- block[[paste0(metric, "_", formats_here[2])]]
        d <- first - second
        complete[i] <- sum(is.finite(d))
        if (complete[i]) means[i] <- mean(d[is.finite(d)])
        z <- matrix(d[picks[[i]]], nrow=8L)
        counts <- colSums(is.finite(z))
        z[!is.finite(z)] <- 0
        draws[,i] <- ifelse(counts > 0, colSums(z)/counts, NA_real_)
        if (metric == "taae_mean") {
          bounds[i,] <- c(mean(ifelse(is.finite(first), first, 0) -
            ifelse(is.finite(second), second, 1)),
            mean(ifelse(is.finite(first), first, 1) -
              ifelse(is.finite(second), second, 0)))
        }
      }
      for (scope in names(scopes)) {
        use <- match(scopes[[scope]], scenarios)
        estimate <- if (all(is.finite(means[use]))) mean(means[use]) else NA_real_
        x <- rowMeans(draws[,use,drop=FALSE])
        undefined <- sum(!is.finite(x))
        degenerate <- !undefined && length(unique(x)) == 1L
        status <- if (!is.finite(estimate)) "undefined_estimate" else if (undefined)
          "undefined_resamples" else if (degenerate) "degenerate_resamples" else "available"
        available <- status == "available"
        result[[length(result)+1L]] <- data.frame(contrast=contrast, metric=metric, scope=scope,
          estimate=estimate, complete_pairs=sum(complete[use]),
          assigned_pairs=8L*length(use),
          lower=if (available) unname(quantile(x, .025, type=1)) else NA_real_,
          upper=if (available) unname(quantile(x, .975, type=1)) else NA_real_,
          standard_error=if (available) sd(x) else NA_real_,
          undefined_resamples=undefined, status=status,
          full_assigned_lower=if (metric == "taae_mean") mean(bounds[use,1]) else NA_real_,
          full_assigned_upper=if (metric == "taae_mean") mean(bounds[use,2]) else NA_real_,
          B=B, seed=seed)
      }
    }
  }
  do.call(rbind, result)
}
# Post-result 16-dataset, two-format template robustness aggregation.
# Source reference.R and pdf_png300_scoring.R first.

pp_sv_pairs <- function(new_scores, old_scores, selected) {
  metrics <- c("taae_mean", "taae_A", "taae_B", "maximum_A", "maximum_B",
    "duration_above_002_A", "duration_above_002_B", "absolute_rmst_error_A",
    "absolute_rmst_error_B", "absolute_rmst_difference_error", "absolute_log_hr_error")
  formats <- c("PDF", "PNG300")
  if (length(selected) != 16L || anyDuplicated(selected) ||
      !all(metrics %in% names(new_scores)) || !all(metrics %in% names(old_scores)))
    stop("Expected 16 unique selected datasets and unchanged score metrics")
  new <- new_scores[new_scores$primary_attempt, ]
  old <- old_scores[old_scores$primary_attempt & old_scores$dataset %in% selected, ]
  if (nrow(new) != 32L || nrow(old) != 32L ||
      anyDuplicated(new[c("dataset", "representation")]) ||
      anyDuplicated(old[c("dataset", "representation")]) ||
      !setequal(new$dataset, selected) || !setequal(old$dataset, selected) ||
      !all(table(factor(new$representation, levels = formats)) == 16L) ||
      !all(table(factor(old$representation, levels = formats)) == 16L))
    stop("Expected one original and one new primary score per format and dataset")
  if (!all(table(factor(new$scenario[!duplicated(new$dataset)],
                        levels = paste0("S", 1:4))) == 4L))
    stop("Expected four independent selected datasets per scenario")
  pairs <- unique(new[c("dataset", "scenario")])
  for (cohort in c("new", "old")) {
    source <- if (cohort == "new") new else old
    for (format in formats) {
      part <- source[source$representation == format,
                     c("dataset", "scenario", "status", metrics)]
      names(part)[-(1:2)] <- paste0(names(part)[-(1:2)], "_", cohort, "_", format)
      pairs <- merge(pairs, part, by = c("dataset", "scenario"), sort = TRUE)
    }
  }
  if (nrow(pairs) != 16L) stop("A historical or new pair is missing")
  list(pairs = pairs, metrics = metrics)
}

pp_sv_bootstrap <- function(pairs, metrics, B = 10000L, seed = 2026092404L) {
  scenarios <- paste0("S", 1:4)
  scopes <- list(all = scenarios, n100 = c("S1", "S2"),
    n1000 = c("S3", "S4"), noncrossing = c("S1", "S3"),
    crossing = c("S2", "S4"), S1 = "S1", S2 = "S2", S3 = "S3", S4 = "S4")
  if (nrow(pairs) != 16L || anyDuplicated(pairs$dataset) ||
      !all(table(factor(pairs$scenario, levels = scenarios)) == 4L))
    stop("Expected four assigned datasets per scenario")
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  picks <- lapply(scenarios, function(s)
    matrix(sample.int(4L, 4L * B, replace = TRUE), nrow = 4L))
  result <- list()
  for (contrast in c("new_PDF_minus_PNG300", "new_minus_old_PDF_PNG_difference")) {
    for (metric in metrics) {
      means <- rep(NA_real_, 4L)
      complete <- integer(4L)
      draws <- matrix(NA_real_, B, 4L)
      bounds <- matrix(NA_real_, 4L, 2L)
      for (i in 1:4) {
        block <- pairs[pairs$scenario == scenarios[i], ]
        pdf <- block[[paste0(metric, "_new_PDF")]]
        png <- block[[paste0(metric, "_new_PNG300")]]
        d <- pdf - png
        if (contrast == "new_minus_old_PDF_PNG_difference")
          d <- d - (block[[paste0(metric, "_old_PDF")]] -
                    block[[paste0(metric, "_old_PNG300")]])
        complete[i] <- sum(is.finite(d))
        if (complete[i]) means[i] <- mean(d[is.finite(d)])
        sampled <- matrix(d[picks[[i]]], nrow = 4L)
        counts <- colSums(is.finite(sampled))
        sampled[!is.finite(sampled)] <- 0
        draws[, i] <- ifelse(counts > 0, colSums(sampled) / counts, NA_real_)
        if (metric == "taae_mean" && contrast == "new_PDF_minus_PNG300") {
          bounds[i, ] <- c(
            mean(ifelse(is.finite(pdf), pdf, 0) -
                   ifelse(is.finite(png), png, 1)),
            mean(ifelse(is.finite(pdf), pdf, 1) -
                   ifelse(is.finite(png), png, 0)))
        }
      }
      for (scope in names(scopes)) {
        use <- match(scopes[[scope]], scenarios)
        estimate <- if (all(is.finite(means[use]))) mean(means[use]) else NA_real_
        x <- rowMeans(draws[, use, drop = FALSE])
        undefined <- sum(!is.finite(x))
        degenerate <- !undefined && length(unique(x)) == 1L
        status <- if (!is.finite(estimate)) "undefined_estimate" else if (undefined)
          "undefined_resamples" else if (degenerate) "degenerate_resamples" else "available"
        available <- status == "available"
        result[[length(result) + 1L]] <- data.frame(
          contrast = contrast, metric = metric, scope = scope,
          estimate = estimate, complete_pairs = sum(complete[use]),
          assigned_pairs = 4L * length(use),
          lower = if (available) unname(quantile(x, .025, type = 1)) else NA_real_,
          upper = if (available) unname(quantile(x, .975, type = 1)) else NA_real_,
          standard_error = if (available) sd(x) else NA_real_,
          undefined_resamples = undefined, status = status,
          full_assigned_lower = if (metric == "taae_mean" &&
                                    contrast == "new_PDF_minus_PNG300")
            mean(bounds[use, 1]) else NA_real_,
          full_assigned_upper = if (metric == "taae_mean" &&
                                    contrast == "new_PDF_minus_PNG300")
            mean(bounds[use, 2]) else NA_real_, B = B, seed = seed)
      }
    }
  }
  do.call(rbind, result)
}

B <- 10000L
seed <- 2026092707L
bootstrap <- function(p, delta, scenario_only) {
  groups <- if (scenario_only) p$scenario else interaction(p$cohort, p$scenario, drop = TRUE)
  strata <- split(seq_len(nrow(p)), groups)
  size <- if (scenario_only) 16L else 8L
  stopifnot(nrow(p) == 64L, !anyDuplicated(p$dataset),
            all(lengths(strata) == size), all(is.finite(delta)))
  RNGkind('Mersenne-Twister', 'Inversion', 'Rejection')
  set.seed(seed)
  picks <- lapply(strata, function(ids) matrix(sample.int(size, size * B, replace = TRUE), nrow = size))
  draws <- sapply(seq_along(strata), function(i) {
    d <- delta[strata[[i]]]
    colMeans(matrix(d[picks[[i]]], nrow = size))
  })
  x <- rowMeans(draws)
  available <- length(unique(x)) > 1L
  c(estimate = mean(vapply(strata, function(ids) mean(delta[ids]), numeric(1))),
    lower = if (available) unname(quantile(x, 0.025, type = 1)) else NA_real_,
    upper = if (available) unname(quantile(x, 0.975, type = 1)) else NA_real_,
    standard_error = if (available) sd(x) else NA_real_)
}

main <- function() {
  scores <- read.csv("build/scores.csv",check.names=FALSE)
  scores$dataset <- scores$dataset_id
  scores$representation <- scores$condition
  scores$primary_attempt <- scores$is_first_run
  primary <- scores[scores$primary64, ]
  original <- pp_pairs(primary[primary$batch=="first", ])
  second <- pp_pairs(primary[primary$batch=="second", ])
  metrics <- original$metrics
  original$pairs$cohort <- "original"
  second$pairs$cohort <- "replication"
  pairs <- rbind(original$pairs,second$pairs)
  write.csv(pairs,"build/primary_pairs.csv",row.names=FALSE,na="")
  out <- list()
  record <- function(family, metric, contrast, left, right) {
    result <- bootstrap(pairs,left-right,TRUE)
    scale <- if(grepl("^(taae|maximum)",metric)) 100 else 1
    out[[length(out)+1L]] <<- data.frame(family=family,metric=metric,contrast=contrast,
      scale=scale,mean_left=mean(left)*scale,mean_right=mean(right)*scale,
      pairs=nrow(pairs),B=B,seed=seed,
      new_estimate=result[["estimate"]]*scale,new_lower=result[["lower"]]*scale,
      new_upper=result[["upper"]]*scale,new_standard_error=result[["standard_error"]]*scale,
      new_status=if(is.finite(result[["lower"]])) "available" else "degenerate_resamples")
  }
  for (m in metrics) record("primary",m,"PDF_minus_PNG300",pairs[[paste0(m,"_PDF")]],pairs[[paste0(m,"_PNG300")]])
  guyot <- scores[scores$method=="Guyot", ]
  guyot <- guyot[match(pairs$dataset,guyot$dataset), ]
  stopifnot(identical(guyot$dataset,pairs$dataset))
  contrasts <- list(Guyot_absolute=c("Guyot","zero"),LLM_PNG300_absolute=c("PNG300","zero"),
                    LLM_PDF_absolute=c("PDF","zero"),Guyot_minus_LLM_PNG300=c("Guyot","PNG300"),
                    LLM_PDF_minus_Guyot=c("PDF","Guyot"))
  for (name in names(contrasts)) for (m in metrics) {
    values <- list(Guyot=guyot[[m]], PDF=pairs[[paste0(m,"_PDF")]], PNG300=pairs[[paste0(m,"_PNG300")]],zero=rep(0,nrow(pairs)))
    record("guyot",m,name,values[[contrasts[[name]][1]]],values[[contrasts[[name]][2]]])
  }
  write.csv(do.call(rbind,out),"build/contrasts_current64.csv",row.names=FALSE,na="")
  cohort_results <- list(
    cbind(analysis="first_batch",pp_bootstrap(original$pairs,metrics,seed=2026092306L)),
    cbind(analysis="second_completed",pp_bootstrap(second$pairs,metrics,seed=2026092706L)))
  # Reconstruct the original assigned tables: missing sessions retain their slots.
  for (population in c("complete_case62","original61")) {
    historic <- primary
    # Remove only the unavailable original runs, retaining their observed partners.
    # Pair inclusion flags mark both partners and cannot be used to mask individual runs.
    unavailable <- historic$run_id %in% c("RP32_R014R1", "RP32_R064R1")
    if (population=="original61") unavailable <- unavailable | historic$run_id=="RP32_R003R1"
    for (m in metrics) historic[[m]][unavailable] <- NA_real_
    historic$status[unavailable] <- "started_unscorable"
    hp <- pp_pairs(historic[historic$batch=="second", ])$pairs
    result <- pooled_bootstrap(original$pairs,hp,metrics)
    name <- if(population=="complete_case62") "historical62" else "historical61"
    cohort_results[[length(cohort_results)+1L]] <- cbind(analysis=name,pp_bootstrap(hp,metrics,seed=2026092706L))
    write.csv(result,paste0("build/",name,".csv"),row.names=FALSE,na="")
    write.csv(rbind(transform(original$pairs,cohort=NULL),hp),paste0("build/",name,"_pairs.csv"),row.names=FALSE,na="")
  }
  write.csv(do.call(rbind,cohort_results),"build/cohort_sensitivity.csv",row.names=FALSE,na="")
  first <- primary[primary$batch=="first", ]
  raster <- pp_three_format_pairs(rbind(first,scores[scores$experiment=="image_only", ]))
  write.csv(raster$pairs,"build/image_only_pairs.csv",row.names=FALSE,na="")
  write.csv(pp_three_format_bootstrap(raster$pairs,raster$metrics),"build/image_only_bootstrap.csv",row.names=FALSE,na="")
  sv <- scores[scores$experiment=="survminer", ]
  template <- pp_sv_pairs(sv,first,unique(sv$dataset))
  write.csv(template$pairs,"build/survminer_pairs.csv",row.names=FALSE,na="")
  write.csv(pp_sv_bootstrap(template$pairs,template$metrics),"build/survminer_bootstrap.csv",row.names=FALSE,na="")
  cat("Computed 64-pair current, 62/61-pair historical, image-only and survminer contrasts.\n")
}
if(sys.nframe()==0L) main()
