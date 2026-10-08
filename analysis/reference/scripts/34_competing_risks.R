source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))
suppressPackageStartupMessages(library(parallel))

log_file <- file.path(CIPDS_LOG, "34_competing_risks.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Replicates:", CIPDS_REPLICATES, " Workers:", CIPDS_THREADS, "\n")

x <- cipds_load_data(FALSE)
dat <- x$dat
common <- cipds_complete_mask(dat, unname(CIPDS_SCORE_VARS))
d <- droplevels(copy(dat[common]))
expected <- fread(file.path(CIPDS_ROOT,"qa/cohort_registry.csv"))[cohort=="common"]
stopifnot(nrow(d)==expected$n, sum(d$Death_AllCause)==expected$events)

causes <- c(CANCER = "competing_cancer", CVD = "competing_cvd")
horizons <- c(5, 10, 15)

weighted_cif_curve <- function(time, status, weight, max_time = Inf) {
  keep <- is.finite(time) & !is.na(status) & is.finite(weight) & weight > 0
  time <- time[keep]; status <- status[keep]; weight <- weight[keep]
  event_times <- sort(unique(time[status > 0 & time <= max_time]))
  if (!length(event_times)) return(data.table(time = 0, cif = 0, survival = 1))
  survival <- 1
  cif <- 0
  rows <- vector("list", length(event_times) + 1L)
  rows[[1]] <- data.table(time = 0, cif = 0, survival = 1)
  for (j in seq_along(event_times)) {
    tt <- event_times[j]
    risk <- sum(weight[time >= tt])
    d_target <- sum(weight[time == tt & status == 1])
    d_all <- sum(weight[time == tt & status > 0])
    cif <- cif + survival * d_target / risk
    survival <- survival * (1 - d_all / risk)
    rows[[j + 1L]] <- data.table(time = tt, cif = cif, survival = survival)
  }
  rbindlist(rows)
}

weighted_cif_at <- function(time, status, weight, horizons) {
  curve <- weighted_cif_curve(time, status, weight, max(horizons))
  vapply(horizons, function(h) {
    j <- findInterval(h, curve$time)
    if (j == 0L) 0 else curve$cif[j]
  }, numeric(1))
}

make_finegray <- function(data_obj, status_var, score_var) {
  dd <- copy(data_obj[!is.na(get(status_var))])
  dd[, competing_factor := factor(
    get(status_var), levels = c(0, 1, 2), labels = c("censor", "target", "competing")
  )]
  carry <- c(score_var, CIPDS_CLINICAL_TERMS, "pooled_mec_weight", "survey_psu", "survey_strata")
  fg_formula <- as.formula(paste(
    "Surv(Follow_Up_Years, competing_factor) ~", paste(carry, collapse = " + ")
  ))
  as.data.table(finegray(fg_formula, data = dd, etype = "target"))
}

fit_finegray <- function(data_obj, status_var, score_var) {
  fg <- make_finegray(data_obj, status_var, score_var)
  fg[, analysis_weight := pooled_mec_weight * fgwt]
  fg[, analysis_weight := analysis_weight / mean(analysis_weight)]
  formula <- as.formula(paste(
    "Surv(fgstart, fgstop, fgstatus) ~",
    paste(c(score_var, CIPDS_CLINICAL_TERMS), collapse = " + ")
  ))
  fit <- coxph(
    formula, data = fg, weights = analysis_weight,
    cluster = survey_psu, robust = TRUE, x = TRUE, y = TRUE, model = TRUE
  )
  list(fit = fit, finegray_n = nrow(fg))
}

finegray_rows <- list()
for (cause in names(causes)) {
  cause_data <- d[!is.na(get(causes[[cause]]))]
  for (score in names(CIPDS_SCORE_VARS)) {
    score_var <- CIPDS_SCORE_VARS[[score]]
    fg <- fit_finegray(cause_data, causes[[cause]], score_var)
    b <- coef(fg$fit)[[score_var]]
    se <- sqrt(vcov(fg$fit)[score_var, score_var])
    finegray_rows[[length(finegray_rows) + 1L]] <- data.table(
      cause = cause, status_variable = causes[[cause]], score = score,
      score_label = CIPDS_SCORE_LABELS[[score]], n = nrow(cause_data),
      target_events = sum(cause_data[[causes[[cause]]]] == 1),
      competing_events = sum(cause_data[[causes[[cause]]]] == 2),
      finegray_expanded_rows = fg$finegray_n,
      subdistribution_hazard_ratio = exp(b),
      ci_lower = exp(b - 1.96 * se), ci_upper = exp(b + 1.96 * se),
      p_value = 2 * pnorm(-abs(b / se)),
      variance_method = "NHANES sampling weights multiplied by Fine-Gray weights; PSU-cluster robust covariance"
    )
  }
}
finegray_effects <- rbindlist(finegray_rows)
finegray_effects[, p_holm := p.adjust(p_value, method = "holm"), by = cause]
fwrite(finegray_effects, file.path(CIPDS_OUT, "competing_risk_finegray_effects.csv"), bom = TRUE)

# Survey-weighted nonparametric cumulative incidence by score quartile.
quartile_registry <- list()
curve_rows <- list()
for (score in names(CIPDS_SCORE_VARS)) {
  score_var <- CIPDS_SCORE_VARS[[score]]
  cuts <- unique(c(-Inf, cipds_weighted_quantile(d[[score_var]], d$pooled_mec_weight, 1:3 / 4), Inf))
  if (length(cuts) != 5L) stop("Quartile collapse: ", score)
  quartile_registry[[score]] <- cut(d[[score_var]], breaks = cuts, labels = FALSE, include.lowest = TRUE)
  for (cause in names(causes)) {
    for (q in 1:4) {
      mask_q <- quartile_registry[[score]] == q & !is.na(d[[causes[[cause]]]])
      curve <- weighted_cif_curve(
        d$Follow_Up_Years[mask_q], d[[causes[[cause]]]][mask_q], d$pooled_mec_weight[mask_q], 15
      )
      curve[, `:=`(
        score = score, score_label = CIPDS_SCORE_LABELS[[score]], cause = cause,
        quartile = q, n = sum(mask_q), target_events = sum(d[[causes[[cause]]]][mask_q] == 1)
      )]
      curve_rows[[length(curve_rows) + 1L]] <- curve
    }
  }
}
cif_curves <- rbindlist(curve_rows)
fwrite(cif_curves, file.path(CIPDS_OUT, "competing_risk_cif_curves.csv"), bom = TRUE)

evaluate_cif_horizons <- function(weight) {
  ans <- c()
  for (score in names(CIPDS_SCORE_VARS)) {
    groups <- quartile_registry[[score]]
    for (cause in names(causes)) {
      status <- d[[causes[[cause]]]]
      for (q in 1:4) {
        in_q <- groups == q & !is.na(status)
        values <- weighted_cif_at(
          d$Follow_Up_Years[in_q], status[in_q], weight[in_q], horizons
        )
        for (j in seq_along(horizons)) {
          ans[paste(score, cause, paste0("q", q), paste0("h", horizons[j]), sep = "__")] <- values[j]
        }
      }
    }
  }
  ans
}

point <- evaluate_cif_horizons(d$pooled_mec_weight)
set.seed(CIPDS_SEED + 134L)
rep_full <- as.svrepdesign(
  x$design_full, type = "bootstrap", replicates = CIPDS_REPLICATES, mse = TRUE
)
rep_d <- rep_full[common, ]
rep_weights <- weights(rep_d, type = "analysis")
cluster <- makeCluster(min(CIPDS_THREADS, CIPDS_REPLICATES))
on.exit(try(stopCluster(cluster), silent = TRUE), add = TRUE)
clusterEvalQ(cluster, {suppressPackageStartupMessages(library(data.table)); NULL})
clusterExport(cluster, c(
  "d", "causes", "horizons", "quartile_registry", "rep_weights", "point", "CIPDS_SCORE_VARS",
  "weighted_cif_curve", "weighted_cif_at", "evaluate_cif_horizons"
), envir = environment())
rep_list <- parLapplyLB(cluster, seq_len(CIPDS_REPLICATES), function(b) {
  tryCatch(evaluate_cif_horizons(rep_weights[, b]), error = function(e) {
    z <- rep(NA_real_, length(point)); names(z) <- names(point); z
  })
})
stopCluster(cluster)
cluster <- NULL
rep_matrix <- do.call(rbind, rep_list)
colnames(rep_matrix) <- names(point)
finite_fraction <- colMeans(is.finite(rep_matrix))
if (any(finite_fraction < 0.98)) stop("CIF bootstrap failure rate exceeded 2%")
for (j in seq_along(point)) rep_matrix[!is.finite(rep_matrix[, j]), j] <- point[j]
variance <- svrVar(rep_matrix, scale = rep_d$scale, rscales = rep_d$rscales,
                   mse = rep_d$mse, coef = point)
se <- sqrt(diag(variance))
df <- degf(x$design_full[common, ])
crit <- qt(0.975, df)

parse_cif_key <- function(z) {
  p <- strsplit(z, "__", fixed = TRUE)[[1]]
  list(score = p[1], cause = p[2], quartile = as.integer(sub("q", "", p[3])),
       horizon = as.integer(sub("h", "", p[4])))
}
cif_horizon_rows <- rbindlist(lapply(seq_along(point), function(j) {
  p <- parse_cif_key(names(point)[j])
  data.table(
    score = p$score, score_label = CIPDS_SCORE_LABELS[[p$score]], cause = p$cause,
    quartile = p$quartile, horizon_years = p$horizon,
    cumulative_incidence = point[j], standard_error = se[j],
    ci_lower = max(0, point[j] - crit * se[j]), ci_upper = min(1, point[j] + crit * se[j]),
    finite_replicate_fraction = finite_fraction[j]
  )
}))
fwrite(cif_horizon_rows, file.path(CIPDS_OUT, "competing_risk_cif_horizon_estimates.csv"), bom = TRUE)

contrast_rows <- list()
for (score in names(CIPDS_SCORE_VARS)) {
  for (cause in names(causes)) {
    for (h in horizons) {
      k4 <- paste(score, cause, "q4", paste0("h", h), sep = "__")
      k1 <- paste(score, cause, "q1", paste0("h", h), sep = "__")
      delta <- point[[k4]] - point[[k1]]
      rep_delta <- rep_matrix[, k4] - rep_matrix[, k1]
      v <- as.numeric(svrVar(matrix(rep_delta, ncol = 1), scale = rep_d$scale,
                             rscales = rep_d$rscales, mse = rep_d$mse, coef = delta))
      s <- sqrt(v)
      contrast_rows[[length(contrast_rows) + 1L]] <- data.table(
        score = score, cause = cause, horizon_years = h,
        q4_minus_q1_absolute_risk_difference = delta,
        ci_lower = delta - crit * s, ci_upper = delta + crit * s,
        p_value = 2 * pt(-abs(delta / s), df = df)
      )
    }
  }
}
cif_contrasts <- rbindlist(contrast_rows)
cif_contrasts[, p_holm := p.adjust(p_value, method = "holm"), by = cause]
fwrite(cif_contrasts, file.path(CIPDS_OUT, "competing_risk_cif_q4_q1_contrasts.csv"), bom = TRUE)

# Temporal absolute-risk calibration of subdistribution-risk mappings.
train <- droplevels(copy(d[CYCLE %in% c("2003-2004", "2005-2006")]))
valid <- droplevels(copy(d[CYCLE %in% c("2007-2008", "2009-2010")]))
calibration_rows <- list()
for (cause in names(causes)) {
  for (score in c("OVERALL", "PHENO")) {
    score_var <- CIPDS_SCORE_VARS[[score]]
    fitted <- fit_finegray(train, causes[[cause]], score_var)$fit
    bh <- basehaz(fitted, centered = TRUE)
    lp <- as.numeric(predict(fitted, newdata = valid, type = "lp", reference = "sample"))
    for (h in c(5, 10)) {
      j <- findInterval(h, bh$time)
      h0 <- if (j == 0L) 0 else bh$hazard[j]
      predicted <- 1 - exp(-h0 * exp(lp))
      cuts <- unique(c(-Inf, cipds_weighted_quantile(predicted, valid$pooled_mec_weight, 1:4 / 5), Inf))
      groups <- cut(predicted, breaks = cuts, labels = FALSE, include.lowest = TRUE)
      for (g in 1:5) {
        in_g <- groups == g
        observed <- weighted_cif_at(
          valid$Follow_Up_Years[in_g], valid[[causes[[cause]]]][in_g], valid$pooled_mec_weight[in_g], h
        )
        expected <- weighted.mean(predicted[in_g], valid$pooled_mec_weight[in_g])
        calibration_rows[[length(calibration_rows) + 1L]] <- data.table(
          cause = cause, score = score, horizon_years = h, risk_group = g,
          n = sum(in_g), target_events = sum(valid[[causes[[cause]]]][in_g] == 1),
          predicted_cumulative_incidence = expected,
          observed_weighted_cumulative_incidence = observed,
          observed_minus_predicted = observed - expected,
          training_cycles = "2003-2006", validation_cycles = "2007-2010"
        )
      }
    }
  }
}
calibration <- rbindlist(calibration_rows)
fwrite(calibration, file.path(CIPDS_OUT, "competing_risk_temporal_absolute_calibration.csv"), bom = TRUE)
saveRDS(list(point = point, replicates = rep_matrix, df = df),
        file.path(CIPDS_OUT, "competing_risk_cif_replicates.rds"), compress = "gzip")

cipds_write_manifest(file.path(CIPDS_OUT, "competing_risk_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  outcomes = list(cancer = "UCOD 2", cardiovascular = "UCOD 1 or 5; 1999-2014"),
  primary_all_cause_outcome_unchanged = TRUE,
  association_sensitivity = "Fine-Gray subdistribution Cox via finegray expansion; NHANES weights times finegray weights; PSU-cluster robust covariance",
  absolute_risk = "survey-weighted nonparametric cumulative incidence by score quartile",
  uncertainty = "1000 survey bootstrap replicate weights for horizon cumulative incidence and Q4-Q1 risk differences",
  temporal_calibration = "early-cycle subdistribution risk mapping applied without refitting to 2007-2010; observed weighted CIF versus predicted CIF quintiles",
  limitation = "Fine-Gray inference does not implement a full stratified multistage survey variance estimator and remains a sensitivity analysis",
  threads = CIPDS_THREADS, seed = CIPDS_SEED + 134L
))
cat("\nFine-Gray effects:\n")
print(finegray_effects)
cat("\nQ4-Q1 cumulative-incidence contrasts:\n")
print(cif_contrasts)
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
