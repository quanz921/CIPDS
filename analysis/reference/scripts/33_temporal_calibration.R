source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))
suppressPackageStartupMessages(library(parallel))

log_file <- file.path(CIPDS_LOG, "33_temporal_calibration.log")
error_audit_file <- file.path(CIPDS_OUT, "temporal_calibration_bootstrap_errors.csv")
if (file.exists(error_audit_file)) unlink(error_audit_file)
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Replicates:", CIPDS_REPLICATES, " Workers:", CIPDS_THREADS, "\n")

x <- cipds_load_data(FALSE)
dat <- x$dat
common <- cipds_complete_mask(dat, unname(CIPDS_SCORE_VARS))
train_mask <- common & dat$CYCLE %in% c("2003-2004", "2005-2006")
valid_mask <- common & dat$CYCLE %in% c("2007-2008", "2009-2010")
train <- droplevels(copy(dat[train_mask]))
valid <- droplevels(copy(dat[valid_mask]))
stopifnot(nrow(train) > 2500L, nrow(valid) > 3000L)

model_sets <- list(
  M0_CLINICAL = character(),
  M1_CLINICAL_OVERALL = "OVERALL_z_age60",
  M2_CLINICAL_PHENO = "PHENO_z_age60",
  M3_CLINICAL_OVERALL_PHENO = c("OVERALL_z_age60", "PHENO_z_age60"),
  M4_CLINICAL_COMPONENTS = c("NM_z_age60", "TB_z_age60", "TC_z_age60"),
  M5_CLINICAL_COMPONENTS_PHENO = c("NM_z_age60", "TB_z_age60", "TC_z_age60", "PHENO_z_age60")
)
model_labels <- c(
  M0_CLINICAL = "Clinical base",
  M1_CLINICAL_OVERALL = "Clinical base + Overall",
  M2_CLINICAL_PHENO = "Clinical base + PhenoAge Acceleration",
  M3_CLINICAL_OVERALL_PHENO = "Clinical base + Overall + PhenoAge Acceleration",
  M4_CLINICAL_COMPONENTS = "Clinical base + NM + TB + TC",
  M5_CLINICAL_COMPONENTS_PHENO = "Clinical base + NM + TB + TC + PhenoAge Acceleration"
)
horizons <- c(5, 10)
grid <- seq(1, 10, by = 1)

train[, normalized_weight := pooled_mec_weight / mean(pooled_mec_weight)]
predictions <- list()
fits <- list()
for (model in names(model_sets)) {
  fit <- coxph(
    cipds_make_formula(model_sets[[model]]), data = train,
    weights = normalized_weight, ties = "efron", x = TRUE, y = TRUE, model = TRUE
  )
  predictions[[model]] <- cipds_predict_risk(fit, valid, grid)
  fits[[model]] <- fit
}

ipcw_binary_data <- function(time, event, predicted_risk, weight, horizon) {
  base_keep <- is.finite(time) & !is.na(event) & is.finite(predicted_risk) &
    is.finite(weight) & weight > 0
  time <- time[base_keep]
  event <- event[base_keep]
  predicted_risk <- predicted_risk[base_keep]
  weight <- weight[base_keep]
  g <- cipds_weighted_km_censoring(time, event, weight)
  cases <- event == 1 & time <= horizon
  controls <- time > horizon
  keep <- cases | controls
  y <- as.integer(cases[keep])
  w <- rep(NA_real_, sum(keep))
  case_index <- which(cases[keep])
  control_index <- which(controls[keep])
  if (length(case_index)) {
    original_case <- which(keep)[case_index]
    w[case_index] <- weight[original_case] /
      pmax(g$before[match(time[original_case], g$time)], 1e-8)
  }
  h_idx <- findInterval(horizon, g$time)
  g_h <- max(if (h_idx == 0L) 1 else g$after[h_idx], 1e-8)
  if (length(control_index)) {
    original_control <- which(keep)[control_index]
    w[control_index] <- weight[original_control] / g_h
  }
  data.table(
    y = y, predicted_risk = pmin(pmax(predicted_risk[keep], 1e-6), 1 - 1e-6),
    weight = w / mean(w)
  )
}

calibration_intercept_slope <- function(time, event, risk, weight, horizon) {
  dd <- ipcw_binary_data(time, event, risk, weight, horizon)
  dd[, cloglog_pred := log(-log(1 - predicted_risk))]
  nll <- function(par, fixed_slope = FALSE) {
    eta <- if (fixed_slope) par[1] + dd$cloglog_pred else par[1] + par[2] * dd$cloglog_pred
    eta <- pmin(pmax(eta, -20), 10)
    probability <- pmin(pmax(1 - exp(-exp(eta)), 1e-10), 1 - 1e-10)
    -sum(dd$weight * (dd$y * log(probability) + (1 - dd$y) * log(1 - probability)))
  }
  intercept_fit <- optim(0, nll, fixed_slope = TRUE, method = "BFGS")
  slope_fit <- optim(c(0, 1), nll, fixed_slope = FALSE, method = "BFGS")
  c(intercept = intercept_fit$par[1], slope = slope_fit$par[2])
}

weighted_km_risk <- function(time, event, weight, horizon) {
  keep <- is.finite(time) & !is.na(event) & is.finite(weight) & weight > 0
  time <- time[keep]; event <- event[keep]; weight <- weight[keep]
  event_times <- sort(unique(time[event == 1 & time <= horizon]))
  if (!length(event_times)) return(0)
  increments <- vapply(event_times, function(tt) {
    sum(weight[event == 1 & time == tt]) / sum(weight[time >= tt])
  }, numeric(1))
  1 - prod(1 - pmin(pmax(increments, 0), 1))
}

group_registry <- list()
for (model in names(predictions)) {
  for (h in horizons) {
    risk <- predictions[[model]]$risk[, match(h, grid)]
    breaks <- unique(c(-Inf, cipds_weighted_quantile(risk, valid$pooled_mec_weight, 1:4 / 5), Inf))
    if (length(breaks) != 6L) stop("Calibration grouping collapsed: ", model, "/", h)
    group_registry[[paste(model, h, sep = "__")]] <- cut(
      risk, breaks = breaks, labels = FALSE, include.lowest = TRUE
    )
  }
}

evaluate_calibration <- function(weight) {
  ans <- c()
  for (model in names(predictions)) {
    pred <- predictions[[model]]
    for (h in horizons) {
      j <- match(h, grid)
      cal <- calibration_intercept_slope(
        valid$Follow_Up_Years, valid$Death_AllCause, pred$risk[, j], weight, h
      )
      ans[paste0("intercept", h, "__", model)] <- cal[["intercept"]]
      ans[paste0("slope", h, "__", model)] <- cal[["slope"]]
      ans[paste0("brier", h, "__", model)] <- cipds_weighted_brier(
        valid$Follow_Up_Years, valid$Death_AllCause, pred$risk[, j], weight, h
      )
      groups <- group_registry[[paste(model, h, sep = "__")]]
      for (g in seq_len(5)) {
        in_g <- groups == g
        ans[paste0("observed", h, "_g", g, "__", model)] <- weighted_km_risk(
          valid$Follow_Up_Years[in_g], valid$Death_AllCause[in_g], weight[in_g], h
        )
        ans[paste0("expected", h, "_g", g, "__", model)] <- weighted.mean(
          pred$risk[in_g, j], weight[in_g]
        )
      }
    }
  }
  ans
}

point <- evaluate_calibration(valid$pooled_mec_weight)
set.seed(CIPDS_SEED + 133L)
rep_full <- as.svrepdesign(
  x$design_full, type = "bootstrap", replicates = CIPDS_REPLICATES, mse = TRUE
)
rep_valid <- rep_full[valid_mask, ]
rep_weights <- weights(rep_valid, type = "analysis")
stopifnot(nrow(rep_weights) == nrow(valid), ncol(rep_weights) == CIPDS_REPLICATES)

run_one_replicate <- function(b) {
  tryCatch(evaluate_calibration(rep_weights[, b]), error = function(e) {
    z <- rep(NA_real_, length(point)); names(z) <- names(point)
    attr(z, "error") <- conditionMessage(e)
    z
  })
}
n_workers <- min(CIPDS_THREADS, CIPDS_REPLICATES)
if (n_workers == 1L) {
  rep_list <- vector("list", CIPDS_REPLICATES)
  for (b in seq_len(CIPDS_REPLICATES)) {
    cat("Calibration replicate", b, "\n"); flush.console()
    rep_list[[b]] <- run_one_replicate(b)
  }
} else {
  cluster <- makeCluster(n_workers)
  on.exit(try(stopCluster(cluster), silent = TRUE), add = TRUE)
  clusterEvalQ(cluster, {
    suppressPackageStartupMessages({library(data.table); library(survival)})
    NULL
  })
  clusterExport(cluster, c(
    "valid", "predictions", "group_registry", "rep_weights", "horizons", "grid", "point",
    "cipds_weighted_km_censoring", "cipds_weighted_brier", "ipcw_binary_data",
    "calibration_intercept_slope", "weighted_km_risk", "evaluate_calibration",
    "run_one_replicate"
  ), envir = environment())
  rep_list <- parLapplyLB(cluster, seq_len(CIPDS_REPLICATES), run_one_replicate)
  stopCluster(cluster)
  cluster <- NULL
}
rep_matrix <- do.call(rbind, rep_list)
colnames(rep_matrix) <- names(point)
cat("Replicate matrix assembled.\n"); flush.console()
error_messages <- vapply(rep_list, function(z) {
  e <- attr(z, "error"); if (is.null(e)) "" else e
}, character(1))
if (any(nzchar(error_messages))) {
  error_audit <- data.table(replicate = which(nzchar(error_messages)), error = error_messages[nzchar(error_messages)])
  fwrite(error_audit, error_audit_file, bom = TRUE)
  cat("Bootstrap errors:", paste(unique(error_audit$error), collapse = " | "), "\n")
}
finite_fraction <- colMeans(is.finite(rep_matrix))
if (any(finite_fraction < 0.98)) stop("Calibration bootstrap failure rate exceeded 2%")
for (j in seq_along(point)) rep_matrix[!is.finite(rep_matrix[, j]), j] <- point[j]
variance <- svrVar(rep_matrix, scale = rep_valid$scale, rscales = rep_valid$rscales,
                   mse = rep_valid$mse, coef = point)
cat("Replicate variance calculated.\n"); flush.console()
se <- sqrt(diag(variance))
df <- degf(x$design_full[valid_mask, ])
crit <- qt(0.975, df)

parse_key <- function(z) {
  p <- strsplit(z, "__", fixed = TRUE)[[1]]
  list(metric = p[1], model = p[2])
}
summary_rows <- rbindlist(lapply(seq_along(point), function(j) {
  p <- parse_key(names(point)[j])
  data.table(
    metric = p$metric, model = p$model, model_label = model_labels[[p$model]],
    n_validation = nrow(valid), deaths_validation = sum(valid$Death_AllCause),
    estimate = point[j], standard_error = se[j],
    ci_lower = point[j] - crit * se[j], ci_upper = point[j] + crit * se[j],
    finite_replicate_fraction = finite_fraction[j]
  )
}))

global <- summary_rows[grepl("^(intercept|slope|brier)", metric)]
global[, horizon_years := as.integer(gsub("[^0-9]", "", metric))]
global[, metric_type := gsub("[0-9]", "", metric)]
grouped <- summary_rows[grepl("^(observed|expected)", metric)]
grouped[, horizon_years := as.integer(sub("^(observed|expected)([0-9]+)_.*$", "\\2", metric))]
grouped[, group := as.integer(sub("^.*_g([0-9]+)$", "\\1", metric))]
grouped[, value_type := sub("^([a-z]+).*$", "\\1", metric)]
cat("Calibration summaries assembled.\n"); flush.console()

fwrite(global, file.path(CIPDS_OUT, "temporal_calibration_metrics.csv"), bom = TRUE)
fwrite(grouped, file.path(CIPDS_OUT, "temporal_calibration_groups.csv"), bom = TRUE)
cat("CSV outputs written.\n"); flush.console()
# The fitted models are reproducible from the registered script and are not
# serialized because coxph formula environments can retain the full session.
saveRDS(list(point = point, replicates = rep_matrix, df = df),
        file.path(CIPDS_OUT, "temporal_calibration_replicates.rds"), compress = "gzip")
cat("Replicate registry written.\n"); flush.console()

cohort_audit <- rbind(
  data.table(role = "TRAIN", cycles = "2003-2006", n = nrow(train), deaths = sum(train$Death_AllCause),
             maximum_follow_up = max(train$Follow_Up_Years)),
  data.table(role = "TEMPORAL_VALIDATION", cycles = "2007-2010", n = nrow(valid), deaths = sum(valid$Death_AllCause),
             maximum_follow_up = max(valid$Follow_Up_Years))
)
fwrite(cohort_audit, file.path(CIPDS_OUT, "temporal_calibration_cohort_audit.csv"), bom = TRUE)

cipds_write_manifest(file.path(CIPDS_OUT, "temporal_calibration_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  training = "NHANES 2003-2006", validation = "NHANES 2007-2010",
  cohort = "same complete-case participants for all six models",
  prediction_mapping = "sampling-weighted Cox model trained only in early cycles",
  time_points = horizons,
  calibration_in_the_large = "IPCW complementary-log-log calibration intercept; ideal 0",
  calibration_slope = "IPCW complementary-log-log slope; ideal 1",
  grouped_calibration = "survey-weighted Kaplan-Meier observed risk versus weighted mean prediction in validation quintiles",
  uncertainty = "1000 survey bootstrap replicate weights in temporal validation; early-cycle models fixed",
  not_estimable = "15-year calibration because validation maximum follow-up is 13.25 years",
  threads = CIPDS_THREADS, seed = CIPDS_SEED + 133L
))
cat("\nCohort audit:\n")
print(cohort_audit)
cat("\nCalibration metrics:\n")
print(global)
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
