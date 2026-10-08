source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))
suppressPackageStartupMessages(library(parallel))

log_file <- file.path(CIPDS_LOG, "30_incremental_prediction.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Replicates:", CIPDS_REPLICATES, " Workers:", CIPDS_THREADS, "\n")

x <- cipds_load_data(FALSE)
dat <- x$dat
common <- cipds_complete_mask(dat, unname(CIPDS_SCORE_VARS))
d <- copy(dat[common])
expected <- fread(file.path(CIPDS_ROOT,"qa/cohort_registry.csv"))[cohort=="common"]
stopifnot(nrow(d)==expected$n, sum(d$Death_AllCause)==expected$events)

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
metric_grid <- seq(1, 15, by = 1)

fit_one_prediction_model <- function(data, weight, predictors, metric_grid, horizons) {
  keep <- is.finite(weight) & weight > 0
  dd <- copy(data[keep])
  ww <- weight[keep] / mean(weight[keep])
  dd[, .analysis_weight := ww]
  formula <- cipds_make_formula(predictors)
  fit <- coxph(
    formula, data = dd, weights = .analysis_weight, ties = "efron",
    x = TRUE, y = TRUE, model = TRUE,
    control = coxph.control(iter.max = 40)
  )
  pred <- cipds_predict_risk(fit, data, metric_grid)
  risk_h <- pred$risk[, match(horizons, metric_grid), drop = FALSE]
  colnames(risk_h) <- paste0("risk", horizons)
  list(fit = fit, lp = pred$lp, risk_grid = pred$risk, risk_h = risk_h)
}

evaluate_one_model <- function(data, weight, predictors, metric_grid, horizons) {
  pred <- fit_one_prediction_model(data, weight, predictors, metric_grid, horizons)
  time <- data$Follow_Up_Years
  event <- data$Death_AllCause
  ans <- c(uno_c15 = cipds_weighted_uno(time, event, pred$lp, weight, 15))
  for (h in horizons) {
    j <- match(h, metric_grid)
    ans[paste0("auc", h)] <- cipds_weighted_td_auc(time, event, pred$lp, weight, h)
    ans[paste0("brier", h)] <- cipds_weighted_brier(time, event, pred$risk_grid[, j], weight, h)
  }
  brier_grid <- vapply(metric_grid, function(h) {
    cipds_weighted_brier(time, event, pred$risk_grid[, match(h, metric_grid)], weight, h)
  }, numeric(1))
  ans["ibs15"] <- sum(diff(metric_grid) * (head(brier_grid, -1) + tail(brier_grid, -1)) / 2) /
    (max(metric_grid) - min(metric_grid))
  ans
}

evaluate_all_models <- function(weight) {
  values <- c()
  for (model in names(model_sets)) {
    z <- evaluate_one_model(d, weight, model_sets[[model]], metric_grid, CIPDS_HORIZONS)
    names(z) <- paste0(names(z), "__", model)
    values <- c(values, z)
  }
  values
}

point_weight <- d$pooled_mec_weight
point <- evaluate_all_models(point_weight)
cat("Point estimates calculated.\n")

set.seed(CIPDS_SEED + 30L)
rep_full <- as.svrepdesign(
  x$design_full, type = "bootstrap", replicates = CIPDS_REPLICATES, mse = TRUE
)
rep_sub <- rep_full[common, ]
rep_weights <- weights(rep_sub, type = "analysis")
if (nrow(rep_weights) != nrow(d) || ncol(rep_weights) != CIPDS_REPLICATES) {
  stop("Replicate-weight dimension drift")
}

n_workers <- min(CIPDS_THREADS, CIPDS_REPLICATES)
cluster <- makeCluster(n_workers)
on.exit(try(stopCluster(cluster), silent = TRUE), add = TRUE)
clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })
  NULL
})
clusterExport(
  cluster,
  c(
    "d", "rep_weights", "point", "model_sets", "metric_grid", "CIPDS_HORIZONS",
    "CIPDS_CLINICAL_TERMS", "cipds_make_formula", "cipds_weighted_km_censoring",
    "cipds_weighted_td_auc", "cipds_weighted_uno", "cipds_weighted_brier",
    "cipds_baseline_hazard", "cipds_predict_risk", "fit_one_prediction_model",
    "evaluate_one_model", "evaluate_all_models"
  ),
  envir = environment()
)
replicate_list <- parLapplyLB(cluster, seq_len(CIPDS_REPLICATES), function(b) {
  tryCatch(
    evaluate_all_models(rep_weights[, b]),
    error = function(e) {
      z <- rep(NA_real_, length(point))
      names(z) <- names(point)
      attr(z, "error") <- conditionMessage(e)
      z
    }
  )
})
stopCluster(cluster)
cluster <- NULL
replicate_matrix <- do.call(rbind, replicate_list)
colnames(replicate_matrix) <- names(point)
finite_fraction <- colMeans(is.finite(replicate_matrix))
if (any(finite_fraction < 0.98)) {
  bad <- names(finite_fraction)[finite_fraction < 0.98]
  stop("Too many failed bootstrap estimates: ", paste(bad, collapse = ", "))
}
# Extremely rare failed replicates are conservatively centered at the point
# estimate so the fixed replicate design and its variance scaling are retained.
for (j in seq_along(point)) {
  replicate_matrix[!is.finite(replicate_matrix[, j]), j] <- point[j]
}

variance <- svrVar(
  replicate_matrix, scale = rep_sub$scale, rscales = rep_sub$rscales,
  mse = rep_sub$mse, coef = point
)
se <- sqrt(diag(variance))
df <- degf(x$design_full[common, ])
crit <- qt(0.975, df)

parse_key <- function(key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  list(metric = parts[1], model = parts[2])
}
estimate_rows <- rbindlist(lapply(seq_along(point), function(j) {
  meta <- parse_key(names(point)[j])
  bounded <- grepl("^(uno|auc|brier|ibs)", meta$metric)
  data.table(
    metric = meta$metric, model = meta$model, model_label = model_labels[[meta$model]],
    n = nrow(d), events = sum(d$Death_AllCause), estimate = point[j],
    standard_error = se[j],
    ci_lower = if (bounded) max(0, point[j] - crit * se[j]) else point[j] - crit * se[j],
    ci_upper = if (bounded) min(1, point[j] + crit * se[j]) else point[j] + crit * se[j],
    survey_bootstrap_replicates = CIPDS_REPLICATES,
    finite_replicate_fraction = finite_fraction[j]
  )
}))

comparison_pairs <- data.table(
  model_a = c(
    "M1_CLINICAL_OVERALL", "M2_CLINICAL_PHENO", "M3_CLINICAL_OVERALL_PHENO",
    "M4_CLINICAL_COMPONENTS", "M5_CLINICAL_COMPONENTS_PHENO",
    "M3_CLINICAL_OVERALL_PHENO", "M3_CLINICAL_OVERALL_PHENO",
    "M4_CLINICAL_COMPONENTS", "M5_CLINICAL_COMPONENTS_PHENO", "M5_CLINICAL_COMPONENTS_PHENO"
  ),
  model_b = c(
    rep("M0_CLINICAL", 5),
    "M2_CLINICAL_PHENO", "M1_CLINICAL_OVERALL",
    "M1_CLINICAL_OVERALL", "M4_CLINICAL_COMPONENTS", "M3_CLINICAL_OVERALL_PHENO"
  ),
  comparison = c(
    "Overall beyond clinical base", "PhenoAge beyond clinical base",
    "Overall plus PhenoAge beyond clinical base", "Three components beyond clinical base",
    "Three components plus PhenoAge beyond clinical base",
    "Overall beyond PhenoAge", "PhenoAge beyond Overall",
    "Three components versus Overall", "PhenoAge beyond three components",
    "Three components versus Overall when both include PhenoAge"
  )
)

comparison_pairs <- rbind(comparison_pairs, data.table(
  model_a="M5_CLINICAL_COMPONENTS_PHENO", model_b="M2_CLINICAL_PHENO",
  comparison="Three components beyond PhenoAge"))
stopifnot(nrow(comparison_pairs)==11L)
metrics <- unique(vapply(names(point), function(z) parse_key(z)$metric, character(1)))
comparison_rows <- list()
for (metric in metrics) {
  higher_better <- grepl("^(uno|auc)", metric)
  for (i in seq_len(nrow(comparison_pairs))) {
    a <- comparison_pairs$model_a[i]
    b <- comparison_pairs$model_b[i]
    ka <- paste0(metric, "__", a)
    kb <- paste0(metric, "__", b)
    delta <- point[[ka]] - point[[kb]]
    rep_delta <- replicate_matrix[, ka] - replicate_matrix[, kb]
    delta_var <- as.numeric(svrVar(
      matrix(rep_delta, ncol = 1), scale = rep_sub$scale, rscales = rep_sub$rscales,
      mse = rep_sub$mse, coef = delta
    ))
    delta_se <- sqrt(delta_var)
    p <- 2 * pt(-abs(delta / delta_se), df = df)
    comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
      metric = metric, higher_is_better = higher_better,
      model_a = a, model_a_label = model_labels[[a]],
      model_b = b, model_b_label = model_labels[[b]],
      comparison = comparison_pairs$comparison[i],
      estimate_a = point[[ka]], estimate_b = point[[kb]],
      paired_difference_a_minus_b = delta,
      benefit_oriented_difference = ifelse(higher_better, delta, -delta),
      difference_standard_error = delta_se,
      difference_ci_lower = delta - crit * delta_se,
      difference_ci_upper = delta + crit * delta_se,
      paired_p_value = p, n = nrow(d), events = sum(d$Death_AllCause),
      survey_bootstrap_replicates = CIPDS_REPLICATES
    )
  }
}
comparisons <- rbindlist(comparison_rows)
comparisons[, paired_p_holm := p.adjust(paired_p_value, method = "holm"), by = metric]
comparisons[, conclusion := fifelse(
  difference_ci_lower > 0,
  ifelse(higher_is_better, paste0(model_a, " better"), paste0(model_b, " better")),
  fifelse(
    difference_ci_upper < 0,
    ifelse(higher_is_better, paste0(model_b, " better"), paste0(model_a, " better")),
    "No statistically resolved difference"
  )
)]

fwrite(estimate_rows, file.path(CIPDS_OUT, "incremental_prediction_estimates.csv"), bom = TRUE)
fwrite(comparisons, file.path(CIPDS_OUT, "incremental_prediction_paired_comparisons.csv"), bom = TRUE)
saveRDS(
  list(point = point, replicates = replicate_matrix, model_sets = model_sets,
       scale = rep_sub$scale, rscales = rep_sub$rscales, mse = rep_sub$mse, df = df),
  file.path(CIPDS_OUT, "incremental_prediction_replicates.rds"), compress = "gzip"
)

cipds_write_manifest(file.path(CIPDS_OUT, "incremental_prediction_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  cohort_n = nrow(d), deaths = sum(d$Death_AllCause),
  cohort = "single common complete-case age >=60 survey domain",
  models = as.list(model_labels),
  primary_metric = "survey-weighted Uno C-index truncated at 15 years",
  secondary_metrics = c("survey-weighted cumulative/dynamic AUC at 5, 10, 15 years",
                        "IPCW Brier score at 5, 10, 15 years", "integrated Brier score from 1 to 15 years"),
  bootstrap = "1000 bootstrap replicate survey weights; every model refitted within every replicate",
  multiplicity = "Holm correction across 11 contrasts including the additional three-component comparison within each metric",
  threads = CIPDS_THREADS, seed = CIPDS_SEED + 30L
))
cat("\nPoint estimates:\n")
print(estimate_rows)
cat("\nKey comparisons:\n")
print(comparisons[comparison %chin% c(
  "Overall beyond clinical base", "PhenoAge beyond clinical base",
  "Overall beyond PhenoAge", "PhenoAge beyond Overall",
  "Three components versus Overall when both include PhenoAge"
), .(metric, comparison, paired_difference_a_minus_b, difference_ci_lower,
     difference_ci_upper, paired_p_value, paired_p_holm, conclusion)])
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
