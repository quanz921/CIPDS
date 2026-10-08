source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))
suppressPackageStartupMessages({
  library(glmnet)
  library(parallel)
})

log_file <- file.path(CIPDS_LOG, "32_matched_target_fair_experiment.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Replicates:", CIPDS_REPLICATES, " Workers:", CIPDS_THREADS, "\n")

x <- cipds_load_data(TRUE)
dat <- x$dat
feature_registry <- fread(file.path(CIPDS_ROOT, "outputs", "nested_a_primary_final_features.csv"))
labs26 <- unique(feature_registry$feature)
labs26 <- labs26[!is.na(labs26) & labs26 != ""]
pheno9 <- c(
  "lab_alb", "lab_crea", "lab_glu", "lab_hcrp", "lab_lym_pct",
  "lab_mcv", "lab_rdw", "lab_alp", "lab_wbc"
)
if (length(labs26) != 26L) stop("Expected 26 unique frozen laboratory features; observed ", length(labs26))
if (length(setdiff(pheno9, names(dat)))) stop("Missing PhenoAge laboratory fields")
labs_union <- unique(c(labs26, pheno9))

required_all <- unique(c(
  "Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", "survey_psu", "survey_strata",
  CIPDS_CLINICAL_TERMS, labs26, pheno9,
  "NM_component", "TB_component", "TC_component", "Overall_expected_burden"
))
base_eligible <- !is.na(dat$Age) & dat$Age >= 60 &
  complete.cases(dat[, ..required_all]) & dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
train_mask <- base_eligible & dat$CYCLE %in% c("2003-2004", "2005-2006")
valid_mask <- base_eligible & dat$CYCLE %in% c("2007-2008", "2009-2010")
train <- droplevels(copy(dat[train_mask]))
valid <- droplevels(copy(dat[valid_mask]))
if (nrow(train) < 2000L || nrow(valid) < 2500L) stop("Temporal fair-experiment cohort unexpectedly small")

cohort_audit <- rbind(
  data.table(role = "TRAIN", cycles = "2003-2004 + 2005-2006", n = nrow(train),
             deaths = sum(train$Death_AllCause), median_follow_up = median(train$Follow_Up_Years),
             maximum_follow_up = max(train$Follow_Up_Years)),
  data.table(role = "TEMPORAL_VALIDATION", cycles = "2007-2008 + 2009-2010", n = nrow(valid),
             deaths = sum(valid$Death_AllCause), median_follow_up = median(valid$Follow_Up_Years),
             maximum_follow_up = max(valid$Follow_Up_Years))
)
for (h in c(5, 10, 15)) {
  cohort_audit[[paste0("cases_", h, "y")]] <- c(
    sum(train$Death_AllCause == 1 & train$Follow_Up_Years < h),
    sum(valid$Death_AllCause == 1 & valid$Follow_Up_Years < h)
  )
  cohort_audit[[paste0("controls_", h, "y")]] <- c(
    sum(train$Follow_Up_Years > h), sum(valid$Follow_Up_Years > h)
  )
}
fwrite(cohort_audit, file.path(CIPDS_OUT, "fair_target_cohort_audit.csv"), bom = TRUE)

# Training-only winsorization and standardization prevents information from the
# temporal validation period entering preprocessing.
transform_registry <- list()
for (v in labs_union) {
  q <- quantile(train[[v]], c(0.01, 0.99), na.rm = TRUE, names = FALSE)
  clipped <- pmin(pmax(train[[v]], q[1]), q[2])
  mu <- mean(clipped)
  sigma <- sd(clipped)
  if (!is.finite(sigma) || sigma <= 0) stop("Invalid training SD: ", v)
  train[[paste0("z_", v)]] <- (clipped - mu) / sigma
  valid[[paste0("z_", v)]] <- (pmin(pmax(valid[[v]], q[1]), q[2]) - mu) / sigma
  transform_registry[[length(transform_registry) + 1L]] <- data.table(
    variable = v, winsor_lower = q[1], winsor_upper = q[2],
    training_mean_after_winsorization = mu, training_sd_after_winsorization = sigma
  )
}
for (v in c("NM_component", "TB_component", "TC_component")) {
  mu <- weighted.mean(train[[v]], train$pooled_mec_weight)
  sigma <- sqrt(sum(train$pooled_mec_weight * (train[[v]] - mu)^2) / sum(train$pooled_mec_weight))
  train[[paste0("zfair_", v)]] <- (train[[v]] - mu) / sigma
  valid[[paste0("zfair_", v)]] <- (valid[[v]] - mu) / sigma
  transform_registry[[length(transform_registry) + 1L]] <- data.table(
    variable = v, winsor_lower = NA_real_, winsor_upper = NA_real_,
    training_mean_after_winsorization = mu, training_sd_after_winsorization = sigma
  )
}
fwrite(rbindlist(transform_registry), file.path(CIPDS_OUT, "fair_target_training_transform_registry.csv"), bom = TRUE)

combined <- rbindlist(list(train, valid), use.names = TRUE, fill = TRUE)
clinical_formula <- as.formula(paste("~", paste(CIPDS_CLINICAL_TERMS, collapse = " + ")))
clinical_matrix <- model.matrix(clinical_formula, data = combined)[, -1, drop = FALSE]
train_index <- seq_len(nrow(train))
valid_index <- nrow(train) + seq_len(nrow(valid))
clinical_train <- clinical_matrix[train_index, , drop = FALSE]
clinical_valid <- clinical_matrix[valid_index, , drop = FALSE]

feature_sets <- list(
  LAB26 = paste0("z_", labs26),
  COMPONENT3 = paste0("zfair_", c("NM_component", "TB_component", "TC_component")),
  PHENO9 = paste0("z_", pheno9)
)
model_labels <- c(
  CLINICAL = "Clinical base",
  LAB26 = "Clinical base + 26 frozen-model laboratory variables",
  COMPONENT3 = "Clinical base + NM + TB + TC",
  PHENO9 = "Clinical base + 9 PhenoAge laboratory variables"
)

set.seed(CIPDS_SEED + 32L)
psu_table <- unique(train[, .(survey_psu, CYCLE)])
psu_table[, fold := sample(rep(seq_len(5), length.out = .N)), by = CYCLE]
fold_id <- psu_table$fold[match(train$survey_psu, psu_table$survey_psu)]
if (anyNA(fold_id) || length(unique(fold_id)) != 5L) stop("PSU-level fold assignment failed")

fit_registry <- list()
prediction_registry <- list()
coefficient_rows <- list()
grid <- seq(1, 10, by = 1)

weighted_breslow_h0 <- function(time, event, lp, weight, horizons) {
  event_times <- sort(unique(time[event == 1 & time <= max(horizons)]))
  increment <- vapply(event_times, function(tt) {
    numerator <- sum(weight[event == 1 & time == tt])
    denominator <- sum(weight[time >= tt] * exp(lp[time >= tt]))
    ifelse(denominator > 0, numerator / denominator, 0)
  }, numeric(1))
  cumulative <- cumsum(increment)
  vapply(horizons, function(h) {
    j <- findInterval(h, event_times)
    if (j == 0L) 0 else cumulative[j]
  }, numeric(1))
}

train[, normalized_weight := pooled_mec_weight / mean(pooled_mec_weight)]
base_formula <- cipds_make_formula(character())
base_fit <- coxph(
  base_formula, data = train, weights = normalized_weight,
  x = TRUE, y = TRUE, model = TRUE, ties = "efron"
)
base_pred_train <- cipds_predict_risk(base_fit, train, grid)
base_pred_valid <- cipds_predict_risk(base_fit, valid, grid)
fit_registry$CLINICAL <- base_fit
prediction_registry$CLINICAL <- base_pred_valid

for (model in names(feature_sets)) {
  features <- feature_sets[[model]]
  x_train <- cbind(clinical_train, as.matrix(train[, ..features]))
  x_valid <- cbind(clinical_valid, as.matrix(valid[, ..features]))
  penalty <- c(rep(0, ncol(clinical_train)), rep(1, length(features)))
  cvfit <- cv.glmnet(
    x = x_train, y = Surv(train$Follow_Up_Years, train$Death_AllCause),
    family = "cox", weights = train$normalized_weight, foldid = fold_id,
    alpha = 0.5, nfolds = 5, type.measure = "deviance",
    standardize = TRUE, penalty.factor = penalty, parallel = FALSE
  )
  lp_train <- as.numeric(predict(cvfit, newx = x_train, s = "lambda.min", type = "link"))
  lp_valid <- as.numeric(predict(cvfit, newx = x_valid, s = "lambda.min", type = "link"))
  h0 <- weighted_breslow_h0(
    train$Follow_Up_Years, train$Death_AllCause, lp_train, train$normalized_weight, grid
  )
  risk_valid <- sapply(h0, function(h) 1 - exp(-h * exp(lp_valid)))
  colnames(risk_valid) <- paste0("risk", grid)
  prediction_registry[[model]] <- list(lp = lp_valid, risk = risk_valid)
  fit_registry[[model]] <- list(cvfit = cvfit, baseline_hazard = h0)
  coefs <- as.matrix(coef(cvfit, s = "lambda.min"))[, 1]
  coefficient_rows[[length(coefficient_rows) + 1L]] <- data.table(
    model = model, variable = names(coefs), coefficient = as.numeric(coefs),
    selected = abs(coefs) > 0, lambda_min = cvfit$lambda.min, lambda_1se = cvfit$lambda.1se
  )
}
fwrite(rbindlist(coefficient_rows), file.path(CIPDS_OUT, "fair_target_model_coefficients.csv"), bom = TRUE)

evaluate_fixed_predictions <- function(weight) {
  ans <- c()
  for (model in names(prediction_registry)) {
    pred <- prediction_registry[[model]]
    z <- c(uno_c10 = cipds_weighted_uno(
      valid$Follow_Up_Years, valid$Death_AllCause, pred$lp, weight, 10
    ))
    for (h in c(5, 10)) {
      z[paste0("auc", h)] <- cipds_weighted_td_auc(
        valid$Follow_Up_Years, valid$Death_AllCause, pred$lp, weight, h
      )
      z[paste0("brier", h)] <- cipds_weighted_brier(
        valid$Follow_Up_Years, valid$Death_AllCause, pred$risk[, match(h, grid)], weight, h
      )
    }
    brier_grid <- vapply(grid, function(h) cipds_weighted_brier(
      valid$Follow_Up_Years, valid$Death_AllCause, pred$risk[, match(h, grid)], weight, h
    ), numeric(1))
    z["ibs10"] <- sum(diff(grid) * (head(brier_grid, -1) + tail(brier_grid, -1)) / 2) /
      (max(grid) - min(grid))
    names(z) <- paste0(names(z), "__", model)
    ans <- c(ans, z)
  }
  ans
}

point <- evaluate_fixed_predictions(valid$pooled_mec_weight)
set.seed(CIPDS_SEED + 132L)
rep_full <- as.svrepdesign(
  x$design_full, type = "bootstrap", replicates = CIPDS_REPLICATES, mse = TRUE
)
rep_valid <- rep_full[valid_mask, ]
rep_weights <- weights(rep_valid, type = "analysis")
stopifnot(nrow(rep_weights) == nrow(valid), ncol(rep_weights) == CIPDS_REPLICATES)

cluster <- makeCluster(min(CIPDS_THREADS, CIPDS_REPLICATES))
on.exit(try(stopCluster(cluster), silent = TRUE), add = TRUE)
clusterEvalQ(cluster, {
  suppressPackageStartupMessages({library(data.table); library(survival)})
  NULL
})
clusterExport(cluster, c(
  "valid", "prediction_registry", "rep_weights", "grid", "point",
  "cipds_weighted_km_censoring", "cipds_weighted_td_auc", "cipds_weighted_uno",
  "cipds_weighted_brier", "evaluate_fixed_predictions"
), envir = environment())
rep_list <- parLapplyLB(cluster, seq_len(CIPDS_REPLICATES), function(b) {
  tryCatch(evaluate_fixed_predictions(rep_weights[, b]), error = function(e) {
    z <- rep(NA_real_, length(point)); names(z) <- names(point); z
  })
})
stopCluster(cluster)
cluster <- NULL
rep_matrix <- do.call(rbind, rep_list)
colnames(rep_matrix) <- names(point)
finite_fraction <- colMeans(is.finite(rep_matrix))
if (any(finite_fraction < 0.98)) stop("Fair-experiment bootstrap failure rate exceeded 2%")
for (j in seq_along(point)) rep_matrix[!is.finite(rep_matrix[, j]), j] <- point[j]

variance <- svrVar(rep_matrix, scale = rep_valid$scale, rscales = rep_valid$rscales,
                   mse = rep_valid$mse, coef = point)
se <- sqrt(diag(variance))
df <- degf(x$design_full[valid_mask, ])
crit <- qt(0.975, df)
parse_key <- function(z) strsplit(z, "__", fixed = TRUE)[[1]]
estimates <- rbindlist(lapply(seq_along(point), function(j) {
  parts <- parse_key(names(point)[j])
  data.table(
    metric = parts[1], model = parts[2], model_label = model_labels[[parts[2]]],
    temporal_training_cycles = "2003-2006", temporal_validation_cycles = "2007-2010",
    n_validation = nrow(valid), validation_deaths = sum(valid$Death_AllCause),
    estimate = point[j], standard_error = se[j],
    ci_lower = max(0, point[j] - crit * se[j]), ci_upper = min(1, point[j] + crit * se[j]),
    survey_bootstrap_replicates = CIPDS_REPLICATES,
    finite_replicate_fraction = finite_fraction[j]
  )
}))

pairs <- data.table(
  model_a = c("LAB26", "COMPONENT3", "PHENO9", "LAB26", "COMPONENT3", "LAB26"),
  model_b = c("CLINICAL", "CLINICAL", "CLINICAL", "PHENO9", "PHENO9", "COMPONENT3"),
  comparison = c(
    "26 labs beyond clinical base", "Three components beyond clinical base",
    "Nine PhenoAge labs beyond clinical base", "26 labs versus nine PhenoAge labs",
    "Three components versus nine PhenoAge labs", "26 labs versus three components"
  )
)
metrics <- unique(vapply(names(point), function(z) parse_key(z)[1], character(1)))
comparison_rows <- list()
for (metric in metrics) {
  higher_better <- grepl("^(uno|auc)", metric)
  for (i in seq_len(nrow(pairs))) {
    ka <- paste0(metric, "__", pairs$model_a[i])
    kb <- paste0(metric, "__", pairs$model_b[i])
    delta <- point[[ka]] - point[[kb]]
    rep_delta <- rep_matrix[, ka] - rep_matrix[, kb]
    v <- as.numeric(svrVar(matrix(rep_delta, ncol = 1), scale = rep_valid$scale,
                           rscales = rep_valid$rscales, mse = rep_valid$mse, coef = delta))
    delta_se <- sqrt(v)
    comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
      metric = metric, higher_is_better = higher_better,
      model_a = pairs$model_a[i], model_b = pairs$model_b[i], comparison = pairs$comparison[i],
      estimate_a = point[[ka]], estimate_b = point[[kb]], paired_difference_a_minus_b = delta,
      benefit_oriented_difference = ifelse(higher_better, delta, -delta),
      difference_ci_lower = delta - crit * delta_se,
      difference_ci_upper = delta + crit * delta_se,
      paired_p_value = 2 * pt(-abs(delta / delta_se), df = df)
    )
  }
}
comparisons <- rbindlist(comparison_rows)
comparisons[, paired_p_holm := p.adjust(paired_p_value, method = "holm"), by = metric]

fwrite(estimates, file.path(CIPDS_OUT, "fair_target_temporal_validation_estimates.csv"), bom = TRUE)
fwrite(comparisons, file.path(CIPDS_OUT, "fair_target_paired_comparisons.csv"), bom = TRUE)
saveRDS(fit_registry, file.path(CIPDS_OUT, "fair_target_fitted_models.rds"), compress = "gzip")
saveRDS(list(point = point, replicates = rep_matrix, df = df),
        file.path(CIPDS_OUT, "fair_target_validation_replicates.rds"), compress = "gzip")

cipds_write_manifest(file.path(CIPDS_OUT, "fair_target_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  estimand = "fair comparison after matching mortality target, algorithm, training period, validation period and clinical base",
  training = "NHANES 2003-2006", validation = "NHANES 2007-2010",
  algorithm = "elastic-net Cox regression, alpha 0.5, fivefold PSU-grouped cross-validation, lambda.min selected for discrimination; lambda.1se value retained in the coefficient registry as a parsimony reference",
  preprocessing = "training-only 1st/99th percentile winsorization and standardization",
  external_validation_uncertainty = "1000 survey bootstrap replicate weights in the temporal validation sample; training models held fixed",
  primary_metric = "Uno C-index truncated at 10 years",
  secondary_metrics = c("AUC at 5 and 10 years", "Brier at 5 and 10 years", "IBS from 1 to 10 years"),
  not_estimable = "15-year temporal-validation performance because 2007-2010 maximum follow-up is 13.17 years",
  threads = CIPDS_THREADS, seed = CIPDS_SEED + 32L
))
cat("\nCohort audit:\n")
print(cohort_audit)
cat("\nValidation estimates:\n")
print(estimates)
cat("\nPaired comparisons:\n")
print(comparisons)
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
