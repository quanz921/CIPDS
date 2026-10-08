source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))
suppressPackageStartupMessages({
  library(glmnet)
  library(mgcv)
  library(parallel)
})

log_file <- file.path(CIPDS_LOG, "40_residual_fidelity_extension.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Replicates:", CIPDS_REPLICATES, " Workers:", CIPDS_THREADS, "\n")

OUT3 <- file.path(CIPDS_ROOT, "outputs", "supplement_v3")
dir.create(OUT3, recursive = TRUE, showWarnings = FALSE)

x <- cipds_load_data(TRUE)
dat <- x$dat

feature_registry <- fread(file.path(CIPDS_ROOT, "outputs", "nested_a_primary_final_features.csv"))
labs26 <- unique(feature_registry$feature)
labs26 <- labs26[!is.na(labs26) & labs26 != ""]
pheno9 <- c(
  "lab_alb", "lab_crea", "lab_glu", "lab_hcrp", "lab_lym_pct",
  "lab_mcv", "lab_rdw", "lab_alp", "lab_wbc"
)
if (length(labs26) != 26L) stop("Expected 26 laboratory variables")

# Retain the exact cohort used by the prespecified matched-target experiment,
# including PhenoAge-laboratory completeness, so every comparison is paired.
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
expected_cohorts <- fread(file.path(CIPDS_OUT,"fair_target_cohort_audit.csv"))
stopifnot(nrow(train) == expected_cohorts[role=="TRAIN",n], sum(train$Death_AllCause) == expected_cohorts[role=="TRAIN",deaths])
stopifnot(nrow(valid) == expected_cohorts[role=="TEMPORAL_VALIDATION",n], sum(valid$Death_AllCause) == expected_cohorts[role=="TEMPORAL_VALIDATION",deaths])

transform_registry <- fread(file.path(CIPDS_OUT, "fair_target_training_transform_registry.csv"))
needed_transform <- unique(c(labs26, pheno9, "NM_component", "TB_component", "TC_component"))
if (length(setdiff(needed_transform, transform_registry$variable))) {
  stop("Stored fair-target transform registry is incomplete")
}

for (v in unique(c(labs26, pheno9))) {
  row <- transform_registry[variable == v]
  lo <- row$winsor_lower[[1]]
  hi <- row$winsor_upper[[1]]
  mu <- row$training_mean_after_winsorization[[1]]
  sigma <- row$training_sd_after_winsorization[[1]]
  train[[paste0("z_", v)]] <- (pmin(pmax(train[[v]], lo), hi) - mu) / sigma
  valid[[paste0("z_", v)]] <- (pmin(pmax(valid[[v]], lo), hi) - mu) / sigma
}
for (v in c("NM_component", "TB_component", "TC_component")) {
  row <- transform_registry[variable == v]
  mu <- row$training_mean_after_winsorization[[1]]
  sigma <- row$training_sd_after_winsorization[[1]]
  train[[paste0("zfair_", v)]] <- (train[[v]] - mu) / sigma
  valid[[paste0("zfair_", v)]] <- (valid[[v]] - mu) / sigma
}

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
fit_registry <- readRDS(file.path(CIPDS_OUT, "fair_target_fitted_models.rds"))

reconstruct_glmnet <- function(model, d, clinical_x) {
  features <- feature_sets[[model]]
  xmat <- cbind(clinical_x, as.matrix(d[, ..features]))
  cvfit <- fit_registry[[model]]$cvfit
  coefficients <- as.matrix(coef(cvfit, s = "lambda.min"))[, 1]
  if (length(setdiff(names(coefficients), colnames(xmat)))) {
    stop("Model matrix coefficient mismatch for ", model)
  }
  xmat <- xmat[, names(coefficients), drop = FALSE]
  lp_matrix <- as.numeric(xmat %*% coefficients)
  lp_predict <- as.numeric(predict(cvfit, newx = xmat, s = "lambda.min", type = "link"))
  if (max(abs(lp_matrix - lp_predict)) > 1e-10) stop("Linear-predictor reconstruction failure: ", model)
  list(
    x = xmat, coefficients = coefficients, lp = lp_matrix,
    baseline_hazard = fit_registry[[model]]$baseline_hazard
  )
}

lab26_train <- reconstruct_glmnet("LAB26", train, clinical_train)
lab26_valid <- reconstruct_glmnet("LAB26", valid, clinical_valid)
component3_valid <- reconstruct_glmnet("COMPONENT3", valid, clinical_valid)
pheno9_valid <- reconstruct_glmnet("PHENO9", valid, clinical_valid)

clinical_names <- colnames(clinical_train)
lab_names <- feature_sets$LAB26
beta_lab26 <- lab26_train$coefficients
eta_clinical_train <- as.numeric(lab26_train$x[, clinical_names, drop = FALSE] %*%
                                   beta_lab26[clinical_names])
eta_clinical_valid <- as.numeric(lab26_valid$x[, clinical_names, drop = FALSE] %*%
                                   beta_lab26[clinical_names])
eta_lab_train <- as.numeric(lab26_train$x[, lab_names, drop = FALSE] %*% beta_lab26[lab_names])
eta_lab_valid <- as.numeric(lab26_valid$x[, lab_names, drop = FALSE] %*% beta_lab26[lab_names])

component_vars <- feature_sets$COMPONENT3
train[, normalized_weight := pooled_mec_weight / mean(pooled_mec_weight)]
valid[, normalized_weight := pooled_mec_weight / mean(pooled_mec_weight)]
train[, eta_lab26 := eta_lab_train]
valid[, eta_lab26 := eta_lab_valid]

# PSU-grouped fivefold cross-fitting estimates how much of the frozen 26-lab
# mortality linear predictor can be represented by the three frozen domains.
set.seed(CIPDS_SEED + 32L)
psu_table <- unique(train[, .(survey_psu, CYCLE)])
psu_table[, fold := sample(rep(seq_len(5), length.out = .N)), by = CYCLE]
fold_id <- psu_table$fold[match(train$survey_psu, psu_table$survey_psu)]
if (anyNA(fold_id) || length(unique(fold_id)) != 5L) stop("PSU fold assignment failure")

mapping_formula <- as.formula(paste(
  "eta_lab26 ~",
  paste(paste0("s(", component_vars, ", bs='cr', k=4)"), collapse = " + ")
))

fit_mapping_fold <- function(k) {
  fit <- gam(
    mapping_formula, data = train[fold_id != k],
    weights = normalized_weight, method = "REML", gamma = 1.2
  )
  list(index = which(fold_id == k), prediction = as.numeric(predict(fit, newdata = train[fold_id == k])))
}
# Five mapping fits are inexpensive and remain sequential to avoid copying the
# wide NHANES table to workers. The 1,000 validation replicates below use all
# requested workers.
fold_predictions <- lapply(seq_len(5), fit_mapping_fold)
g_oof <- rep(NA_real_, nrow(train))
for (z in fold_predictions) g_oof[z$index] <- z$prediction
if (anyNA(g_oof)) stop("Cross-fitted mapping predictions incomplete")

g_full <- gam(
  mapping_formula, data = train, weights = normalized_weight,
  method = "REML", gamma = 1.2
)
g_train <- as.numeric(predict(g_full, newdata = train))
g_valid <- as.numeric(predict(g_full, newdata = valid))
residual_train_oof <- eta_lab_train - g_oof
residual_valid <- eta_lab_valid - g_valid

weighted_summary <- function(observed, fitted, weight, label, role) {
  mu <- weighted.mean(observed, weight)
  sse <- sum(weight * (observed - fitted)^2)
  sst <- sum(weight * (observed - mu)^2)
  corr <- cov.wt(cbind(observed, fitted), wt = weight, cor = TRUE)$cor[1, 2]
  data.table(
    role = role, representation = label, n = length(observed),
    weighted_r_squared = 1 - sse / sst,
    weighted_rmse = sqrt(sse / sum(weight)),
    weighted_correlation = corr,
    observed_weighted_sd = sqrt(sum(weight * (observed - mu)^2) / sum(weight)),
    residual_weighted_sd = sqrt(sum(weight * ((observed - fitted) -
                                      weighted.mean(observed - fitted, weight))^2) / sum(weight))
  )
}

mapping_quality <- rbind(
  weighted_summary(eta_lab_train, g_oof, train$pooled_mec_weight,
                   "Cross-fitted three-domain GAM", "TRAIN_OOF"),
  weighted_summary(eta_lab_train, g_train, train$pooled_mec_weight,
                   "Full-training three-domain GAM", "TRAIN_APPARENT"),
  weighted_summary(eta_lab_valid, g_valid, valid$pooled_mec_weight,
                   "Frozen three-domain GAM", "TEMPORAL_VALIDATION")
)
fwrite(mapping_quality, file.path(OUT3, "residual_mapping_quality.csv"), bom = TRUE)

# Exact decomposition uses the LAB26 model's own clinical contribution plus
# the projected laboratory contribution and its residual.
lp_projected <- eta_clinical_valid + g_valid
lp_hybrid <- eta_clinical_valid + g_valid + residual_valid
lp_exact <- lab26_valid$lp
grid <- seq_len(10)
h0_lab26 <- lab26_valid$baseline_hazard
make_risk <- function(lp, h0) {
  risk <- sapply(h0, function(h) 1 - exp(-h * exp(lp)))
  colnames(risk) <- paste0("risk", seq_along(h0))
  risk
}

predictions <- list(
  PROJECTED3 = list(lp = lp_projected, risk = make_risk(lp_projected, h0_lab26)),
  HYBRID4 = list(lp = lp_hybrid, risk = make_risk(lp_hybrid, h0_lab26)),
  LAB26 = list(lp = lp_exact, risk = make_risk(lp_exact, h0_lab26)),
  COMPONENT3 = list(
    lp = component3_valid$lp,
    risk = make_risk(component3_valid$lp, component3_valid$baseline_hazard)
  ),
  PHENO9 = list(
    lp = pheno9_valid$lp,
    risk = make_risk(pheno9_valid$lp, pheno9_valid$baseline_hazard)
  )
)

identity_audit <- data.table(
  check = c(
    "LAB26 total LP equals reconstructed clinical plus laboratory contributions",
    "Hybrid4 LP equals LAB26 LP",
    "Hybrid4 5-year risk equals LAB26 risk",
    "Hybrid4 10-year risk equals LAB26 risk"
  ),
  maximum_absolute_error = c(
    max(abs(lp_exact - (eta_clinical_valid + eta_lab_valid))),
    max(abs(lp_hybrid - lp_exact)),
    max(abs(predictions$HYBRID4$risk[, 5] - predictions$LAB26$risk[, 5])),
    max(abs(predictions$HYBRID4$risk[, 10] - predictions$LAB26$risk[, 10]))
  ),
  tolerance = c(1e-10, 1e-10, 1e-8, 1e-8)
)
identity_audit[, passed := maximum_absolute_error <= tolerance]
if (!all(identity_audit$passed)) stop("Exact residual decomposition identity failed")
fwrite(identity_audit, file.path(OUT3, "residual_fidelity_identity_audit.csv"), bom = TRUE)

evaluate_fixed_predictions <- function(weight) {
  ans <- c()
  for (model in names(predictions)) {
    pred <- predictions[[model]]
    z <- c(uno_c10 = cipds_weighted_uno(
      valid$Follow_Up_Years, valid$Death_AllCause, pred$lp, weight, 10
    ))
    for (h in c(5, 10)) {
      z[paste0("auc", h)] <- cipds_weighted_td_auc(
        valid$Follow_Up_Years, valid$Death_AllCause, pred$lp, weight, h
      )
      z[paste0("brier", h)] <- cipds_weighted_brier(
        valid$Follow_Up_Years, valid$Death_AllCause, pred$risk[, h], weight, h
      )
    }
    brier_grid <- vapply(grid, function(h) cipds_weighted_brier(
      valid$Follow_Up_Years, valid$Death_AllCause, pred$risk[, h], weight, h
    ), numeric(1))
    z["ibs10"] <- sum(diff(grid) * (head(brier_grid, -1) + tail(brier_grid, -1)) / 2) /
      (max(grid) - min(grid))
    names(z) <- paste0(names(z), "__", model)
    ans <- c(ans, z)
  }
  ans
}

point <- evaluate_fixed_predictions(valid$pooled_mec_weight)
set.seed(CIPDS_SEED + 140L)
rep_full <- as.svrepdesign(x$design_full, type = "bootstrap", replicates = CIPDS_REPLICATES, mse = TRUE)
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
  "valid", "predictions", "rep_weights", "grid", "point",
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
if (any(finite_fraction < 0.98)) stop("Residual-fidelity bootstrap failure rate exceeded 2%")
for (j in seq_along(point)) rep_matrix[!is.finite(rep_matrix[, j]), j] <- point[j]

variance <- svrVar(
  rep_matrix, scale = rep_valid$scale, rscales = rep_valid$rscales,
  mse = rep_valid$mse, coef = point
)
se <- sqrt(diag(variance))
df <- degf(x$design_full[valid_mask, ])
crit <- qt(0.975, df)
parse_key <- function(z) strsplit(z, "__", fixed = TRUE)[[1]]
model_labels <- c(
  PROJECTED3 = "LAB26 clinical contribution + three-domain projected laboratory contribution",
  HYBRID4 = "NM + TB + TC representation plus residual laboratory mortality dimension",
  LAB26 = "Full 26-laboratory mortality model",
  COMPONENT3 = "Clinical base + NM + TB + TC mortality model",
  PHENO9 = "Clinical base + nine PhenoAge laboratory mortality model"
)

estimates <- rbindlist(lapply(seq_along(point), function(j) {
  parts <- parse_key(names(point)[j])
  data.table(
    metric = parts[1], model = parts[2], model_label = model_labels[[parts[2]]],
    temporal_training_cycles = "2003-2006", temporal_validation_cycles = "2007-2010",
    n_validation = nrow(valid), validation_deaths = sum(valid$Death_AllCause),
    estimate = point[j], standard_error = se[j],
    ci_lower = max(0, point[j] - crit * se[j]),
    ci_upper = min(1, point[j] + crit * se[j]),
    survey_bootstrap_replicates = CIPDS_REPLICATES,
    finite_replicate_fraction = finite_fraction[j]
  )
}))

pairs <- data.table(
  model_a = c("LAB26", "HYBRID4", "LAB26", "LAB26", "LAB26", "PROJECTED3"),
  model_b = c("PROJECTED3", "LAB26", "COMPONENT3", "PHENO9", "HYBRID4", "COMPONENT3"),
  comparison = c(
    "Residual information recovered: LAB26 versus projected three-domain representation",
    "Exact-fidelity audit: Hybrid4 versus LAB26",
    "LAB26 versus independently fitted three-component mortality model",
    "LAB26 versus nine-PhenoAge-laboratory mortality model",
    "Exact-fidelity audit: LAB26 versus Hybrid4",
    "Projected three-domain representation versus independently fitted three-component mortality model"
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
    v <- as.numeric(svrVar(
      matrix(rep_delta, ncol = 1), scale = rep_valid$scale,
      rscales = rep_valid$rscales, mse = rep_valid$mse, coef = delta
    ))
    delta_se <- sqrt(max(v, 0))
    p <- if (delta_se <= .Machine$double.eps) {
      ifelse(abs(delta) <= 1e-12, 1, 0)
    } else {
      2 * pt(-abs(delta / delta_se), df = df)
    }
    comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
      metric = metric, higher_is_better = higher_better,
      model_a = pairs$model_a[i], model_b = pairs$model_b[i], comparison = pairs$comparison[i],
      estimate_a = point[[ka]], estimate_b = point[[kb]], paired_difference_a_minus_b = delta,
      benefit_oriented_difference = ifelse(higher_better, delta, -delta),
      difference_standard_error = delta_se,
      difference_ci_lower = delta - crit * delta_se,
      difference_ci_upper = delta + crit * delta_se,
      paired_p_value = p
    )
  }
}
comparisons <- rbindlist(comparison_rows)
comparisons[, paired_p_holm := p.adjust(paired_p_value, method = "holm"), by = metric]

fwrite(estimates, file.path(OUT3, "residual_fidelity_validation_estimates.csv"), bom = TRUE)
fwrite(comparisons, file.path(OUT3, "residual_fidelity_paired_comparisons.csv"), bom = TRUE)
saveRDS(
  list(point = point, replicates = rep_matrix, df = df),
  file.path(OUT3, "residual_fidelity_validation_replicates.rds"), compress = "gzip"
)
saveRDS(
  list(mapping_model = g_full, formula = mapping_formula, component_variables = component_vars),
  file.path(OUT3, "residual_fidelity_mapping_model.rds"), compress = "gzip"
)

weighted_cor <- function(a, b, w) {
  keep <- is.finite(a) & is.finite(b) & is.finite(w) & w > 0
  cov.wt(cbind(a[keep], b[keep]), wt = w[keep], cor = TRUE)$cor[1, 2]
}
residual_correlations <- rbindlist(lapply(labs26, function(v) {
  z <- paste0("z_", v)
  data.table(
    laboratory = v,
    glmnet_coefficient_in_lab26_model = beta_lab26[[z]],
    weighted_correlation_with_residual = weighted_cor(valid[[z]], residual_valid, valid$pooled_mec_weight)
  )
}))
residual_correlations[, absolute_weighted_correlation := abs(weighted_correlation_with_residual)]
setorder(residual_correlations, -absolute_weighted_correlation)
residual_correlations[, descriptive_rank := seq_len(.N)]
fwrite(residual_correlations, file.path(OUT3, "residual_dimension_lab_correlations.csv"), bom = TRUE)

cipds_write_manifest(file.path(OUT3, "residual_fidelity_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  estimand = "target-specific fidelity decomposition of the laboratory contribution to the frozen matched-target LAB26 mortality model",
  training = "NHANES 2003-2006", validation = "NHANES 2007-2010",
  cohort = list(training_n = nrow(train), training_deaths = sum(train$Death_AllCause),
                validation_n = nrow(valid), validation_deaths = sum(valid$Death_AllCause)),
  lab26_model = "prespecified elastic-net Cox LAB26 model from the matched-target experiment",
  mapping = "weighted GAM of the LAB26 laboratory-only linear predictor on NM, TB, and TC; fivefold PSU-grouped cross-fitting in training; full-training mapping frozen before validation",
  residual = "LAB26 laboratory-only linear predictor minus the frozen three-domain GAM projection",
  exact_identity = "LAB26 total linear predictor = LAB26 clinical contribution + three-domain projection + residual",
  uncertainty = "1000 NHANES survey bootstrap replicate weights; all models and mappings held fixed",
  nonclaim = c(
    "This is not lossless compression of the 26 raw laboratory values.",
    "The residual dimension is not a fourth biological domain or aging archetype.",
    "This NHANES mortality-targeted analysis does not modify the hospital-trained primary scores."
  ),
  threads = CIPDS_THREADS,
  seed = CIPDS_SEED + 140L
))

cat("\nMapping quality:\n")
print(mapping_quality)
cat("\nIdentity audit:\n")
print(identity_audit)
cat("\nValidation estimates:\n")
print(estimates)
cat("\nKey paired comparisons:\n")
print(comparisons[grepl("Residual information recovered|Exact-fidelity", comparison)])
cat("\nTop residual laboratory correlates:\n")
print(head(residual_correlations, 10))
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
