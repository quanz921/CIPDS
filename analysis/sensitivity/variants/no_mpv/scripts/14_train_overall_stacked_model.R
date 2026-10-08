options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(catboost)
  library(data.table)
  library(pROC)
  library(MASS)
  library(jsonlite)
})

package_dir <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
out_dir <- file.path(package_dir, "outputs")
log_dir <- file.path(package_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "14_train_overall_stacked_model.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

THREADS <- 23L
OUTER_K <- 5L
BASE_SEED <- 20260831L
EPS <- 1e-6

OUTCOMES <- c(
  "Outcome_NutriMetab",
  "Outcome_TumorBurden",
  "Outcome_TreatComp"
)
SHORT <- c(
  Outcome_NutriMetab = "NM",
  Outcome_TumorBurden = "TB",
  Outcome_TreatComp = "TC"
)

cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Purpose: cross-fitted stacked multidomain laboratory vulnerability model\n")
cat("Base models: frozen A-primary outcome-specific CatBoost models\n")
cat("PhenoAge and NHANES outcomes are not model inputs\n")
cat("CatBoost threads:", THREADS, "\n")

data_file <- file.path(out_dir, "hospital_patient_model_ready_129.csv")
features_file <- file.path(out_dir, "nested_a_primary_final_features.csv")
boruta_file <- file.path(out_dir, "nested_a_primary_boruta_importance.csv")
outer_file <- file.path(out_dir, "nested_a_primary_outer_performance.csv")
registry_file <- file.path(out_dir, "nested_a_primary_calibration_and_model_registry.csv")
required_files <- c(data_file, features_file, boruta_file, outer_file, registry_file)
if (any(!file.exists(required_files))) {
  stop("Missing prerequisite A-primary nested outputs: ",
       paste(basename(required_files[!file.exists(required_files)]), collapse = ", "))
}

dat <- fread(data_file, na.strings = c("", "NA"), showProgress = TRUE)
features_dt <- fread(features_file)
boruta_dt <- fread(boruta_file)
outer_dt <- fread(outer_file)
registry_dt <- fread(registry_file)

if (nrow(dat) != 75248L || uniqueN(dat$patient_key) != 75248L) {
  stop("Hospital patient grain drift")
}
if (!all(OUTCOMES %in% names(dat))) stop("Missing component outcomes")
if (any(registry_dt$scenario != "A_PRIMARY")) stop("Only A-primary registry is allowed")

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold_id <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold_id
}

make_pool <- function(x, y = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(y)) {
    catboost.load_pool(data = x, feature_names = as.list(colnames(x)))
  } else {
    catboost.load_pool(
      data = x,
      label = as.numeric(y),
      feature_names = as.list(colnames(x))
    )
  }
}

class_weights <- function(y) c(1, sum(y == 0) / sum(y == 1))

fit_base <- function(x, y, depth, learning_rate, l2_leaf_reg, seed) {
  pool <- make_pool(x, y)
  model <- catboost.train(pool, params = list(
    loss_function = "Logloss",
    eval_metric = "AUC",
    iterations = 500L,
    random_seed = as.integer(seed),
    thread_count = THREADS,
    class_weights = class_weights(y),
    depth = as.integer(depth),
    learning_rate = as.numeric(learning_rate),
    l2_leaf_reg = as.numeric(l2_leaf_reg),
    logging_level = "Silent",
    allow_writing_files = FALSE
  ))
  model
}

predict_base <- function(model, x) {
  as.numeric(catboost.predict(model, make_pool(x), prediction_type = "Probability"))
}

apply_platt <- function(raw_probability, intercept, slope) {
  lp <- qlogis(pmin(pmax(raw_probability, EPS), 1 - EPS))
  plogis(as.numeric(intercept) + as.numeric(slope) * lp)
}

auc_value <- function(y, p) {
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))
}

average_precision <- function(y, p) {
  ok <- !is.na(y) & is.finite(p)
  y <- as.integer(y[ok]); p <- p[ok]
  if (!sum(y == 1L)) return(NA_real_)
  ord <- order(p, decreasing = TRUE)
  yy <- as.integer(y[ord] == 1L); pp <- p[ord]
  last <- which(c(diff(pp) != 0, TRUE))
  tp <- cumsum(yy)[last]
  sum(diff(c(0, tp / sum(yy))) * tp / last)
}

binary_auc_ci <- function(y, p) {
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  c(auc = as.numeric(pROC::auc(roc_obj)), lower = ci[1], upper = ci[3])
}

calibration_terms <- function(y, p) {
  lp <- qlogis(pmin(pmax(p, EPS), 1 - EPS))
  fit <- tryCatch(glm(y ~ lp, family = binomial()), error = function(e) NULL)
  if (is.null(fit)) return(c(intercept = NA_real_, slope = NA_real_))
  c(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]))
}

# Refit the already selected/tuned outer-fold models solely to recover strict
# patient-level OOF probabilities. No feature selection or tuning is repeated.
all_oof <- list()
oof_checks <- list()

for (outcome_idx in seq_along(OUTCOMES)) {
  outcome <- OUTCOMES[outcome_idx]
  outcome_name <- outcome
  short <- SHORT[[outcome]]
  d <- dat[split_calendar_entry == "train" & !is.na(get(outcome))]
  setorder(d, patient_key)
  fold_id <- make_stratified_folds(d[[outcome]], OUTER_K, BASE_SEED + outcome_idx)
  oof <- rep(NA_real_, nrow(d))

  cat("\nRecovering OOF predictions for", outcome, "N=", nrow(d), "\n")
  for (outer_fold in seq_len(OUTER_K)) {
    fold_number <- outer_fold
    tr <- which(fold_id != outer_fold)
    ho <- which(fold_id == outer_fold)
    selected <- boruta_dt[
      outcome == outcome_name & outer_fold == fold_number & selected == TRUE,
      feature
    ]
    if (length(selected) < 2L) stop("Stored outer features missing for ", outcome)
    pars <- outer_dt[outcome == outcome_name & outer_fold == fold_number]
    if (nrow(pars) != 1L) stop("Stored outer parameters missing for ", outcome)
    fold_seed <- BASE_SEED + outcome_idx * 10000L + outer_fold * 1000L
    model <- fit_base(
      d[tr, ..selected], d[[outcome]][tr],
      pars$depth, pars$learning_rate, pars$l2_leaf_reg,
      fold_seed + 900L
    )
    oof[ho] <- predict_base(model, d[ho, ..selected])
    reproduced_auc <- auc_value(d[[outcome]][ho], oof[ho])
    stored_auc <- as.numeric(pars$outer_auc)
    oof_checks[[length(oof_checks) + 1L]] <- data.table(
      outcome = outcome,
      outer_fold = outer_fold,
      n = length(ho),
      positives = sum(d[[outcome]][ho] == 1),
      stored_outer_auc = stored_auc,
      reproduced_outer_auc = reproduced_auc,
      absolute_auc_difference = abs(reproduced_auc - stored_auc),
      reproduced_within_tolerance = abs(reproduced_auc - stored_auc) < 1e-10
    )
    cat("  fold", outer_fold, "stored/reproduced AUC:",
        round(stored_auc, 6), "/", round(reproduced_auc, 6), "\n")
  }
  if (anyNA(oof)) stop("OOF coverage failure for ", outcome)
  reg <- registry_dt[outcome == outcome_name]
  if (nrow(reg) != 1L) stop("Calibration registry mismatch for ", outcome)
  calibrated <- apply_platt(oof, reg$intercept, reg$slope)
  all_oof[[outcome]] <- data.table(
    patient_key = d$patient_key,
    observed = d[[outcome]],
    raw_probability = oof,
    calibrated_probability = calibrated
  )
  setnames(
    all_oof[[outcome]],
    c("observed", "raw_probability", "calibrated_probability"),
    paste0(c("observed_", "raw_", "p_"), short)
  )
}

oof_checks_dt <- rbindlist(oof_checks)
if (!all(oof_checks_dt$reproduced_within_tolerance)) {
  stop("At least one stored outer model could not be reproduced exactly")
}

meta_train <- Reduce(
  function(x, y) merge(x, y, by = "patient_key", all = FALSE),
  all_oof
)
train_labels <- dat[
  split_calendar_entry == "train",
  c("patient_key", OUTCOMES),
  with = FALSE
]
meta_train <- merge(meta_train, train_labels, by = "patient_key", all.x = TRUE)
if (anyNA(meta_train[, ..OUTCOMES])) stop("Unexpected missing labels after OOF join")
meta_train[, burden := as.integer(Outcome_NutriMetab + Outcome_TumorBurden + Outcome_TreatComp)]
if (!all(meta_train$burden %in% 0:3)) stop("Invalid multidomain burden")

for (short in c("NM", "TB", "TC")) {
  meta_train[[paste0("logit_", short)]] <- qlogis(
    pmin(pmax(meta_train[[paste0("p_", short)]], EPS), 1 - EPS)
  )
}

scaling_rows <- lapply(c("NM", "TB", "TC"), function(short) {
  x <- meta_train[[paste0("logit_", short)]]
  data.table(component = short, mean = mean(x), sd = sd(x))
})
scaling_dt <- rbindlist(scaling_rows)
if (any(!is.finite(scaling_dt$sd) | scaling_dt$sd <= 0)) stop("Invalid OOF logit scaling")
for (short in c("NM", "TB", "TC")) {
  row <- scaling_dt[component == short]
  meta_train[[paste0("z_", short)]] <-
    (meta_train[[paste0("logit_", short)]] - row$mean) / row$sd
}
meta_train[, burden_ordered := ordered(burden, levels = 0:3)]

main_formula <- burden_ordered ~ z_NM + z_TB + z_TC
interaction_formula <- burden_ordered ~ z_NM + z_TB + z_TC +
  z_NM:z_TB + z_NM:z_TC + z_TB:z_TC

main_model <- polr(main_formula, data = meta_train, method = "logistic", Hess = TRUE)
interaction_model <- tryCatch(
  polr(interaction_formula, data = meta_train, method = "logistic", Hess = TRUE),
  error = function(e) NULL
)
if (main_model$convergence != 0L) stop("Main-effects ordinal stack did not converge")
if (!is.null(interaction_model) && interaction_model$convergence != 0L) interaction_model <- NULL

score_final_components <- function(split_name) {
  d <- dat[
    split_calendar_entry == split_name &
      !is.na(Outcome_NutriMetab) &
      !is.na(Outcome_TumorBurden) &
      !is.na(Outcome_TreatComp)
  ]
  setorder(d, patient_key)
  result <- d[, .(
    patient_key,
    first_observed_test_date,
    Outcome_NutriMetab,
    Outcome_TumorBurden,
    Outcome_TreatComp
  )]
  for (outcome in OUTCOMES) {
    outcome_name <- outcome
    short <- SHORT[[outcome]]
    selected <- features_dt[outcome == outcome_name, feature]
    reg <- registry_dt[outcome == outcome_name]
    if (!length(selected) || nrow(reg) != 1L) stop("Final model registry error")
    model <- catboost.load_model(file.path(out_dir, reg$model_file))
    raw <- predict_base(model, d[, ..selected])
    result[[paste0("raw_", short)]] <- raw
    result[[paste0("p_", short)]] <- apply_platt(raw, reg$intercept, reg$slope)
  }
  result[, burden := as.integer(Outcome_NutriMetab + Outcome_TumorBurden + Outcome_TreatComp)]
  result
}

validation_scores <- score_final_components("validation")
test_scores <- score_final_components("test")

# Confirm that recomputed test component probabilities are byte-level consistent
# up to floating-point tolerance with the previously frozen outputs.
stored_test <- fread(file.path(out_dir, "nested_a_primary_calendar_test_predictions.csv"))
test_component_checks <- list()
for (outcome in OUTCOMES) {
  outcome_name <- outcome
  short <- SHORT[[outcome]]
  stored <- stored_test[outcome == outcome_name, .(patient_key, calibrated_probability)]
  current <- test_scores[, c("patient_key", paste0("p_", short)), with = FALSE]
  chk <- merge(stored, current, by = "patient_key", all = FALSE)
  diff <- abs(chk$calibrated_probability - chk[[paste0("p_", short)]])
  test_component_checks[[length(test_component_checks) + 1L]] <- data.table(
    outcome = outcome,
    n_compared = nrow(chk),
    max_absolute_probability_difference = max(diff),
    within_tolerance = max(diff) < 1e-12
  )
}
test_component_checks_dt <- rbindlist(test_component_checks)
if (!all(test_component_checks_dt$within_tolerance)) {
  stop("Final test component predictions do not reproduce frozen outputs")
}

prepare_meta_features <- function(d) {
  ans <- copy(d)
  for (short in c("NM", "TB", "TC")) {
    logit_name <- paste0("logit_", short)
    z_name <- paste0("z_", short)
    ans[[logit_name]] <- qlogis(pmin(pmax(ans[[paste0("p_", short)]], EPS), 1 - EPS))
    row <- scaling_dt[component == short]
    ans[[z_name]] <- (ans[[logit_name]] - row$mean) / row$sd
  }
  ans
}

validation_meta <- prepare_meta_features(validation_scores)
test_meta <- prepare_meta_features(test_scores)

predict_ordinal <- function(model, d) {
  prob <- as.matrix(predict(model, newdata = d, type = "probs"))
  expected_cols <- as.character(0:3)
  if (!all(expected_cols %in% colnames(prob))) stop("Ordinal probability columns missing")
  prob <- prob[, expected_cols, drop = FALSE]
  expected <- as.numeric(prob %*% 0:3)
  list(
    prob = prob,
    expected = expected,
    probability_any = 1 - prob[, "0"],
    probability_multidomain = prob[, "2"] + prob[, "3"]
  )
}

ranked_probability_score <- function(y, prob) {
  obs <- model.matrix(~ factor(y, levels = 0:3) - 1)
  cum_prob <- t(apply(prob, 1, cumsum))[, 1:3, drop = FALSE]
  cum_obs <- t(apply(obs, 1, cumsum))[, 1:3, drop = FALSE]
  mean(rowSums((cum_prob - cum_obs)^2) / 3)
}

multiclass_logloss <- function(y, prob) {
  idx <- cbind(seq_along(y), y + 1L)
  -mean(log(pmax(prob[idx], 1e-15)))
}

ordinal_performance <- function(model_name, model, d, split_name) {
  pred <- predict_ordinal(model, d)
  y <- d$burden
  any_y <- as.integer(y >= 1L)
  multi_y <- as.integer(y >= 2L)
  any_auc <- binary_auc_ci(any_y, pred$probability_any)
  multi_auc <- binary_auc_ci(multi_y, pred$probability_multidomain)
  data.table(
    candidate_model = model_name,
    split = split_name,
    n = length(y),
    burden_0 = sum(y == 0L),
    burden_1 = sum(y == 1L),
    burden_2 = sum(y == 2L),
    burden_3 = sum(y == 3L),
    ordinal_rps = ranked_probability_score(y, pred$prob),
    multiclass_logloss = multiclass_logloss(y, pred$prob),
    expected_burden_mae = mean(abs(pred$expected - y)),
    expected_burden_rmse = sqrt(mean((pred$expected - y)^2)),
    expected_burden_spearman = suppressWarnings(cor(pred$expected, y, method = "spearman")),
    any_auc = any_auc["auc"],
    any_auc_ci_lower = any_auc["lower"],
    any_auc_ci_upper = any_auc["upper"],
    any_average_precision = average_precision(any_y, pred$probability_any),
    multidomain_auc = multi_auc["auc"],
    multidomain_auc_ci_lower = multi_auc["lower"],
    multidomain_auc_ci_upper = multi_auc["upper"],
    multidomain_average_precision = average_precision(multi_y, pred$probability_multidomain)
  )
}

candidates <- list(main_effects = main_model)
if (!is.null(interaction_model)) candidates$pairwise_interactions <- interaction_model

candidate_validation <- rbindlist(lapply(names(candidates), function(nm) {
  ordinal_performance(nm, candidates[[nm]], validation_meta, "calendar_validation")
}))

chosen_name <- "main_effects"
selection_reason <- "Pre-specified parsimonious main-effects ordinal stack"
if ("pairwise_interactions" %in% candidate_validation$candidate_model) {
  main_rps <- candidate_validation[candidate_model == "main_effects", ordinal_rps]
  int_rps <- candidate_validation[candidate_model == "pairwise_interactions", ordinal_rps]
  main_multi_auc <- candidate_validation[candidate_model == "main_effects", multidomain_auc]
  int_multi_auc <- candidate_validation[candidate_model == "pairwise_interactions", multidomain_auc]
  if (int_rps <= main_rps * 0.995 && int_multi_auc >= main_multi_auc - 0.002) {
    chosen_name <- "pairwise_interactions"
    selection_reason <- paste(
      "Interaction stack selected prospectively because validation ranked-probability",
      "score improved by at least 0.5% without >0.002 multidomain-AUC loss"
    )
  }
}
chosen_model <- candidates[[chosen_name]]

chosen_validation <- ordinal_performance(
  chosen_name, chosen_model, validation_meta, "calendar_validation"
)
chosen_test <- ordinal_performance(
  chosen_name, chosen_model, test_meta, "calendar_test"
)
performance_dt <- rbindlist(list(candidate_validation, chosen_test), fill = TRUE)

make_prediction_table <- function(d, split_name) {
  pred <- predict_ordinal(chosen_model, d)
  data.table(
    split = split_name,
    patient_key = d$patient_key,
    first_observed_test_date = d$first_observed_test_date,
    observed_NM = d$Outcome_NutriMetab,
    observed_TB = d$Outcome_TumorBurden,
    observed_TC = d$Outcome_TreatComp,
    observed_burden = d$burden,
    p_NM = d$p_NM,
    p_TB = d$p_TB,
    p_TC = d$p_TC,
    probability_burden_0 = pred$prob[, "0"],
    probability_burden_1 = pred$prob[, "1"],
    probability_burden_2 = pred$prob[, "2"],
    probability_burden_3 = pred$prob[, "3"],
    overall_expected_burden = pred$expected,
    probability_any_domain = pred$probability_any,
    probability_multidomain = pred$probability_multidomain,
    equal_sum_component_probability = d$p_NM + d$p_TB + d$p_TC
  )
}

validation_predictions <- make_prediction_table(validation_meta, "calendar_validation")
test_predictions <- make_prediction_table(test_meta, "calendar_test")
predictions_dt <- rbindlist(list(validation_predictions, test_predictions))

# Paired comparison against the transparent equal-sum benchmark. This benchmark
# is not promoted as the final score; it tests whether learned stacking adds value.
paired_rows <- list()
for (target in c("any", "multidomain")) {
  y <- if (target == "any") {
    as.integer(test_predictions$observed_burden >= 1L)
  } else {
    as.integer(test_predictions$observed_burden >= 2L)
  }
  learned <- if (target == "any") {
    test_predictions$probability_any_domain
  } else {
    test_predictions$probability_multidomain
  }
  baseline <- test_predictions$equal_sum_component_probability
  learned_roc <- pROC::roc(y, learned, quiet = TRUE, direction = "<")
  baseline_roc <- pROC::roc(y, baseline, quiet = TRUE, direction = "<")
  test_result <- pROC::roc.test(learned_roc, baseline_roc, paired = TRUE, method = "delong")
  paired_rows[[length(paired_rows) + 1L]] <- data.table(
    target = target,
    n = length(y),
    positives = sum(y == 1L),
    learned_stack_auc = as.numeric(pROC::auc(learned_roc)),
    equal_sum_auc = as.numeric(pROC::auc(baseline_roc)),
    auc_difference_learned_minus_equal_sum =
      as.numeric(pROC::auc(learned_roc) - pROC::auc(baseline_roc)),
    paired_delong_p = as.numeric(test_result$p.value),
    learned_stack_average_precision = average_precision(y, learned),
    equal_sum_average_precision = average_precision(y, baseline)
  )
}
paired_dt <- rbindlist(paired_rows)

calibration_rows <- list()
for (split_name in c("calendar_validation", "calendar_test")) {
  d <- predictions_dt[split == split_name]
  for (target in c("any", "multidomain")) {
    y <- if (target == "any") as.integer(d$observed_burden >= 1L) else as.integer(d$observed_burden >= 2L)
    p <- if (target == "any") d$probability_any_domain else d$probability_multidomain
    cal <- calibration_terms(y, p)
    calibration_rows[[length(calibration_rows) + 1L]] <- data.table(
      split = split_name,
      target = target,
      n = length(y),
      positives = sum(y == 1L),
      calibration_intercept = cal["intercept"],
      calibration_slope = cal["slope"],
      brier = mean((p - y)^2)
    )
  }
}
calibration_dt <- rbindlist(calibration_rows)

state_counts <- dat[
  !is.na(Outcome_NutriMetab) & !is.na(Outcome_TumorBurden) & !is.na(Outcome_TreatComp),
  .(
    n = .N,
    first_date = min(first_observed_test_date),
    last_date = max(first_observed_test_date)
  ),
  by = .(
    split_calendar_entry,
    Outcome_NutriMetab,
    Outcome_TumorBurden,
    Outcome_TreatComp
  )
]
state_counts[, burden := Outcome_NutriMetab + Outcome_TumorBurden + Outcome_TreatComp]
state_counts[, state := paste0(Outcome_NutriMetab, Outcome_TumorBurden, Outcome_TreatComp)]
setorder(state_counts, split_calendar_entry, state)

model_bundle <- list(
  model = chosen_model,
  chosen_model_name = chosen_name,
  selection_reason = selection_reason,
  component_order = c("NM", "TB", "TC"),
  component_input = "Platt-calibrated probabilities from frozen A-primary base models",
  logit_scaling = as.data.frame(scaling_dt),
  output_definition = paste(
    "Ordinal probabilities for 0,1,2,3 affected domains; overall score is",
    "the predicted expected number of affected domains"
  ),
  training_population = "Hospital calendar-entry training period, strict base-model OOF predictions",
  phenoage_used = FALSE,
  nhanes_outcomes_used = FALSE,
  analysis_status = "AUDIT_ONLY_LEGACY_UNDATED_OUTCOMES"
)
saveRDS(model_bundle, file.path(out_dir, "overall_stacked_ordinal_model.rds"))

model_coefficients <- data.table(
  parameter = c(names(coef(chosen_model)), names(chosen_model$zeta)),
  parameter_type = c(
    rep("coefficient", length(coef(chosen_model))),
    rep("threshold", length(chosen_model$zeta))
  ),
  estimate = c(unname(coef(chosen_model)), unname(chosen_model$zeta))
)

fwrite(oof_checks_dt, file.path(out_dir, "overall_oof_reproduction_checks.csv"), bom = TRUE)
fwrite(test_component_checks_dt, file.path(out_dir, "overall_frozen_component_reproduction_checks.csv"), bom = TRUE)
fwrite(meta_train[, .(
  patient_key, observed_NM, observed_TB, observed_TC,
  p_NM, p_TB, p_TC, burden
)], file.path(out_dir, "overall_meta_training_oof.csv"), bom = TRUE)
fwrite(scaling_dt, file.path(out_dir, "overall_meta_logit_scaling.csv"), bom = TRUE)
fwrite(state_counts, file.path(out_dir, "overall_hospital_outcome_state_counts.csv"), bom = TRUE)
fwrite(performance_dt, file.path(out_dir, "overall_calendar_performance.csv"), bom = TRUE)
fwrite(predictions_dt, file.path(out_dir, "overall_calendar_predictions.csv"), bom = TRUE)
fwrite(paired_dt, file.path(out_dir, "overall_calendar_paired_comparison.csv"), bom = TRUE)
fwrite(calibration_dt, file.path(out_dir, "overall_calendar_calibration.csv"), bom = TRUE)
fwrite(model_coefficients, file.path(out_dir, "overall_stacked_model_coefficients.csv"), bom = TRUE)

manifest <- list(
  completed_at = as.character(Sys.time()),
  model_name = "Multidomain laboratory vulnerability overall model",
  architecture = "cross-fitted outcome-specific CatBoost base models plus ordinal logistic stack",
  base_models_retrained = FALSE,
  base_model_oof_recovery = paste(
    "Stored fold-specific features, hyperparameters, folds, and seeds were used",
    "to reproduce the original outer-fold predictions"
  ),
  chosen_stack = chosen_name,
  stack_selection_reason = selection_reason,
  training_n = nrow(meta_train),
  training_burden_counts = as.list(table(meta_train$burden)),
  validation_n = nrow(validation_meta),
  test_n = nrow(test_meta),
  test_burden_counts = as.list(table(test_meta$burden)),
  phenoage_used_as_input = FALSE,
  nhanes_used_for_training = FALSE,
  component_models_remain_separately_reported = TRUE,
  no_auc_weighted_composite = TRUE,
  thread_count = THREADS,
  analysis_status = "AUDIT_ONLY_LEGACY_UNDATED_OUTCOMES"
)
write_json(
  manifest,
  file.path(out_dir, "overall_stacked_run_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\nChosen stack:", chosen_name, "\n")
cat("Selection reason:", selection_reason, "\n")
cat("\nCalendar performance:\n")
print(performance_dt)
cat("\nPaired test comparison:\n")
print(paired_dt)
cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
