suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

package_dir <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
out_dir <- file.path(package_dir, "outputs")
outcomes <- c("Outcome_NutriMetab", "Outcome_TumorBurden", "Outcome_TreatComp")

pred_a <- fread(file.path(out_dir, "nested_a_primary_calendar_test_predictions.csv"))
pred_ab <- fread(file.path(out_dir, "nested_ab_sensitivity_calendar_test_predictions.csv"))
outer_a <- fread(file.path(out_dir, "nested_a_primary_outer_performance.csv"))
outer_ab <- fread(file.path(out_dir, "nested_ab_sensitivity_outer_performance.csv"))
features_ab <- fread(file.path(out_dir, "nested_ab_sensitivity_final_features.csv"))

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

comparison_rows <- list()
for (outcome in outcomes) {
  outcome_name <- outcome
  a <- pred_a[outcome == outcome_name, .(
    patient_key, observed_a = observed, probability_a = calibrated_probability
  )]
  ab <- pred_ab[outcome == outcome_name, .(
    patient_key, observed_ab = observed, probability_ab = calibrated_probability
  )]
  d <- merge(a, ab, by = "patient_key", all = FALSE)
  if (nrow(d) != 11125L || any(d$observed_a != d$observed_ab)) {
    stop("Paired test set mismatch for ", outcome)
  }
  roc_a <- pROC::roc(d$observed_a, d$probability_a, quiet = TRUE, direction = "<")
  roc_ab <- pROC::roc(d$observed_a, d$probability_ab, quiet = TRUE, direction = "<")
  paired <- pROC::roc.test(
    roc_a, roc_ab, paired = TRUE, method = "delong", conf.level = 0.95
  )
  ci <- paired$conf.int
  comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
    outcome = outcome,
    n = nrow(d),
    positives = sum(d$observed_a == 1),
    auc_a_primary = as.numeric(pROC::auc(roc_a)),
    auc_ab_sensitivity = as.numeric(pROC::auc(roc_ab)),
    auc_difference_ab_minus_a = as.numeric(pROC::auc(roc_ab) - pROC::auc(roc_a)),
    # pROC reports the interval in roc1-roc2 orientation (A-AB). Convert it
    # explicitly to the reported AB-A orientation.
    auc_difference_ci_lower = if (length(ci) == 2L) -as.numeric(ci[2]) else NA_real_,
    auc_difference_ci_upper = if (length(ci) == 2L) -as.numeric(ci[1]) else NA_real_,
    auc_difference_ci_orientation = "AB_SENSITIVITY minus A_PRIMARY",
    paired_delong_p = as.numeric(paired$p.value),
    average_precision_a = average_precision(d$observed_a, d$probability_a),
    average_precision_ab = average_precision(d$observed_a, d$probability_ab),
    average_precision_difference = average_precision(d$observed_a, d$probability_ab) -
      average_precision(d$observed_a, d$probability_a),
    brier_a = mean((d$probability_a - d$observed_a)^2),
    brier_ab = mean((d$probability_ab - d$observed_a)^2),
    brier_difference_ab_minus_a = mean((d$probability_ab - d$observed_a)^2) -
      mean((d$probability_a - d$observed_a)^2)
  )
}
comparison <- rbindlist(comparison_rows)
fwrite(comparison, file.path(out_dir, "nested_a_vs_ab_paired_test_comparison.csv"), bom = TRUE)

outer_pair <- merge(
  outer_a[, .(outcome, outer_fold, outer_auc_a = outer_auc)],
  outer_ab[, .(outcome, outer_fold, outer_auc_ab = outer_auc)],
  by = c("outcome", "outer_fold"), all = FALSE
)
outer_pair[, outer_auc_difference_ab_minus_a := outer_auc_ab - outer_auc_a]
outer_summary <- outer_pair[, .(
  mean_outer_auc_a = mean(outer_auc_a),
  mean_outer_auc_ab = mean(outer_auc_ab),
  mean_paired_difference = mean(outer_auc_difference_ab_minus_a),
  sd_paired_difference = sd(outer_auc_difference_ab_minus_a),
  min_paired_difference = min(outer_auc_difference_ab_minus_a),
  max_paired_difference = max(outer_auc_difference_ab_minus_a)
), by = outcome]
fwrite(outer_pair, file.path(out_dir, "nested_a_vs_ab_outer_fold_pairs.csv"), bom = TRUE)
fwrite(outer_summary, file.path(out_dir, "nested_a_vs_ab_outer_fold_summary.csv"), bom = TRUE)

b_features <- c("lab_k", "lab_na", "lab_cl", "lab_glu")
b_selection <- features_ab[feature %in% b_features, .(outcome, feature, votes, status)]
fwrite(b_selection, file.path(out_dir, "nested_ab_selected_b_tier_features.csv"), bom = TRUE)

# Diagnostic: can the mere availability pattern of B-tier tests classify each
# legacy outcome? This is not a competing model; it audits workflow signal.
dat <- fread(file.path(out_dir, "hospital_patient_model_ready_129.csv"),
             select = c("patient_key", "split_calendar_entry", outcomes, b_features))
missingness_rows <- list()
for (outcome in outcomes) {
  train <- dat[split_calendar_entry == "train" & !is.na(get(outcome))]
  validation <- dat[split_calendar_entry == "validation" & !is.na(get(outcome))]
  test <- dat[split_calendar_entry == "test" & !is.na(get(outcome))]
  for (f in b_features) {
    train[[paste0("miss_", f)]] <- as.integer(is.na(train[[f]]))
    validation[[paste0("miss_", f)]] <- as.integer(is.na(validation[[f]]))
    test[[paste0("miss_", f)]] <- as.integer(is.na(test[[f]]))
  }
  miss_terms <- paste0("miss_", b_features)
  formula <- as.formula(paste(outcome, "~", paste(miss_terms, collapse = " + ")))
  fit <- suppressWarnings(glm(formula, data = train, family = binomial()))
  for (split_name in c("validation", "test")) {
    target <- if (split_name == "validation") validation else test
    p <- as.numeric(predict(fit, newdata = target, type = "response"))
    roc_obj <- pROC::roc(target[[outcome]], p, quiet = TRUE, direction = "<")
    missingness_rows[[length(missingness_rows) + 1L]] <- data.table(
      outcome = outcome,
      split = paste0("calendar_", split_name),
      n = nrow(target),
      positives = sum(target[[outcome]] == 1),
      missingness_only_auc = as.numeric(pROC::auc(roc_obj)),
      missingness_only_average_precision = average_precision(target[[outcome]], p)
    )
  }
}
missingness_diagnostic <- rbindlist(missingness_rows)
fwrite(missingness_diagnostic,
       file.path(out_dir, "nested_b_tier_missingness_only_diagnostic.csv"), bom = TRUE)

cat("Paired test comparison:\n")
print(comparison)
cat("\nB-tier selected features:\n")
print(b_selection)
cat("\nMissingness-only diagnostic:\n")
print(missingness_diagnostic)
