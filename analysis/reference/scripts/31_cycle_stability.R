source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))

log_file <- file.path(CIPDS_LOG, "31_cycle_stability.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

x <- cipds_load_data(FALSE)
dat <- x$dat

fit_effect <- function(data_obj, score_var, score, cycle, standardization) {
  tryCatch({
    data_obj <- copy(data_obj)
    data_obj[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]
    fit <- coxph(
      cipds_make_formula(score_var), data = data_obj,
      weights = normalized_sampling_weight, cluster = survey_psu, robust = TRUE
    )
    b <- coef(fit)[[score_var]]
    se <- sqrt(vcov(fit)[score_var, score_var])
    data.table(
      score = score, score_label = CIPDS_SCORE_LABELS[[score]], cycle = cycle,
      standardization = standardization, n = nrow(data_obj), events = sum(data_obj$Death_AllCause),
      log_hazard_ratio = b, standard_error = se,
      hazard_ratio = exp(b), ci_lower = exp(b - 1.96 * se), ci_upper = exp(b + 1.96 * se),
      p_value = 2 * pnorm(-abs(b / se)),
      design_degrees_of_freedom = uniqueN(data_obj$survey_psu) - uniqueN(data_obj$survey_strata),
      estimable = TRUE, error = NA_character_
    )
  }, error = function(e) data.table(
    score = score, score_label = CIPDS_SCORE_LABELS[[score]], cycle = cycle,
    standardization = standardization, n = nrow(data_obj), events = sum(data_obj$Death_AllCause),
    log_hazard_ratio = NA_real_, standard_error = NA_real_, hazard_ratio = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
    design_degrees_of_freedom = uniqueN(data_obj$survey_psu) - uniqueN(data_obj$survey_strata),
    estimable = FALSE, error = conditionMessage(e)
  ))
}

robust_block_test <- function(fit, terms) {
  b <- coef(fit)[terms]
  v <- vcov(fit)[terms, terms, drop = FALSE]
  keep <- is.finite(b) & apply(v, 1, function(z) all(is.finite(z)))
  b <- b[keep]
  v <- v[keep, keep, drop = FALSE]
  if (!length(b)) return(c(statistic = NA_real_, df = 0, p = NA_real_))
  rank <- qr(v)$rank
  statistic <- as.numeric(t(b) %*% MASS::ginv(v) %*% b)
  c(statistic = statistic, df = rank, p = pchisq(statistic, df = rank, lower.tail = FALSE))
}

cycle_rows <- list()
interaction_rows <- list()
audit_rows <- list()

for (score in names(CIPDS_SCORE_VARS)) {
  score_code <- score
  score_var <- CIPDS_SCORE_VARS[[score]]
  mask <- cipds_complete_mask(dat, score_var)
  ds <- droplevels(copy(dat[mask]))
  ds[, CYCLE := droplevels(CYCLE)]
  ds[, era3 := factor(fcase(
    CYCLE %in% c("1999-2000", "2001-2002", "2003-2004"), "Early: 1999-2004",
    CYCLE %in% c("2005-2006", "2007-2008", "2009-2010"), "Middle: 2005-2010",
    default = "Late: 2011-2016"
  ), levels = c("Early: 1999-2004", "Middle: 2005-2010", "Late: 2011-2016"))]
  design_score <- svydesign(
    ids = ~survey_psu, strata = ~survey_strata, weights = ~pooled_mec_weight,
    data = ds, nest = TRUE
  )

  for (cycle in levels(ds$CYCLE)) {
    mask_cycle <- ds$CYCLE == cycle
    dc <- ds[mask_cycle]
    cycle_design <- design_score[mask_cycle, ]
    cycle_rows[[length(cycle_rows) + 1L]] <- fit_effect(
      dc, score_var, score, cycle, "global age>=60 survey-weighted SD"
    )
    mu <- as.numeric(coef(svymean(as.formula(paste0("~", score_var)), cycle_design, na.rm = TRUE)))
    sigma <- sqrt(as.numeric(svyvar(as.formula(paste0("~", score_var)), cycle_design, na.rm = TRUE)))
    temp <- paste0("cycle_z_", tolower(score))
    dc[[temp]] <- (dc[[score_var]] - mu) / sigma
    cycle_rows[[length(cycle_rows) + 1L]] <- fit_effect(
      dc, temp, score, cycle, "cycle-specific survey-weighted SD"
    )
    audit_rows[[length(audit_rows) + 1L]] <- data.table(
      score = score, score_label = CIPDS_SCORE_LABELS[[score]], cycle = cycle,
      n = nrow(dc), events = sum(dc$Death_AllCause),
      weighted_population = sum(dc$pooled_mec_weight),
      median_follow_up = median(dc$Follow_Up_Years), maximum_follow_up = max(dc$Follow_Up_Years),
      cases_5y = sum(dc$Death_AllCause == 1 & dc$Follow_Up_Years < 5),
      controls_5y = sum(dc$Follow_Up_Years > 5),
      cases_10y = sum(dc$Death_AllCause == 1 & dc$Follow_Up_Years < 10),
      controls_10y = sum(dc$Follow_Up_Years > 10),
      cases_15y = sum(dc$Death_AllCause == 1 & dc$Follow_Up_Years < 15),
      controls_15y = sum(dc$Follow_Up_Years > 15)
    )
  }

  ds[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]
  full_formula <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    paste(c(paste0(score_var, " * CYCLE"), CIPDS_CLINICAL_TERMS), collapse = " + ")
  ))
  fit_cycle <- coxph(
    full_formula, data = ds, weights = normalized_sampling_weight,
    cluster = survey_psu, robust = TRUE, x = TRUE
  )
  interaction_terms <- names(coef(fit_cycle))[
    grepl(paste0("^", score_var, ":CYCLE"), names(coef(fit_cycle)))
  ]
  cycle_test <- robust_block_test(fit_cycle, interaction_terms)

  d_el <- droplevels(ds[era3 != "Middle: 2005-2010"])
  d_el[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]
  el_formula <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~", score_var,
    "* era3 +", paste(CIPDS_CLINICAL_TERMS, collapse = " + ")
  ))
  fit_el <- coxph(
    el_formula, data = d_el, weights = normalized_sampling_weight,
    cluster = survey_psu, robust = TRUE, x = TRUE
  )
  el_term <- names(coef(fit_el))[
    grepl(paste0("^", score_var, ":era3"), names(coef(fit_el)))
  ]
  el_test <- robust_block_test(fit_el, el_term)
  ratio_beta <- coef(fit_el)[[el_term]]
  ratio_se <- sqrt(vcov(fit_el)[el_term, el_term])

  fixed <- rbindlist(cycle_rows)[
    score == score_code & standardization == "global age>=60 survey-weighted SD" & estimable == TRUE
  ]
  precision <- 1 / fixed$standard_error^2
  pooled_beta <- sum(precision * fixed$log_hazard_ratio) / sum(precision)
  q <- sum(precision * (fixed$log_hazard_ratio - pooled_beta)^2)
  q_df <- nrow(fixed) - 1L
  q_p <- pchisq(q, df = q_df, lower.tail = FALSE)
  i2 <- ifelse(q > 0, max(0, (q - q_df) / q) * 100, 0)

  interaction_rows[[length(interaction_rows) + 1L]] <- data.table(
    score = score, score_label = CIPDS_SCORE_LABELS[[score]], n = nrow(ds),
    events = sum(ds$Death_AllCause), cycles_observed = nlevels(ds$CYCLE),
    cycle_interaction_df = cycle_test[["df"]], cycle_interaction_p = cycle_test[["p"]],
    early_late_ratio_of_hrs = exp(ratio_beta),
    early_late_ratio_ci_lower = exp(ratio_beta - 1.96 * ratio_se),
    early_late_ratio_ci_upper = exp(ratio_beta + 1.96 * ratio_se),
    early_late_interaction_p = el_test[["p"]],
    descriptive_cochran_q = q, descriptive_q_df = q_df,
    descriptive_q_p = q_p, descriptive_i2_percent = i2
  )
}

cycle_effects <- rbindlist(cycle_rows, fill = TRUE)
interaction <- rbindlist(interaction_rows)
cycle_audit <- rbindlist(audit_rows)
interaction[, cycle_interaction_p_holm := p.adjust(cycle_interaction_p, method = "holm")]
interaction[, early_late_interaction_p_holm := p.adjust(early_late_interaction_p, method = "holm")]

# Demographic-adjusted sensitivity retains early cycles that lack some lifestyle
# covariates and therefore distinguishes temporal instability from covariate
# availability. PhenoAge remains limited to cycles with all formula inputs.
demographic_terms <- c("Age", "Sex_F_design", "Race_F_design", "Cancer_F_design")
make_demographic_formula <- function(score_term, interaction_term = NULL) {
  predictors <- if (is.null(interaction_term)) score_term else interaction_term
  as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    paste(c(predictors, demographic_terms), collapse = " + ")
  ))
}
demographic_effect_rows <- list()
demographic_interaction_rows <- list()
for (score in names(CIPDS_SCORE_VARS)) {
  score_var <- CIPDS_SCORE_VARS[[score]]
  required <- c("Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", score_var, demographic_terms)
  mask <- !is.na(dat$Age) & dat$Age >= 60 & complete.cases(dat[, ..required]) &
    dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
  ds <- droplevels(copy(dat[mask]))
  ds[, CYCLE := droplevels(CYCLE)]
  ds[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]
  for (cycle in levels(ds$CYCLE)) {
    dc <- droplevels(ds[CYCLE == cycle])
    dc[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]
    fit <- coxph(
      make_demographic_formula(score_var), data = dc,
      weights = normalized_sampling_weight, cluster = survey_psu, robust = TRUE
    )
    b <- coef(fit)[[score_var]]
    se <- sqrt(vcov(fit)[score_var, score_var])
    demographic_effect_rows[[length(demographic_effect_rows) + 1L]] <- data.table(
      score = score, score_label = CIPDS_SCORE_LABELS[[score]], cycle = cycle,
      n = nrow(dc), events = sum(dc$Death_AllCause),
      hazard_ratio = exp(b), ci_lower = exp(b - 1.96 * se), ci_upper = exp(b + 1.96 * se),
      p_value = 2 * pnorm(-abs(b / se)), adjustment = "age + sex + race/ethnicity + cancer history"
    )
  }
  fit_interaction <- coxph(
    make_demographic_formula(score_var, paste0(score_var, " * CYCLE")),
    data = ds, weights = normalized_sampling_weight, cluster = survey_psu, robust = TRUE
  )
  terms <- names(coef(fit_interaction))[
    grepl(paste0("^", score_var, ":CYCLE"), names(coef(fit_interaction)))
  ]
  test <- robust_block_test(fit_interaction, terms)
  demographic_interaction_rows[[length(demographic_interaction_rows) + 1L]] <- data.table(
    score = score, score_label = CIPDS_SCORE_LABELS[[score]], n = nrow(ds),
    events = sum(ds$Death_AllCause), cycles_observed = nlevels(ds$CYCLE),
    cycle_interaction_df = test[["df"]], cycle_interaction_p = test[["p"]],
    adjustment = "age + sex + race/ethnicity + cancer history"
  )
}
demographic_effects <- rbindlist(demographic_effect_rows)
demographic_interactions <- rbindlist(demographic_interaction_rows)
demographic_interactions[, cycle_interaction_p_holm := p.adjust(cycle_interaction_p, method = "holm")]

fwrite(cycle_effects, file.path(CIPDS_OUT, "cycle_specific_score_effects.csv"), bom = TRUE)
fwrite(interaction, file.path(CIPDS_OUT, "cycle_stability_interactions.csv"), bom = TRUE)
fwrite(cycle_audit, file.path(CIPDS_OUT, "cycle_stability_cohort_audit.csv"), bom = TRUE)
fwrite(demographic_effects, file.path(CIPDS_OUT, "cycle_specific_demographic_adjusted_effects.csv"), bom = TRUE)
fwrite(demographic_interactions, file.path(CIPDS_OUT, "cycle_stability_demographic_adjusted_interactions.csv"), bom = TRUE)

cipds_write_manifest(file.path(CIPDS_OUT, "cycle_stability_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  cohort_strategy = "score-specific fully adjusted complete-case cohort so components and Overall retain all available NHANES cycles",
  main_standardization = "global survey-weighted SD in the age >=60 domain",
  sensitivity_standardization = "within-cycle survey-weighted SD",
  early_cycle_availability_sensitivity = "demographic adjustment (age, sex, race/ethnicity, cancer history) retains all cycles in which each score is computable",
  primary_heterogeneity_test = "sampling-weighted score by cycle interaction with PSU-cluster robust covariance",
  secondary_heterogeneity_test = "sampling-weighted early 1999-2004 versus late 2011-2016 score interaction with PSU-cluster robust covariance",
  descriptive_only = "Cochran Q and I2 calculated from cycle-specific log hazard ratios",
  multiplicity = "Holm correction across five scores"
))
cat("\nCycle audit:\n")
print(cycle_audit)
cat("\nInteraction tests:\n")
print(interaction)
cat("\nDemographic-adjusted interaction sensitivity:\n")
print(demographic_interactions)
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
