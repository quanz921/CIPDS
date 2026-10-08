source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))

log_file <- file.path(CIPDS_LOG, "39_time_reverse_age_extension.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

OUT3 <- file.path(CIPDS_ROOT, "outputs", "supplement_v3")
dir.create(OUT3, recursive = TRUE, showWarnings = FALSE)

x <- cipds_load_data(FALSE)
dat <- x$dat
common <- cipds_complete_mask(dat, unname(CIPDS_SCORE_VARS))
expected <- fread(file.path(CIPDS_ROOT,"qa/cohort_registry.csv"))[cohort=="common"]
stopifnot(sum(common)==expected$n, sum(dat$Death_AllCause[common])==expected$events)

clinical_no_age <- setdiff(CIPDS_CLINICAL_TERMS, "Age")

make_design <- function(d) {
  svydesign(
    ids = ~survey_psu, strata = ~survey_strata, weights = ~pooled_mec_weight,
    data = d, nest = TRUE
  )
}

extract_svy_term <- function(fit, term) {
  beta <- coef(fit)[[term]]
  se <- sqrt(vcov(fit)[term, term])
  data.table(
    beta = beta, standard_error = se, hazard_ratio = exp(beta),
    ci_lower = exp(beta - qnorm(0.975) * se),
    ci_upper = exp(beta + qnorm(0.975) * se),
    p_value = 2 * pnorm(-abs(beta / se))
  )
}

# -------------------------------------------------------------------------
# 1) Baseline and correctly defined 1- and 2-year landmark analyses.
# Only participants still under observation and alive at the landmark enter;
# follow-up is restarted at the landmark to avoid immortal-time distortion.
# -------------------------------------------------------------------------
landmark_rows <- list()
for (landmark in c(0, 1, 2)) {
  eligible <- common & dat$Follow_Up_Years > landmark
  dat[, time_from_landmark := Follow_Up_Years - landmark]
  d <- copy(dat[eligible])
  # Build the multistage design before domain restriction, as required for
  # NHANES subpopulation inference.
  design <- make_design(dat)[eligible, ]
  for (score in names(CIPDS_SCORE_VARS)) {
    score_var <- CIPDS_SCORE_VARS[[score]]
    f <- as.formula(paste(
      "Surv(time_from_landmark, Death_AllCause) ~",
      paste(c(score_var, CIPDS_CLINICAL_TERMS), collapse = " + ")
    ))
    fit <- svycoxph(f, design = design)
    e <- extract_svy_term(fit, score_var)
    e[, `:=`(
      score = score,
      score_label = CIPDS_SCORE_LABELS[[score]],
      landmark_years = landmark,
      analysis = ifelse(
        landmark == 0, "Baseline follow-up-time analysis",
        paste0(landmark, "-year landmark; early deaths and early censoring excluded")
      ),
      n = nrow(d),
      events_after_landmark = sum(d$Death_AllCause),
      early_deaths_excluded = sum(dat$Death_AllCause[common] == 1 &
                                    dat$Follow_Up_Years[common] <= landmark),
      early_censoring_excluded = sum(dat$Death_AllCause[common] == 0 &
                                      dat$Follow_Up_Years[common] <= landmark),
      time_origin = ifelse(landmark == 0, "NHANES examination", paste0(landmark, "-year landmark")),
      variance = "NHANES stratified multistage survey design"
    )]
    landmark_rows[[length(landmark_rows) + 1L]] <- e
  }
}
landmark_results <- rbindlist(landmark_rows, fill = TRUE)
landmark_results[, p_holm_within_landmark := p.adjust(p_value, method = "holm"), by = landmark_years]
setcolorder(landmark_results, c(
  "score", "score_label", "landmark_years", "analysis", "n", "events_after_landmark",
  "early_deaths_excluded", "early_censoring_excluded", "hazard_ratio", "ci_lower",
  "ci_upper", "standard_error", "p_value", "p_holm_within_landmark", "beta",
  "time_origin", "variance"
))
fwrite(landmark_results, file.path(OUT3, "early_death_landmark_sensitivity.csv"), bom = TRUE)

# -------------------------------------------------------------------------
# 2) Attained age as the Cox time scale. Baseline age is deliberately omitted
# from the right-hand side because it defines entry time.
# -------------------------------------------------------------------------
d_common <- copy(dat[common])
dat[, attained_age_exit := Age + Follow_Up_Years]
d_common[, attained_age_exit := Age + Follow_Up_Years]
design_common <- make_design(dat)[common, ]

age_rows <- list()
for (score in names(CIPDS_SCORE_VARS)) {
  score_var <- CIPDS_SCORE_VARS[[score]]
  follow_formula <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    paste(c(score_var, CIPDS_CLINICAL_TERMS), collapse = " + ")
  ))
  age_formula <- as.formula(paste(
    "Surv(Age, attained_age_exit, Death_AllCause) ~",
    paste(c(score_var, clinical_no_age), collapse = " + ")
  ))
  for (definition in c("FOLLOW_UP_TIME", "ATTAINED_AGE")) {
    fit <- if (definition == "FOLLOW_UP_TIME") {
      svycoxph(follow_formula, design = design_common)
    } else {
      svycoxph(age_formula, design = design_common)
    }
    e <- extract_svy_term(fit, score_var)
    e[, `:=`(
      score = score,
      score_label = CIPDS_SCORE_LABELS[[score]],
      time_scale = definition,
      n = nrow(d_common),
      events = sum(d_common$Death_AllCause),
      baseline_age_rhs = definition == "FOLLOW_UP_TIME",
      time_definition = ifelse(
        definition == "FOLLOW_UP_TIME",
        "time since NHANES examination, with baseline age adjusted",
        "attained age from examination age to age at event/censoring; baseline age omitted from RHS"
      ),
      variance = "NHANES stratified multistage survey design"
    )]
    age_rows[[length(age_rows) + 1L]] <- e
  }
}
age_results <- rbindlist(age_rows, fill = TRUE)
age_results[, p_holm_within_time_scale := p.adjust(p_value, method = "holm"), by = time_scale]
setcolorder(age_results, c(
  "score", "score_label", "time_scale", "n", "events", "hazard_ratio",
  "ci_lower", "ci_upper", "standard_error", "p_value", "p_holm_within_time_scale",
  "beta", "baseline_age_rhs", "time_definition", "variance"
))
fwrite(age_results, file.path(OUT3, "age_time_scale_sensitivity.csv"), bom = TRUE)

# -------------------------------------------------------------------------
# 3) Piecewise score effects in 0-5, 5-10, and >10 years. A weighted robust
# counting-process Cox model permits separate baseline hazards by interval.
# -------------------------------------------------------------------------
d_split <- survSplit(
  Surv(Follow_Up_Years, Death_AllCause) ~ ., data = as.data.frame(d_common),
  cut = c(5, 10), start = "tstart", end = "tstop", event = "event_interval",
  episode = "interval_id"
)
d_split <- as.data.table(d_split)
d_split[, interval := factor(
  interval_id, levels = 1:3, labels = c("0-5 years", "5-10 years", ">10 years")
)]
d_split[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]

piecewise_rows <- list()
heterogeneity_rows <- list()
for (score in names(CIPDS_SCORE_VARS)) {
  score_var <- CIPDS_SCORE_VARS[[score]]
  f <- as.formula(paste(
    "Surv(tstart, tstop, event_interval) ~ 0 +",
    paste0(score_var, ":interval"), "+ strata(interval) +",
    paste(CIPDS_CLINICAL_TERMS, collapse = " + ")
  ))
  fit <- coxph(
    f, data = d_split, weights = normalized_sampling_weight,
    cluster = survey_psu, robust = TRUE, ties = "efron", x = TRUE, y = TRUE
  )
  term_names <- names(coef(fit))[grepl(paste0("^", score_var, ":interval"), names(coef(fit)))]
  if (length(term_names) != 3L) stop("Expected three interval coefficients for ", score)
  beta <- coef(fit)[term_names]
  variance <- vcov(fit)[term_names, term_names, drop = FALSE]
  se <- sqrt(diag(variance))
  interval_labels <- levels(d_split$interval)
  piecewise_rows[[length(piecewise_rows) + 1L]] <- data.table(
    score = score,
    score_label = CIPDS_SCORE_LABELS[[score]],
    interval = interval_labels,
    interval_coefficient = term_names,
    n_persons = nrow(d_common),
    person_intervals = nrow(d_split),
    events_in_interval = vapply(interval_labels, function(z) {
      sum(d_split$event_interval[d_split$interval == z])
    }, numeric(1)),
    beta = as.numeric(beta),
    standard_error = as.numeric(se),
    hazard_ratio = exp(as.numeric(beta)),
    ci_lower = exp(as.numeric(beta) - qnorm(0.975) * as.numeric(se)),
    ci_upper = exp(as.numeric(beta) + qnorm(0.975) * as.numeric(se)),
    p_value = 2 * pnorm(-abs(as.numeric(beta) / as.numeric(se))),
    model = "sampling-weighted Cox counting-process model, PSU-cluster robust variance, interval-stratified baseline hazard"
  )
  contrast <- rbind(c(1, -1, 0), c(1, 0, -1))
  difference <- as.numeric(contrast %*% beta)
  contrast_var <- contrast %*% variance %*% t(contrast)
  statistic <- as.numeric(t(difference) %*% solve(contrast_var, difference))
  heterogeneity_rows[[length(heterogeneity_rows) + 1L]] <- data.table(
    score = score,
    score_label = CIPDS_SCORE_LABELS[[score]],
    wald_chisq = statistic,
    df = 2L,
    interval_heterogeneity_p = pchisq(statistic, df = 2, lower.tail = FALSE)
  )
}
piecewise_results <- rbindlist(piecewise_rows)
piecewise_results[, p_holm_within_interval := p.adjust(p_value, method = "holm"), by = interval]
heterogeneity_results <- rbindlist(heterogeneity_rows)
heterogeneity_results[, interval_heterogeneity_p_holm := p.adjust(
  interval_heterogeneity_p, method = "holm"
)]
fwrite(piecewise_results, file.path(OUT3, "time_varying_piecewise_hr.csv"), bom = TRUE)
fwrite(heterogeneity_results, file.path(OUT3, "time_varying_interval_heterogeneity.csv"), bom = TRUE)

# -------------------------------------------------------------------------
# 4) Smooth HR(t) from the prespecified score x log(time) interaction.
# This is descriptive of the fitted time trend and complements the intervals.
# -------------------------------------------------------------------------
time_grid <- seq(0.5, 15, by = 0.25)
curve_rows <- list()
trend_rows <- list()
d_common[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]
for (score in names(CIPDS_SCORE_VARS)) {
  score_var <- CIPDS_SCORE_VARS[[score]]
  f <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    paste(c(score_var, CIPDS_CLINICAL_TERMS, paste0("tt(", score_var, ")")), collapse = " + ")
  ))
  fit <- coxph(
    f, data = d_common, weights = normalized_sampling_weight,
    cluster = survey_psu, robust = TRUE,
    tt = function(z, t, ...) z * log(pmax(t, 0.25)), ties = "efron"
  )
  score_name <- score_var
  tt_name <- names(coef(fit))[grepl("^tt\\(", names(coef(fit)))]
  if (length(tt_name) != 1L) stop("Time interaction coefficient missing for ", score)
  beta <- coef(fit)[c(score_name, tt_name)]
  variance <- vcov(fit)[c(score_name, tt_name), c(score_name, tt_name), drop = FALSE]
  design_matrix <- cbind(1, log(time_grid))
  log_hr <- as.numeric(design_matrix %*% beta)
  se <- sqrt(rowSums((design_matrix %*% variance) * design_matrix))
  curve_rows[[length(curve_rows) + 1L]] <- data.table(
    score = score,
    score_label = CIPDS_SCORE_LABELS[[score]],
    follow_up_years = time_grid,
    hazard_ratio_per_weighted_sd = exp(log_hr),
    ci_lower = exp(log_hr - qnorm(0.975) * se),
    ci_upper = exp(log_hr + qnorm(0.975) * se),
    model = "sampling-weighted Cox with PSU-cluster robust score by log(time) interaction"
  )
  gamma <- beta[[tt_name]]
  gamma_se <- sqrt(variance[tt_name, tt_name])
  trend_rows[[length(trend_rows) + 1L]] <- data.table(
    score = score,
    score_label = CIPDS_SCORE_LABELS[[score]],
    score_beta_at_one_year = beta[[score_name]],
    score_log_time_beta = gamma,
    score_log_time_se = gamma_se,
    score_log_time_p = 2 * pnorm(-abs(gamma / gamma_se))
  )
}
curve_results <- rbindlist(curve_rows)
trend_results <- rbindlist(trend_rows)
trend_results[, score_log_time_p_holm := p.adjust(score_log_time_p, method = "holm")]
fwrite(curve_results, file.path(OUT3, "time_varying_continuous_hr_curves.csv"), bom = TRUE)
fwrite(trend_results, file.path(OUT3, "time_varying_logtime_tests.csv"), bom = TRUE)

cipds_write_manifest(file.path(OUT3, "time_reverse_age_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  source_cohort_n = nrow(d_common),
  source_deaths = sum(d_common$Death_AllCause),
  scores = names(CIPDS_SCORE_VARS),
  time_varying_effects = c(
    "0-5, 5-10, and >10 year interval-specific HRs",
    "continuous HR(t) from score by log(time) interaction"
  ),
  reverse_causation = "correct 1-year and 2-year landmark analyses; early deaths and early censoring excluded, follow-up restarted",
  age_time_scale = "left-truncated attained-age Cox model; baseline age omitted from the right-hand side",
  multiplicity = "Holm correction across five scores within each analysis family",
  caution = "attained-age sensitivity inherits NHANES public-release age top-coding"
))

cat("\nLandmark results:\n")
print(landmark_results[, .(score, landmark_years, n, events_after_landmark, hazard_ratio, ci_lower, ci_upper, p_holm_within_landmark)])
cat("\nAge time-scale results:\n")
print(age_results[, .(score, time_scale, hazard_ratio, ci_lower, ci_upper, p_holm_within_time_scale)])
cat("\nPiecewise results:\n")
print(piecewise_results[, .(score, interval, events_in_interval, hazard_ratio, ci_lower, ci_upper)])
cat("\nInterval heterogeneity:\n")
print(heterogeneity_results)
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
