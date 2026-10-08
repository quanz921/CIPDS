source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))
suppressPackageStartupMessages({
  library(Hmisc)
  library(splines)
})

log_file <- file.path(CIPDS_LOG, "29_ph_nonlinearity.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

x <- cipds_load_data(FALSE)
dat <- x$dat
design_full <- x$design_full
common <- cipds_complete_mask(dat, unname(CIPDS_SCORE_VARS))
d <- copy(dat[common])
d[, normalized_sampling_weight := pooled_mec_weight / mean(pooled_mec_weight)]
design <- design_full[common, ]
expected <- fread(file.path(CIPDS_ROOT,"qa/cohort_registry.csv"))[cohort=="common"]
stopifnot(nrow(d)==expected$n, sum(d$Death_AllCause)==expected$events)

ph_rows <- list()
spline_rows <- list()
curve_rows <- list()

for (score in names(CIPDS_SCORE_VARS)) {
  score_var <- CIPDS_SCORE_VARS[[score]]
  formula_svy <- cipds_make_formula(score_var)
  linear_fit <- svycoxph(formula_svy, design = design)
  b <- coef(linear_fit)[[score_var]]
  ci <- confint(linear_fit)[score_var, ]
  p <- summary(linear_fit)$coefficients[score_var, "Pr(>|z|)"]

  formula_robust <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    paste(c(score_var, CIPDS_CLINICAL_TERMS), collapse = " + ")
  ))
  robust_fit <- coxph(
    formula_robust, data = d, weights = normalized_sampling_weight,
    cluster = survey_psu, robust = TRUE, x = TRUE, y = TRUE
  )
  zph <- cox.zph(robust_fit, transform = "km", global = TRUE)
  zph_table <- as.data.frame(zph$table)
  zph_score_p <- zph_table[score_var, "p"]
  zph_global_p <- zph_table["GLOBAL", "p"]

  formula_tt <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    paste(c(score_var, CIPDS_CLINICAL_TERMS, paste0("tt(", score_var, ")")), collapse = " + ")
  ))
  tt_fit <- coxph(
    formula_tt, data = d, weights = normalized_sampling_weight,
    cluster = survey_psu, robust = TRUE,
    tt = function(x, t, ...) x * log(pmax(t, 0.25))
  )
  tt_name <- names(coef(tt_fit))[grepl("^tt\\(", names(coef(tt_fit)))]
  tt_beta <- coef(tt_fit)[[tt_name]]
  tt_se <- sqrt(vcov(tt_fit)[tt_name, tt_name])
  tt_p <- 2 * pnorm(-abs(tt_beta / tt_se))

  ph_rows[[length(ph_rows) + 1L]] <- data.table(
    score = score, score_label = CIPDS_SCORE_LABELS[[score]], n = nrow(d),
    events = sum(d$Death_AllCause), hazard_ratio_per_weighted_sd = exp(b),
    ci_lower = exp(ci[1]), ci_upper = exp(ci[2]), linear_association_p = p,
    schoenfeld_score_p = zph_score_p, schoenfeld_global_p = zph_global_p,
    score_log_time_interaction_beta = tt_beta,
    score_log_time_interaction_se = tt_se,
    score_log_time_interaction_p = tt_p,
    ph_test_method = "sampling-weighted Schoenfeld score test; separate PSU-cluster-robust score-by-log-time interaction"
  )

  knots <- cipds_weighted_quantile(d[[score_var]], d$pooled_mec_weight, c(0.05, 0.35, 0.65, 0.95))
  basis <- rcspline.eval(d[[score_var]], knots = knots, inclx = TRUE)
  basis_names <- paste0("rcs_", tolower(score), "_", seq_len(ncol(basis)))
  for (j in seq_len(ncol(basis))) d[[basis_names[j]]] <- basis[, j]
  design_spline <- svydesign(
    ids = ~survey_psu, strata = ~survey_strata, weights = ~pooled_mec_weight,
    data = d, nest = TRUE
  )
  spline_formula <- cipds_make_formula(basis_names)
  spline_fit <- svycoxph(spline_formula, design = design_spline)
  nonlinear_terms <- basis_names[-1]
  p_overall <- as.numeric(regTermTest(
    spline_fit, as.formula(paste("~", paste(basis_names, collapse = " + ")))
  )$p)
  p_nonlinear <- as.numeric(regTermTest(
    spline_fit, as.formula(paste("~", paste(nonlinear_terms, collapse = " + ")))
  )$p)
  spline_rows[[length(spline_rows) + 1L]] <- data.table(
    score = score, score_label = CIPDS_SCORE_LABELS[[score]], n = nrow(d),
    events = sum(d$Death_AllCause), knot_1 = knots[1], knot_2 = knots[2],
    knot_3 = knots[3], knot_4 = knots[4], spline_overall_p = p_overall,
    spline_nonlinearity_p = p_nonlinear,
    spline_definition = "restricted cubic spline with survey-weighted 5th, 35th, 65th, and 95th percentile knots"
  )

  beta <- coef(spline_fit)[basis_names]
  variance <- vcov(spline_fit)[basis_names, basis_names, drop = FALSE]
  grid_range <- cipds_weighted_quantile(d[[score_var]], d$pooled_mec_weight, c(0.01, 0.99))
  grid <- seq(grid_range[1], grid_range[2], length.out = 121)
  ref <- cipds_weighted_quantile(d[[score_var]], d$pooled_mec_weight, 0.50)
  grid_basis <- rcspline.eval(grid, knots = knots, inclx = TRUE)
  ref_basis <- as.numeric(rcspline.eval(ref, knots = knots, inclx = TRUE))
  contrast <- sweep(grid_basis, 2, ref_basis, "-")
  log_hr <- as.numeric(contrast %*% beta)
  se_hr <- sqrt(rowSums((contrast %*% variance) * contrast))
  density <- density(d[[score_var]], weights = d$pooled_mec_weight / sum(d$pooled_mec_weight),
                     from = grid_range[1], to = grid_range[2], n = 512)
  curve_rows[[length(curve_rows) + 1L]] <- data.table(
    score = score, score_label = CIPDS_SCORE_LABELS[[score]], score_value = grid,
    reference_value = ref, hazard_ratio = exp(log_hr),
    ci_lower = exp(log_hr - 1.96 * se_hr), ci_upper = exp(log_hr + 1.96 * se_hr),
    weighted_density = approx(density$x, density$y, xout = grid, rule = 2)$y
  )
}

ph <- rbindlist(ph_rows)
spline <- rbindlist(spline_rows)
curves <- rbindlist(curve_rows)
ph[, schoenfeld_score_p_holm := p.adjust(schoenfeld_score_p, method = "holm")]
ph[, score_log_time_interaction_p_holm := p.adjust(score_log_time_interaction_p, method = "holm")]
spline[, spline_nonlinearity_p_holm := p.adjust(spline_nonlinearity_p, method = "holm")]
spline[, spline_overall_p_holm := p.adjust(spline_overall_p, method = "holm")]

fwrite(ph, file.path(CIPDS_OUT, "ph_assumption_tests.csv"), bom = TRUE)
fwrite(spline, file.path(CIPDS_OUT, "restricted_cubic_spline_tests.csv"), bom = TRUE)
fwrite(curves, file.path(CIPDS_OUT, "restricted_cubic_spline_curves.csv"), bom = TRUE)

cipds_write_manifest(file.path(CIPDS_OUT, "ph_nonlinearity_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  cohort_n = nrow(d), events = sum(d$Death_AllCause),
  common_complete_case = TRUE,
  scores = names(CIPDS_SCORE_VARS),
  proportional_hazards = c(
    "sampling-weighted Schoenfeld score test (model-based diagnostic)",
    "score by log(time) interaction sensitivity"
  ),
  nonlinearity = "survey Cox model with four-knot restricted cubic spline",
  multiplicity = "Holm correction across five scores within each test family"
))
cat("\nPH tests:\n")
print(ph[, .(score, schoenfeld_score_p, schoenfeld_score_p_holm,
             score_log_time_interaction_p, score_log_time_interaction_p_holm)])
cat("\nSpline tests:\n")
print(spline[, .(score, spline_overall_p, spline_nonlinearity_p, spline_nonlinearity_p_holm)])
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
