source(file.path(
  Sys.getenv("CIPDS_PACKAGE_DIR", unset = "cipds_rebuild_20260831"),
  "scripts", "27_supplement_common.R"
))
suppressPackageStartupMessages(library(splines))

log_file <- file.path(CIPDS_LOG, "28_selection_bias_ipw.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

x <- cipds_load_data(FALSE)
dat <- x$dat
design_full <- x$design_full
source_mask <- x$age60 & is.finite(dat$pooled_mec_weight) & dat$pooled_mec_weight > 0 &
  is.finite(dat$Follow_Up_Years) & dat$Follow_Up_Years > 0 & !is.na(dat$Death_AllCause)
overall_mask <- source_mask & !is.na(dat$OVERALL_z_age60)
component_mask <- cipds_complete_mask(dat, c("NM_z_age60", "TB_z_age60", "TC_z_age60"))
paired_mask <- source_mask & complete.cases(dat[, .(
  NM_z_age60, TB_z_age60, TC_z_age60, OVERALL_z_age60, PHENO_z_age60
)])
incremental_mask <- cipds_complete_mask(dat, unname(CIPDS_SCORE_VARS))

flow <- data.table(
  node = c("source_age60", "overall_available", "component_fully_adjusted",
           "paired_score_only", "incremental_common_complete"),
  parent = c(NA_character_, "source_age60", "source_age60", "source_age60", "source_age60"),
  analysis_role = c("SOURCE", "OVERALL_ASSOCIATION", "COMPONENT_ASSOCIATION",
                    "DIRECT_DISCRIMINATION", "INCREMENTAL_PREDICTION"),
  n = c(sum(source_mask), sum(overall_mask), sum(component_mask), sum(paired_mask), sum(incremental_mask)),
  deaths = c(sum(dat$Death_AllCause[source_mask]), sum(dat$Death_AllCause[overall_mask]),
             sum(dat$Death_AllCause[component_mask]), sum(dat$Death_AllCause[paired_mask]),
             sum(dat$Death_AllCause[incremental_mask]))
)
flow[, excluded_from_source_n := n[1] - n]
flow[, retained_percent := 100 * n / n[1]]
fwrite(flow, file.path(CIPDS_OUT, "selection_participant_flow.csv"), bom = TRUE)

weighted_mean_sd <- function(v, w) {
  keep <- is.finite(v) & is.finite(w) & w > 0
  if (sum(keep) < 2L) return(c(mean = NA_real_, sd = NA_real_))
  mu <- weighted.mean(v[keep], w[keep])
  variance <- sum(w[keep] * (v[keep] - mu)^2) / sum(w[keep])
  c(mean = mu, sd = sqrt(variance))
}

baseline_rows <- list()
continuous_vars <- c(
  Age = "Age", PIR = "Poverty-income ratio", BMI = "Body mass index",
  Follow_Up_Years = "Follow-up, years", Overall_expected_burden = "Overall score",
  PhenoAge_acceleration = "PhenoAge Acceleration, years"
)
for (v in names(continuous_vars)) {
  inc <- incremental_mask
  exc <- source_mask & !incremental_mask
  a <- weighted_mean_sd(dat[[v]][inc], dat$pooled_mec_weight[inc])
  b <- weighted_mean_sd(dat[[v]][exc], dat$pooled_mec_weight[exc])
  pooled_sd <- sqrt((a[["sd"]]^2 + b[["sd"]]^2) / 2)
  baseline_rows[[length(baseline_rows) + 1L]] <- data.table(
    variable = v, variable_label = continuous_vars[[v]], type = "continuous", level = NA_character_,
    included_unweighted_n = sum(inc & !is.na(dat[[v]])),
    excluded_unweighted_n = sum(exc & !is.na(dat[[v]])),
    included_weighted_value = a[["mean"]], excluded_weighted_value = b[["mean"]],
    included_weighted_sd = a[["sd"]], excluded_weighted_sd = b[["sd"]],
    standardized_mean_difference = ifelse(is.finite(pooled_sd) && pooled_sd > 0,
                                          (a[["mean"]] - b[["mean"]]) / pooled_sd, NA_real_),
    included_missing_percent = 100 * sum(inc & is.na(dat[[v]])) / sum(inc),
    excluded_missing_percent = 100 * sum(exc & is.na(dat[[v]])) / sum(exc)
  )
}

categorical_vars <- c(
  Sex = "Sex", RIDRETH1 = "Race/ethnicity", Education_clean = "Education",
  Smoking_F = "Smoking", Alcohol_F = "Alcohol", Diabetes = "Diabetes",
  Hypertension = "Hypertension", CVD = "Cardiovascular disease",
  Cancer_Diagnosed = "Cancer history", Death_AllCause = "All-cause death", CYCLE = "NHANES cycle"
)
for (v in names(categorical_vars)) {
  z <- as.character(dat[[v]])
  z[is.na(z) | z == ""] <- "Missing"
  levels_v <- sort(unique(z[source_mask]))
  for (lev in levels_v) {
    inc <- incremental_mask
    exc <- source_mask & !incremental_mask
    p_inc <- sum(dat$pooled_mec_weight[inc] * (z[inc] == lev)) / sum(dat$pooled_mec_weight[inc])
    p_exc <- sum(dat$pooled_mec_weight[exc] * (z[exc] == lev)) / sum(dat$pooled_mec_weight[exc])
    p_pool <- (p_inc + p_exc) / 2
    denom <- sqrt(p_pool * (1 - p_pool))
    baseline_rows[[length(baseline_rows) + 1L]] <- data.table(
      variable = v, variable_label = categorical_vars[[v]], type = "categorical", level = lev,
      included_unweighted_n = sum(inc & z == lev), excluded_unweighted_n = sum(exc & z == lev),
      included_weighted_value = p_inc, excluded_weighted_value = p_exc,
      included_weighted_sd = NA_real_, excluded_weighted_sd = NA_real_,
      standardized_mean_difference = ifelse(denom > 0, (p_inc - p_exc) / denom, 0),
      included_missing_percent = 100 * sum(inc & z == "Missing") / sum(inc),
      excluded_missing_percent = 100 * sum(exc & z == "Missing") / sum(exc)
    )
  }
}
baseline <- rbindlist(baseline_rows, fill = TRUE)
baseline[, absolute_smd := abs(standardized_mean_difference)]
setorder(baseline, -absolute_smd, variable, level)
fwrite(baseline, file.path(CIPDS_OUT, "selection_weighted_baseline_smd.csv"), bom = TRUE)

# Complete-case inverse-probability weighting. Selection is modeled in the full
# age >=60 survey domain using variables observed before defining complete cases.
dat[, selected_incremental := as.integer(incremental_mask)]
dat[, Sex_sel := addNA(factor(Sex), ifany = TRUE)]
dat[, Race_sel := addNA(factor(RIDRETH1), ifany = TRUE)]
dat[, Cancer_sel := factor(Cancer_Diagnosed)]
supported_cycles <- unique(as.character(dat$CYCLE[incremental_mask]))
selection_source <- source_mask & as.character(dat$CYCLE) %in% supported_cycles & !is.na(dat$Cancer_Diagnosed)
support_audit <- dat[source_mask, .(source_n=.N, selected_n=sum(selected_incremental)), by=CYCLE]
support_audit[, supported := selected_n>0]
fwrite(support_audit,file.path(CIPDS_OUT,"selection_cycle_support.csv"))
selection_design <- design_full[selection_source, ]
selection_fit <- svyglm(
  selected_incremental ~ ns(Age, df = 4) + Sex_sel + Race_sel + Cancer_sel +
    CYCLE + Death_AllCause + ns(Follow_Up_Years, df = 4),
  design = selection_design,
  family = quasibinomial()
)
dat[selection_source, selection_probability := as.numeric(predict(selection_fit, type = "response"))]
selection_prevalence <- sum(dat$pooled_mec_weight[selection_source] * dat$selected_incremental[selection_source]) /
  sum(dat$pooled_mec_weight[selection_source])
dat[incremental_mask, ipw_raw := selection_prevalence / pmax(selection_probability, 0.01)]
trim <- quantile(dat$ipw_raw[incremental_mask], c(0.01, 0.99), na.rm = TRUE, names = FALSE)
dat[incremental_mask, ipw_trimmed := pmin(pmax(ipw_raw, trim[1]), trim[2])]
dat[incremental_mask, combined_weight_ipw := pooled_mec_weight * ipw_trimmed]

ipw_audit <- data.table(
  source_n = sum(selection_source), full_source_n=sum(source_mask), structurally_excluded_n=sum(source_mask & !selection_source), excluded_cycle_n=sum(source_mask & !(as.character(dat$CYCLE) %in% supported_cycles)), excluded_unknown_cancer_n=sum(source_mask & as.character(dat$CYCLE) %in% supported_cycles & is.na(dat$Cancer_Diagnosed)), selected_n = sum(incremental_mask),
  source_deaths = sum(dat$Death_AllCause[selection_source]),
  selected_deaths = sum(dat$Death_AllCause[incremental_mask]),
  survey_weighted_selection_prevalence = selection_prevalence,
  selection_probability_min = min(dat$selection_probability[selection_source]),
  selection_probability_p01 = quantile(dat$selection_probability[incremental_mask], 0.01),
  selection_probability_median = median(dat$selection_probability[incremental_mask]),
  selection_probability_p99 = quantile(dat$selection_probability[incremental_mask], 0.99),
  selection_probability_max = max(dat$selection_probability[selection_source]),
  ipw_trim_lower = trim[1], ipw_trim_upper = trim[2],
  ipw_effective_sample_size = with(dat[incremental_mask], sum(combined_weight_ipw)^2 / sum(combined_weight_ipw^2))
)
fwrite(ipw_audit, file.path(CIPDS_OUT, "selection_ipw_audit.csv"), bom = TRUE)

model_sets <- list(
  M0_CLINICAL = character(),
  M1_CLINICAL_OVERALL = "OVERALL_z_age60",
  M2_CLINICAL_PHENO = "PHENO_z_age60",
  M3_CLINICAL_OVERALL_PHENO = c("OVERALL_z_age60", "PHENO_z_age60"),
  M4_CLINICAL_COMPONENTS = c("NM_z_age60", "TB_z_age60", "TC_z_age60"),
  M5_CLINICAL_COMPONENTS_PHENO = c("NM_z_age60", "TB_z_age60", "TC_z_age60", "PHENO_z_age60")
)
design_standard <- design_full[incremental_mask, ]
design_ipw <- svydesign(
  ids = ~survey_psu, strata = ~survey_strata, weights = ~combined_weight_ipw,
  data = dat[incremental_mask], nest = TRUE
)

ipw_rows <- list()
for (model_name in names(model_sets)) {
  predictors <- model_sets[[model_name]]
  formula <- cipds_make_formula(predictors)
  for (weighting in c("ORIGINAL_SURVEY", "COMPLETE_CASE_IPW")) {
    design <- if (weighting == "ORIGINAL_SURVEY") design_standard else design_ipw
    fit <- svycoxph(formula, design = design)
    b <- coef(fit)
    ci <- confint(fit)
    p <- summary(fit)$coefficients[, "Pr(>|z|)"]
    terms <- intersect(predictors, names(b))
    if (!length(terms)) terms <- "__MODEL__"
    for (term in terms) {
      if (term == "__MODEL__") {
        ipw_rows[[length(ipw_rows) + 1L]] <- data.table(
          model = model_name, weighting = weighting, term = term, score = NA_character_,
          n = sum(incremental_mask), events = sum(dat$Death_AllCause[incremental_mask]),
          hazard_ratio = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_
        )
      } else {
        score <- names(CIPDS_SCORE_VARS)[match(term, CIPDS_SCORE_VARS)]
        ipw_rows[[length(ipw_rows) + 1L]] <- data.table(
          model = model_name, weighting = weighting, term = term, score = score,
          n = sum(incremental_mask), events = sum(dat$Death_AllCause[incremental_mask]),
          hazard_ratio = exp(b[[term]]), ci_lower = exp(ci[term, 1]), ci_upper = exp(ci[term, 2]),
          p_value = p[[term]]
        )
      }
    }
  }
}
ipw_effects <- rbindlist(ipw_rows, fill = TRUE)
fwrite(ipw_effects, file.path(CIPDS_OUT, "selection_ipw_association_sensitivity.csv"), bom = TRUE)
saveRDS(selection_fit, file.path(CIPDS_OUT, "selection_ipw_model.rds"), compress = "xz")

stopifnot(flow[node == "source_age60", n] == 15048L)
stopifnot(flow[node == "overall_available", n] == 14715L)
stopifnot(flow[node == "component_fully_adjusted", n] == fread(file.path(CIPDS_ROOT,"qa/cohort_registry.csv"))[cohort=="components", n])
stopifnot(flow[node == "paired_score_only", n] == fread(file.path(CIPDS_ROOT,"qa/cohort_registry.csv"))[cohort=="paired", n])
stopifnot(flow[node == "incremental_common_complete", n] == fread(file.path(CIPDS_ROOT,"qa/cohort_registry.csv"))[cohort=="common", n])

cipds_write_manifest(file.path(CIPDS_OUT, "selection_manifest.json"), list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  source_population = "NHANES 1999-2016 survey domain age >=60 years",
  full_design_before_domain_restriction = TRUE,
  primary_selection_endpoint = "complete clinical covariates plus all five score variables",
  sensitivity = "stabilized complete-case inverse-probability weights trimmed at the 1st and 99th percentiles",
  note = "Analysis cohorts are parallel branches. IPW targets older adults with known cancer history in cycles with observable complete cases; structurally unsupported cycles are excluded from the selection model."
))
cat("\nParticipant flow:\n")
print(flow)
cat("\nIPW audit:\n")
print(ipw_audit)
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
