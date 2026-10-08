source(file.path(Sys.getenv("CIPDS_PACKAGE_DIR"), "scripts", "weight_contract.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survey)
  library(survival)
  library(MASS)
})

options(survey.lonely.psu = "adjust")

package_dir <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
out_dir <- file.path(package_dir, "outputs")
log_dir <- file.path(package_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, "15_nhanes_overall_model_primary.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

SCENARIO <- "A_PRIMARY_FROZEN_PANCANCER_OVERALL_STACK"
N_CYCLES <- 9L
EPS <- 1e-6
RAW_SCORES <- c(
  OVERALL = "Overall_expected_burden",
  NM = "NM_component",
  TB = "TB_component",
  TC = "TC_component"
)
LABELS <- c(
  OVERALL = "Multidomain laboratory vulnerability overall score",
  NM = "Nutritional-Metabolic laboratory pattern",
  TB = "Tumor-burden-related laboratory pattern",
  TC = "Treatment-complication-related laboratory pattern",
  PHENO = "PhenoAge Acceleration"
)

cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Overall model is frozen from hospital OOF component predictions.\n")
cat("PhenoAge and NHANES mortality were not used to train the overall model.\n")

load(file.path(package_dir, "inputs", "nhanes_phase7_batch1.RData"))
if (!exists("nhanes_valid")) stop("nhanes_valid missing")
dat <- as.data.table(nhanes_valid)
if (nrow(dat) != 44772L || uniqueN(dat$SEQN) != 44772L) stop("NHANES grain drift")
if (uniqueN(dat$CYCLE) != N_CYCLES) stop("NHANES cycle count drift")

scores <- fread(file.path(out_dir, "nested_a_primary_nhanes_component_scores.csv"))
if (nrow(scores) != 44772L || uniqueN(scores$SEQN) != 44772L) stop("Component score grain drift")
score_idx <- match(dat$SEQN, scores$SEQN)
if (anyNA(score_idx)) stop("Component score join failure")
for (v in c("NM_component", "TB_component", "TC_component")) {
  dat[[v]] <- scores[[v]][score_idx]
  if (anyNA(dat[[v]]) || any(dat[[v]] < 0 | dat[[v]] > 1)) stop("Invalid frozen score: ", v)
}

# Preserve the component-specific observation gate used by the older-adult analysis.
mapped <- readRDS(file.path(out_dir, "nhanes_expanded_candidate_matrix.rds"))
mapped_idx <- match(dat$SEQN, mapped$SEQN)
if (anyNA(mapped_idx)) stop("Mapped feature join failure")
features_dt <- fread(file.path(out_dir, "nested_a_primary_final_features.csv"))
outcome_to_component <- c(
  Outcome_NutriMetab = "NM",
  Outcome_TumorBurden = "TB",
  Outcome_TreatComp = "TC"
)
for (outcome_name in names(outcome_to_component)) {
  component <- outcome_to_component[[outcome_name]]
  selected <- features_dt[outcome == outcome_name, feature]
  if (!length(selected) || any(!selected %in% names(mapped))) {
    stop("Frozen final feature mapping failure for ", outcome_name)
  }
  observed_count <- rowSums(!is.na(as.matrix(mapped[mapped_idx, ..selected])))
  dat[[paste0(component, "_input_observed_n")]] <- observed_count
  component_col <- paste0(component, "_component")
  dat[observed_count == 0, (component_col) := NA_real_]
}

bundle <- readRDS(file.path(out_dir, "overall_stacked_ordinal_model.rds"))
if (isTRUE(bundle$phenoage_used) || isTRUE(bundle$nhanes_outcomes_used)) {
  stop("Overall model provenance violation")
}
if (!identical(bundle$component_order, c("NM", "TB", "TC"))) {
  stop("Overall model component order drift")
}
overall_eligible <- complete.cases(dat[, .(NM_component, TB_component, TC_component)])
meta_new <- data.frame(
  p_NM = dat$NM_component[overall_eligible],
  p_TB = dat$TB_component[overall_eligible],
  p_TC = dat$TC_component[overall_eligible]
)
scaling <- as.data.table(bundle$logit_scaling)
for (short in c("NM", "TB", "TC")) {
  logit_value <- qlogis(pmin(pmax(meta_new[[paste0("p_", short)]], EPS), 1 - EPS))
  scale_row <- scaling[component == short]
  if (nrow(scale_row) != 1L || !is.finite(scale_row$sd) || scale_row$sd <= 0) {
    stop("Overall model scaling registry error for ", short)
  }
  meta_new[[paste0("z_", short)]] <- (logit_value - scale_row$mean) / scale_row$sd
}
overall_prob <- as.matrix(predict(bundle$model, newdata = meta_new, type = "probs"))
if (!all(as.character(0:3) %in% colnames(overall_prob))) stop("Overall probability output drift")
overall_prob <- overall_prob[, as.character(0:3), drop = FALSE]
dat[, Overall_expected_burden := NA_real_]
dat[, Overall_probability_any_domain := NA_real_]
dat[, Overall_probability_multidomain := NA_real_]
dat[overall_eligible, Overall_expected_burden := as.numeric(overall_prob %*% 0:3)]
dat[overall_eligible, Overall_probability_any_domain := 1 - overall_prob[, "0"]]
dat[overall_eligible, Overall_probability_multidomain := overall_prob[, "2"] + overall_prob[, "3"]]
if (any(dat$Overall_expected_burden[overall_eligible] < 0 |
        dat$Overall_expected_burden[overall_eligible] > 3)) {
  stop("Overall score outside 0-3 range")
}

fwrite(
  dat[, .(
    SEQN, CYCLE, Age, Cancer_Diagnosed, Dead, Follow_Up_Years,
    NM_component, TB_component, TC_component,
    NM_input_observed_n, TB_input_observed_n, TC_input_observed_n,
    Overall_expected_burden,
    Overall_probability_any_domain,
    Overall_probability_multidomain
  )],
  file.path(out_dir, "overall_nhanes_scores.csv"),
  bom = TRUE
)

required <- c(
  "CYCLE", "WTMEC2YR", "SDMVSTRA", "SDMVPSU", "Follow_Up_Years", "Dead",
  "UCOD_LEADING", "Age", "Sex", "RIDRETH1", "Cancer_Diagnosed", "Education",
  "PIR", "BMI", "Smoking_F", "Alcohol_F", "Diabetes", "Hypertension", "CVD",
  unname(RAW_SCORES), "PhenoAge"
)
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required fields: ", paste(missing_required, collapse = ", "))

dat[, Death_AllCause := as.integer(Dead == 1)]
dat[, Death_Cancer := as.integer(Dead == 1 & UCOD_LEADING == 2)]
dat[, Death_CVD := fifelse(CYCLE == "2015-2016", NA_integer_, as.integer(Dead == 1 & UCOD_LEADING %in% c(1, 5)))]
dat[, Death_CLRD := fifelse(CYCLE == "2015-2016", NA_integer_, as.integer(Dead == 1 & UCOD_LEADING == 3))]
dat[, Death_Alzheimer := fifelse(CYCLE == "2015-2016", NA_integer_, as.integer(Dead == 1 & UCOD_LEADING == 6))]
dat[, Death_Diabetes_UCOD := fifelse(CYCLE == "2015-2016", NA_integer_, as.integer(Dead == 1 & UCOD_LEADING == 7))]
dat[, Death_Kidney := fifelse(CYCLE == "2015-2016", NA_integer_, as.integer(Dead == 1 & UCOD_LEADING == 9))]
if (sum(dat$Death_CVD, na.rm=TRUE) != sum(dat$Dead == 1 & dat$UCOD_LEADING %in% c(1, 5) & dat$CYCLE != "2015-2016")) {
  stop("CVD mortality definition drift")
}
if (any(dat$Dead == 1 & is.na(dat$UCOD_LEADING))) stop("Deceased participant missing UCOD")

dat[, pooled_mec_weight := cipds_pooled_mec_weights(dat, N_CYCLES)]
dat[, survey_strata := interaction(CYCLE, SDMVSTRA, drop = TRUE, lex.order = TRUE)]
dat[, survey_psu := interaction(CYCLE, SDMVSTRA, SDMVPSU, drop = TRUE, lex.order = TRUE)]
dat[, Sex_F_design := factor(Sex)]
dat[, Race_F_design := factor(RIDRETH1)]
dat[, Education_clean := ifelse(Education %in% 1:5, Education, NA_real_)]
dat[, Education_F_design := factor(Education_clean, levels = 1:5)]
dat[, Smoking_F_design := factor(Smoking_F)]
dat[, Alcohol_F_design := factor(Alcohol_F)]
dat[, Cancer_F_design := factor(Cancer_Diagnosed, levels = c(0, 1))]

# The complete nine-cycle survey design is created before age-domain restriction.
design_raw <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~pooled_mec_weight,
  data = dat,
  nest = TRUE
)
domains <- list(
  age60 = list(
    label = "Age >= 60 years", role = "PRIMARY",
    mask = !is.na(dat$Age) & dat$Age >= 60
  ),
  age65 = list(
    label = "Age >= 65 years", role = "SENSITIVITY",
    mask = !is.na(dat$Age) & dat$Age >= 65
  )
)

# Recompute survey-weighted means and SDs within each age domain. PhenoAge
# Acceleration remains a survey-weighted residual from PhenoAge ~ chronological age.
scaling_rows <- list()
for (domain_name in names(domains)) {
  domain <- domains[[domain_name]]
  domain_design <- design_raw[domain$mask, ]
  for (short_name in names(RAW_SCORES)) {
    raw_name <- RAW_SCORES[[short_name]]
    mean_value <- as.numeric(coef(svymean(as.formula(paste0("~", raw_name)), domain_design, na.rm = TRUE)))
    sd_value <- sqrt(as.numeric(svyvar(as.formula(paste0("~", raw_name)), domain_design, na.rm = TRUE)))
    if (!is.finite(sd_value) || sd_value <= 0) stop("Invalid survey SD for ", raw_name)
    z_name <- paste0(short_name, "_z_", domain_name)
    dat[[z_name]] <- (dat[[raw_name]] - mean_value) / sd_value
    scaling_rows[[length(scaling_rows) + 1L]] <- data.table(
      scenario = SCENARIO,
      domain = domain_name,
      domain_label = domain$label,
      analysis_role = domain$role,
      standardized_variable = z_name,
      source_variable = raw_name,
      survey_weighted_mean = mean_value,
      survey_weighted_sd = sd_value,
      age_regression_intercept = NA_real_,
      age_regression_slope = NA_real_,
      standardization_definition = "survey-weighted mean and SD within age domain"
    )
  }
  pheno_fit <- svyglm(PhenoAge ~ Age, design = domain_design)
  pheno_coef <- coef(pheno_fit)
  residual_name <- paste0("PhenoAge_residual_", domain_name)
  dat[[residual_name]] <- dat$PhenoAge -
    (pheno_coef[["(Intercept)"]] + pheno_coef[["Age"]] * dat$Age)
  design_residual <- svydesign(
    ids = ~survey_psu, strata = ~survey_strata,
    weights = ~pooled_mec_weight, data = dat, nest = TRUE
  )[domain$mask, ]
  residual_mean <- as.numeric(coef(svymean(as.formula(paste0("~", residual_name)), design_residual, na.rm = TRUE)))
  residual_sd <- sqrt(as.numeric(svyvar(as.formula(paste0("~", residual_name)), design_residual, na.rm = TRUE)))
  if (!is.finite(residual_sd) || residual_sd <= 0) stop("Invalid PhenoAge residual SD")
  pheno_z <- paste0("PHENO_z_", domain_name)
  dat[[pheno_z]] <- (dat[[residual_name]] - residual_mean) / residual_sd
  scaling_rows[[length(scaling_rows) + 1L]] <- data.table(
    scenario = SCENARIO,
    domain = domain_name,
    domain_label = domain$label,
    analysis_role = domain$role,
    standardized_variable = pheno_z,
    source_variable = "PhenoAge",
    survey_weighted_mean = residual_mean,
    survey_weighted_sd = residual_sd,
    age_regression_intercept = pheno_coef[["(Intercept)"]],
    age_regression_slope = pheno_coef[["Age"]],
    standardization_definition = paste(
      "survey-weighted residual from PhenoAge ~ chronological age,",
      "then survey-weighted SD within age domain"
    )
  )
}
scaling_dt <- rbindlist(scaling_rows)
fwrite(scaling_dt, file.path(out_dir, "overall_nhanes_domain_scaling.csv"), bom = TRUE)

design_all <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~pooled_mec_weight,
  data = dat,
  nest = TRUE
)

population_masks <- function(domain_mask) list(
  Overall = domain_mask,
  Cancer = domain_mask & !is.na(dat$Cancer_Diagnosed) & dat$Cancer_Diagnosed == 1,
  Noncancer = domain_mask & !is.na(dat$Cancer_Diagnosed) & dat$Cancer_Diagnosed == 0
)

adjustment_terms <- function(population, level = c("demographics", "fully_adjusted")) {
  level <- match.arg(level)
  terms <- if (level == "demographics") {
    c("Age", "Sex_F_design", "Race_F_design")
  } else {
    c(
      "Age", "Sex_F_design", "Race_F_design", "Education_F_design", "PIR", "BMI",
      "Smoking_F_design", "Alcohol_F_design", "Diabetes", "Hypertension", "CVD"
    )
  }
  if (population == "Overall") terms <- c(terms, "Cancer_F_design")
  terms
}

safe_term_p <- function(fit, terms) {
  tryCatch(
    as.numeric(regTermTest(fit, as.formula(paste("~", paste(terms, collapse = " + "))))$p),
    error = function(e) NA_real_
  )
}

fit_cox <- function(domain_name, population, population_mask, event_var,
                    predictor_terms, adjustment_level, model_type, result_role) {
  adjustments <- adjustment_terms(population, adjustment_level)
  required_vars <- unique(c(
    "Follow_Up_Years", event_var, "pooled_mec_weight",
    predictor_terms, adjustments
  ))
  complete <- population_mask & complete.cases(dat[, ..required_vars]) &
    dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
  if (sum(dat[[event_var]][complete] == 1) < 20L) stop("Insufficient events")
  sub_design <- design_all[complete, ]
  formula <- as.formula(paste(
    "Surv(Follow_Up_Years,", event_var, ") ~",
    paste(c(predictor_terms, adjustments), collapse = " + ")
  ))
  fit <- svycoxph(formula, design = sub_design)
  b <- coef(fit)
  conf <- confint(fit)
  coef_table <- summary(fit)$coefficients
  rows <- rbindlist(lapply(predictor_terms, function(term) {
    short <- sub("_z_(age60|age65)$", "", term)
    data.table(
      scenario = SCENARIO,
      domain = domain_name,
      domain_label = domains[[domain_name]]$label,
      analysis_role = result_role,
      population = population,
      mortality_outcome = event_var,
      model_type = model_type,
      adjustment = adjustment_level,
      term = term,
      component = short,
      component_label = LABELS[[short]],
      n = sum(complete),
      events = sum(dat[[event_var]][complete] == 1),
      hazard_ratio_per_domain_weighted_sd = exp(b[term]),
      ci_lower = exp(conf[term, 1]),
      ci_upper = exp(conf[term, 2]),
      p_value = coef_table[term, "Pr(>|z|)"],
      design_degrees_of_freedom = degf(sub_design)
    )
  }))
  list(fit = fit, rows = rows, complete = complete,
       joint_p = if (length(predictor_terms) > 1L) safe_term_p(fit, predictor_terms) else NA_real_)
}

# Overall-score all-cause association.
allcause_rows <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  score_term <- paste0("OVERALL_z_", domain_name)
  for (pop_name in names(pops)) {
    for (adjustment_level in c("demographics", "fully_adjusted")) {
      role <- ifelse(domain_name == "age60" && adjustment_level == "fully_adjusted",
                     "PRIMARY", "SENSITIVITY")
      fitted <- fit_cox(
        domain_name, pop_name, pops[[pop_name]], "Death_AllCause", score_term,
        adjustment_level, "overall_score_separate", role
      )
      allcause_rows[[length(allcause_rows) + 1L]] <- fitted$rows
    }
  }
}
allcause_dt <- rbindlist(allcause_rows)
fwrite(allcause_dt, file.path(out_dir, "overall_nhanes_allcause.csv"), bom = TRUE)

# Cancer-history interaction for the overall score.
interaction_rows <- list()
for (domain_name in names(domains)) {
  score_term <- paste0("OVERALL_z_", domain_name)
  domain_mask <- domains[[domain_name]]$mask & !is.na(dat$Cancer_Diagnosed)
  adjustments <- adjustment_terms("Overall", "fully_adjusted")
  required_vars <- unique(c("Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", score_term, adjustments))
  complete <- domain_mask & complete.cases(dat[, ..required_vars]) &
    dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
  sub_design <- design_all[complete, ]
  formula <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    score_term, "* Cancer_F_design +",
    paste(setdiff(adjustments, "Cancer_F_design"), collapse = " + ")
  ))
  fit <- svycoxph(formula, design = sub_design)
  b <- coef(fit)
  v <- vcov(fit)
  int_term <- names(b)[grepl(paste0("^", score_term, ":Cancer_F_design"), names(b))]
  if (length(int_term) != 1L) stop("Overall interaction term not found")
  beta_non <- b[[score_term]]
  beta_int <- b[[int_term]]
  beta_can <- beta_non + beta_int
  se_non <- sqrt(v[score_term, score_term])
  se_int <- sqrt(v[int_term, int_term])
  se_can <- sqrt(v[score_term, score_term] + v[int_term, int_term] + 2 * v[score_term, int_term])
  interaction_rows[[length(interaction_rows) + 1L]] <- data.table(
    scenario = SCENARIO,
    domain = domain_name,
    analysis_role = ifelse(domain_name == "age60", "PRIMARY", "SENSITIVITY"),
    component = "OVERALL",
    component_label = LABELS[["OVERALL"]],
    n = sum(complete),
    events = sum(dat$Death_AllCause[complete] == 1),
    hr_non_cancer_history = exp(beta_non),
    non_cancer_ci_lower = exp(beta_non - 1.96 * se_non),
    non_cancer_ci_upper = exp(beta_non + 1.96 * se_non),
    hr_cancer_history = exp(beta_can),
    cancer_ci_lower = exp(beta_can - 1.96 * se_can),
    cancer_ci_upper = exp(beta_can + 1.96 * se_can),
    interaction_ratio_of_hrs = exp(beta_int),
    interaction_ci_lower = exp(beta_int - 1.96 * se_int),
    interaction_ci_upper = exp(beta_int + 1.96 * se_int),
    interaction_p = summary(fit)$coefficients[int_term, "Pr(>|z|)"],
    design_degrees_of_freedom = degf(sub_design)
  )
}
interaction_dt <- rbindlist(interaction_rows)
fwrite(interaction_dt, file.path(out_dir, "overall_nhanes_cancer_history_interaction.csv"), bom = TRUE)

# Formal common-complete-case comparison of overall score and PhenoAge.
pheno_tests <- list()
pheno_effects <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  overall <- paste0("OVERALL_z_", domain_name)
  pheno <- paste0("PHENO_z_", domain_name)
  for (pop_name in names(pops)) {
    adjustments <- adjustment_terms(pop_name, "fully_adjusted")
    required_vars <- unique(c("Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", overall, pheno, adjustments))
    complete <- pops[[pop_name]] & complete.cases(dat[, ..required_vars]) &
      dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
    sub_design <- design_all[complete, ]
    make_formula <- function(predictors) as.formula(paste(
      "Surv(Follow_Up_Years, Death_AllCause) ~",
      paste(c(predictors, adjustments), collapse = " + ")
    ))
    fit_m0 <- svycoxph(make_formula(character()), design = sub_design)
    fit_m1 <- svycoxph(make_formula(pheno), design = sub_design)
    fit_m2 <- svycoxph(make_formula(overall), design = sub_design)
    fit_m3 <- svycoxph(make_formula(c(pheno, overall)), design = sub_design)
    tests <- list(
      list("M0_to_M1", "PhenoAge beyond clinical base", fit_m1, pheno),
      list("M0_to_M2", "Overall score beyond clinical base", fit_m2, overall),
      list("M1_to_M3", "Overall score beyond PhenoAge and clinical base", fit_m3, overall),
      list("M2_to_M3", "PhenoAge beyond overall score and clinical base", fit_m3, pheno)
    )
    for (test in tests) {
      pheno_tests[[length(pheno_tests) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        analysis_role = ifelse(domain_name == "age60", "PRIMARY", "SENSITIVITY"),
        population = pop_name,
        comparison = test[[1]],
        interpretation = test[[2]],
        n = sum(complete),
        events = sum(dat$Death_AllCause[complete] == 1),
        design_adjusted_block_wald_p = safe_term_p(test[[3]], test[[4]]),
        common_complete_case_cohort = TRUE
      )
    }
    b <- coef(fit_m3)
    conf <- confint(fit_m3)
    coef_table <- summary(fit_m3)$coefficients
    for (term in c(pheno, overall)) {
      short <- sub("_z_(age60|age65)$", "", term)
      pheno_effects[[length(pheno_effects) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        analysis_role = ifelse(domain_name == "age60", "PRIMARY", "SENSITIVITY"),
        population = pop_name,
        model = "M3_clinical_plus_phenoage_plus_overall",
        term = term,
        component = short,
        n = sum(complete),
        events = sum(dat$Death_AllCause[complete] == 1),
        hazard_ratio_per_domain_weighted_sd = exp(b[term]),
        ci_lower = exp(conf[term, 1]),
        ci_upper = exp(conf[term, 2]),
        p_value = coef_table[term, "Pr(>|z|)"]
      )
    }
  }
}
pheno_tests_dt <- rbindlist(pheno_tests)
pheno_effects_dt <- rbindlist(pheno_effects)
fwrite(pheno_tests_dt, file.path(out_dir, "overall_nhanes_phenoage_incremental_tests.csv"), bom = TRUE)
fwrite(pheno_effects_dt, file.path(out_dir, "overall_nhanes_phenoage_incremental_effects.csv"), bom = TRUE)

# Determine how much information is retained or lost by collapsing the three
# component scores into one scalar overall score.
component_tests <- list()
component_effects <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  overall <- paste0("OVERALL_z_", domain_name)
  components <- paste0(c("NM", "TB", "TC"), "_z_", domain_name)
  for (pop_name in names(pops)) {
    adjustments <- adjustment_terms(pop_name, "fully_adjusted")
    required_vars <- unique(c("Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", overall, components, adjustments))
    complete <- pops[[pop_name]] & complete.cases(dat[, ..required_vars]) &
      dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
    sub_design <- design_all[complete, ]
    make_formula <- function(predictors) as.formula(paste(
      "Surv(Follow_Up_Years, Death_AllCause) ~",
      paste(c(predictors, adjustments), collapse = " + ")
    ))
    fit_base <- svycoxph(make_formula(character()), design = sub_design)
    fit_overall <- svycoxph(make_formula(overall), design = sub_design)
    fit_components <- svycoxph(make_formula(components), design = sub_design)
    fit_both <- svycoxph(make_formula(c(overall, components)), design = sub_design)
    comparisons <- list(
      list("BASE_to_OVERALL", "Overall score beyond clinical base", fit_overall, overall),
      list("BASE_to_COMPONENTS", "Three components beyond clinical base", fit_components, components),
      list("COMPONENTS_to_BOTH", "Nonlinear overall summary beyond three linear component terms", fit_both, overall),
      list("OVERALL_to_BOTH", "Three component terms beyond overall summary", fit_both, components)
    )
    for (comparison in comparisons) {
      component_tests[[length(component_tests) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        analysis_role = ifelse(domain_name == "age60", "PRIMARY", "SENSITIVITY"),
        population = pop_name,
        comparison = comparison[[1]],
        interpretation = comparison[[2]],
        n = sum(complete),
        events = sum(dat$Death_AllCause[complete] == 1),
        design_adjusted_block_wald_p = safe_term_p(comparison[[3]], comparison[[4]]),
        common_complete_case_cohort = TRUE
      )
    }
    b <- coef(fit_both)
    conf <- confint(fit_both)
    coef_table <- summary(fit_both)$coefficients
    for (term in c(overall, components)) {
      short <- sub("_z_(age60|age65)$", "", term)
      component_effects[[length(component_effects) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        analysis_role = ifelse(domain_name == "age60", "PRIMARY", "SENSITIVITY"),
        population = pop_name,
        model = "clinical_plus_overall_plus_three_components",
        term = term,
        component = short,
        n = sum(complete),
        events = sum(dat$Death_AllCause[complete] == 1),
        hazard_ratio_per_domain_weighted_sd = exp(b[term]),
        ci_lower = exp(conf[term, 1]),
        ci_upper = exp(conf[term, 2]),
        p_value = coef_table[term, "Pr(>|z|)"]
      )
    }
  }
}
component_tests_dt <- rbindlist(component_tests)
component_effects_dt <- rbindlist(component_effects)
fwrite(component_tests_dt, file.path(out_dir, "overall_nhanes_component_incremental_tests.csv"), bom = TRUE)
fwrite(component_effects_dt, file.path(out_dir, "overall_nhanes_component_incremental_effects.csv"), bom = TRUE)

# Cause-specific mortality for the overall score; other deaths are censored.
cause_definitions <- list(
  Death_Cancer = list(role = "KEY_SECONDARY", label = "Cancer mortality"),
  Death_CVD = list(role = "KEY_SECONDARY", label = "Cardiovascular mortality; UCOD 1 + 5"),
  Death_CLRD = list(role = "EXPLORATORY", label = "Chronic lower respiratory disease mortality"),
  Death_Alzheimer = list(role = "EXPLORATORY", label = "Alzheimer disease mortality"),
  Death_Diabetes_UCOD = list(role = "EXPLORATORY", label = "Diabetes mortality"),
  Death_Kidney = list(role = "EXPLORATORY", label = "Kidney disease mortality")
)
cause_rows <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  overall <- paste0("OVERALL_z_", domain_name)
  for (event_var in names(cause_definitions)) {
    info <- cause_definitions[[event_var]]
    if (domain_name == "age65" && info$role == "EXPLORATORY") next
    eligible_pops <- if (info$role == "EXPLORATORY") "Overall" else names(pops)
    for (pop_name in eligible_pops) {
      fitted <- fit_cox(
        domain_name, pop_name, pops[[pop_name]], event_var, overall,
        "fully_adjusted", "cause_specific_overall_score", info$role
      )
      fitted$rows[, `:=`(
        cause_label = info$label,
        competing_deaths_handling = "censored_at_death_time"
      )]
      cause_rows[[length(cause_rows) + 1L]] <- fitted$rows
    }
  }
}
cause_dt <- rbindlist(cause_rows, fill = TRUE)
cause_dt[, p_fdr_bh := p.adjust(p_value, method = "BH"), by = .(domain, population)]
fwrite(cause_dt, file.path(out_dir, "overall_nhanes_cause_specific.csv"), bom = TRUE)

coverage_dt <- rbindlist(lapply(names(domains), function(domain_name) {
  mask <- domains[[domain_name]]$mask
  data.table(
    scenario = SCENARIO,
    domain = domain_name,
    source_n = sum(mask),
    overall_score_nonmissing_n = sum(mask & !is.na(dat$Overall_expected_burden)),
    overall_score_missing_n = sum(mask & is.na(dat$Overall_expected_burden)),
    missing_reason = "at least one frozen component had zero observed selected inputs"
  )
}))
fwrite(coverage_dt, file.path(out_dir, "overall_nhanes_coverage.csv"), bom = TRUE)

design_audit <- data.table(
  item = c(
    "source_scope", "full_design_before_age_domain", "pooled_weight_formula",
    "primary_domain", "age65_sensitivity", "overall_training_source",
    "overall_inputs", "phenoage_used_as_overall_input", "nhanes_outcomes_used_for_training",
    "overall_observation_gate", "primary_outcome", "key_secondary_outcomes",
    "cvd_definition", "phenoage_acceleration_definition", "component_reporting"
  ),
  value = c(
    "1999-2016; nine cycles; 44,772 participants with valid follow-up",
    "YES", CIPDS_WEIGHT_DESCRIPTION, "Age >= 60 years", "Age >= 65 years",
    "hospital calendar-training strict OOF component predictions",
    "NM, TB and TC calibrated probabilities only", "FALSE", "FALSE",
    "all three frozen components require at least one observed selected input",
    "all-cause mortality", "cancer mortality; cardiovascular mortality",
    "UCOD 1 + 5", "survey-weighted residual of PhenoAge ~ chronological age within domain",
    "Overall plus NM plus TB plus TC"
  )
)
fwrite(design_audit, file.path(out_dir, "overall_nhanes_design_audit.csv"), bom = TRUE)

cat("\nNHANES overall-score analysis complete.\n")
cat("Coverage:\n")
print(coverage_dt)
cat("\nPrimary >=60 fully adjusted overall-score results:\n")
print(allcause_dt[domain == "age60" & analysis_role == "PRIMARY"])
cat("\nPrimary >=60 PhenoAge comparisons:\n")
print(pheno_tests_dt[domain == "age60"])
cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
