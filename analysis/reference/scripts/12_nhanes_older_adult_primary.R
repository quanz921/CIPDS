source(file.path(Sys.getenv("CIPDS_PACKAGE_DIR"), "scripts", "weight_contract.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survey)
  library(survival)
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
log_file <- file.path(log_dir, "12_nhanes_older_adult_primary.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

SCENARIO <- "A_PRIMARY_FROZEN_PANCANCER_MODEL"
N_CYCLES <- 9L
COMPONENT_RAW <- c(
  NM = "NM_component",
  TB = "TB_component",
  TC = "TC_component"
)
COMPONENT_LABELS <- c(
  NM = "Nutritional-Metabolic laboratory pattern",
  TB = "Tumor-burden-related laboratory pattern",
  TC = "Treatment-complication-related laboratory pattern",
  PHENO = "PhenoAge Acceleration"
)

load(file.path(package_dir, "inputs", "nhanes_phase7_batch1.RData"))
if (!exists("nhanes_valid")) stop("nhanes_valid missing")
dat <- as.data.table(nhanes_valid)
if (nrow(dat) != 44772L || uniqueN(dat$SEQN) != 44772L) stop("NHANES grain drift")
if (uniqueN(dat$CYCLE) != N_CYCLES) stop("NHANES cycle count drift")

scores <- fread(file.path(out_dir, "nested_a_primary_nhanes_component_scores.csv"))
if (nrow(scores) != 44772L || uniqueN(scores$SEQN) != 44772L) stop("Frozen score grain drift")
score_idx <- match(dat$SEQN, scores$SEQN)
if (anyNA(score_idx)) stop("Frozen score join failure")
for (v in c("NM_component", "TB_component", "TC_component")) {
  dat[[v]] <- scores[[v]][score_idx]
  if (anyNA(dat[[v]]) || any(dat[[v]] < 0 | dat[[v]] > 1)) stop("Invalid frozen score: ", v)
}

# A mathematical CatBoost prediction with every selected input missing is not treated as
# observed physiology. Preserve the frozen model but gate each component to at least one
# actually observed selected laboratory value in that participant.
mapped <- readRDS(file.path(out_dir, "nhanes_expanded_candidate_matrix.rds"))
mapped_idx <- match(dat$SEQN, mapped$SEQN)
if (anyNA(mapped_idx)) stop("Mapped feature join failure")
features_dt <- fread(file.path(out_dir, "nested_a_primary_final_features.csv"))
final_features <- unique(features_dt$feature)
if (any(!final_features %in% names(mapped))) stop("Frozen final feature missing from NHANES matrix")
outcome_to_component <- c(
  Outcome_NutriMetab = "NM",
  Outcome_TumorBurden = "TB",
  Outcome_TreatComp = "TC"
)
for (outcome_name in names(outcome_to_component)) {
  component <- outcome_to_component[[outcome_name]]
  selected <- features_dt[outcome == outcome_name, feature]
  x <- as.matrix(mapped[mapped_idx, ..selected])
  observed_count <- rowSums(!is.na(x))
  dat[[paste0(component, "_input_observed_n")]] <- observed_count
  component_col <- paste0(component, "_component")
  dat[observed_count == 0, (component_col) := NA_real_]
}

required <- c(
  "CYCLE", "WTMEC2YR", "SDMVSTRA", "SDMVPSU", "Follow_Up_Years", "Dead",
  "UCOD_LEADING", "Age", "Sex", "RIDRETH1", "Cancer_Diagnosed", "Education",
  "PIR", "BMI", "Smoking_F", "Alcohol_F", "Diabetes", "Hypertension", "CVD",
  unname(COMPONENT_RAW), "PhenoAge"
)
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required fields: ", paste(missing_required, collapse = ", "))

# Mortality definitions are locked to the NHANES linked-mortality recode.
# Critical correction: code 3 is chronic lower respiratory disease; cerebrovascular disease is code 5.
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
if (any(dat$Dead == 1 & is.na(dat$UCOD_LEADING))) stop("Deceased participant missing UCOD_LEADING")

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

# The full 44,772-person survey design is deliberately created before any age-domain restriction.
design_raw <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~pooled_mec_weight,
  data = dat,
  nest = TRUE
)

domains <- list(
  age60 = list(
    label = "Age >= 60 years",
    role = "PRIMARY",
    threshold = 60,
    mask = !is.na(dat$Age) & dat$Age >= 60
  ),
  age65 = list(
    label = "Age >= 65 years",
    role = "SENSITIVITY",
    threshold = 65,
    mask = !is.na(dat$Age) & dat$Age >= 65
  )
)

# Survey-weighted standardization is recomputed inside each older-adult domain.
scaling_rows <- list()
for (domain_name in names(domains)) {
  domain <- domains[[domain_name]]
  domain_design <- design_raw[domain$mask, ]
  for (short_name in names(COMPONENT_RAW)) {
    raw_name <- COMPONENT_RAW[[short_name]]
    mean_value <- as.numeric(coef(svymean(as.formula(paste0("~", raw_name)), domain_design, na.rm = TRUE)))
    sd_value <- sqrt(as.numeric(svyvar(as.formula(paste0("~", raw_name)), domain_design, na.rm = TRUE)))
    if (!is.finite(sd_value) || sd_value <= 0) stop("Invalid survey SD for ", raw_name, " in ", domain_name)
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
  # PhenoAge Acceleration is the survey-weighted residual of PhenoAge on
  # chronological age within the target age domain, not the simple PhenoAge-Age difference.
  pheno_fit <- svyglm(PhenoAge ~ Age, design = domain_design)
  pheno_coef <- coef(pheno_fit)
  pheno_residual_name <- paste0("PhenoAge_residual_", domain_name)
  dat[[pheno_residual_name]] <- dat$PhenoAge -
    (pheno_coef[["(Intercept)"]] + pheno_coef[["Age"]] * dat$Age)
  design_with_residual <- svydesign(
    ids = ~survey_psu,
    strata = ~survey_strata,
    weights = ~pooled_mec_weight,
    data = dat,
    nest = TRUE
  )
  residual_design <- design_with_residual[domain$mask, ]
  residual_mean <- as.numeric(coef(svymean(
    as.formula(paste0("~", pheno_residual_name)), residual_design, na.rm = TRUE
  )))
  residual_sd <- sqrt(as.numeric(svyvar(
    as.formula(paste0("~", pheno_residual_name)), residual_design, na.rm = TRUE
  )))
  if (!is.finite(residual_sd) || residual_sd <= 0) stop("Invalid PhenoAge residual SD")
  pheno_z_name <- paste0("PHENO_z_", domain_name)
  dat[[pheno_z_name]] <- (dat[[pheno_residual_name]] - residual_mean) / residual_sd
  scaling_rows[[length(scaling_rows) + 1L]] <- data.table(
    scenario = SCENARIO,
    domain = domain_name,
    domain_label = domain$label,
    analysis_role = domain$role,
    standardized_variable = pheno_z_name,
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
fwrite(scaling_dt, file.path(out_dir, "older_nhanes_domain_scaling.csv"), bom = TRUE)

# Rebuild the same full survey design so domain-specific standardized variables are available.
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

fit_svycox_terms <- function(domain_name, population, population_mask, event_var,
                             predictor_terms, adjustment_level, model_type,
                             result_role) {
  adjustments <- adjustment_terms(population, adjustment_level)
  required_vars <- unique(c(
    "Follow_Up_Years", event_var, "pooled_mec_weight", "survey_strata", "survey_psu",
    predictor_terms, adjustments
  ))
  complete <- population_mask & complete.cases(dat[, ..required_vars]) &
    dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
  if (sum(dat[[event_var]][complete] == 1) < 20L) {
    stop("Insufficient events for ", domain_name, "/", population, "/", event_var)
  }
  sub_design <- design_all[complete, ]
  formula <- as.formula(paste(
    "Surv(Follow_Up_Years,", event_var, ") ~",
    paste(c(predictor_terms, adjustments), collapse = " + ")
  ))
  fit <- svycoxph(formula, design = sub_design)
  coef_table <- summary(fit)$coefficients
  conf <- confint(fit)
  rows <- lapply(predictor_terms, function(term) data.table(
    scenario = SCENARIO,
    domain = domain_name,
    domain_label = domains[[domain_name]]$label,
    analysis_role = result_role,
    population = population,
    mortality_outcome = event_var,
    model_type = model_type,
    adjustment = adjustment_level,
    term = term,
    component = sub("_z_(age60|age65)$", "", term),
    component_label = COMPONENT_LABELS[[sub("_z_(age60|age65)$", "", term)]],
    n = sum(complete),
    events = sum(dat[[event_var]][complete] == 1),
    hazard_ratio_per_domain_weighted_sd = exp(coef(fit)[term]),
    ci_lower = exp(conf[term, 1]),
    ci_upper = exp(conf[term, 2]),
    p_value = coef_table[term, "Pr(>|z|)"],
    design_degrees_of_freedom = degf(sub_design)
  ))
  list(
    fit = fit,
    rows = rbindlist(rows),
    complete = complete,
    joint_p = if (length(predictor_terms) > 1L) safe_term_p(fit, predictor_terms) else NA_real_
  )
}

# Cohort/event audit is done before complete-case modeling.
event_vars <- c(
  "Death_AllCause", "Death_Cancer", "Death_CVD", "Death_CLRD",
  "Death_Alzheimer", "Death_Diabetes_UCOD", "Death_Kidney"
)
event_audit <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  for (pop_name in names(pops)) {
    mask <- pops[[pop_name]]
    row <- data.table(
      scenario = SCENARIO,
      domain = domain_name,
      domain_label = domains[[domain_name]]$label,
      analysis_role = domains[[domain_name]]$role,
      population = pop_name,
      n = sum(mask),
      cycles = uniqueN(dat$CYCLE[mask]),
      median_follow_up_years = median(dat$Follow_Up_Years[mask], na.rm = TRUE)
    )
    for (event_var in event_vars) row[[event_var]] <- sum(dat[[event_var]][mask] == 1, na.rm = TRUE)
    event_audit[[length(event_audit) + 1L]] <- row
  }
}
event_audit_dt <- rbindlist(event_audit, fill = TRUE)
fwrite(event_audit_dt, file.path(out_dir, "older_nhanes_cohort_event_audit.csv"), bom = TRUE)

# Audit the exact frozen A-primary inputs in each older-adult domain and cycle.
coverage_rows <- list()
for (domain_name in names(domains)) {
  domain_mask <- domains[[domain_name]]$mask
  cycle_levels <- c("ALL", sort(unique(dat$CYCLE[domain_mask])))
  for (cycle_value in cycle_levels) {
    mask <- domain_mask
    if (cycle_value != "ALL") mask <- mask & dat$CYCLE == cycle_value
    for (feature in final_features) {
      values <- mapped[[feature]][mapped_idx][mask]
      coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        cycle = cycle_value,
        feature = feature,
        n = length(values),
        nonmissing_n = sum(!is.na(values)),
        missing_pct = 100 * mean(is.na(values))
      )
    }
  }
}
coverage_dt <- rbindlist(coverage_rows)
fwrite(coverage_dt, file.path(out_dir, "older_nhanes_final_feature_coverage.csv"), bom = TRUE)

component_coverage_rows <- list()
for (domain_name in names(domains)) {
  mask <- domains[[domain_name]]$mask
  for (outcome_name in unique(features_dt$outcome)) {
    selected <- features_dt[outcome == outcome_name, feature]
    x <- as.matrix(mapped[mapped_idx[mask], ..selected])
    observed_count <- rowSums(!is.na(x))
    component_coverage_rows[[length(component_coverage_rows) + 1L]] <- data.table(
      scenario = SCENARIO,
      domain = domain_name,
      outcome = outcome_name,
      n_features = length(selected),
      n = length(observed_count),
      median_observed_features = median(observed_count),
      p10_observed_features = as.numeric(quantile(observed_count, 0.10)),
      min_observed_features = min(observed_count),
      zero_observed_n = sum(observed_count == 0)
    )
  }
}
component_coverage_dt <- rbindlist(component_coverage_rows)
fwrite(component_coverage_dt, file.path(out_dir, "older_nhanes_component_input_coverage.csv"), bom = TRUE)

# Main all-cause analyses: components and PhenoAge separately, then components jointly.
allcause_separate <- list()
allcause_joint <- list()
allcause_joint_tests <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  component_terms <- paste0(c("NM", "TB", "TC"), "_z_", domain_name)
  pheno_term <- paste0("PHENO_z_", domain_name)
  for (pop_name in names(pops)) {
    for (adjustment_level in c("demographics", "fully_adjusted")) {
      for (term in c(component_terms, pheno_term)) {
        result_role <- if (
          domain_name == "age60" && adjustment_level == "fully_adjusted"
        ) "PRIMARY" else "SENSITIVITY"
        fitted <- fit_svycox_terms(
          domain_name, pop_name, pops[[pop_name]], "Death_AllCause", term,
          adjustment_level,
          ifelse(grepl("PHENO", term), "phenoage_benchmark_separate", "component_separate"),
          result_role
        )
        allcause_separate[[length(allcause_separate) + 1L]] <- fitted$rows
      }
      result_role <- if (
        domain_name == "age60" && adjustment_level == "fully_adjusted"
      ) "PRIMARY" else "SENSITIVITY"
      joint <- fit_svycox_terms(
        domain_name, pop_name, pops[[pop_name]], "Death_AllCause", component_terms,
        adjustment_level, "three_components_joint_no_composite", result_role
      )
      allcause_joint[[length(allcause_joint) + 1L]] <- joint$rows
      allcause_joint_tests[[length(allcause_joint_tests) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        analysis_role = result_role,
        population = pop_name,
        adjustment = adjustment_level,
        n = sum(joint$complete),
        events = sum(dat$Death_AllCause[joint$complete] == 1),
        component_block_wald_p = joint$joint_p
      )
    }
  }
}
allcause_separate_dt <- rbindlist(allcause_separate)
allcause_joint_dt <- rbindlist(allcause_joint)
allcause_joint_tests_dt <- rbindlist(allcause_joint_tests)
fwrite(allcause_separate_dt, file.path(out_dir, "older_nhanes_allcause_component_separate.csv"), bom = TRUE)
fwrite(allcause_joint_dt, file.path(out_dir, "older_nhanes_allcause_components_joint.csv"), bom = TRUE)
fwrite(allcause_joint_tests_dt, file.path(out_dir, "older_nhanes_allcause_joint_wald.csv"), bom = TRUE)

# Cancer-history interactions: the interaction ratio compares component HRs in cancer vs noncancer.
interaction_rows <- list()
for (domain_name in names(domains)) {
  domain_mask <- domains[[domain_name]]$mask & !is.na(dat$Cancer_Diagnosed)
  for (component in c("NM", "TB", "TC")) {
    score_term <- paste0(component, "_z_", domain_name)
    adjustments <- adjustment_terms("Overall", "fully_adjusted")
    required_vars <- unique(c(
      "Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", score_term, adjustments
    ))
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
    interaction_term <- names(b)[grepl(paste0("^", score_term, ":Cancer_F_design"), names(b))]
    if (length(interaction_term) != 1L) stop("Interaction term not found for ", score_term)
    beta_non <- b[[score_term]]
    beta_int <- b[[interaction_term]]
    beta_can <- beta_non + beta_int
    se_non <- sqrt(v[score_term, score_term])
    se_int <- sqrt(v[interaction_term, interaction_term])
    se_can <- sqrt(
      v[score_term, score_term] + v[interaction_term, interaction_term] +
        2 * v[score_term, interaction_term]
    )
    interaction_rows[[length(interaction_rows) + 1L]] <- data.table(
      scenario = SCENARIO,
      domain = domain_name,
      analysis_role = ifelse(domain_name == "age60", "PRIMARY", "SENSITIVITY"),
      component = component,
      component_label = COMPONENT_LABELS[[component]],
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
      interaction_p = summary(fit)$coefficients[interaction_term, "Pr(>|z|)"],
      design_degrees_of_freedom = degf(sub_design)
    )
  }
}
interaction_dt <- rbindlist(interaction_rows)
interaction_dt[, interaction_fdr_bh := p.adjust(interaction_p, method = "BH"), by = domain]
fwrite(interaction_dt, file.path(out_dir, "older_nhanes_cancer_history_interactions.csv"), bom = TRUE)

# Formal incremental comparison with a common complete-case cohort.
incremental_tests <- list()
incremental_effects <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  components <- paste0(c("NM", "TB", "TC"), "_z_", domain_name)
  pheno <- paste0("PHENO_z_", domain_name)
  for (pop_name in names(pops)) {
    adjustments <- adjustment_terms(pop_name, "fully_adjusted")
    required_vars <- unique(c(
      "Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", components, pheno, adjustments
    ))
    complete <- pops[[pop_name]] & complete.cases(dat[, ..required_vars]) &
      dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
    sub_design <- design_all[complete, ]
    make_formula <- function(predictors) as.formula(paste(
      "Surv(Follow_Up_Years, Death_AllCause) ~",
      paste(c(predictors, adjustments), collapse = " + ")
    ))
    fit_m0 <- svycoxph(make_formula(character()), design = sub_design)
    fit_m1 <- svycoxph(make_formula(pheno), design = sub_design)
    fit_m2 <- svycoxph(make_formula(components), design = sub_design)
    fit_m3 <- svycoxph(make_formula(c(pheno, components)), design = sub_design)
    tests <- list(
      list("M0_to_M1", "PhenoAge beyond clinical base", fit_m1, pheno),
      list("M0_to_M2", "Three components beyond clinical base", fit_m2, components),
      list("M1_to_M3", "Three components beyond PhenoAge and clinical base", fit_m3, components),
      list("M2_to_M3", "PhenoAge beyond three components and clinical base", fit_m3, pheno)
    )
    for (test in tests) {
      incremental_tests[[length(incremental_tests) + 1L]] <- data.table(
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
    for (term in c(pheno, components)) {
      incremental_effects[[length(incremental_effects) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        analysis_role = ifelse(domain_name == "age60", "PRIMARY", "SENSITIVITY"),
        population = pop_name,
        model = "M3_clinical_plus_phenoage_plus_three_components",
        term = term,
        component = sub("_z_(age60|age65)$", "", term),
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
incremental_tests_dt <- rbindlist(incremental_tests)
incremental_effects_dt <- rbindlist(incremental_effects)
fwrite(incremental_tests_dt, file.path(out_dir, "older_nhanes_phenoage_incremental_tests.csv"), bom = TRUE)
fwrite(incremental_effects_dt, file.path(out_dir, "older_nhanes_phenoage_incremental_effects.csv"), bom = TRUE)

# Key secondary cause-specific mortality analyses; competing deaths are censored.
cause_definitions <- list(
  Death_Cancer = list(role = "KEY_SECONDARY", label = "Cancer mortality"),
  Death_CVD = list(role = "KEY_SECONDARY", label = "Cardiovascular mortality; UCOD 1 + 5"),
  Death_CLRD = list(role = "EXPLORATORY", label = "Chronic lower respiratory disease mortality"),
  Death_Alzheimer = list(role = "EXPLORATORY", label = "Alzheimer disease mortality"),
  Death_Diabetes_UCOD = list(role = "EXPLORATORY", label = "Diabetes mortality"),
  Death_Kidney = list(role = "EXPLORATORY", label = "Kidney disease mortality")
)
cause_separate <- list()
cause_joint <- list()
cause_joint_tests <- list()
for (domain_name in names(domains)) {
  pops <- population_masks(domains[[domain_name]]$mask)
  component_terms <- paste0(c("NM", "TB", "TC"), "_z_", domain_name)
  for (event_var in names(cause_definitions)) {
    cause_info <- cause_definitions[[event_var]]
    # Exploratory causes are kept to the primary >=60 overall domain only.
    eligible_populations <- if (cause_info$role == "EXPLORATORY") "Overall" else names(pops)
    if (domain_name == "age65" && cause_info$role == "EXPLORATORY") next
    for (pop_name in eligible_populations) {
      for (term in component_terms) {
        fitted <- fit_svycox_terms(
          domain_name, pop_name, pops[[pop_name]], event_var, term,
          "fully_adjusted", "cause_specific_component_separate", cause_info$role
        )
        fitted$rows[, cause_label := cause_info$label]
        cause_separate[[length(cause_separate) + 1L]] <- fitted$rows
      }
      joint <- fit_svycox_terms(
        domain_name, pop_name, pops[[pop_name]], event_var, component_terms,
        "fully_adjusted", "cause_specific_three_components_joint", cause_info$role
      )
      joint$rows[, cause_label := cause_info$label]
      cause_joint[[length(cause_joint) + 1L]] <- joint$rows
      cause_joint_tests[[length(cause_joint_tests) + 1L]] <- data.table(
        scenario = SCENARIO,
        domain = domain_name,
        analysis_role = cause_info$role,
        population = pop_name,
        mortality_outcome = event_var,
        cause_label = cause_info$label,
        n = sum(joint$complete),
        events = sum(dat[[event_var]][joint$complete] == 1),
        component_block_wald_p = joint$joint_p,
        competing_deaths_handling = "censored_at_death_time"
      )
    }
  }
}
cause_separate_dt <- rbindlist(cause_separate, fill = TRUE)
cause_joint_dt <- rbindlist(cause_joint, fill = TRUE)
cause_joint_tests_dt <- rbindlist(cause_joint_tests, fill = TRUE)
cause_separate_dt[, p_fdr_bh := p.adjust(p_value, method = "BH"),
                  by = .(domain, population, mortality_outcome)]
cause_joint_dt[, p_fdr_bh := p.adjust(p_value, method = "BH"),
               by = .(domain, population, mortality_outcome)]
fwrite(cause_separate_dt, file.path(out_dir, "older_nhanes_cause_specific_component_separate.csv"), bom = TRUE)
fwrite(cause_joint_dt, file.path(out_dir, "older_nhanes_cause_specific_components_joint.csv"), bom = TRUE)
fwrite(cause_joint_tests_dt, file.path(out_dir, "older_nhanes_cause_specific_joint_wald.csv"), bom = TRUE)

design_audit <- data.table(
  scenario = SCENARIO,
  item = c(
    "raw_source_cycle_scope", "analytic_cycle_scope",
    "excluded_2017_2018_reason", "full_design_rows_before_domain_subset",
    "full_design_cycles", "weight_variable",
    "pooled_weight_formula", "cycle_specific_strata", "cycle_specific_psu",
    "full_design_degrees_of_freedom", "primary_domain", "age65_sensitivity_domain",
    "primary_outcome", "key_secondary_outcomes", "cvd_ucod_definition",
    "ucod_3_definition", "domain_standardization", "phenoage_acceleration_definition",
    "individual_component_observation_gate",
    "cancer_history_adjustment_overall",
    "weighted_composite_generated", "hospital_model_retrained", "all_age_nhanes_role"
  ),
  value = c(
    "1999-2018; 10 cycles; 49,774 raw participants",
    "1999-2016; 9 cycles; 44,772 participants with valid follow-up time",
    "5,002 participants and zero nonmissing Follow_Up_Years in locked raw source",
    nrow(dat), N_CYCLES, "WTMEC4YR for 1999-2002; WTMEC2YR for 2003-2016", CIPDS_WEIGHT_DESCRIPTION, uniqueN(dat$survey_strata),
    uniqueN(dat$survey_psu), degf(design_all), "Age >= 60 years", "Age >= 65 years",
    "All-cause mortality", "Cancer mortality; cardiovascular mortality",
    "1 (heart disease) + 5 (cerebrovascular disease)",
    "Chronic lower respiratory disease; never treated as cardiovascular",
    "survey-weighted mean and SD within each age domain",
    "survey-weighted residual from PhenoAge regressed on chronological age within each domain",
    "at least one observed selected laboratory input per component", "YES", "FALSE", "FALSE",
    "SENSITIVITY"
  )
)
fwrite(design_audit, file.path(out_dir, "older_nhanes_design_audit.csv"), bom = TRUE)

cat("Older-adult NHANES analysis complete.\n")
cat("Full survey design rows:", nrow(dat), "\n")
print(event_audit_dt)
cat("\nPrimary >=60 fully adjusted separate all-cause results:\n")
print(allcause_separate_dt[domain == "age60" & analysis_role == "PRIMARY"])
cat("\nPrimary >=60 PhenoAge incremental tests:\n")
print(incremental_tests_dt[domain == "age60"])
cat("\nNo weighted composite generated; hospital models remained frozen.\n")
