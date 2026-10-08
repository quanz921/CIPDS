source(file.path(Sys.getenv("CIPDS_PACKAGE_DIR"), "scripts", "weight_contract.R"))
suppressPackageStartupMessages({
  library(catboost)
  library(data.table)
  library(survey)
  library(survival)
})

options(survey.lonely.psu = "adjust")

package_dir <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
scenario <- Sys.getenv("CIPDS_CANDIDATE_SCENARIO", unset = "A_PRIMARY")
if (!scenario %in% c("A_PRIMARY", "AB_SENSITIVITY")) stop("Invalid scenario")
scenario_slug <- tolower(scenario)
out_dir <- file.path(package_dir, "outputs")
log_dir <- file.path(package_dir, "logs")
log_file <- file.path(log_dir, paste0("10_score_nhanes_nested_survey_", scenario_slug, ".log"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

load(file.path(package_dir, "inputs", "nhanes_phase7_batch1.RData"))
if (!exists("nhanes_valid")) stop("nhanes_valid missing")
dat <- as.data.table(nhanes_valid)
mapped <- readRDS(file.path(out_dir, "nhanes_expanded_candidate_matrix.rds"))
if (nrow(dat) != 44772L || uniqueN(dat$SEQN) != 44772L) stop("NHANES grain drift")
if (nrow(mapped) != 44772L || uniqueN(mapped$SEQN) != 44772L) stop("Mapped NHANES grain drift")
idx <- match(dat$SEQN, mapped$SEQN)
if (anyNA(idx)) stop("NHANES feature matrix join failure")
if (any(dat$CYCLE != mapped$CYCLE[idx])) stop("NHANES cycle mismatch after feature join")

features_file <- file.path(out_dir, paste0("nested_", scenario_slug, "_final_features.csv"))
model_registry_file <- file.path(
  out_dir, paste0("nested_", scenario_slug, "_calibration_and_model_registry.csv")
)
features_dt <- fread(features_file)
model_registry <- fread(model_registry_file)
outcomes <- c("Outcome_NutriMetab", "Outcome_TumorBurden", "Outcome_TreatComp")
component_names <- c(
  Outcome_NutriMetab = "NM_component",
  Outcome_TumorBurden = "TB_component",
  Outcome_TreatComp = "TC_component"
)

eps <- 1e-6
score_registry <- list()
for (outcome in outcomes) {
  outcome_name <- outcome
  selected <- features_dt[outcome == outcome_name, feature]
  if (!length(selected)) stop("No selected features for ", outcome)
  if (any(!selected %in% names(mapped))) stop("NHANES mapping missing final features for ", outcome)
  x <- as.matrix(mapped[idx, ..selected])
  storage.mode(x) <- "double"
  pool <- catboost.load_pool(data = x, feature_names = as.list(selected))
  reg <- model_registry[outcome == outcome_name]
  if (nrow(reg) != 1L) stop("Model registry mismatch for ", outcome)
  model_path <- file.path(out_dir, reg$model_file)
  model <- catboost.load_model(model_path)
  raw <- as.numeric(catboost.predict(model, pool, prediction_type = "Probability"))
  lp <- qlogis(pmin(pmax(raw, eps), 1 - eps))
  calibrated <- plogis(reg$intercept + reg$slope * lp)
  dat[[component_names[[outcome]]]] <- calibrated
  score_registry[[length(score_registry) + 1L]] <- data.table(
    scenario = scenario,
    outcome = outcome,
    component_column = component_names[[outcome]],
    n_features = length(selected),
    features = paste(selected, collapse = ";"),
    model_file = reg$model_file,
    calibration_intercept = reg$intercept,
    calibration_slope = reg$slope,
    score_min = min(calibrated),
    score_max = max(calibrated),
    score_missing_n = sum(is.na(calibrated))
  )
}
score_registry_dt <- rbindlist(score_registry)
fwrite(
  score_registry_dt,
  file.path(out_dir, paste0("nested_", scenario_slug, "_nhanes_score_registry.csv")),
  bom = TRUE
)
fwrite(
  dat[, .(
    SEQN, CYCLE, Cancer_Diagnosed, Dead, Follow_Up_Years,
    NM_component, TB_component, TC_component
  )],
  file.path(out_dir, paste0("nested_", scenario_slug, "_nhanes_component_scores.csv")),
  bom = TRUE
)

required <- c(
  "CYCLE", "WTMEC2YR", "SDMVSTRA", "SDMVPSU", "Follow_Up_Years", "Dead",
  "Age", "Sex", "RIDRETH1", "Cancer_Diagnosed", "Education", "PIR", "BMI",
  "Smoking_F", "Alcohol_F", "Diabetes", "Hypertension", "CVD",
  "NM_component", "TB_component", "TC_component", "PhenoAge_Accel"
)
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing survey variables: ", paste(missing_required, collapse = ", "))

n_cycles <- uniqueN(dat$CYCLE)
if (n_cycles != 9L) stop("Expected nine NHANES cycles")
dat[, pooled_mec_weight := cipds_pooled_mec_weights(dat, n_cycles)]
dat[, survey_strata := interaction(CYCLE, SDMVSTRA, drop = TRUE, lex.order = TRUE)]
dat[, survey_psu := interaction(CYCLE, SDMVSTRA, SDMVPSU, drop = TRUE, lex.order = TRUE)]
dat[, Sex_F_design := factor(Sex)]
dat[, Race_F_design := factor(RIDRETH1)]
dat[, Education_clean := ifelse(Education %in% 1:5, Education, NA_real_)]
dat[, Education_F_design := factor(Education_clean, levels = 1:5)]
dat[, Smoking_F_design := factor(Smoking_F)]
dat[, Alcohol_F_design := factor(Alcohol_F)]

design_base <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~pooled_mec_weight,
  data = dat,
  nest = TRUE
)

# Recompute the all-age PhenoAge residual using the corrected pooled MEC weights.
pheno_fit_corrected <- svyglm(PhenoAge ~ Age, design = design_base)
dat[, PhenoAge_Accel := PhenoAge - as.numeric(predict(pheno_fit_corrected, newdata = dat))]
design_base$variables$PhenoAge_Accel <- dat$PhenoAge_Accel

raw_scores <- c(
  NM_z = "NM_component",
  TB_z = "TB_component",
  TC_z = "TC_component",
  PhenoAge_Accel_z = "PhenoAge_Accel"
)
scaling_rows <- list()
for (z_name in names(raw_scores)) {
  raw_name <- raw_scores[[z_name]]
  mean_value <- as.numeric(coef(svymean(as.formula(paste0("~", raw_name)), design_base, na.rm = TRUE)))
  sd_value <- sqrt(as.numeric(svyvar(as.formula(paste0("~", raw_name)), design_base, na.rm = TRUE)))
  dat[[z_name]] <- (dat[[raw_name]] - mean_value) / sd_value
  scaling_rows[[length(scaling_rows) + 1L]] <- data.table(
    scenario = scenario,
    standardized_variable = z_name,
    source_variable = raw_name,
    survey_weighted_mean = mean_value,
    survey_weighted_sd = sd_value
  )
}
scaling_dt <- rbindlist(scaling_rows)
fwrite(
  scaling_dt,
  file.path(out_dir, paste0("nested_", scenario_slug, "_nhanes_survey_scaling.csv")),
  bom = TRUE
)

design_all <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~pooled_mec_weight,
  data = dat,
  nest = TRUE
)
populations <- list(
  Overall = rep(TRUE, nrow(dat)),
  Cancer = !is.na(dat$Cancer_Diagnosed) & dat$Cancer_Diagnosed == 1,
  Noncancer = !is.na(dat$Cancer_Diagnosed) & dat$Cancer_Diagnosed == 0
)
adjustments <- list(
  demographics_adjusted = c("Age", "Sex_F_design", "Race_F_design"),
  fully_adjusted = c(
    "Age", "Sex_F_design", "Race_F_design", "Education_F_design", "PIR", "BMI",
    "Smoking_F_design", "Alcohol_F_design", "Diabetes", "Hypertension", "CVD"
  )
)
component_labels <- c(
  NM_z = "Nutritional-Metabolic",
  TB_z = "Tumor Burden",
  TC_z = "Treatment Complications",
  PhenoAge_Accel_z = "PhenoAge Acceleration benchmark"
)

fit_svycox <- function(pop_name, pop_mask, score_terms, adjustment_name,
                       adjustment_terms, model_type) {
  required_vars <- unique(c(
    "Follow_Up_Years", "Dead", "pooled_mec_weight", "survey_strata", "survey_psu",
    score_terms, adjustment_terms
  ))
  complete <- complete.cases(dat[, ..required_vars]) & dat$pooled_mec_weight > 0 &
    dat$Follow_Up_Years > 0 & pop_mask
  sub_design <- design_all[complete, ]
  formula <- as.formula(paste(
    "Surv(Follow_Up_Years, Dead) ~",
    paste(c(score_terms, adjustment_terms), collapse = " + ")
  ))
  fit <- svycoxph(formula, design = sub_design)
  coef_table <- summary(fit)$coefficients
  conf <- confint(fit)
  rows <- lapply(score_terms, function(term) data.table(
    scenario = scenario,
    analysis_role = ifelse(adjustment_name == "fully_adjusted", "PRIMARY", "SENSITIVITY"),
    weighting = "NHANES complex survey design",
    population = pop_name,
    model_type = model_type,
    adjustment = adjustment_name,
    term = term,
    component_label = component_labels[[term]],
    n = sum(complete),
    events = sum(dat$Dead[complete] == 1),
    hazard_ratio_per_survey_weighted_sd = exp(coef(fit)[term]),
    ci_lower = exp(conf[term, 1]),
    ci_upper = exp(conf[term, 2]),
    p_value = coef_table[term, "Pr(>|z|)"],
    design_degrees_of_freedom = degf(sub_design)
  ))
  joint_p <- if (length(score_terms) > 1L) {
    tryCatch(
      as.numeric(regTermTest(
        fit, as.formula(paste("~", paste(score_terms, collapse = " + ")))
      )$p),
      error = function(e) NA_real_
    )
  } else NA_real_
  list(rows = rbindlist(rows), joint_p = joint_p, complete = complete)
}

separate_results <- list()
joint_results <- list()
joint_tests <- list()
for (pop_name in names(populations)) {
  for (adjustment_name in names(adjustments)) {
    pop_mask <- populations[[pop_name]]
    adjustment_terms <- adjustments[[adjustment_name]]
    for (score in c("NM_z", "TB_z", "TC_z", "PhenoAge_Accel_z")) {
      fitted <- fit_svycox(
        pop_name, pop_mask, score, adjustment_name, adjustment_terms,
        ifelse(score == "PhenoAge_Accel_z", "benchmark_separate", "component_separate")
      )
      separate_results[[length(separate_results) + 1L]] <- fitted$rows
    }
    joint <- fit_svycox(
      pop_name, pop_mask, c("NM_z", "TB_z", "TC_z"), adjustment_name,
      adjustment_terms, "three_components_joint_no_weighted_composite"
    )
    joint_results[[length(joint_results) + 1L]] <- joint$rows
    joint_tests[[length(joint_tests) + 1L]] <- data.table(
      scenario = scenario,
      analysis_role = ifelse(adjustment_name == "fully_adjusted", "PRIMARY", "SENSITIVITY"),
      population = pop_name,
      adjustment = adjustment_name,
      model_type = "three_components_joint_no_weighted_composite",
      n = sum(joint$complete),
      events = sum(dat$Dead[joint$complete] == 1),
      joint_wald_p = joint$joint_p
    )
  }
}

separate_dt <- rbindlist(separate_results)
joint_dt <- rbindlist(joint_results)
joint_test_dt <- rbindlist(joint_tests)
fwrite(
  separate_dt,
  file.path(out_dir, paste0("nested_", scenario_slug, "_nhanes_survey_component_separate.csv")),
  bom = TRUE
)
fwrite(
  joint_dt,
  file.path(out_dir, paste0("nested_", scenario_slug, "_nhanes_survey_components_joint.csv")),
  bom = TRUE
)
fwrite(
  joint_test_dt,
  file.path(out_dir, paste0("nested_", scenario_slug, "_nhanes_survey_joint_wald.csv")),
  bom = TRUE
)

design_audit <- data.table(
  scenario = scenario,
  item = c(
    "analysis_rows", "cycles", "raw_strata_codes", "cycle_specific_strata",
    "cycle_specific_psu", "design_degrees_of_freedom", "weight_variable",
    "pooled_weight_formula", "lonely_psu_handling", "primary_analysis",
    "weighted_composite_generated"
  ),
  value = c(
    nrow(dat), n_cycles, uniqueN(dat$SDMVSTRA), uniqueN(dat$survey_strata),
    uniqueN(dat$survey_psu), degf(design_all), "WTMEC2YR",
    CIPDS_WEIGHT_DESCRIPTION, "adjust",
    "fully adjusted survey-weighted Cox; components separate and joint", "FALSE"
  )
)
fwrite(
  design_audit,
  file.path(out_dir, paste0("nested_", scenario_slug, "_nhanes_survey_design_audit.csv")),
  bom = TRUE
)

cat("Scenario:", scenario, "\n")
cat("NHANES complex-survey analysis complete; no weighted composite generated.\n")
print(separate_dt[analysis_role == "PRIMARY" & model_type == "component_separate"])
cat("\nJoint model:\n")
print(joint_dt[analysis_role == "PRIMARY"])
