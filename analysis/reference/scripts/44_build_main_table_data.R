suppressPackageStartupMessages({
  library(data.table)
  library(survey)
  library(survival)
})

options(survey.lonely.psu = "adjust")
setDTthreads(min(23L, max(1L, parallel::detectCores(logical = TRUE))))

ROOT <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
OUT <- file.path(ROOT, "outputs", "main_tables")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

HOSPITAL_RELEASE <- Sys.getenv(
  "CIPDS_HOSPITAL_RELEASE",
  unset = "hospital_source"
)

fmt_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_smd <- function(x) {
  if (!is.finite(x)) return(NA_character_)
  sprintf("%.3f", abs(x))
}

fmt_mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  sprintf("%.2f (%.2f)", mean(x), sd(x))
}

fmt_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE, type = 2)
  sprintf("%.2f [%.2f, %.2f]", q[2], q[1], q[3])
}

fmt_n_pct <- function(x, level = 1) {
  ok <- !is.na(x)
  if (!any(ok)) return(NA_character_)
  n <- sum(x[ok] == level)
  sprintf("%s (%.1f%%)", format(n, big.mark = ",", scientific = FALSE), 100 * n / sum(ok))
}

continuous_smd <- function(a, b) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (length(a) < 2 || length(b) < 2) return(NA_real_)
  denom <- sqrt((var(a) + var(b)) / 2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  (mean(b) - mean(a)) / denom
}

binary_smd <- function(p1, p0) {
  denom <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  (p1 - p0) / denom
}

safe_wilcox_p <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  if (sum(ok) < 4 || length(unique(g[ok])) < 2) return(NA_real_)
  tryCatch(wilcox.test(x[ok] ~ g[ok], exact = FALSE)$p.value, error = function(e) NA_real_)
}

safe_t_p <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  if (sum(ok) < 4 || length(unique(g[ok])) < 2) return(NA_real_)
  tryCatch(t.test(x[ok] ~ g[ok], var.equal = FALSE)$p.value, error = function(e) NA_real_)
}

safe_chisq_p <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  if (sum(ok) < 4 || length(unique(x[ok])) < 2 || length(unique(g[ok])) < 2) return(NA_real_)
  tab <- table(x[ok], g[ok])
  tryCatch(suppressWarnings(chisq.test(tab, correct = FALSE)$p.value), error = function(e) NA_real_)
}

hospital_cont_row <- function(dat, var, label, section, display = c("median", "mean")) {
  display <- match.arg(display)
  a <- as.numeric(dat[cohort == "Development"][[var]])
  b <- as.numeric(dat[cohort == "Calendar test"][[var]])
  stat <- if (display == "median") fmt_median_iqr else fmt_mean_sd
  p <- if (display == "median") safe_wilcox_p(as.numeric(dat[[var]]), dat$cohort) else safe_t_p(as.numeric(dat[[var]]), dat$cohort)
  data.table(
    section = section,
    characteristic = label,
    development = stat(a),
    development_missing_pct = sprintf("%.1f", 100 * mean(!is.finite(a))),
    calendar_test = stat(b),
    calendar_test_missing_pct = sprintf("%.1f", 100 * mean(!is.finite(b))),
    p_value = fmt_p(p),
    smd = fmt_smd(continuous_smd(a, b)),
    method = if (display == "median") "Wilcoxon rank-sum" else "Welch t test"
  )
}

hospital_binary_row <- function(dat, var, label, section, level = 1) {
  a <- dat[cohort == "Development"][[var]]
  b <- dat[cohort == "Calendar test"][[var]]
  p1 <- mean(b == level, na.rm = TRUE)
  p0 <- mean(a == level, na.rm = TRUE)
  data.table(
    section = section,
    characteristic = label,
    development = fmt_n_pct(a, level),
    development_missing_pct = sprintf("%.1f", 100 * mean(is.na(a))),
    calendar_test = fmt_n_pct(b, level),
    calendar_test_missing_pct = sprintf("%.1f", 100 * mean(is.na(b))),
    p_value = fmt_p(safe_chisq_p(dat[[var]], dat$cohort)),
    smd = fmt_smd(binary_smd(p1, p0)),
    method = "Pearson chi-square"
  )
}

# -----------------------------------------------------------------------------
# Table 1, Panel A: hospital development period vs independent calendar test
# -----------------------------------------------------------------------------
model <- fread(file.path(ROOT, "outputs", "hospital_patient_model_ready_129.csv"))

core_cols <- c(
  "patient_key", "age_at_first_observed_check", "sex_status",
  "current_primary_site_systems_documented_ever",
  "flag_semantic__multiple_current_primary_documented_ever",
  "flag_patient_source_outpatient_ever", "flag_patient_source_inpatient_ever",
  "flag_care_icu_ever", "flag_care_emergency_ever",
  "flag_oncology_department_context_ever", "flag_radiotherapy_ever_final",
  "flag_chemotherapy_ever_final", "flag_intervention_ever_final",
  "flag_surgery_ever_final"
)
core <- fread(file.path(HOSPITAL_RELEASE, "142_patient_nonlab_core_01na_v1.csv"), select = core_cols)
stopifnot(nrow(model) == 75248L, uniqueN(model$patient_key) == 75248L)
stopifnot(nrow(core) == 75248L, uniqueN(core$patient_key) == 75248L)
hospital <- merge(model, core, by = "patient_key", all.x = TRUE, sort = FALSE)
stopifnot(nrow(hospital) == 75248L)
hospital[, cohort := fifelse(split_calendar_entry == "test", "Calendar test", "Development")]
hospital[, cohort := factor(cohort, levels = c("Development", "Calendar test"))]
hospital[, age_num := suppressWarnings(as.numeric(age_at_first_observed_check))]
hospital[, male := fifelse(sex_status == "MALE", 1, fifelse(sex_status == "FEMALE", 0, NA_real_))]
hospital[, age_ge60 := fifelse(is.finite(age_num), as.numeric(age_num >= 60), NA_real_)]
hospital[, age_ge65 := fifelse(is.finite(age_num), as.numeric(age_num >= 65), NA_real_)]
hospital[, age_ge75 := fifelse(is.finite(age_num), as.numeric(age_num >= 75), NA_real_)]

site <- tolower(fifelse(
  is.na(hospital$current_primary_site_systems_documented_ever),
  "", hospital$current_primary_site_systems_documented_ever
))
site_documented <- nzchar(site) & !site %chin% c("unknown", "not_applicable")
site_patterns <- list(
  site_lung = "lung|respiratory",
  site_gastrointestinal = "colorectal|stomach|esophagus|gastro|small_intestine",
  site_hepatobiliary_pancreatic = "liver|biliary|bile_duct|pancreas",
  site_breast = "breast",
  site_gynecologic = "cervix|endometrium|uterus|ovary|fallopian|vulva|gynecologic",
  site_genitourinary = "kidney|renal|bladder|prostate|testis|ureter|urethra|urinary|penis|urachal",
  site_head_neck = "nasopharynx|oropharynx|hypopharynx|larynx|oral_cavity|salivary|sinonasal|head_face|neck|thyroid",
  site_hematologic = "hematologic|hematolymphoid|lymph",
  site_cns = "central_nervous_system|cranial_nerve|neurologic",
  site_bone_soft_tissue = "bone|connective_tissue|musculoskeletal"
)
site_any <- rep(FALSE, nrow(hospital))
for (nm in names(site_patterns)) {
  hospital[[nm]] <- as.integer(site_documented & grepl(site_patterns[[nm]], site))
  site_any <- site_any | hospital[[nm]] == 1
}
hospital[, site_other_documented := as.integer(site_documented & !site_any)]
hospital[, site_undocumented := as.integer(!site_documented)]

panel_a <- rbindlist(list(
  data.table(
    section = "Cohort definition",
    characteristic = "Participants, n",
    development = format(sum(hospital$cohort == "Development"), big.mark = ","),
    development_missing_pct = "0.0",
    calendar_test = format(sum(hospital$cohort == "Calendar test"), big.mark = ","),
    calendar_test_missing_pct = "0.0",
    p_value = NA_character_, smd = NA_character_, method = "Not applicable"
  ),
  data.table(
    section = "Cohort definition",
    characteristic = "First test observation date range",
    development = paste(min(hospital[cohort == "Development"]$first_observed_test_date),
                        max(hospital[cohort == "Development"]$first_observed_test_date), sep = " to "),
    development_missing_pct = sprintf("%.1f", 100 * mean(is.na(hospital[cohort == "Development"]$first_observed_test_date))),
    calendar_test = paste(min(hospital[cohort == "Calendar test"]$first_observed_test_date),
                          max(hospital[cohort == "Calendar test"]$first_observed_test_date), sep = " to "),
    calendar_test_missing_pct = sprintf("%.1f", 100 * mean(is.na(hospital[cohort == "Calendar test"]$first_observed_test_date))),
    p_value = NA_character_, smd = NA_character_, method = "Not applicable"
  ),
  hospital_cont_row(hospital, "age_num", "Age at first test observation, years", "Demographics", "mean"),
  hospital_binary_row(hospital, "male", "Male sex", "Demographics"),
  hospital_binary_row(hospital, "age_ge60", "Age >=60 years", "Demographics"),
  hospital_binary_row(hospital, "age_ge65", "Age >=65 years", "Demographics"),
  hospital_binary_row(hospital, "age_ge75", "Age >=75 years", "Demographics"),
  hospital_cont_row(hospital, "observation_span_days", "Test observation span, days", "Observation structure", "median"),
  hospital_cont_row(hospital, "encounter_count", "Encounters per patient", "Observation structure", "median"),
  hospital_cont_row(hospital, "test_record_count", "Laboratory records per patient", "Observation structure", "median"),
  hospital_binary_row(hospital, "flag_patient_source_outpatient_ever", "Outpatient source documented ever", "Care context"),
  hospital_binary_row(hospital, "flag_patient_source_inpatient_ever", "Inpatient source documented ever", "Care context"),
  hospital_binary_row(hospital, "flag_oncology_department_context_ever", "Oncology department context documented ever", "Care context"),
  hospital_binary_row(hospital, "flag_care_icu_ever", "Intensive care context documented ever", "Care context"),
  hospital_binary_row(hospital, "flag_care_emergency_ever", "Emergency context documented ever", "Care context"),
  hospital_binary_row(hospital, "flag_radiotherapy_ever_final", "Radiotherapy documented ever", "Treatment context"),
  hospital_binary_row(hospital, "flag_chemotherapy_ever_final", "Chemotherapy documented ever", "Treatment context"),
  hospital_binary_row(hospital, "flag_intervention_ever_final", "Interventional treatment documented ever", "Treatment context"),
  hospital_binary_row(hospital, "flag_surgery_ever_final", "Surgery documented ever", "Treatment context")
), use.names = TRUE, fill = TRUE)

site_labels <- c(
  site_lung = "Lung or respiratory",
  site_gastrointestinal = "Gastrointestinal",
  site_hepatobiliary_pancreatic = "Hepatobiliary or pancreatic",
  site_breast = "Breast",
  site_gynecologic = "Gynecologic",
  site_genitourinary = "Genitourinary",
  site_head_neck = "Head and neck or thyroid",
  site_hematologic = "Hematologic or lymphatic",
  site_cns = "Central nervous system",
  site_bone_soft_tissue = "Bone or soft tissue",
  site_other_documented = "Other documented primary site",
  site_undocumented = "Current primary site not documented"
)
panel_a <- rbind(
  panel_a,
  rbindlist(lapply(names(site_labels), function(v) {
    hospital_binary_row(hospital, v, site_labels[[v]], "Documented current primary site systems, non-mutually exclusive")
  })),
  hospital_binary_row(hospital, "flag_semantic__multiple_current_primary_documented_ever",
                      "Multiple current primary sites documented ever",
                      "Documented current primary site systems, non-mutually exclusive")
)

outcome_labels <- c(
  Outcome_NutriMetab = "Nutritional-metabolic domain positive",
  Outcome_TumorBurden = "Tumor-burden-related domain positive",
  Outcome_TreatComp = "Treatment-complication-related domain positive"
)
for (v in names(outcome_labels)) {
  panel_a <- rbind(panel_a, hospital_binary_row(hospital, v, outcome_labels[[v]], "Legacy patient-ever domain states"))
}

features <- unique(fread(file.path(ROOT, "outputs", "nested_a_primary_final_features.csv"))$feature)
elig <- fread(file.path(ROOT, "outputs", "lab_129_crosscohort_eligibility.csv"))
lab_meta <- elig[final_variable_id %chin% features,
                 .(final_variable_id, lab_order, canonical_name_english, analysis_unit)]
setorder(lab_meta, lab_order)
stopifnot(nrow(lab_meta) == 26L)

# Add the two PhenoAge-specific laboratory inputs not retained by the frozen
# component-model feature union, so Table 1 covers both compared constructs.
benchmark_extra <- elig[final_variable_id %chin% c("lab_glu", "lab_hcrp"),
                        .(final_variable_id, lab_order, canonical_name_english, analysis_unit)]
lab_meta_table1 <- unique(rbind(lab_meta, benchmark_extra), by = "final_variable_id")

for (i in seq_len(nrow(lab_meta_table1))) {
  m <- lab_meta_table1[i]
  label <- sprintf("%s, %s", m$canonical_name_english, m$analysis_unit)
  section <- if (m$final_variable_id %chin% features) {
    "Laboratory variables retained by at least one frozen component model"
  } else {
    "Additional laboratory variables required for the PhenoAge benchmark"
  }
  panel_a <- rbind(panel_a, hospital_cont_row(hospital, m$final_variable_id, label, section, "median"))
}

fwrite(panel_a, file.path(OUT, "table1_panel_a_hospital.csv"), bom = TRUE)

# -----------------------------------------------------------------------------
# Table 1, Panel B: NHANES age >=60 survey domain, by cancer history
# -----------------------------------------------------------------------------
source(file.path(ROOT, "scripts", "27_supplement_common.R"))
loaded <- cipds_load_data(include_labs = TRUE)
nh <- loaded$dat
design_full <- loaded$design_full
age60_mask <- loaded$age60
design_age60 <- design_full[age60_mask, ]
design_compare <- design_full[age60_mask & !is.na(nh$Cancer_Diagnosed), ]
design_cancer <- design_full[age60_mask & !is.na(nh$Cancer_Diagnosed) & nh$Cancer_Diagnosed == 1, ]
design_noncancer <- design_full[age60_mask & !is.na(nh$Cancer_Diagnosed) & nh$Cancer_Diagnosed == 0, ]

weighted_missing_pct <- function(des, var) {
  x <- des$variables[[var]]
  w <- weights(des, type = "sampling")
  if (!length(w) || sum(w, na.rm = TRUE) <= 0) return(NA_real_)
  100 * weighted.mean(is.na(x), w, na.rm = TRUE)
}

weighted_mean_sd <- function(des, var) {
  f <- as.formula(paste0("~", var))
  mu <- tryCatch(as.numeric(coef(svymean(f, des, na.rm = TRUE)))[1], error = function(e) NA_real_)
  vv <- tryCatch(as.numeric(svyvar(f, des, na.rm = TRUE))[1], error = function(e) NA_real_)
  c(mean = mu, sd = sqrt(vv))
}

weighted_prop <- function(des, var, level) {
  x <- des$variables[[var]]
  w <- weights(des, type = "sampling")
  weighted.mean(x == level, w, na.rm = TRUE)
}

weighted_n_pct <- function(des, var, level) {
  x <- des$variables[[var]]
  n <- sum(x == level, na.rm = TRUE)
  p <- weighted_prop(des, var, level)
  sprintf("%s (%.1f%%)", format(n, big.mark = ",", scientific = FALSE), 100 * p)
}

survey_cont_p <- function(var) {
  f <- as.formula(paste(var, "~ Cancer_F_design"))
  fit <- tryCatch(svyglm(f, design = design_compare), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  tryCatch(regTermTest(fit, ~Cancer_F_design)$p, error = function(e) NA_real_)
}

survey_cat_p <- function(var) {
  f <- as.formula(paste("~", var, "+ Cancer_F_design"))
  tryCatch(as.numeric(svychisq(f, design_compare, statistic = "F")$p.value), error = function(e) NA_real_)
}

nhanes_cont_row <- function(var, label, section) {
  all <- weighted_mean_sd(design_age60, var)
  ca <- weighted_mean_sd(design_cancer, var)
  no <- weighted_mean_sd(design_noncancer, var)
  smd <- (ca[["mean"]] - no[["mean"]]) / sqrt((ca[["sd"]]^2 + no[["sd"]]^2) / 2)
  data.table(
    section = section,
    characteristic = label,
    all_age60 = sprintf("%.2f (%.2f)", all[["mean"]], all[["sd"]]),
    all_missing_pct = sprintf("%.1f", weighted_missing_pct(design_age60, var)),
    cancer_history = sprintf("%.2f (%.2f)", ca[["mean"]], ca[["sd"]]),
    cancer_missing_pct = sprintf("%.1f", weighted_missing_pct(design_cancer, var)),
    no_cancer_history = sprintf("%.2f (%.2f)", no[["mean"]], no[["sd"]]),
    no_cancer_missing_pct = sprintf("%.1f", weighted_missing_pct(design_noncancer, var)),
    p_value = fmt_p(survey_cont_p(var)),
    smd = fmt_smd(smd),
    method = "Survey-weighted Wald test"
  )
}

nhanes_cat_rows <- function(var, label, section, levels, level_labels) {
  p <- survey_cat_p(var)
  rbindlist(lapply(seq_along(levels), function(i) {
    lv <- levels[[i]]
    pca <- weighted_prop(design_cancer, var, lv)
    pno <- weighted_prop(design_noncancer, var, lv)
    data.table(
      section = section,
      characteristic = if (length(levels) == 1L) label else paste0(label, ": ", level_labels[[i]]),
      all_age60 = weighted_n_pct(design_age60, var, lv),
      all_missing_pct = sprintf("%.1f", weighted_missing_pct(design_age60, var)),
      cancer_history = weighted_n_pct(design_cancer, var, lv),
      cancer_missing_pct = sprintf("%.1f", weighted_missing_pct(design_cancer, var)),
      no_cancer_history = weighted_n_pct(design_noncancer, var, lv),
      no_cancer_missing_pct = sprintf("%.1f", weighted_missing_pct(design_noncancer, var)),
      p_value = if (i == 1L) fmt_p(p) else NA_character_,
      smd = fmt_smd(binary_smd(pca, pno)),
      method = if (i == 1L) "Rao-Scott design-adjusted chi-square" else NA_character_
    )
  }))
}

n_all <- nrow(design_age60$variables)
n_ca <- nrow(design_cancer$variables)
n_no <- nrow(design_noncancer$variables)
panel_b <- data.table(
  section = "Cohort definition",
  characteristic = "Participants, unweighted n",
  all_age60 = format(n_all, big.mark = ","), all_missing_pct = "0.0",
  cancer_history = format(n_ca, big.mark = ","), cancer_missing_pct = "0.0",
  no_cancer_history = format(n_no, big.mark = ","), no_cancer_missing_pct = "0.0",
  p_value = NA_character_, smd = NA_character_, method = "Not applicable"
)

nh[, age_ge75 := fifelse(is.finite(Age), as.integer(Age >= 75), NA_integer_)]
nh[, poverty_pir_lt1 := fifelse(is.finite(PIR), as.integer(PIR < 1), NA_integer_)]
# Rebuild the designs after adding derived Table 1 variables.
design_full_t1 <- svydesign(ids = ~survey_psu, strata = ~survey_strata,
                            weights = ~pooled_mec_weight, data = nh, nest = TRUE)
design_age60 <- design_full_t1[age60_mask, ]
design_compare <- design_full_t1[age60_mask & !is.na(nh$Cancer_Diagnosed), ]
design_cancer <- design_full_t1[age60_mask & !is.na(nh$Cancer_Diagnosed) & nh$Cancer_Diagnosed == 1, ]
design_noncancer <- design_full_t1[age60_mask & !is.na(nh$Cancer_Diagnosed) & nh$Cancer_Diagnosed == 0, ]

panel_b <- rbind(
  panel_b,
  nhanes_cont_row("Age", "Age, years", "Demographics and socioeconomic factors"),
  nhanes_cat_rows("age_ge75", "Age >=75 years", "Demographics and socioeconomic factors", list(1L), list("Yes")),
  nhanes_cat_rows("Sex", "Male sex", "Demographics and socioeconomic factors", list(1L), list("Male")),
  nhanes_cat_rows("RIDRETH1", "Race or ethnicity", "Demographics and socioeconomic factors",
                  as.list(1:5), as.list(c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Other or multiracial"))),
  nhanes_cat_rows("Education_clean", "Education", "Demographics and socioeconomic factors",
                  as.list(1:5), as.list(c("Less than 9th grade", "9th-11th grade", "High school or GED", "Some college or associate degree", "College graduate or above"))),
  nhanes_cont_row("PIR", "Family poverty-income ratio", "Demographics and socioeconomic factors"),
  nhanes_cat_rows("poverty_pir_lt1", "Family poverty-income ratio <1", "Demographics and socioeconomic factors", list(1L), list("Yes")),
  nhanes_cat_rows("CYCLE", "NHANES cycle", "Survey cycle",
                  as.list(levels(nh$CYCLE)), as.list(levels(nh$CYCLE))),
  nhanes_cont_row("BMI", "Body mass index, kg/m2", "Lifestyle and clinical characteristics"),
  nhanes_cat_rows("Smoking_F", "Smoking status", "Lifestyle and clinical characteristics",
                  as.list(c("Never", "Former", "Current")), as.list(c("Never", "Former", "Current"))),
  nhanes_cat_rows("Alcohol_F", "Alcohol history", "Lifestyle and clinical characteristics",
                  as.list(c("Non-drinker", "Drinker")), as.list(c("Below threshold", "At least 12 drinks in any year"))),
  nhanes_cat_rows("Diabetes", "Diabetes", "Lifestyle and clinical characteristics", list(1L), list("Yes")),
  nhanes_cat_rows("Hypertension", "Hypertension", "Lifestyle and clinical characteristics", list(1L), list("Yes")),
  nhanes_cat_rows("CVD", "Cardiovascular disease", "Lifestyle and clinical characteristics", list(1L), list("Yes")),
  nhanes_cont_row("Follow_Up_Years", "Follow-up, years", "Mortality follow-up"),
  nhanes_cat_rows("Death_AllCause", "All-cause death", "Mortality follow-up", list(1L), list("Yes")),
  nhanes_cat_rows("Death_Cancer", "Cancer death", "Mortality follow-up", list(1L), list("Yes")),
  nhanes_cat_rows("Death_CVD", "Cardiovascular death", "Mortality follow-up", list(1L), list("Yes")),
  nhanes_cont_row("NM_component", "Nutritional-metabolic component probability", "Transferred laboratory scores"),
  nhanes_cont_row("TB_component", "Tumor-burden-related component probability", "Transferred laboratory scores"),
  nhanes_cont_row("TC_component", "Treatment-complication-related component probability", "Transferred laboratory scores"),
  nhanes_cont_row("Overall_expected_burden", "Overall expected burden, range 0-3", "Transferred laboratory scores"),
  nhanes_cont_row("PhenoAge_acceleration", "PhenoAge Acceleration, years", "Transferred laboratory scores")
)

for (i in seq_len(nrow(lab_meta_table1))) {
  m <- lab_meta_table1[i]
  label <- sprintf("%s, %s", m$canonical_name_english, m$analysis_unit)
  section <- if (m$final_variable_id %chin% features) {
    "Laboratory variables retained by at least one frozen component model"
  } else {
    "Additional laboratory variables required for the PhenoAge benchmark"
  }
  panel_b <- rbind(panel_b, nhanes_cont_row(m$final_variable_id, label, section))
}

fwrite(panel_b, file.path(OUT, "table1_panel_b_nhanes.csv"), bom = TRUE)

# -----------------------------------------------------------------------------
# Table 2 supporting calculation: paired cancer-history interaction for all
# five scores, including the PhenoAge benchmark on its own complete-case cohort.
# -----------------------------------------------------------------------------
clinical_no_cancer <- c(
  "Age", "Sex_F_design", "Race_F_design", "Education_F_design", "PIR", "BMI",
  "Smoking_F_design", "Alcohol_F_design", "Diabetes", "Hypertension", "CVD"
)
score_vars <- c(
  OVERALL = "OVERALL_z_age60", NM = "NM_z_age60", TB = "TB_z_age60",
  TC = "TC_z_age60", PHENO = "PHENO_z_age60"
)
score_labels <- c(
  OVERALL = "Multidomain laboratory vulnerability overall score",
  NM = "Nutritional-metabolic laboratory pattern",
  TB = "Tumor-burden-related laboratory pattern",
  TC = "Treatment-complication-related laboratory pattern",
  PHENO = "PhenoAge Acceleration"
)

interaction_rows <- rbindlist(lapply(names(score_vars), function(short) {
  score <- score_vars[[short]]
  required <- c("Follow_Up_Years", "Death_AllCause", score, "Cancer_F_design", clinical_no_cancer)
  mask <- age60_mask & complete.cases(nh[, ..required]) &
    is.finite(nh$Follow_Up_Years) & nh$Follow_Up_Years > 0 &
    is.finite(nh$pooled_mec_weight) & nh$pooled_mec_weight > 0
  des <- design_full_t1[mask, ]
  f <- as.formula(paste(
    "Surv(Follow_Up_Years, Death_AllCause) ~",
    paste0(score, " * Cancer_F_design + ", paste(clinical_no_cancer, collapse = " + "))
  ))
  fit <- svycoxph(f, design = des)
  b <- coef(fit)
  V <- vcov(fit)
  main <- score
  interaction <- grep(paste0("^", score, ":Cancer_F_design1$|^Cancer_F_design1:", score, "$"), names(b), value = TRUE)
  if (length(interaction) != 1L) stop("Interaction coefficient not uniquely identified for ", short)
  bi <- interaction[[1]]
  z <- qnorm(0.975)
  beta0 <- b[[main]]
  se0 <- sqrt(V[main, main])
  beta1 <- b[[main]] + b[[bi]]
  se1 <- sqrt(V[main, main] + V[bi, bi] + 2 * V[main, bi])
  int_se <- sqrt(V[bi, bi])
  int_p <- 2 * pnorm(-abs(b[[bi]] / int_se))
  data.table(
    component = short,
    component_label = score_labels[[short]],
    n = nrow(des$variables),
    events = sum(des$variables$Death_AllCause),
    hr_non_cancer_history = exp(beta0),
    non_cancer_ci_lower = exp(beta0 - z * se0),
    non_cancer_ci_upper = exp(beta0 + z * se0),
    hr_cancer_history = exp(beta1),
    cancer_ci_lower = exp(beta1 - z * se1),
    cancer_ci_upper = exp(beta1 + z * se1),
    interaction_ratio_of_hrs = exp(b[[bi]]),
    interaction_ci_lower = exp(b[[bi]] - z * int_se),
    interaction_ci_upper = exp(b[[bi]] + z * int_se),
    interaction_p = int_p,
    design_degrees_of_freedom = degf(des)
  )
}))
interaction_rows[, interaction_p_holm := p.adjust(interaction_p, method = "holm")]
fwrite(interaction_rows, file.path(OUT, "table2_panel_d_cancer_history_interactions_all_scores.csv"), bom = TRUE)

# -----------------------------------------------------------------------------
# Table 2 supporting calculation: bootstrap confidence intervals for AP.
# No model is refitted; only frozen calendar-test predictions are resampled.
# -----------------------------------------------------------------------------
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

ap_boot_ci <- function(y, p, B = 1000L, seed = 20260904L) {
  ok <- !is.na(y) & is.finite(p)
  y <- as.integer(y[ok])
  p <- p[ok]
  pos <- which(y == 1L)
  neg <- which(y == 0L)
  set.seed(seed)
  vals <- vapply(seq_len(B), function(i) {
    idx <- c(sample(pos, length(pos), replace = TRUE), sample(neg, length(neg), replace = TRUE))
    average_precision(y[idx], p[idx])
  }, numeric(1))
  c(point = average_precision(y, p), lower = quantile(vals, 0.025, na.rm = TRUE, names = FALSE),
    upper = quantile(vals, 0.975, na.rm = TRUE, names = FALSE), replicates = B)
}

pred_component <- fread(file.path(ROOT, "outputs", "nested_a_primary_calendar_test_predictions.csv"))
ap_rows <- rbindlist(lapply(unique(pred_component$outcome), function(outcome_id) {
  z <- pred_component[outcome == outcome_id]
  ci <- ap_boot_ci(z$observed, z$calibrated_probability,
                   seed = 20260904L + match(outcome_id, unique(pred_component$outcome)))
  data.table(model = outcome_id, target = outcome_id, n = nrow(z), positives = sum(z$observed),
             average_precision = ci[["point"]], ap_ci_lower = ci[["lower"]],
             ap_ci_upper = ci[["upper"]], bootstrap_replicates = as.integer(ci[["replicates"]]))
}))

pred_overall <- fread(file.path(ROOT, "outputs", "overall_calendar_predictions.csv"))[split == "calendar_test"]
for (target in c("any", "multidomain")) {
  y <- if (target == "any") as.integer(pred_overall$observed_burden >= 1L) else as.integer(pred_overall$observed_burden >= 2L)
  p <- if (target == "any") pred_overall$probability_any_domain else pred_overall$probability_multidomain
  ci <- ap_boot_ci(y, p, seed = 20260920L + match(target, c("any", "multidomain")))
  ap_rows <- rbind(ap_rows, data.table(
    model = "Overall", target = target, n = length(y), positives = sum(y),
    average_precision = ci[["point"]], ap_ci_lower = ci[["lower"]],
    ap_ci_upper = ci[["upper"]], bootstrap_replicates = as.integer(ci[["replicates"]])
  ))
}
fwrite(ap_rows, file.path(OUT, "table2_average_precision_bootstrap_ci.csv"), bom = TRUE)

# -----------------------------------------------------------------------------
# Multiplicity correction for the five fully adjusted all-cause score models.
# -----------------------------------------------------------------------------
overall_assoc <- fread(file.path(ROOT, "outputs", "overall_nhanes_allcause.csv"))[
  domain == "age60" & population == "Overall"
]
component_assoc <- fread(file.path(ROOT, "outputs", "older_nhanes_allcause_component_separate.csv"))[
  domain == "age60" & population == "Overall"
]
allcause_full <- rbind(
  overall_assoc[adjustment == "fully_adjusted",
                .(component, component_label, n, events,
                  hr = hazard_ratio_per_domain_weighted_sd, ci_lower, ci_upper, p_value)],
  component_assoc[adjustment == "fully_adjusted",
                  .(component, component_label, n, events,
                    hr = hazard_ratio_per_domain_weighted_sd, ci_lower, ci_upper, p_value)]
)
allcause_full <- unique(allcause_full, by = "component")
allcause_full <- allcause_full[match(c("OVERALL", "NM", "TB", "TC", "PHENO"), component)]
allcause_full[, p_holm := p.adjust(p_value, method = "holm")]
fwrite(allcause_full, file.path(OUT, "table2_panel_c_allcause_holm.csv"), bom = TRUE)

# QA checks for cohort denominators and frozen performance reproduction.
metric <- fread(file.path(ROOT, "outputs", "nested_a_primary_calendar_metrics.csv"))[split == "calendar_test"]
overall_perf <- fread(file.path(ROOT, "outputs", "overall_calendar_performance.csv"))[
  split == "calendar_test" & candidate_model == "main_effects"
]
qa <- rbindlist(list(
  data.table(check = "hospital_total_n", passed = nrow(hospital) == 75248L,
             observed = as.character(nrow(hospital)), expected = "75248"),
  data.table(check = "hospital_development_n", passed = sum(hospital$cohort == "Development") == 63983L,
             observed = as.character(sum(hospital$cohort == "Development")), expected = "63983"),
  data.table(check = "hospital_calendar_test_n", passed = sum(hospital$cohort == "Calendar test") == 11265L,
             observed = as.character(sum(hospital$cohort == "Calendar test")), expected = "11265"),
  data.table(check = "nhanes_age60_n", passed = n_all == 15048L,
             observed = as.character(n_all), expected = "15048"),
  data.table(check = "nhanes_cancer_n", passed = n_ca == 2897L,
             observed = as.character(n_ca), expected = "2897"),
  data.table(check = "nhanes_noncancer_n", passed = n_no == 12128L,
             observed = as.character(n_no), expected = "12128"),
  data.table(check = "component_ap_points_match", passed = all(abs(
    ap_rows[model != "Overall"]$average_precision -
      metric[match(ap_rows[model != "Overall"]$model, outcome)]$average_precision
  ) < 1e-12), observed = "recomputed", expected = "absolute difference <1e-12"),
  data.table(check = "overall_ap_any_matches", passed = abs(
    ap_rows[model == "Overall" & target == "any"]$average_precision - overall_perf$any_average_precision
  ) < 1e-12, observed = "recomputed", expected = "absolute difference <1e-12"),
  data.table(check = "overall_ap_multidomain_matches", passed = abs(
    ap_rows[model == "Overall" & target == "multidomain"]$average_precision - overall_perf$multidomain_average_precision
  ) < 1e-12, observed = "recomputed", expected = "absolute difference <1e-12"),
  data.table(check = "interaction_rows_complete", passed = nrow(interaction_rows) == 5L,
             observed = as.character(nrow(interaction_rows)), expected = "5"),
  data.table(check = "allcause_rows_complete", passed = nrow(allcause_full) == 5L,
             observed = as.character(nrow(allcause_full)), expected = "5")
))
fwrite(qa, file.path(OUT, "main_table_build_qa.csv"), bom = TRUE)
if (!all(qa$passed)) stop("Main-table data QA failed")

cat("Main-table data created in", OUT, "\n")
print(qa)
