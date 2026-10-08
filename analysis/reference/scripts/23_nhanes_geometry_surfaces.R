#!/usr/bin/env Rscript

# Formal, design-based comparison of the frozen Overall score geometry with
# PhenoAge Acceleration in the NHANES age >=60 common-complete-case cohort.

suppressPackageStartupMessages({
  library(data.table)
  library(survey)
  library(jsonlite)
})

options(survey.lonely.psu = "adjust")
set.seed(20260831)
package_dir <- Sys.getenv("CIPDS_PACKAGE_DIR", ".")
out_dir <- file.path(package_dir, "outputs")
geo_dir <- file.path(out_dir, "geometry_v1")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

input_path <- file.path(geo_dir, "nhanes_geometry_analysis.csv")
coef_path <- file.path(out_dir, "overall_stacked_model_coefficients.csv")
scale_path <- file.path(out_dir, "overall_meta_logit_scaling.csv")
stopifnot(file.exists(input_path), file.exists(coef_path), file.exists(scale_path))

dat <- fread(input_path)
if (nrow(dat) != as.integer(fread(file.path(package_dir,"qa/cohort_registry.csv"))[cohort=="paired",n]) || uniqueN(dat$SEQN) != as.integer(fread(file.path(package_dir,"qa/cohort_registry.csv"))[cohort=="paired",n])) stop("Primary paired cohort grain drift")
required <- c(
  "SEQN", "survey_strata", "survey_psu", "pooled_mec_weight",
  "Overall_expected_burden", "PhenoAge_acceleration",
  "NM_component", "TB_component", "TC_component"
)
if (!all(required %in% names(dat))) stop("Geometry surface input schema drift")
if (any(!is.finite(dat$pooled_mec_weight) | dat$pooled_mec_weight <= 0)) stop("Invalid survey weight")
dat[, `:=`(
  survey_strata = factor(survey_strata),
  survey_psu = factor(survey_psu)
)]

meta_scale <- fread(scale_path)
meta_coef <- fread(coef_path)
component_order <- c("NM", "TB", "TC")
epsilon <- 1e-12
for (component_name in component_order) {
  scaling <- meta_scale[component == component_name]
  if (nrow(scaling) != 1L || !is.finite(scaling$sd) || scaling$sd <= 0) {
    stop("Invalid frozen scaling for ", component_name)
  }
  probability <- pmin(pmax(dat[[paste0(component_name, "_component")]], epsilon), 1 - epsilon)
  dat[[paste0("z_", component_name)]] <- (qlogis(probability) - scaling$mean) / scaling$sd
}

beta_overall <- setNames(
  meta_coef[parameter_type == "coefficient"][match(paste0("z_", component_order), parameter), estimate],
  paste0("z_", component_order)
)
threshold <- setNames(
  meta_coef[parameter_type == "threshold"][match(c("0|1", "1|2", "2|3"), parameter), estimate],
  c("0|1", "1|2", "2|3")
)
if (any(!is.finite(beta_overall)) || any(!is.finite(threshold))) stop("Frozen ordinal parameters unavailable")

dat[, Overall_eta :=
  beta_overall[["z_NM"]] * z_NM +
  beta_overall[["z_TB"]] * z_TB +
  beta_overall[["z_TC"]] * z_TC]
cum0 <- plogis(threshold[["0|1"]] - dat$Overall_eta)
cum1 <- plogis(threshold[["1|2"]] - dat$Overall_eta)
cum2 <- plogis(threshold[["2|3"]] - dat$Overall_eta)
expected_reproduced <- (cum1 - cum0) + 2 * (cum2 - cum1) + 3 * (1 - cum2)
max_reproduction_error <- max(abs(expected_reproduced - dat$Overall_expected_burden))
if (!is.finite(max_reproduction_error) || max_reproduction_error > 1e-10) {
  stop("Frozen Overall score did not reproduce: max error=", max_reproduction_error)
}

design <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~pooled_mec_weight,
  data = dat,
  nest = TRUE
)

linear_formula <- PhenoAge_acceleration ~ z_NM + z_TB + z_TC
quadratic_formula <- PhenoAge_acceleration ~ z_NM + z_TB + z_TC +
  I(z_NM^2) + I(z_TB^2) + I(z_TC^2) +
  z_NM:z_TB + z_NM:z_TC + z_TB:z_TC
linear_fit <- svyglm(linear_formula, design = design)
quadratic_fit <- svyglm(quadratic_formula, design = design)
nonlinear_test <- regTermTest(
  quadratic_fit,
  ~I(z_NM^2) + I(z_TB^2) + I(z_TC^2) + z_NM:z_TB + z_NM:z_TC + z_TB:z_TC,
  method = "Wald"
)
nonlinear_p <- as.numeric(nonlinear_test$p)
selected_model <- if (is.finite(nonlinear_p) && nonlinear_p < 0.05) "quadratic" else "linear"
selected_fit <- if (selected_model == "quadratic") quadratic_fit else linear_fit
selected_formula <- if (selected_model == "quadratic") quadratic_formula else linear_formula

weighted_mean <- function(x, w) sum(w * x) / sum(w)
weighted_correlation <- function(x, y, w) {
  mx <- weighted_mean(x, w)
  my <- weighted_mean(y, w)
  covariance <- weighted_mean((x - mx) * (y - my), w)
  covariance / sqrt(weighted_mean((x - mx)^2, w) * weighted_mean((y - my)^2, w))
}
weighted_median <- function(x, w) {
  ordering <- order(x)
  x <- x[ordering]
  w <- w[ordering]
  x[which(cumsum(w) >= 0.5 * sum(w))[1L]]
}

gradient_matrix <- function(coefficients, frame, model_type) {
  coefficient <- function(name) if (name %in% names(coefficients)) unname(coefficients[[name]]) else 0
  gradient <- matrix(NA_real_, nrow(frame), 3L, dimnames = list(NULL, c("z_NM", "z_TB", "z_TC")))
  gradient[, "z_NM"] <- coefficient("z_NM")
  gradient[, "z_TB"] <- coefficient("z_TB")
  gradient[, "z_TC"] <- coefficient("z_TC")
  if (model_type == "quadratic") {
    gradient[, "z_NM"] <- gradient[, "z_NM"] +
      2 * coefficient("I(z_NM^2)") * frame$z_NM +
      coefficient("z_NM:z_TB") * frame$z_TB +
      coefficient("z_NM:z_TC") * frame$z_TC
    gradient[, "z_TB"] <- gradient[, "z_TB"] +
      2 * coefficient("I(z_TB^2)") * frame$z_TB +
      coefficient("z_NM:z_TB") * frame$z_NM +
      coefficient("z_TB:z_TC") * frame$z_TC
    gradient[, "z_TC"] <- gradient[, "z_TC"] +
      2 * coefficient("I(z_TC^2)") * frame$z_TC +
      coefficient("z_NM:z_TC") * frame$z_NM +
      coefficient("z_TB:z_TC") * frame$z_TB
  }
  gradient
}

direction_metrics <- function(coefficients, frame, weights, model_type) {
  gradient <- gradient_matrix(coefficients, frame, model_type)
  overall_norm <- sqrt(sum(beta_overall^2))
  origin_gradient <- unname(gradient_matrix(
    coefficients,
    data.frame(z_NM = 0, z_TB = 0, z_TC = 0),
    model_type
  )[1, ])
  origin_gradient_unit <- origin_gradient / sqrt(sum(origin_gradient^2))
  gradient_norm <- sqrt(rowSums(gradient^2))
  cosine <- as.numeric(gradient %*% beta_overall) / (gradient_norm * overall_norm)
  cosine <- pmin(pmax(cosine, -1), 1)
  angle <- acos(cosine) * 180 / pi
  predictions <- as.numeric(model.matrix(selected_formula, data = frame) %*% coefficients)
  y_mean <- weighted_mean(frame$PhenoAge_acceleration, weights)
  r_squared <- 1 - sum(weights * (frame$PhenoAge_acceleration - predictions)^2) /
    sum(weights * (frame$PhenoAge_acceleration - y_mean)^2)
  c(
    weighted_mean_cosine = weighted_mean(cosine, weights),
    weighted_mean_angle_degrees = weighted_mean(angle, weights),
    weighted_median_angle_degrees = weighted_median(angle, weights),
    weighted_overall_pheno_correlation = weighted_correlation(
      frame$Overall_expected_burden, frame$PhenoAge_acceleration, weights
    ),
    weighted_model_r_squared = r_squared,
    gradient_NM_at_origin = origin_gradient[1],
    gradient_TB_at_origin = origin_gradient[2],
    gradient_TC_at_origin = origin_gradient[3],
    unit_gradient_NM_at_origin = origin_gradient_unit[1],
    unit_gradient_TB_at_origin = origin_gradient_unit[2],
    unit_gradient_TC_at_origin = origin_gradient_unit[3]
  )
}

# Rao-Wu bootstrap replicate weights preserve the complex survey design. The
# custom statistic refits the selected surface in every replicate and therefore
# propagates both coefficient and population-distribution uncertainty.
replicate_design <- as.svrepdesign(
  design,
  type = "bootstrap",
  replicates = 1000,
  mse = TRUE
)
replicate_result <- withReplicates(
  replicate_design,
  theta = function(weights, data) {
    data$.replicate_weight <- weights
    fit <- lm(selected_formula, data = data, weights = .replicate_weight)
    direction_metrics(coef(fit), data, weights, selected_model)
  },
  return.replicates = TRUE
)
metric_estimate <- coef(replicate_result)
metric_se <- SE(replicate_result)
critical <- qt(0.975, df = degf(design))
metric_table <- data.table(
  metric = names(metric_estimate),
  estimate = as.numeric(metric_estimate),
  standard_error = as.numeric(metric_se),
  ci_lower = as.numeric(metric_estimate - critical * metric_se),
  ci_upper = as.numeric(metric_estimate + critical * metric_se),
  survey_bootstrap_replicates = 1000L,
  design_degrees_of_freedom = degf(design)
)
metric_table[metric == "weighted_mean_cosine", `:=`(
  ci_lower = pmax(ci_lower, -1), ci_upper = pmin(ci_upper, 1)
)]
metric_table[metric == "weighted_overall_pheno_correlation", `:=`(
  ci_lower = pmax(ci_lower, -1), ci_upper = pmin(ci_upper, 1)
)]
metric_table[grepl("^unit_gradient_", metric), `:=`(
  ci_lower = pmax(ci_lower, -1), ci_upper = pmin(ci_upper, 1)
)]
fwrite(metric_table, file.path(geo_dir, "geometry_direction_metrics.csv"), bom = TRUE)

coefficient_vector <- coef(selected_fit)
coefficient_vcov <- vcov(selected_fit)
coefficient_table <- data.table(
  term = names(coefficient_vector),
  estimate = as.numeric(coefficient_vector),
  standard_error = sqrt(diag(coefficient_vcov)),
  ci_lower = as.numeric(coefficient_vector - critical * sqrt(diag(coefficient_vcov))),
  ci_upper = as.numeric(coefficient_vector + critical * sqrt(diag(coefficient_vcov))),
  selected_model = selected_model,
  nonlinear_block_wald_p = nonlinear_p
)
fwrite(coefficient_table, file.path(geo_dir, "geometry_pheno_surface_coefficients.csv"), bom = TRUE)

gradient <- gradient_matrix(coefficient_vector, dat, selected_model)
gradient_norm <- sqrt(rowSums(gradient^2))
cosine <- as.numeric(gradient %*% beta_overall) / (gradient_norm * sqrt(sum(beta_overall^2)))
dat[, `:=`(
  Pheno_surface_prediction = as.numeric(predict(selected_fit, newdata = dat)),
  Pheno_gradient_cosine = pmin(pmax(cosine, -1), 1),
  Pheno_gradient_angle_degrees = acos(pmin(pmax(cosine, -1), 1)) * 180 / pi
)]
fwrite(
  dat[, .(
    SEQN, CYCLE, survey_strata, survey_psu, pooled_mec_weight,
    Follow_Up_Years, Death_AllCause, Cancer_Diagnosed, Age, Sex, Race,
    z_NM, z_TB, z_TC, Overall_eta, Overall_expected_burden,
    PhenoAge_acceleration, Pheno_surface_prediction,
    Pheno_gradient_cosine, Pheno_gradient_angle_degrees
  )],
  file.path(geo_dir, "geometry_model_space_data.csv"),
  bom = TRUE
)

axis_limits <- rbindlist(lapply(paste0("z_", component_order), function(variable) {
  values <- dat[[variable]]
  data.table(
    axis = variable,
    lower_01 = unname(quantile(values, 0.01)),
    upper_99 = unname(quantile(values, 0.99)),
    lower_005 = unname(quantile(values, 0.005)),
    upper_995 = unname(quantile(values, 0.995))
  )
}))
fwrite(axis_limits, file.path(geo_dir, "geometry_axis_limits.csv"), bom = TRUE)

frozen_planes <- data.table(
  boundary = c("P(any domain) = 0.50", "P(at least two domains) = 0.50", "P(all three domains) = 0.50"),
  ordinal_threshold = c("0|1", "1|2", "2|3"),
  eta = unname(threshold),
  beta_z_NM = beta_overall[["z_NM"]],
  beta_z_TB = beta_overall[["z_TB"]],
  beta_z_TC = beta_overall[["z_TC"]]
)
fwrite(frozen_planes, file.path(geo_dir, "geometry_overall_iso_probability_planes.csv"), bom = TRUE)

qa <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  cohort_n = nrow(dat),
  unique_seqn_n = uniqueN(dat$SEQN),
  deaths_n = sum(dat$Death_AllCause == 1),
  design_degrees_of_freedom = degf(design),
  survey_bootstrap_replicates = 1000L,
  selected_pheno_surface = selected_model,
  nonlinear_block_wald_p = nonlinear_p,
  overall_score_max_reproduction_error = max_reproduction_error,
  overall_gradient = as.list(beta_overall),
  ordinal_thresholds = as.list(threshold),
  phenoage_used_in_overall_training = FALSE,
  nhanes_mortality_used_in_overall_training = FALSE,
  checks_passed = TRUE
)
write_json(qa, file.path(geo_dir, "geometry_surface_qa.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)

cat("Selected PhenoAge surface:", selected_model, "\n")
cat("Nonlinear block Wald P:", format(nonlin_p <- nonlinear_p, digits = 5), "\n")
cat("Overall reproduction max error:", format(max_reproduction_error, scientific = TRUE), "\n")
print(metric_table)
