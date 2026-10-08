source(file.path(Sys.getenv("CIPDS_PACKAGE_DIR"), "scripts", "weight_contract.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(parallel)
  library(survey)
  library(survival)
  library(timeROC)
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

n_replicates <- as.integer(Sys.getenv("CIPDS_DISCRIMINATION_REPLICATES", unset = "1000"))
n_threads <- as.integer(Sys.getenv("CIPDS_THREADS", unset = "23"))
seed <- as.integer(Sys.getenv("CIPDS_DISCRIMINATION_SEED", unset = "20260831"))
horizons <- c(5, 10, 15)
stopifnot(n_replicates >= 200L, n_threads >= 1L)

log_file <- file.path(log_dir, "19_nhanes_paired_discrimination.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Survey bootstrap replicates:", n_replicates, "\n")
cat("Parallel workers:", n_threads, "\n")
cat("Pre-specified horizons:", paste(horizons, collapse = ", "), "years\n")

load(file.path(package_dir, "inputs", "nhanes_phase7_batch1.RData"))
if (!exists("nhanes_valid")) stop("nhanes_valid missing")
dat <- as.data.table(nhanes_valid)
if (nrow(dat) != 44772L || uniqueN(dat$SEQN) != 44772L) stop("NHANES grain drift")

scores <- fread(file.path(out_dir, "overall_nhanes_scores.csv"))
score_cols <- c(
  "NM_component", "TB_component", "TC_component", "Overall_expected_burden"
)
if (nrow(scores) != 44772L || uniqueN(scores$SEQN) != 44772L) stop("Score grain drift")
score_idx <- match(dat$SEQN, scores$SEQN)
if (anyNA(score_idx)) stop("Score join failure")
for (v in score_cols) dat[[v]] <- scores[[v]][score_idx]

required <- c(
  "SEQN", "CYCLE", "WTMEC2YR", "SDMVSTRA", "SDMVPSU", "Follow_Up_Years",
  "Dead", "Age", "Cancer_Diagnosed", "PhenoAge", score_cols
)
missing_required <- setdiff(required, names(dat))
if (length(missing_required)) stop("Missing required variables: ", paste(missing_required, collapse = ", "))

dat[, pooled_mec_weight := cipds_pooled_mec_weights(dat, 9)]
dat[, survey_strata := interaction(CYCLE, SDMVSTRA, drop = TRUE, lex.order = TRUE)]
dat[, survey_psu := interaction(CYCLE, SDMVSTRA, SDMVPSU, drop = TRUE, lex.order = TRUE)]
dat[, Death_AllCause := as.integer(Dead == 1)]

# The complete nine-cycle design and its bootstrap replicate design are created
# before restricting to an age domain or cancer-history subgroup.
design_full <- svydesign(
  ids = ~survey_psu,
  strata = ~survey_strata,
  weights = ~pooled_mec_weight,
  data = dat,
  nest = TRUE
)
set.seed(seed)
rep_design_full <- as.svrepdesign(
  design_full,
  type = "bootstrap",
  replicates = n_replicates,
  mse = TRUE
)

domains <- list(
  age60 = list(age = 60, label = "Age >=60 years", role = "PRIMARY"),
  age65 = list(age = 65, label = "Age >=65 years", role = "SENSITIVITY")
)
models <- c(
  OVERALL = "Overall_expected_burden",
  PHENO = "PhenoAge_acceleration",
  NM = "NM_component",
  TB = "TB_component",
  TC = "TC_component"
)
model_labels <- c(
  OVERALL = "Multidomain laboratory vulnerability overall score",
  PHENO = "PhenoAge Acceleration",
  NM = "Nutritional-Metabolic laboratory pattern",
  TB = "Tumor-burden-related laboratory pattern",
  TC = "Treatment-complication-related laboratory pattern"
)

weighted_km_censoring <- function(time, event, weight) {
  keep <- is.finite(time) & !is.na(event) & is.finite(weight) & weight > 0
  time <- time[keep]
  event <- event[keep]
  weight <- weight[keep]
  ord <- order(time)
  time <- time[ord]
  event <- event[ord]
  weight <- weight[ord]
  unique_time <- sort(unique(time))
  group_index <- match(time, unique_time)
  total_weight <- as.numeric(rowsum(weight, group_index, reorder = FALSE))
  censor_weight <- as.numeric(rowsum(weight * (event == 0), group_index, reorder = FALSE))
  risk_weight <- rev(cumsum(rev(total_weight)))
  hazard <- ifelse(risk_weight > 0, censor_weight / risk_weight, 0)
  hazard <- pmin(pmax(hazard, 0), 1)
  survival_before <- c(1, head(cumprod(1 - hazard), -1))
  survival_after <- cumprod(1 - hazard)
  list(
    time = unique_time,
    before = survival_before,
    after = survival_after
  )
}

weighted_td_auc <- function(time, event, score, weight, horizon) {
  keep <- is.finite(time) & !is.na(event) & is.finite(score) & is.finite(weight) & weight > 0
  time <- time[keep]
  event <- event[keep]
  score <- score[keep]
  weight <- weight[keep]
  # Match timeROC's cumulative/dynamic convention: cases fail strictly before
  # the horizon and controls remain event-free strictly beyond the horizon.
  cases <- event == 1 & time < horizon
  controls <- time > horizon
  if (sum(cases) < 2L || sum(controls) < 2L) return(NA_real_)

  g <- weighted_km_censoring(time, event, weight)
  case_time_index <- match(time[cases], g$time)
  g_case <- pmax(g$before[case_time_index], 1e-8)
  horizon_index <- findInterval(horizon, g$time)
  g_horizon <- if (horizon_index == 0L) 1 else g$after[horizon_index]
  g_horizon <- max(g_horizon, 1e-8)

  case_weight <- weight[cases] / g_case
  control_weight <- weight[controls] / g_horizon
  case_score <- score[cases]
  control_score <- score[controls]

  control_groups <- data.table(score = control_score, weight = control_weight)[
    , .(weight = sum(weight)), by = score
  ][order(score)]
  cumulative_control_weight <- cumsum(control_groups$weight)
  less_or_equal_index <- findInterval(case_score, control_groups$score)
  less_or_equal_weight <- ifelse(
    less_or_equal_index > 0L,
    cumulative_control_weight[pmax(less_or_equal_index, 1L)],
    0
  )
  equal_index <- match(case_score, control_groups$score, nomatch = 0L)
  equal_weight <- ifelse(equal_index > 0L, control_groups$weight[pmax(equal_index, 1L)], 0)
  concordant_control_weight <- less_or_equal_weight - 0.5 * equal_weight
  denominator <- sum(case_weight) * sum(control_weight)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(case_weight * concordant_control_weight) / denominator
}

metric_bundle <- function(time, event, score_matrix, weight, horizons) {
  keep <- is.finite(time) & !is.na(event) & is.finite(weight) & weight > 0 &
    rowSums(is.finite(score_matrix)) == ncol(score_matrix)
  time <- time[keep]
  event <- event[keep]
  score_matrix <- score_matrix[keep, , drop = FALSE]
  weight <- weight[keep]
  if (sum(event == 1) < 20L) stop("Insufficient events in a replicate")

  y <- Surv(time, event)
  uno <- concordancefit(
    y = y,
    x = score_matrix,
    weights = weight,
    ymin = 0,
    ymax = max(horizons),
    timewt = "n/G2",
    reverse = TRUE,
    std.err = TRUE
  )$concordance
  harrell <- concordancefit(
    y = y,
    x = score_matrix,
    weights = weight,
    timewt = "n",
    reverse = TRUE,
    std.err = TRUE
  )$concordance
  names(uno) <- colnames(score_matrix)
  names(harrell) <- colnames(score_matrix)

  values <- c()
  for (model in colnames(score_matrix)) {
    values[paste0("uno_c15__", model)] <- uno[[model]]
    values[paste0("harrell_c__", model)] <- harrell[[model]]
    for (horizon in horizons) {
      values[paste0("auc", horizon, "__", model)] <- weighted_td_auc(
        time, event, score_matrix[, model], weight, horizon
      )
    }
  }
  values
}

metric_metadata <- function(metric_key) {
  parts <- strsplit(metric_key, "__", fixed = TRUE)[[1]]
  metric_code <- parts[1]
  model <- parts[2]
  horizon <- if (grepl("^auc", metric_code)) as.numeric(sub("^auc", "", metric_code)) else {
    if (metric_code == "uno_c15") 15 else NA_real_
  }
  metric_label <- if (metric_code == "uno_c15") {
    "Survey-weighted Uno C-index truncated at 15 years"
  } else if (metric_code == "harrell_c") {
    "Survey-weighted Harrell C-index over observed follow-up"
  } else {
    paste0("Survey-weighted cumulative/dynamic AUC at ", horizon, " years")
  }
  metric_role <- if (metric_code == "uno_c15") {
    "PRIMARY"
  } else if (grepl("^auc", metric_code)) {
    "KEY_SECONDARY"
  } else {
    "SENSITIVITY"
  }
  list(
    metric_code = metric_code,
    metric_label = metric_label,
    metric_role = metric_role,
    horizon_years = horizon,
    model = model
  )
}

estimate_rows <- list()
comparison_rows <- list()
cohort_rows <- list()
replicate_registry <- list()
crosscheck_rows <- list()
primary_cohort_written <- FALSE

n_threads <- min(n_threads, n_replicates)
cluster <- makeCluster(n_threads)
on.exit(stopCluster(cluster), add = TRUE)
clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
  })
  NULL
})
clusterExport(
  cluster,
  c("weighted_km_censoring", "weighted_td_auc", "metric_bundle", "horizons"),
  envir = environment()
)

for (domain_name in names(domains)) {
  domain <- domains[[domain_name]]
  domain_mask <- !is.na(dat$Age) & dat$Age >= domain$age
  domain_design <- design_full[domain_mask, ]

  # Match the locked definition used by the existing NHANES analysis:
  # survey-weighted residual from PhenoAge ~ chronological age within domain.
  pheno_fit <- svyglm(PhenoAge ~ Age, design = domain_design)
  pheno_coef <- coef(pheno_fit)
  dat[, PhenoAge_acceleration := PhenoAge -
        (pheno_coef[["(Intercept)"]] + pheno_coef[["Age"]] * Age)]

  populations <- list(
    Overall = domain_mask,
    Cancer = domain_mask & !is.na(dat$Cancer_Diagnosed) & dat$Cancer_Diagnosed == 1,
    Noncancer = domain_mask & !is.na(dat$Cancer_Diagnosed) & dat$Cancer_Diagnosed == 0
  )

  for (population_name in names(populations)) {
    required_metric_vars <- c(
      "Follow_Up_Years", "Death_AllCause", "pooled_mec_weight", "survey_strata",
      "survey_psu", unname(models)
    )
    complete <- populations[[population_name]] &
      complete.cases(dat[, ..required_metric_vars]) &
      dat$pooled_mec_weight > 0 & dat$Follow_Up_Years > 0
    n <- sum(complete)
    events <- sum(dat$Death_AllCause[complete] == 1)
    if (n < 500L || events < 100L) stop("Insufficient common cohort: ", domain_name, "/", population_name)

    sub_dat <- dat[complete]
    score_matrix <- as.matrix(sub_dat[, ..models])
    colnames(score_matrix) <- names(models)
    time <- sub_dat$Follow_Up_Years
    event <- sub_dat$Death_AllCause
    sampling_weight <- sub_dat$pooled_mec_weight

    design_sub <- design_full[complete, ]
    rep_design_sub <- rep_design_full[complete, ]
    replicate_weight <- weights(rep_design_sub, type = "analysis")
    if (nrow(replicate_weight) != n || ncol(replicate_weight) != n_replicates) {
      stop("Replicate-weight dimension drift")
    }

    cases_by_horizon <- vapply(horizons, function(h) sum(event == 1 & time < h), integer(1))
    controls_by_horizon <- vapply(horizons, function(h) sum(time > h), integer(1))
    cohort_rows[[length(cohort_rows) + 1L]] <- data.table(
      domain = domain_name,
      domain_label = domain$label,
      analysis_role = ifelse(
        domain_name == "age60" && population_name == "Overall",
        "PRIMARY",
        ifelse(domain_name == "age60", "SUBGROUP", "SENSITIVITY")
      ),
      population = population_name,
      n = n,
      events = events,
      weighted_population = sum(sampling_weight),
      design_degrees_of_freedom = degf(design_sub),
      cases_5y = cases_by_horizon[1],
      controls_5y = controls_by_horizon[1],
      cases_10y = cases_by_horizon[2],
      controls_10y = controls_by_horizon[2],
      cases_15y = cases_by_horizon[3],
      controls_15y = controls_by_horizon[3],
      common_complete_case = TRUE
    )

    if (!primary_cohort_written && domain_name == "age60" && population_name == "Overall") {
      fwrite(
        sub_dat[, .(
          SEQN, CYCLE, survey_strata, survey_psu, pooled_mec_weight,
          Follow_Up_Years, Death_AllCause, Cancer_Diagnosed, Age,
          Overall_expected_burden, PhenoAge_acceleration,
          NM_component, TB_component, TC_component
        )],
        file.path(out_dir, "nhanes_paired_discrimination_primary_cohort.csv"),
        bom = TRUE
      )
      primary_cohort_written <- TRUE
    }

    point <- metric_bundle(time, event, score_matrix, sampling_weight, horizons)
    cat("\n", domain_name, population_name, "n/events", n, events, "\n")
    print(round(point, 5))

    # Unweighted marginal-IPCW cross-check against timeROC on the primary cohort.
    if (domain_name == "age60" && population_name == "Overall") {
      for (model in colnames(score_matrix)) {
        tr <- timeROC(
          T = time,
          delta = event,
          marker = score_matrix[, model],
          cause = 1,
          weighting = "marginal",
          times = horizons,
          iid = FALSE
        )
        ours_unweighted <- vapply(
          horizons,
          function(h) weighted_td_auc(time, event, score_matrix[, model], rep(1, n), h),
          numeric(1)
        )
        crosscheck_rows[[length(crosscheck_rows) + 1L]] <- data.table(
          model = model,
          horizon_years = horizons,
          custom_unweighted_auc = ours_unweighted,
          timeROC_unweighted_auc = as.numeric(tr$AUC),
          absolute_difference = abs(ours_unweighted - as.numeric(tr$AUC))
        )
      }
    }

    # Parallel paired survey-bootstrap calculation. Every model in a replicate
    # uses the same participants and the same replicate weights.
    clusterExport(
      cluster,
      c("time", "event", "score_matrix", "replicate_weight"),
      envir = environment()
    )
    replicate_list <- parLapplyLB(cluster, seq_len(n_replicates), function(b) {
      metric_bundle(time, event, score_matrix, replicate_weight[, b], horizons)
    })
    replicate_matrix <- do.call(rbind, replicate_list)
    if (!identical(colnames(replicate_matrix), names(point))) stop("Replicate metric-name drift")
    if (any(!is.finite(replicate_matrix))) stop("Non-finite bootstrap metric")

    variance <- svrVar(
      replicate_matrix,
      scale = rep_design_sub$scale,
      rscales = rep_design_sub$rscales,
      mse = rep_design_sub$mse,
      coef = point
    )
    standard_error <- sqrt(diag(variance))
    names(standard_error) <- names(point)
    df <- degf(design_sub)
    critical <- qt(0.975, df = df)

    for (metric_key in names(point)) {
      meta <- metric_metadata(metric_key)
      estimate_rows[[length(estimate_rows) + 1L]] <- data.table(
        scenario = "A_PRIMARY_FROZEN_PANCANCER_OVERALL_STACK",
        domain = domain_name,
        domain_label = domain$label,
        analysis_role = ifelse(
          domain_name == "age60" && population_name == "Overall",
          meta$metric_role,
          ifelse(domain_name == "age60", "SUBGROUP", "SENSITIVITY")
        ),
        population = population_name,
        metric = meta$metric_code,
        metric_label = meta$metric_label,
        metric_role = meta$metric_role,
        horizon_years = meta$horizon_years,
        model = meta$model,
        model_label = model_labels[[meta$model]],
        n = n,
        events = events,
        estimate = point[[metric_key]],
        standard_error = standard_error[[metric_key]],
        ci_lower = max(0, point[[metric_key]] - critical * standard_error[[metric_key]]),
        ci_upper = min(1, point[[metric_key]] + critical * standard_error[[metric_key]]),
        design_degrees_of_freedom = df,
        survey_bootstrap_replicates = n_replicates,
        common_complete_case = TRUE
      )
    }

    comparison_pairs <- rbind(
      data.table(
        model_a = "OVERALL", model_b = "PHENO",
        comparison_family = "PRIMARY_OVERALL_VS_PHENO"
      ),
      data.table(
        model_a = c("NM", "TB", "TC"), model_b = "PHENO",
        comparison_family = "SECONDARY_COMPONENT_VS_PHENO"
      ),
      data.table(
        model_a = "OVERALL", model_b = c("NM", "TB", "TC"),
        comparison_family = "SECONDARY_OVERALL_VS_COMPONENT"
      )
    )
    metric_codes <- unique(vapply(names(point), function(x) metric_metadata(x)$metric_code, character(1)))
    for (metric_code in metric_codes) {
      metric_meta <- metric_metadata(paste0(metric_code, "__OVERALL"))
      for (pair_index in seq_len(nrow(comparison_pairs))) {
        pair <- comparison_pairs[pair_index]
        key_a <- paste0(metric_code, "__", pair$model_a)
        key_b <- paste0(metric_code, "__", pair$model_b)
        difference <- point[[key_a]] - point[[key_b]]
        difference_replicates <- replicate_matrix[, key_a] - replicate_matrix[, key_b]
        difference_variance <- as.numeric(svrVar(
          matrix(difference_replicates, ncol = 1),
          scale = rep_design_sub$scale,
          rscales = rep_design_sub$rscales,
          mse = rep_design_sub$mse,
          coef = difference
        ))
        difference_se <- sqrt(difference_variance)
        statistic <- difference / difference_se
        p_value <- 2 * pt(-abs(statistic), df = df)
        comparison_rows[[length(comparison_rows) + 1L]] <- data.table(
          scenario = "A_PRIMARY_FROZEN_PANCANCER_OVERALL_STACK",
          domain = domain_name,
          domain_label = domain$label,
          population = population_name,
          metric = metric_code,
          metric_label = metric_meta$metric_label,
          metric_role = metric_meta$metric_role,
          horizon_years = metric_meta$horizon_years,
          comparison_family = pair$comparison_family,
          model_a = pair$model_a,
          model_a_label = model_labels[[pair$model_a]],
          model_b = pair$model_b,
          model_b_label = model_labels[[pair$model_b]],
          difference_orientation = paste0(pair$model_a, " minus ", pair$model_b),
          n = n,
          events = events,
          estimate_a = point[[key_a]],
          estimate_b = point[[key_b]],
          paired_difference = difference,
          difference_standard_error = difference_se,
          difference_ci_lower = difference - critical * difference_se,
          difference_ci_upper = difference + critical * difference_se,
          paired_t_statistic = statistic,
          paired_p_value = p_value,
          design_degrees_of_freedom = df,
          survey_bootstrap_replicates = n_replicates,
          common_complete_case = TRUE
        )
      }
    }

    replicate_registry[[paste(domain_name, population_name, sep = "__")]] <- list(
      domain = domain_name,
      population = population_name,
      point = point,
      replicate_estimates = replicate_matrix,
      scale = rep_design_sub$scale,
      rscales = rep_design_sub$rscales,
      mse = rep_design_sub$mse,
      design_degrees_of_freedom = df
    )
  }
}

estimate_dt <- rbindlist(estimate_rows, fill = TRUE)
comparison_dt <- rbindlist(comparison_rows, fill = TRUE)
cohort_dt <- rbindlist(cohort_rows, fill = TRUE)
crosscheck_dt <- rbindlist(crosscheck_rows, fill = TRUE)

# The Uno C-index is the single primary metric. The three time-dependent AUC
# horizons are multiplicity-controlled with Holm correction within population.
comparison_dt[, paired_p_holm := paired_p_value]
comparison_dt[
  comparison_family == "PRIMARY_OVERALL_VS_PHENO" & metric %in% paste0("auc", horizons),
  paired_p_holm := p.adjust(paired_p_value, method = "holm"),
  by = .(domain, population, comparison_family)
]
comparison_dt[
  comparison_family != "PRIMARY_OVERALL_VS_PHENO",
  paired_p_fdr_bh := p.adjust(paired_p_value, method = "BH"),
  by = .(domain, population, comparison_family)
]

comparison_dt[, conclusion := fifelse(
  difference_ci_lower > 0,
  paste0(model_a, " higher discrimination"),
  fifelse(difference_ci_upper < 0, paste0(model_b, " higher discrimination"), "No statistically resolved difference")
)]

setorder(estimate_dt, domain, population, metric_role, metric, model)
setorder(comparison_dt, domain, population, comparison_family, metric)
setorder(cohort_dt, domain, population)
setorder(crosscheck_dt, model, horizon_years)

fwrite(estimate_dt, file.path(out_dir, "nhanes_paired_discrimination_estimates.csv"), bom = TRUE)
fwrite(comparison_dt, file.path(out_dir, "nhanes_paired_discrimination_comparisons.csv"), bom = TRUE)
fwrite(cohort_dt, file.path(out_dir, "nhanes_paired_discrimination_cohort_audit.csv"), bom = TRUE)
fwrite(crosscheck_dt, file.path(out_dir, "nhanes_paired_discrimination_timeROC_crosscheck.csv"), bom = TRUE)
saveRDS(replicate_registry, file.path(out_dir, "nhanes_paired_discrimination_replicates.rds"), compress = "xz")

manifest <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  scenario = "A_PRIMARY_FROZEN_PANCANCER_OVERALL_STACK",
  source_scope = "NHANES 1999-2016, nine cycles, mortality-linked release",
  primary_population = "Age >=60 years, overall population",
  subgroups = c("cancer history", "no cancer history"),
  sensitivity_domain = "Age >=65 years",
  outcome = "all-cause mortality",
  common_complete_case = TRUE,
  primary_metric = "survey-weighted Uno C-index truncated at 15 years",
  key_secondary_metrics = paste0("survey-weighted cumulative/dynamic AUC at ", horizons, " years"),
  sensitivity_metric = "survey-weighted Harrell C-index over observed follow-up",
  primary_comparison = "OVERALL minus PHENO",
  score_models_retrained = FALSE,
  phenoage_definition = "survey-weighted residual from PhenoAge ~ chronological age within age domain",
  full_survey_design_created_before_domain_restriction = TRUE,
  pooled_weight = CIPDS_WEIGHT_DESCRIPTION,
  bootstrap_type = "survey::as.svrepdesign type=bootstrap, MSE variance",
  bootstrap_replicates = n_replicates,
  threads = n_threads,
  random_seed = seed,
  time_dependent_auc_definition = paste(
    "cumulative cases and dynamic controls with marginal IPCW for censoring;",
    "survey sampling weights enter the weighted case-control concordance probability"
  ),
  multiplicity = list(
    primary_uno_c = "single primary comparison; unadjusted paired P value",
    time_auc = "Holm correction across 5, 10 and 15 years within domain/population",
    component_comparisons = "Benjamini-Hochberg FDR within comparison family and domain/population"
  ),
  limitations = c(
    "Discrimination does not establish calibration or clinical utility",
    "No claim of universal superiority is allowed if metrics or horizons disagree",
    "The direct score comparison evaluates frozen markers; clinical-adjusted association remains reported separately"
  )
)
write_json(
  manifest,
  file.path(out_dir, "nhanes_paired_discrimination_manifest.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("\nPrimary Overall versus PhenoAge results:\n")
print(comparison_dt[
  domain == "age60" & population == "Overall" &
    comparison_family == "PRIMARY_OVERALL_VS_PHENO",
  .(metric, estimate_a, estimate_b, paired_difference, difference_ci_lower,
    difference_ci_upper, paired_p_value, paired_p_holm, conclusion)
])
cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
