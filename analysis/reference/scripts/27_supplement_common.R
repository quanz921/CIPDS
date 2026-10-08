source(file.path(Sys.getenv("CIPDS_PACKAGE_DIR"), "scripts", "weight_contract.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(survey)
  library(survival)
})

options(survey.lonely.psu = "adjust")

CIPDS_ROOT <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
CIPDS_OUT <- file.path(CIPDS_ROOT, "outputs", "supplement_v2")
CIPDS_FIG <- file.path(CIPDS_ROOT, "figures", "supplement_v2")
CIPDS_REPORT <- file.path(CIPDS_ROOT, "reports", "supplement_v2")
CIPDS_LOG <- file.path(CIPDS_ROOT, "logs")
dir.create(CIPDS_OUT, recursive = TRUE, showWarnings = FALSE)
dir.create(CIPDS_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(CIPDS_REPORT, recursive = TRUE, showWarnings = FALSE)
dir.create(CIPDS_LOG, recursive = TRUE, showWarnings = FALSE)

CIPDS_THREADS <- min(
  23L,
  max(1L, as.integer(Sys.getenv("CIPDS_THREADS", unset = "23")))
)
CIPDS_REPLICATES <- max(
  200L,
  as.integer(Sys.getenv("CIPDS_SUPPLEMENT_REPLICATES", unset = "1000"))
)
CIPDS_SEED <- as.integer(Sys.getenv("CIPDS_SUPPLEMENT_SEED", unset = "20260903"))
CIPDS_HORIZONS <- c(5, 10, 15)
CIPDS_CYCLE_ORDER <- c(
  "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
  "2009-2010", "2011-2012", "2013-2014", "2015-2016"
)
CIPDS_SCORE_VARS <- c(
  NM = "NM_z_age60", TB = "TB_z_age60", TC = "TC_z_age60",
  OVERALL = "OVERALL_z_age60", PHENO = "PHENO_z_age60"
)
CIPDS_SCORE_LABELS <- c(
  NM = "Nutritional-metabolic pattern",
  TB = "Tumor-burden-related pattern",
  TC = "Treatment-complication-related pattern",
  OVERALL = "Multidomain overall score",
  PHENO = "PhenoAge Acceleration"
)
CIPDS_CLINICAL_TERMS <- c(
  "Age", "Sex_F_design", "Race_F_design", "Education_F_design", "PIR", "BMI",
  "Smoking_F_design", "Alcohol_F_design", "Diabetes", "Hypertension", "CVD",
  "Cancer_F_design"
)
CIPDS_CLINICAL_RAW <- c(
  "Age", "Sex", "RIDRETH1", "Education_clean", "PIR", "BMI",
  "Smoking_F", "Alcohol_F", "Diabetes", "Hypertension", "CVD", "Cancer_Diagnosed"
)

cipds_load_data <- function(include_labs = FALSE) {
  env <- new.env(parent = emptyenv())
  load(file.path(CIPDS_ROOT, "inputs", "nhanes_phase7_batch1.RData"), envir = env)
  if (!exists("nhanes_valid", envir = env, inherits = FALSE)) stop("nhanes_valid missing")
  dat <- as.data.table(get("nhanes_valid", envir = env))
  if (nrow(dat) != 44772L || uniqueN(dat$SEQN) != 44772L) stop("NHANES grain drift")

  scores <- fread(file.path(CIPDS_ROOT, "outputs", "overall_nhanes_scores.csv"))
  score_raw <- c("NM_component", "TB_component", "TC_component", "Overall_expected_burden")
  if (nrow(scores) != 44772L || uniqueN(scores$SEQN) != 44772L) stop("Score grain drift")
  idx <- match(dat$SEQN, scores$SEQN)
  if (anyNA(idx)) stop("Score join failure")
  for (v in score_raw) dat[[v]] <- scores[[v]][idx]

  if (include_labs) {
    labs <- as.data.table(readRDS(file.path(CIPDS_ROOT, "outputs", "nhanes_expanded_candidate_matrix.rds")))
    if (nrow(labs) != 44772L || uniqueN(labs$SEQN) != 44772L) stop("Laboratory matrix grain drift")
    lab_cols <- setdiff(names(labs), c("SEQN", "CYCLE"))
    lab_idx <- match(dat$SEQN, labs$SEQN)
    if (anyNA(lab_idx)) stop("Laboratory join failure")
    for (v in lab_cols) dat[[v]] <- labs[[v]][lab_idx]
  }

  required <- c(
    "SEQN", "CYCLE", "WTMEC2YR", "SDMVSTRA", "SDMVPSU", "Follow_Up_Years",
    "Dead", "UCOD_LEADING", "Age", "Sex", "RIDRETH1", "Cancer_Diagnosed",
    "Education", "PIR", "BMI", "Smoking_F", "Alcohol_F", "Diabetes",
    "Hypertension", "CVD", "PhenoAge", score_raw
  )
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required)) stop("Missing variables: ", paste(missing_required, collapse = ", "))

  dat[, CYCLE := factor(CYCLE, levels = CIPDS_CYCLE_ORDER, ordered = FALSE)]
  dat[, pooled_mec_weight := cipds_pooled_mec_weights(dat, 9)]
  dat[, survey_strata := interaction(CYCLE, SDMVSTRA, drop = TRUE, lex.order = TRUE)]
  dat[, survey_psu := interaction(CYCLE, SDMVSTRA, SDMVPSU, drop = TRUE, lex.order = TRUE)]
  dat[, Death_AllCause := as.integer(Dead == 1)]
  dat[, Death_Cancer := as.integer(Dead == 1 & UCOD_LEADING == 2)]
  dat[, Death_CVD := fifelse(CYCLE == "2015-2016", NA_integer_, as.integer(Dead == 1 & UCOD_LEADING %in% c(1, 5)))]
  dat[, competing_cancer := fifelse(Dead == 0, 0L, fifelse(UCOD_LEADING == 2, 1L, 2L))]
  dat[, competing_cvd := fifelse(CYCLE == "2015-2016", NA_integer_, fifelse(Dead == 0, 0L, fifelse(UCOD_LEADING %in% c(1, 5), 1L, 2L)))]
  dat[, Sex_F_design := factor(Sex)]
  dat[, Race_F_design := factor(RIDRETH1)]
  dat[, Education_clean := ifelse(Education %in% 1:5, Education, NA_real_)]
  dat[, Education_F_design := factor(Education_clean, levels = 1:5)]
  dat[, Smoking_F_design := factor(Smoking_F)]
  dat[, Alcohol_F_design := factor(Alcohol_F)]
  dat[, Cancer_F_design := factor(Cancer_Diagnosed, levels = c(0, 1))]

  design_full <- svydesign(
    ids = ~survey_psu,
    strata = ~survey_strata,
    weights = ~pooled_mec_weight,
    data = dat,
    nest = TRUE
  )
  age60 <- !is.na(dat$Age) & dat$Age >= 60
  design_age60 <- design_full[age60, ]

  raw_scores <- c(
    NM = "NM_component", TB = "TB_component", TC = "TC_component",
    OVERALL = "Overall_expected_burden"
  )
  for (short in names(raw_scores)) {
    raw <- raw_scores[[short]]
    mu <- as.numeric(coef(svymean(as.formula(paste0("~", raw)), design_age60, na.rm = TRUE)))
    sigma <- sqrt(as.numeric(svyvar(as.formula(paste0("~", raw)), design_age60, na.rm = TRUE)))
    if (!is.finite(sigma) || sigma <= 0) stop("Invalid score SD: ", raw)
    dat[[paste0(short, "_z_age60")]] <- (dat[[raw]] - mu) / sigma
  }
  pheno_fit <- svyglm(PhenoAge ~ Age, design = design_age60)
  pheno_coef <- coef(pheno_fit)
  dat[, PhenoAge_acceleration := PhenoAge -
        (pheno_coef[["(Intercept)"]] + pheno_coef[["Age"]] * Age)]
  design_pheno <- svydesign(
    ids = ~survey_psu, strata = ~survey_strata, weights = ~pooled_mec_weight,
    data = dat, nest = TRUE
  )[age60, ]
  pheno_mu <- as.numeric(coef(svymean(~PhenoAge_acceleration, design_pheno, na.rm = TRUE)))
  pheno_sd <- sqrt(as.numeric(svyvar(~PhenoAge_acceleration, design_pheno, na.rm = TRUE)))
  dat[, PHENO_z_age60 := (PhenoAge_acceleration - pheno_mu) / pheno_sd]

  design_full <- svydesign(
    ids = ~survey_psu, strata = ~survey_strata, weights = ~pooled_mec_weight,
    data = dat, nest = TRUE
  )
  list(dat = dat, design_full = design_full, age60 = age60)
}

cipds_complete_mask <- function(dat, extra = character(), outcome = "Death_AllCause") {
  required <- unique(c(
    "Follow_Up_Years", outcome, "pooled_mec_weight", CIPDS_CLINICAL_TERMS, extra
  ))
  !is.na(dat$Age) & dat$Age >= 60 &
    complete.cases(dat[, ..required]) &
    is.finite(dat$pooled_mec_weight) & dat$pooled_mec_weight > 0 &
    is.finite(dat$Follow_Up_Years) & dat$Follow_Up_Years > 0
}

cipds_make_formula <- function(predictors, outcome = "Death_AllCause") {
  as.formula(paste(
    "Surv(Follow_Up_Years,", outcome, ") ~",
    paste(c(predictors, CIPDS_CLINICAL_TERMS), collapse = " + ")
  ))
}

cipds_weighted_quantile <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]
  w <- w[keep]
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

cipds_weighted_km_censoring <- function(time, event, weight) {
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
  list(
    time = unique_time,
    before = c(1, head(cumprod(1 - hazard), -1)),
    after = cumprod(1 - hazard)
  )
}

cipds_weighted_td_auc <- function(time, event, score, weight, horizon) {
  keep <- is.finite(time) & !is.na(event) & is.finite(score) & is.finite(weight) & weight > 0
  time <- time[keep]
  event <- event[keep]
  score <- score[keep]
  weight <- weight[keep]
  cases <- event == 1 & time < horizon
  controls <- time > horizon
  if (sum(cases) < 2L || sum(controls) < 2L) return(NA_real_)
  g <- cipds_weighted_km_censoring(time, event, weight)
  case_time_index <- match(time[cases], g$time)
  g_case <- pmax(g$before[case_time_index], 1e-8)
  horizon_index <- findInterval(horizon, g$time)
  g_horizon <- if (horizon_index == 0L) 1 else g$after[horizon_index]
  g_horizon <- max(g_horizon, 1e-8)
  case_weight <- weight[cases] / g_case
  control_weight <- weight[controls] / g_horizon
  case_score <- score[cases]
  control_groups <- data.table(score = score[controls], weight = control_weight)[
    , .(weight = sum(weight)), by = score
  ][order(score)]
  cumulative_control_weight <- cumsum(control_groups$weight)
  le_index <- findInterval(case_score, control_groups$score)
  le_weight <- ifelse(le_index > 0L, cumulative_control_weight[pmax(le_index, 1L)], 0)
  equal_index <- match(case_score, control_groups$score, nomatch = 0L)
  equal_weight <- ifelse(equal_index > 0L, control_groups$weight[pmax(equal_index, 1L)], 0)
  denominator <- sum(case_weight) * sum(control_weight)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(case_weight * (le_weight - 0.5 * equal_weight)) / denominator
}

cipds_weighted_uno <- function(time, event, score, weight, max_time = 15) {
  keep <- is.finite(time) & !is.na(event) & is.finite(score) & is.finite(weight) & weight > 0
  if (sum(event[keep] == 1) < 20L) return(NA_real_)
  as.numeric(concordancefit(
    y = Surv(time[keep], event[keep]), x = score[keep], weights = weight[keep],
    ymin = 0, ymax = max_time, timewt = "n/G2", reverse = TRUE,
    std.err = FALSE
  )$concordance)
}

cipds_weighted_brier <- function(time, event, risk, weight, horizon) {
  keep <- is.finite(time) & !is.na(event) & is.finite(risk) & is.finite(weight) & weight > 0
  time <- time[keep]
  event <- event[keep]
  risk <- pmin(pmax(risk[keep], 0), 1)
  weight <- weight[keep]
  g <- cipds_weighted_km_censoring(time, event, weight)
  cases <- event == 1 & time <= horizon
  controls <- time > horizon
  if (sum(cases) < 2L || sum(controls) < 2L) return(NA_real_)
  g_case <- pmax(g$before[match(time[cases], g$time)], 1e-8)
  h_idx <- findInterval(horizon, g$time)
  g_h <- max(if (h_idx == 0L) 1 else g$after[h_idx], 1e-8)
  case_w <- weight[cases] / g_case
  control_w <- weight[controls] / g_h
  sum(case_w * (1 - risk[cases])^2) + sum(control_w * risk[controls]^2) -> numerator
  numerator / (sum(case_w) + sum(control_w))
}

cipds_baseline_hazard <- function(fit, horizons) {
  bh <- basehaz(fit, centered = TRUE)
  vapply(horizons, function(h) {
    idx <- findInterval(h, bh$time)
    if (idx == 0L) 0 else bh$hazard[idx]
  }, numeric(1))
}

cipds_predict_risk <- function(fit, newdata, horizons) {
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "sample"))
  h0 <- cipds_baseline_hazard(fit, horizons)
  risk <- sapply(h0, function(v) 1 - exp(-v * exp(lp)))
  colnames(risk) <- paste0("risk", horizons)
  list(lp = lp, risk = risk)
}

cipds_metric_vector <- function(time, event, lp, risk, weight, horizons = CIPDS_HORIZONS) {
  ans <- c(uno_c15 = cipds_weighted_uno(time, event, lp, weight, max(horizons)))
  for (j in seq_along(horizons)) {
    h <- horizons[j]
    ans[paste0("auc", h)] <- cipds_weighted_td_auc(time, event, lp, weight, h)
    ans[paste0("brier", h)] <- cipds_weighted_brier(time, event, risk[, j], weight, h)
  }
  grid <- seq(1, max(horizons), by = 1)
  if (nrow(risk) && length(grid) > 1L) {
    ans["ibs15"] <- NA_real_
  }
  ans
}

cipds_svr_summary <- function(point, replicate_matrix, rep_design, df, bounds = NULL) {
  v <- svrVar(
    replicate_matrix,
    scale = rep_design$scale,
    rscales = rep_design$rscales,
    mse = rep_design$mse,
    coef = point
  )
  se <- sqrt(diag(v))
  crit <- qt(0.975, df = df)
  out <- data.table(
    key = names(point), estimate = as.numeric(point), standard_error = as.numeric(se),
    ci_lower = as.numeric(point) - crit * as.numeric(se),
    ci_upper = as.numeric(point) + crit * as.numeric(se)
  )
  if (!is.null(bounds)) {
    out[, ci_lower := pmax(bounds[1], ci_lower)]
    out[, ci_upper := pmin(bounds[2], ci_upper)]
  }
  out
}

cipds_design_replicates <- function(design_full, mask, replicates = CIPDS_REPLICATES, seed = CIPDS_SEED) {
  set.seed(seed)
  rep_full <- as.svrepdesign(
    design_full, type = "bootstrap", replicates = replicates, mse = TRUE
  )
  list(
    design = design_full[mask, ],
    replicate_design = rep_full[mask, ],
    weights = weights(rep_full[mask, ], type = "analysis")
  )
}

cipds_write_manifest <- function(path, payload) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite required")
  jsonlite::write_json(
    payload, path, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null"
  )
}
