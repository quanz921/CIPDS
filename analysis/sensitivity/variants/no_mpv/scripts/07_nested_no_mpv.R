options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(catboost)
  library(data.table)
  library(pROC)
  library(jsonlite)
  library(parallel)
})

package_dir <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
out_dir <- file.path(package_dir, "outputs")
log_dir <- file.path(package_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

scenario <- Sys.getenv("CIPDS_CANDIDATE_SCENARIO", unset = "A_PRIMARY")
if (!scenario %in% c("A_PRIMARY", "AB_SENSITIVITY")) {
  stop("CIPDS_CANDIDATE_SCENARIO must be A_PRIMARY or AB_SENSITIVITY")
}
scenario_slug <- tolower(scenario)
log_file <- file.path(log_dir, paste0("07_nested_5x5_boruta_catboost_", scenario_slug, ".log"))
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

THREADS <- 23L
OUTER_K <- 5L
INNER_K <- 5L
BORUTA_ITERATIONS <- 300L
MODEL_ITERATIONS <- 500L
BASE_SEED <- 20260831L

cat("Scenario:", scenario, "\n")
cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("CatBoost threads:", THREADS, "\n")
cat("Nested CV:", OUTER_K, "outer x", INNER_K, "inner\n")

data_file <- file.path(out_dir, "hospital_patient_model_ready_129.csv")
registry_file <- file.path(out_dir, "lab_129_crosscohort_eligibility.csv")
if (!file.exists(data_file) || !file.exists(registry_file)) stop("Missing audit outputs")

dat <- fread(data_file, na.strings = c("", "NA"), showProgress = TRUE)
registry <- fread(registry_file)
if (nrow(dat) != 75248L || uniqueN(dat$patient_key) != 75248L) stop("Patient grain drift")

if (scenario == "A_PRIMARY") {
  base_features <- registry[primary_crosscohort_eligible == TRUE, final_variable_id]
} else {
  base_features <- registry[full_period_sensitivity_eligible == TRUE, final_variable_id]
}
base_features <- setdiff(base_features, "lab_mpv")
if (!length(base_features)) stop("No eligible cross-cohort features")

outcome_spec <- list(
  Outcome_NutriMetab = list(
    label = "Nutritional-Metabolic",
    excluded = c(
      "lab_alb", "lab_a_ratio_g", "lab_glob", "lab_tp",
      "lab_hgb", "lab_hct", "lab_rbc", "lab_rdw", "lab_mcv", "lab_mch", "lab_mchc"
    ),
    rationale = paste(
      "Exclude albumin/protein constituents and direct or closely overlapping",
      "erythrocyte indices for the hypoalbuminemia-or-anemia label"
    )
  ),
  Outcome_TumorBurden = list(
    label = "Tumor Burden",
    excluded = character(),
    rationale = "No laboratory feature directly constitutes the five metastasis flags"
  ),
  Outcome_TreatComp = list(
    label = "Treatment Complications",
    excluded = c(
      "lab_baso", "lab_baso_pct", "lab_eos", "lab_eos_pct", "lab_hct", "lab_hgb",
      "lab_lym", "lab_lym_pct", "lab_mch", "lab_mchc", "lab_mcv", "lab_mono",
      "lab_mono_pct", "lab_neu", "lab_neu_pct", "lab_plt", "lab_rbc", "lab_rdw",
      "lab_wbc", "lab_mpv"
    ),
    rationale = paste(
      "Exclude the complete CBC family because marrow suppression and bleeding",
      "components directly or closely overlap hematologic measurements"
    )
  )
)

param_grid <- CJ(
  depth = c(4L, 6L, 8L),
  learning_rate = c(0.03, 0.10),
  l2_leaf_reg = c(1, 5),
  sorted = FALSE
)

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold_id <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold_id
}

make_pool <- function(x, y = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (is.null(y)) {
    catboost.load_pool(data = x, feature_names = as.list(colnames(x)))
  } else {
    catboost.load_pool(data = x, label = as.numeric(y), feature_names = as.list(colnames(x)))
  }
}

class_weights <- function(y) c(1, sum(y == 0) / sum(y == 1))

fit_catboost <- function(x, y, params_extra, seed, iterations = MODEL_ITERATIONS) {
  pool <- make_pool(x, y)
  params <- c(list(
    loss_function = "Logloss",
    eval_metric = "AUC",
    iterations = as.integer(iterations),
    random_seed = as.integer(seed),
    thread_count = THREADS,
    class_weights = class_weights(y),
    logging_level = "Silent",
    allow_writing_files = FALSE
  ), params_extra)
  list(model = catboost.train(pool, params = params), pool = pool)
}

predict_catboost <- function(model, x) {
  pool <- make_pool(x)
  as.numeric(catboost.predict(model, pool, prediction_type = "Probability"))
}

auc_value <- function(y, p) {
  as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))
}

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

boruta_once <- function(x, y, seed, outcome, outer_fold) {
  set.seed(seed)
  x_real <- as.data.frame(x)
  x_shadow <- as.data.frame(lapply(x_real, function(z) sample(z, length(z), replace = FALSE)))
  names(x_shadow) <- paste0("shadow__", names(x_real))
  x_all <- cbind(x_real, x_shadow)
  fit <- fit_catboost(
    x_all,
    y,
    params_extra = list(depth = 6L, learning_rate = 0.05, l2_leaf_reg = 5),
    seed = seed,
    iterations = BORUTA_ITERATIONS
  )
  imp <- as.numeric(catboost.get_feature_importance(
    fit$model, fit$pool, type = "PredictionValuesChange"
  ))
  imp_dt <- data.table(feature = names(x_all), importance = imp)
  threshold <- max(imp_dt[grepl("^shadow__", feature), importance], na.rm = TRUE)
  real_dt <- imp_dt[!grepl("^shadow__", feature)]
  real_dt[, `:=`(
    selected = importance > threshold,
    shadow_threshold = threshold,
    outcome = outcome,
    outer_fold = outer_fold,
    scenario = scenario
  )]
  selected <- real_dt[selected == TRUE, feature]
  if (length(selected) < 2L) {
    selected <- real_dt[order(-importance)][seq_len(min(2L, .N)), feature]
    real_dt[feature %in% selected, selected := TRUE]
  }
  list(selected = selected, importance = real_dt)
}

TUNE_CLUSTER <- makeCluster(23L)
clusterEvalQ(TUNE_CLUSTER, {suppressPackageStartupMessages({library(catboost);library(data.table);library(pROC)});setDTthreads(1);NULL})
clusterExport(TUNE_CLUSTER,c('make_stratified_folds','make_pool','class_weights','fit_catboost','predict_catboost','auc_value','MODEL_ITERATIONS'))
clusterEvalQ(TUNE_CLUSTER,{THREADS<-1L;NULL})
tune_inner <- function(d, features, outcome, outer_fold, seed) {
 cache_dir<-file.path(out_dir,'tuning_cache');dir.create(cache_dir,showWarnings=FALSE)
 cache_file<-file.path(cache_dir,paste0(outcome,'_outer_',outer_fold,'.rds'))
 if(file.exists(cache_file)) {cat('    Cached tuning:',basename(cache_file),'\n');return(readRDS(cache_file))}
 inner_id<-make_stratified_folds(d[[outcome]],INNER_K,seed)
 tasks<-CJ(g=seq_len(nrow(param_grid)),inner_fold=seq_len(INNER_K))
 clusterExport(TUNE_CLUSTER,c('d','features','outcome','outer_fold','seed','inner_id','param_grid','scenario'),envir=environment())
 rows<-parLapplyLB(TUNE_CLUSTER,seq_len(nrow(tasks)),function(j,tasks){
   g<-tasks$g[j];inner_fold<-tasks$inner_fold[j]
   tr<-which(inner_id!=inner_fold);va<-which(inner_id==inner_fold)
   pars<-list(depth=as.integer(param_grid$depth[g]),learning_rate=as.numeric(param_grid$learning_rate[g]),l2_leaf_reg=as.numeric(param_grid$l2_leaf_reg[g]))
   fit<-fit_catboost(d[tr,..features],d[[outcome]][tr],pars,seed=seed+g*100L+inner_fold)
   p<-predict_catboost(fit$model,d[va,..features])
   data.table(scenario=scenario,outcome=outcome,outer_fold=outer_fold,inner_fold=inner_fold,depth=pars$depth,learning_rate=pars$learning_rate,l2_leaf_reg=pars$l2_leaf_reg,auc=auc_value(d[[outcome]][va],p),n_features=length(features))
 },tasks=tasks)
 detail<-rbindlist(rows);summary<-detail[,.(mean_auc=mean(auc),sd_auc=sd(auc)),by=.(depth,learning_rate,l2_leaf_reg,n_features)]
 setorder(summary,-mean_auc,sd_auc,depth,learning_rate,l2_leaf_reg)
 ans<-list(detail=detail,summary=summary,best=summary[1L]);saveRDS(ans,cache_file);ans
}

metric_row <- function(y, p, threshold, outcome, split_name, n_features) {
  roc_obj <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))
  pred <- as.integer(p >= threshold)
  tp <- sum(pred == 1 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)
  lp <- qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  cal <- tryCatch(glm(y ~ lp, family = binomial()), error = function(e) NULL)
  data.table(
    scenario = scenario,
    outcome = outcome,
    split = split_name,
    n = length(y),
    positives = sum(y == 1),
    prevalence = mean(y),
    n_features = n_features,
    auc = as.numeric(pROC::auc(roc_obj)),
    auc_ci_lower = ci[1],
    auc_ci_upper = ci[3],
    average_precision = average_precision(y, p),
    brier = mean((p - y)^2),
    calibration_intercept = if (is.null(cal)) NA_real_ else unname(coef(cal)[1]),
    calibration_slope = if (is.null(cal)) NA_real_ else unname(coef(cal)[2]),
    threshold_from_validation = threshold,
    sensitivity = if ((tp + fn) == 0) NA_real_ else tp / (tp + fn),
    specificity = if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  )
}

all_outer_perf <- list()
all_boruta <- list()
all_tuning <- list()
all_feature_votes <- list()
all_final_features <- list()
all_final_tuning <- list()
all_metrics <- list()
all_predictions <- list()
all_calibration <- list()
list_index <- 0L

for (outcome_idx in 1:2) {
  outcome <- names(outcome_spec)[outcome_idx]
  spec <- outcome_spec[[outcome_idx]]
  candidates <- setdiff(base_features, spec$excluded)
  candidates <- candidates[candidates %in% names(dat)]
  if (length(candidates) < 2L) stop("Too few candidates for ", outcome)

  train <- dat[split_calendar_entry == "train" & !is.na(get(outcome))]
  validation <- dat[split_calendar_entry == "validation" & !is.na(get(outcome))]
  test <- dat[split_calendar_entry == "test" & !is.na(get(outcome))]
  setorder(train, patient_key)
  setorder(validation, patient_key)
  setorder(test, patient_key)

  cat("\nOutcome:", outcome, "-", spec$label, "\n")
  cat("Candidates after overlap exclusion:", length(candidates), "\n")
  cat("Train/validation/test:", nrow(train), nrow(validation), nrow(test), "\n")
  cat("Train positives:", sum(train[[outcome]] == 1), "\n")

  outer_id <- make_stratified_folds(train[[outcome]], OUTER_K, BASE_SEED + outcome_idx)
  selected_by_outer <- vector("list", OUTER_K)

  for (outer_fold in seq_len(OUTER_K)) {
    cat("  Outer fold", outer_fold, "of", OUTER_K, "\n")
    outer_train <- train[outer_id != outer_fold]
    outer_holdout <- train[outer_id == outer_fold]
    fold_seed <- BASE_SEED + outcome_idx * 10000L + outer_fold * 1000L

    boruta <- boruta_once(
      outer_train[, ..candidates], outer_train[[outcome]], fold_seed,
      outcome, outer_fold
    )
    selected_by_outer[[outer_fold]] <- boruta$selected
    all_boruta[[length(all_boruta) + 1L]] <- boruta$importance
    cat("    Boruta selected:", length(boruta$selected), "\n")

    tuning <- tune_inner(
      outer_train, boruta$selected, outcome, outer_fold,
      seed = fold_seed + 500L
    )
    all_tuning[[length(all_tuning) + 1L]] <- tuning$detail
    best <- tuning$best
    cat(
      "    Best inner parameters: depth", best$depth,
      "lr", best$learning_rate, "l2", best$l2_leaf_reg,
      "mean AUC", round(best$mean_auc, 4), "\n"
    )

    selected_outer <- boruta$selected
    outer_fit <- fit_catboost(
      outer_train[, ..selected_outer], outer_train[[outcome]],
      params_extra = list(
        depth = as.integer(best$depth),
        learning_rate = as.numeric(best$learning_rate),
        l2_leaf_reg = as.numeric(best$l2_leaf_reg)
      ),
      seed = fold_seed + 900L
    )
    p_outer <- predict_catboost(outer_fit$model, outer_holdout[, ..selected_outer])
    outer_auc <- auc_value(outer_holdout[[outcome]], p_outer)
    all_outer_perf[[length(all_outer_perf) + 1L]] <- data.table(
      scenario = scenario,
      outcome = outcome,
      outer_fold = outer_fold,
      n_holdout = nrow(outer_holdout),
      positives_holdout = sum(outer_holdout[[outcome]] == 1),
      n_selected = length(boruta$selected),
      depth = as.integer(best$depth),
      learning_rate = as.numeric(best$learning_rate),
      l2_leaf_reg = as.numeric(best$l2_leaf_reg),
      inner_mean_auc = as.numeric(best$mean_auc),
      outer_auc = outer_auc
    )
    cat("    Outer holdout AUC:", round(outer_auc, 4), "\n")

    checkpoint <- list(
      scenario = scenario,
      outcome = outcome,
      completed_outer_fold = outer_fold,
      selected_by_outer = selected_by_outer[seq_len(outer_fold)],
      updated_at = as.character(Sys.time())
    )
    saveRDS(
      checkpoint,
      file.path(out_dir, paste0("nested_checkpoint_", scenario_slug, "_", outcome, ".rds"))
    )
  }

  votes <- data.table(feature = candidates)
  votes[, votes := vapply(feature, function(f) {
    sum(vapply(selected_by_outer, function(z) f %in% z, logical(1)))
  }, integer(1))]
  votes[, status := fifelse(votes >= 4L, "CONFIRMED",
                            fifelse(votes == 3L, "TENTATIVE", "REJECTED"))]
  votes[, `:=`(scenario = scenario, outcome = outcome)]
  final_features <- votes[votes >= 3L, feature]
  if (length(final_features) < 2L) {
    final_features <- votes[order(-votes)][seq_len(min(2L, .N)), feature]
  }
  all_feature_votes[[length(all_feature_votes) + 1L]] <- votes
  all_final_features[[length(all_final_features) + 1L]] <- data.table(
    scenario = scenario,
    outcome = outcome,
    feature = final_features,
    votes = votes[match(final_features, feature), votes],
    status = votes[match(final_features, feature), status],
    preexcluded_features = paste(intersect(base_features, spec$excluded), collapse = ";"),
    exclusion_rationale = spec$rationale
  )
  cat("  Stable final features (>=3/5 outer votes):", length(final_features), "\n")

  final_tune <- tune_inner(
    train, final_features, outcome, outer_fold = 0L,
    seed = BASE_SEED + outcome_idx * 100000L + 777L
  )
  final_tune$detail[, tuning_scope := "full_calendar_training"]
  all_final_tuning[[length(all_final_tuning) + 1L]] <- final_tune$detail
  final_best <- final_tune$best

  final_fit <- fit_catboost(
    train[, ..final_features], train[[outcome]],
    params_extra = list(
      depth = as.integer(final_best$depth),
      learning_rate = as.numeric(final_best$learning_rate),
      l2_leaf_reg = as.numeric(final_best$l2_leaf_reg)
    ),
    seed = BASE_SEED + outcome_idx * 100000L + 999L
  )
  model_file <- file.path(
    out_dir, paste0("catboost_nested_", scenario_slug, "_", outcome, ".cbm")
  )
  catboost.save_model(final_fit$model, model_file)

  p_val_raw <- predict_catboost(final_fit$model, validation[, ..final_features])
  p_test_raw <- predict_catboost(final_fit$model, test[, ..final_features])
  eps <- 1e-6
  val_lp <- qlogis(pmin(pmax(p_val_raw, eps), 1 - eps))
  calibration <- glm(validation[[outcome]] ~ val_lp, family = binomial())
  p_val <- as.numeric(predict(calibration, type = "response"))
  test_lp <- qlogis(pmin(pmax(p_test_raw, eps), 1 - eps))
  p_test <- plogis(unname(coef(calibration)[1]) + unname(coef(calibration)[2]) * test_lp)
  roc_val <- pROC::roc(validation[[outcome]], p_val, quiet = TRUE, direction = "<")
  threshold <- as.numeric(pROC::coords(
    roc_val, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE
  )[1])

  all_metrics[[length(all_metrics) + 1L]] <- metric_row(
    validation[[outcome]], p_val, threshold, outcome, "calendar_validation", length(final_features)
  )
  all_metrics[[length(all_metrics) + 1L]] <- metric_row(
    test[[outcome]], p_test, threshold, outcome, "calendar_test", length(final_features)
  )
  all_predictions[[length(all_predictions) + 1L]] <- data.table(
    scenario = scenario,
    outcome = outcome,
    patient_key = test$patient_key,
    first_observed_test_date = test$first_observed_test_date,
    observed = test[[outcome]],
    raw_probability = p_test_raw,
    calibrated_probability = p_test
  )
  all_calibration[[length(all_calibration) + 1L]] <- data.table(
    scenario = scenario,
    outcome = outcome,
    intercept = unname(coef(calibration)[1]),
    slope = unname(coef(calibration)[2]),
    threshold_from_calendar_validation = threshold,
    depth = as.integer(final_best$depth),
    learning_rate = as.numeric(final_best$learning_rate),
    l2_leaf_reg = as.numeric(final_best$l2_leaf_reg),
    iterations = MODEL_ITERATIONS,
    thread_count = THREADS,
    model_file = basename(model_file)
  )
  cat("  Calendar test AUC:", round(auc_value(test[[outcome]], p_test), 4), "\n")
}

stopCluster(TUNE_CLUSTER)
all_outer_perf[[length(all_outer_perf)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_outer_performance.csv"))[outcome=="Outcome_TreatComp"]
all_boruta[[length(all_boruta)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_boruta_importance.csv"))[outcome=="Outcome_TreatComp"]
all_tuning[[length(all_tuning)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_inner_tuning.csv"))[outcome=="Outcome_TreatComp"]
all_feature_votes[[length(all_feature_votes)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_feature_votes.csv"))[outcome=="Outcome_TreatComp"]
all_final_features[[length(all_final_features)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_final_features.csv"))[outcome=="Outcome_TreatComp"]
all_final_tuning[[length(all_final_tuning)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_final_tuning.csv"))[outcome=="Outcome_TreatComp"]
all_metrics[[length(all_metrics)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_calendar_metrics.csv"))[outcome=="Outcome_TreatComp"]
all_predictions[[length(all_predictions)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_calendar_test_predictions.csv"))[outcome=="Outcome_TreatComp"]
all_calibration[[length(all_calibration)+1L]] <- fread(file.path(Sys.getenv("CIPDS_REFERENCE_DIR"),"outputs","nested_a_primary_calibration_and_model_registry.csv"))[outcome=="Outcome_TreatComp"]
outer_perf <- rbindlist(all_outer_perf, fill = TRUE)
boruta_detail <- rbindlist(all_boruta, fill = TRUE)
tuning_detail <- rbindlist(all_tuning, fill = TRUE)
feature_votes <- rbindlist(all_feature_votes, fill = TRUE)
final_features_dt <- rbindlist(all_final_features, fill = TRUE)
final_tuning_dt <- rbindlist(all_final_tuning, fill = TRUE)
metrics_dt <- rbindlist(all_metrics, fill = TRUE)
predictions_dt <- rbindlist(all_predictions, fill = TRUE)
calibration_dt <- rbindlist(all_calibration, fill = TRUE)

prefix <- paste0("nested_", scenario_slug, "_")
fwrite(outer_perf, file.path(out_dir, paste0(prefix, "outer_performance.csv")), bom = TRUE)
fwrite(boruta_detail, file.path(out_dir, paste0(prefix, "boruta_importance.csv")), bom = TRUE)
fwrite(tuning_detail, file.path(out_dir, paste0(prefix, "inner_tuning.csv")), bom = TRUE)
fwrite(feature_votes, file.path(out_dir, paste0(prefix, "feature_votes.csv")), bom = TRUE)
fwrite(final_features_dt, file.path(out_dir, paste0(prefix, "final_features.csv")), bom = TRUE)
fwrite(final_tuning_dt, file.path(out_dir, paste0(prefix, "final_tuning.csv")), bom = TRUE)
fwrite(metrics_dt, file.path(out_dir, paste0(prefix, "calendar_metrics.csv")), bom = TRUE)
fwrite(predictions_dt, file.path(out_dir, paste0(prefix, "calendar_test_predictions.csv")), bom = TRUE)
fwrite(calibration_dt, file.path(out_dir, paste0(prefix, "calibration_and_model_registry.csv")), bom = TRUE)

manifest <- list(
  scenario = scenario,
  completed_at = as.character(Sys.time()),
  patient_grain = "one row per pseudonymous patient",
  outer_split = "fixed calendar entry train/validation/test; nested CV uses training period only",
  nested_cv = list(outer_folds = OUTER_K, inner_folds = INNER_K),
  boruta = list(
    style = "CatBoost real features versus within-fold permuted shadow features",
    iterations = BORUTA_ITERATIONS,
    final_rule = "confirmed 4-5/5 or tentative 3/5 outer selections"
  ),
  catboost = list(
    version = as.character(packageVersion("catboost")),
    thread_count = THREADS,
    iterations = MODEL_ITERATIONS,
    grid = as.data.frame(param_grid)
  ),
  candidate_count_before_outcome_exclusion = length(base_features),
  selected_feature_counts = as.list(table(final_features_dt$outcome)),
  test_auc = setNames(
    as.list(metrics_dt[split == "calendar_test", auc]),
    metrics_dt[split == "calendar_test", outcome]
  ),
  analysis_status = "NO_MPV_SENSITIVITY"
)
write_json(
  manifest,
  file.path(out_dir, paste0(prefix, "run_manifest.json")),
  pretty = TRUE, auto_unbox = TRUE
)
cat("\nCompleted:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
print(metrics_dt)
