suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

package_dir <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "."
)
outputs_dir <- file.path(package_dir, "outputs")
out_dir <- file.path(outputs_dir, "geometry_v1")
log_dir <- file.path(package_dir, "logs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "21_prepare_geometry_inputs.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

cat("Start:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

eligibility <- fread(file.path(outputs_dir, "lab_129_crosscohort_eligibility.csv"))
eligible <- eligibility[
  primary_crosscohort_eligible == TRUE & primary_geometry_eligible == TRUE
]
geometry_labs <- eligible$final_variable_id
if (length(geometry_labs) != 22L) {
  stop("Expected 22 primary cross-cohort geometry laboratories, observed ", length(geometry_labs))
}

expanded <- as.data.table(readRDS(file.path(outputs_dir, "nhanes_expanded_candidate_matrix.rds")))
if (nrow(expanded) != 44772L || uniqueN(expanded$SEQN) != 44772L) {
  stop("Expanded NHANES laboratory matrix grain drift")
}
missing_labs <- setdiff(geometry_labs, names(expanded))
if (length(missing_labs)) {
  stop("Eligible laboratory variables missing from NHANES matrix: ", paste(missing_labs, collapse = ", "))
}

paired <- fread(file.path(outputs_dir, "nhanes_paired_discrimination_primary_cohort.csv"))
if (nrow(paired) != as.integer(fread(file.path(package_dir,"qa/cohort_registry.csv"))[cohort=="paired",n]) || uniqueN(paired$SEQN) != as.integer(fread(file.path(package_dir,"qa/cohort_registry.csv"))[cohort=="paired",n])) {
  stop("Primary paired-discrimination cohort grain drift")
}

load(file.path(package_dir, "inputs", "nhanes_phase7_batch1.RData"))
if (!exists("nhanes_valid")) stop("nhanes_valid missing")
nhanes <- as.data.table(nhanes_valid)
if (nrow(nhanes) != 44772L || uniqueN(nhanes$SEQN) != 44772L) {
  stop("NHANES source grain drift")
}

lab_index <- match(paired$SEQN, expanded$SEQN)
source_index <- match(paired$SEQN, nhanes$SEQN)
if (anyNA(lab_index) || anyNA(source_index)) stop("NHANES geometry merge failure")

analysis <- copy(paired)
for (v in geometry_labs) analysis[[v]] <- expanded[[v]][lab_index]
analysis[, Sex := nhanes$Sex[source_index]]
analysis[, Race := nhanes$Race[source_index]]

if (anyDuplicated(analysis$SEQN)) stop("Duplicate SEQN after geometry merge")
if (any(!is.finite(analysis$pooled_mec_weight) | analysis$pooled_mec_weight <= 0)) {
  stop("Invalid MEC weights in geometry cohort")
}

fwrite(
  analysis,
  file.path(out_dir, "nhanes_geometry_analysis.csv"),
  na = ""
)

dictionary <- eligible[, .(
  variable = final_variable_id,
  canonical_name_english,
  canonical_name_chinese,
  clinical_domain,
  analysis_unit,
  hospital_missingness_tier,
  overall_missing_pct,
  nhanes_missing_pct,
  nhanes_source_or_transform
)]
fwrite(dictionary, file.path(out_dir, "geometry_lab_dictionary.csv"))

qa <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  paired_cohort_n = nrow(analysis),
  paired_cohort_unique_seqn = uniqueN(analysis$SEQN),
  death_n = sum(analysis$Death_AllCause == 1L),
  cancer_history_n = sum(analysis$Cancer_Diagnosed == 1L, na.rm = TRUE),
  geometry_lab_n = length(geometry_labs),
  all_weights_positive = all(is.finite(analysis$pooled_mec_weight) & analysis$pooled_mec_weight > 0),
  source_grain_n = nrow(nhanes),
  expanded_grain_n = nrow(expanded),
  checks_passed = TRUE
)
write_json(
  qa,
  file.path(out_dir, "geometry_input_qa.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

cat("Prepared NHANES geometry cohort:", nrow(analysis), "participants\n")
cat("Geometry laboratories:", length(geometry_labs), "\n")
cat("Deaths:", qa$death_n, " Cancer history:", qa$cancer_history_n, "\n")
cat("End:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
