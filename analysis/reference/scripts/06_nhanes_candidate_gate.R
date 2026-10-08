options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

package_dir <- Sys.getenv(
  "CIPDS_PACKAGE_DIR",
  unset = "reference"
)
out_dir <- file.path(package_dir, "outputs")
input_rdata <- file.path(package_dir, "inputs", "nhanes_phase7_batch1.RData")
hospital_audit_file <- file.path(out_dir, "lab_129_hospital_audit.csv")

load(input_rdata)
if (!exists("nhanes_valid") || !is.data.frame(nhanes_valid)) {
  stop("nhanes_valid was not found in the locked NHANES input")
}
nh <- as.data.table(nhanes_valid)
expected_cycles <- c(
  "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
  "2009-2010", "2011-2012", "2013-2014", "2015-2016"
)
if (!setequal(unique(nh$CYCLE), expected_cycles)) stop("NHANES cycle set drift")
if (nrow(nh) != 44772L || uniqueN(nh$SEQN) != 44772L) stop("NHANES row/key drift")

coalesce_numeric <- function(a, b) {
  ans <- a
  ans[is.na(ans)] <- b[is.na(ans)]
  ans
}

# Exact cross-database mappings use the locked NHANES variables already loaded
# in Phase 2. Deterministic derivations are explicit and unit harmonized.
mapped <- list(
  lab_a_ratio_g = nh$Albumin_Globulin_Ratio,
  lab_alb = coalesce_numeric(nh$LBXSAL * 10, nh$LBDSALSI),
  lab_alp = coalesce_numeric(nh$LBXSAPSI, nh$LBDSAPSI),
  lab_alt = nh$Alanine_Aminotransferase,
  lab_ast = nh$Aspartate_Aminotransferase,
  lab_baso = nh$LBDBANO,
  lab_baso_pct = nh$LBXBAPCT,
  lab_ca = coalesce_numeric(nh$LBXSCA / 4, nh$LBDSCASI),
  lab_cl = nh$Chloride,
  lab_crea = nh$Creatinine,
  lab_eos = nh$LBDEONO,
  lab_eos_pct = nh$Eosinophil_Percentage,
  lab_ggt = nh$Gamma_Glutamyl_Transferase,
  lab_glob = coalesce_numeric(nh$LBXSGB * 10, nh$LBDSGBSI),
  lab_glu = coalesce_numeric(nh$LBXSGL / 18, nh$LBDSGLSI),
  lab_hba1c = nh$HbA1c,
  lab_hcrp = coalesce_numeric(nh$LBXHSCRP, nh$LBXCRP * 10),
  lab_hct = nh$Hematocrit,
  lab_hdl_c = nh$LBDHDDSI,
  lab_hgb = nh$Hemoglobin,
  lab_k = nh$Potassium,
  lab_ldh = coalesce_numeric(nh$LBXSLDSI, nh$LBDSLDSI),
  lab_lym = nh$LBDLYMNO,
  lab_lym_pct = nh$Lymphocyte_Percentage,
  lab_mch = nh$Hemoglobin / nh$Red_Blood_Cell_Count,
  lab_mchc = 100 * nh$Hemoglobin / nh$Hematocrit,
  lab_mcv = 10 * nh$Hematocrit / nh$Red_Blood_Cell_Count,
  lab_mono = nh$LBDMONO,
  lab_mono_pct = nh$Monocyte_Percentage,
  lab_mpv = nh$Mean_Platelet_Volume,
  lab_na = nh$Sodium,
  lab_neu = nh$LBDNENO,
  lab_neu_pct = nh$Neutrophil_Percentage,
  lab_phos = coalesce_numeric(nh$LBXSPH / 3.1, nh$LBDSPHSI),
  lab_plt = nh$Platelet_Count,
  lab_rbc = nh$Red_Blood_Cell_Count,
  lab_rdw = nh$Red_Cell_Distribution_Width,
  lab_t_ch = nh$LBDTCSI,
  lab_tbil = coalesce_numeric(nh$LBXSTB * 17.1, nh$LBDSTBSI),
  lab_tco2 = nh$LBXSC3SI,
  lab_tp = nh$LBXSTP * 10,
  lab_ua = coalesce_numeric(nh$LBXSUA * 59.48, nh$LBDSUASI),
  lab_urea = coalesce_numeric(nh$LBXSBU / 2.8, nh$LBDSBUSI),
  lab_wbc = nh$White_Blood_Cell_Count
)

source_note <- c(
  lab_a_ratio_g = "Albumin/Globulin",
  lab_alb = "coalesce(LBXSAL*10,LBDSALSI)",
  lab_alp = "coalesce(LBXSAPSI,LBDSAPSI)", lab_alt = "LBXSATSI",
  lab_ast = "LBXSASSI", lab_baso = "LBDBANO", lab_baso_pct = "LBXBAPCT",
  lab_ca = "coalesce(LBXSCA/4,LBDSCASI)", lab_cl = "LBXSCLSI",
  lab_crea = "coalesce(LBXSCR*88.4,LBDSCRSI)",
  lab_eos = "LBDEONO", lab_eos_pct = "LBXEOPCT", lab_ggt = "LBXSGTSI",
  lab_glob = "coalesce(LBXSGB*10,LBDSGBSI)",
  lab_glu = "coalesce(LBXSGL/18,LBDSGLSI)", lab_hba1c = "LBXGH",
  lab_hcrp = "LBXHSCRP or LBXCRP*10", lab_hct = "LBXHCT",
  lab_hdl_c = "LBDHDDSI", lab_hgb = "LBXHGB*10", lab_k = "LBXSKSI",
  lab_ldh = "LBXSLDSI or LBDSLDSI", lab_lym = "LBDLYMNO",
  lab_lym_pct = "LBXLYPCT", lab_mch = "Hemoglobin/RBC",
  lab_mchc = "100*Hemoglobin/Hematocrit", lab_mcv = "10*Hematocrit/RBC",
  lab_mono = "LBDMONO", lab_mono_pct = "LBXMOPCT", lab_mpv = "LBXMPSI",
  lab_na = "LBXSNASI", lab_neu = "LBDNENO", lab_neu_pct = "LBXNEPCT",
  lab_phos = "coalesce(LBXSPH/3.1,LBDSPHSI)", lab_plt = "LBXPLTSI",
  lab_rbc = "LBXRBCSI", lab_rdw = "LBXRDW", lab_t_ch = "LBDTCSI",
  lab_tbil = "coalesce(LBXSTB*17.1,LBDSTBSI)",
  lab_tco2 = "LBXSC3SI", lab_tp = "LBXSTP*10",
  lab_ua = "coalesce(LBXSUA*59.48,LBDSUASI)",
  lab_urea = "coalesce(LBXSBU/2.8,LBDSBUSI)", lab_wbc = "LBXWBCSI"
)

mapped_dt <- data.table(SEQN = nh$SEQN, CYCLE = nh$CYCLE)
for (id in names(mapped)) mapped_dt[[id]] <- as.numeric(mapped[[id]])
saveRDS(mapped_dt, file.path(out_dir, "nhanes_expanded_candidate_matrix.rds"), compress = "xz")

cycle_rows <- list()
overall_rows <- list()
for (id in names(mapped)) {
  overall_rows[[id]] <- data.table(
    final_variable_id = id,
    nhanes_n = nrow(mapped_dt),
    nhanes_nonmissing_n = sum(!is.na(mapped_dt[[id]])),
    nhanes_missing_pct = 100 * mean(is.na(mapped_dt[[id]])),
    nhanes_source_or_transform = unname(source_note[id])
  )
  cycle_rows[[id]] <- mapped_dt[, .(
    n = .N,
    nonmissing_n = sum(!is.na(get(id))),
    missing_pct = 100 * mean(is.na(get(id)))
  ), by = CYCLE][, final_variable_id := id]
}
overall <- rbindlist(overall_rows, use.names = TRUE)
cycle_long <- rbindlist(cycle_rows, use.names = TRUE)
setcolorder(cycle_long, c("final_variable_id", "CYCLE", "n", "nonmissing_n", "missing_pct"))
fwrite(cycle_long, file.path(out_dir, "nhanes_129_missingness_by_cycle.csv"), bom = TRUE)

cycle_summary <- cycle_long[, .(
  nhanes_cycles_with_any_data = sum(nonmissing_n > 0),
  nhanes_max_cycle_missing_pct = max(missing_pct),
  nhanes_min_cycle_missing_pct = min(missing_pct),
  nhanes_cycle_missingness_range_pp = max(missing_pct) - min(missing_pct)
), by = final_variable_id]

hospital <- fread(hospital_audit_file)
registry <- merge(hospital, overall, by = "final_variable_id", all.x = TRUE)
registry <- merge(registry, cycle_summary, by = "final_variable_id", all.x = TRUE)
registry[, nhanes_mapped := !is.na(nhanes_n)]
registry[, nhanes_all_9_cycles_available := nhanes_cycles_with_any_data == 9L]

registry[, crosscohort_tier := fifelse(
  value_type_is_numeric != TRUE, "C_NONNUMERIC",
  fifelse(
    hospital_missingness_tier == "C", "C_HOSPITAL_MISSING_GT50",
    fifelse(
      nhanes_mapped != TRUE, "C_NO_LOCKED_NHANES_MAPPING",
      fifelse(
        calendar_missingness_drift_gt20pp == TRUE, "C_HOSPITAL_CALENDAR_DRIFT_GT20PP",
        fifelse(
          hospital_missingness_tier == "A" & nhanes_all_9_cycles_available == TRUE &
            nhanes_max_cycle_missing_pct <= 30,
          "A_PRIMARY",
          fifelse(
            hospital_missingness_tier %in% c("A", "B") &
              nhanes_all_9_cycles_available == TRUE & nhanes_max_cycle_missing_pct <= 50,
            "B_FULL_PERIOD_SENSITIVITY",
            "B_RESTRICTED_CYCLE_OR_HIGH_NHANES_MISSING"
          )
        )
      )
    )
  )
)]
registry[, primary_crosscohort_eligible := crosscohort_tier == "A_PRIMARY"]
registry[, full_period_sensitivity_eligible := crosscohort_tier %in% c(
  "A_PRIMARY", "B_FULL_PERIOD_SENSITIVITY"
)]

setorder(registry, crosscohort_tier, train_missing_pct, final_variable_id)
fwrite(registry, file.path(out_dir, "lab_129_crosscohort_eligibility.csv"), bom = TRUE)

crosswalk <- data.table(
  final_variable_id = names(mapped),
  nhanes_source_or_transform = unname(source_note[names(mapped)])
)
fwrite(crosswalk, file.path(out_dir, "nhanes_lab_crosswalk_locked.csv"), bom = TRUE)

summary <- list(
  nhanes_rows = nrow(mapped_dt),
  nhanes_unique_seqn = uniqueN(mapped_dt$SEQN),
  nhanes_cycles = length(unique(mapped_dt$CYCLE)),
  locked_nhanes_mappings = length(mapped),
  crosscohort_tier_counts = as.list(table(registry$crosscohort_tier)),
  primary_crosscohort_candidate_count = sum(registry$primary_crosscohort_eligible),
  full_period_sensitivity_candidate_count = sum(registry$full_period_sensitivity_eligible),
  gate = paste(
    "A: hospital training missing <=30%, no >20pp calendar drift,",
    "NHANES <=30% missing in every cycle; B: hospital <=50% and NHANES <=50%",
    "in every cycle; other mapped features restricted-cycle only"
  )
)
write_json(summary, file.path(out_dir, "lab_129_crosscohort_eligibility_qa.json"),
           pretty = TRUE, auto_unbox = TRUE)
print(summary)
