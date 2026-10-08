from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


RELEASE = Path("hospital_source")
CURRENT_BASELINE = Path("external_inputs/Combined_Baseline_Data.csv")
PACKAGE = Path("reference")
OUT = PACKAGE / "outputs"


# These are exact, locked laboratory identifiers from the formal release.
# They intentionally do not use fuzzy matching.
LAB_MAP = {
    "Albumin_Globulin_Ratio": "lab_a_ratio_g",
    "Alpha_Fetoprotein": "lab_afp",
    "Albumin": "lab_alb",
    "Alkaline_Phosphatase": "lab_alp",
    "Alanine_Aminotransferase": "lab_alt",
    "Alpha_Amylase": "lab_amy",
    "Anti_Thyroglobulin_Antibody": "lab_anti_tg",
    "Anti_Thyroid_Microsomal_Antibody": "lab_anti_tpo",
    "Apolipoprotein_A1": "lab_apoa1",
    "Apolipoprotein_B": "lab_apob",
    "Activated_Partial_Thromboplastin_Time": "lab_aptt",
    "Aspartate_Aminotransferase": "lab_ast",
    "Vitamin_B12": "lab_b12",
    "Basophil_Count": "lab_baso",
    "Basophil_Percentage": "lab_baso_pct",
    "Total_Calcium": "lab_ca",
    "Cancer_Antigen_125": "lab_ca1252",
    "Cancer_Antigen_153": "lab_ca1532",
    "Cancer_Antigen_199": "lab_ca199",
    "Cancer_Antigen_724": "lab_ca724",
    "Carcinoembryonic_Antigen": "lab_cea",
    "Cholinesterase": "lab_che",
    "Creatine_Kinase": "lab_ck",
    "Chloride": "lab_cl",
    "Carbon_Dioxide_Combining_Power": "lab_co2_cp",
    "Cortisol": "lab_cortisol",
    "Creatinine": "lab_crea",
    "Cardiac_Troponin_I": "lab_ctni",
    "Cytokeratin_19_Fragment": "lab_cyfra21_1",
    "Cystatin_C": "lab_cys_c",
    "Direct_Bilirubin": "lab_dbil",
    "D_Dimer": "lab_ddi",
    "Eosinophil_Count": "lab_eos",
    "Eosinophil_Percentage": "lab_eos_pct",
    "Free_Prostate_Specific_Antigen": "lab_f_psa",
    "Free_to_Total_PSA_Ratio": "lab_f_psa_ratio_t_psa",
    # The old English table calls source field FBG "Fibrinogen".  The formal
    # dictionary currently labels it a fasting-glucose field.  It is retained
    # only for deterministic row linkage and is not used as a model predictor.
    "Fibrinogen": "lab_fbg",
    "Serum_Iron": "lab_fe",
    "Ferritin": "lab_ferr",
    "Folate": "lab_folate",
    "Free_Triiodothyronine": "lab_ft3",
    "Free_Thyroxine": "lab_ft4",
    "Gamma_Glutamyl_Transferase": "lab_ggt",
    "Globulin": "lab_glob",
    "Glucose": "lab_glu",
    "Glycated_Hemoglobin": "lab_hba1c",
    "Alpha_Hydroxybutyrate_Dehydrogenase": "lab_hbdh",
    "High_Sensitivity_CRP": "lab_hcrp",
    "Hematocrit": "lab_hct",
    "HDL_Cholesterol": "lab_hdl_c",
    "Human_Epididymis_Protein_4": "lab_he4",
    "Hemoglobin": "lab_hgb",
    "Indirect_Bilirubin": "lab_idbil",
    "International_Normalized_Ratio": "lab_inr",
    "Potassium": "lab_k",
    "Lactate_Dehydrogenase": "lab_ldh",
    "LDL_Cholesterol": "lab_ldl_c",
    "Lipoprotein_A": "lab_lpa",
    "Lymphocyte_Count": "lab_lym",
    "Lymphocyte_Percentage": "lab_lym_pct",
    "Mean_Corpuscular_Hemoglobin": "lab_mch",
    "Mean_Corpuscular_Hemoglobin_Concentration": "lab_mchc",
    "Mean_Corpuscular_Volume": "lab_mcv",
    "Magnesium": "lab_mg",
    "Monocyte_Count": "lab_mono",
    "Monocyte_Percentage": "lab_mono_pct",
    "Mean_Platelet_Volume": "lab_mpv",
    "Sodium": "lab_na",
    "Neutrophil_Count": "lab_neu",
    "Neutrophil_Percentage": "lab_neu_pct",
    "Neuron_Specific_Enolase": "lab_nse",
    "Prealbumin": "lab_palb",
    "Plateletcrit": "lab_plateletcrit",
    "Platelet_Distribution_Width": "lab_pdw",
    "Phosphorus": "lab_phos",
    "Protein_Induced_by_Vitamin_K_Absence": "lab_pivka",
    "Platelet_Count": "lab_plt",
    "Pro_Gastrin_Releasing_Peptide": "lab_progrp",
    "Prothrombin_Time": "lab_pt",
    "Prothrombin_Activity": "lab_pt_pct",
    "Prothrombin_Ratio": "lab_ptr",
    "Red_Blood_Cell_Count": "lab_rbc",
    "Red_Cell_Distribution_Width": "lab_rdw",
    "Reticulocyte_Count": "lab_ret_count",
    "Reticulocyte_Percentage": "lab_ret_pct",
    "Risk_of_Ovarian_Malignancy_Algorithm": "lab_roma",
    "Total_Cholesterol": "lab_t_ch",
    "Total_Prostate_Specific_Antigen": "lab_t_psa",
    "Total_Bile_Acid": "lab_tba",
    "Total_Bilirubin": "lab_tbil",
    "Total_Carbon_Dioxide": "lab_tco2",
    "Triglycerides": "lab_tg",
    "Total_Protein": "lab_tp",
    "Thyroid_Stimulating_Hormone": "lab_tsh",
    "Thrombin_Time": "lab_tt",
    "Uric_Acid": "lab_ua",
    "Blood_Urea_Nitrogen": "lab_urea",
    "White_Blood_Cell_Count": "lab_wbc",
    "Beta_Human_Chorionic_Gonadotropin": "lab_beta_hcg",
    "Beta_2_Microglobulin": "lab_beta_2_mg",
}


MODEL_FEATURES = [
    "Hematocrit", "Basophil_Count", "Eosinophil_Percentage", "Hemoglobin",
    "Lymphocyte_Percentage", "Monocyte_Percentage", "Neutrophil_Percentage",
    "Platelet_Count", "Red_Blood_Cell_Count", "Red_Cell_Distribution_Width",
    "White_Blood_Cell_Count", "Mean_Platelet_Volume", "Albumin",
    "Albumin_Globulin_Ratio", "Alkaline_Phosphatase",
    "Alanine_Aminotransferase", "Aspartate_Aminotransferase",
    "Gamma_Glutamyl_Transferase", "Globulin", "Total_Bilirubin", "Creatinine",
    "Uric_Acid", "Blood_Urea_Nitrogen", "Potassium", "Chloride", "Sodium",
    "Total_Calcium", "Phosphorus", "Glucose",
]


PANELS = {
    "cbc": [
        "Hemoglobin", "Hematocrit", "White_Blood_Cell_Count",
        "Red_Blood_Cell_Count", "Platelet_Count", "Red_Cell_Distribution_Width",
        "Mean_Platelet_Volume", "Lymphocyte_Percentage", "Monocyte_Percentage",
        "Neutrophil_Percentage", "Eosinophil_Percentage", "Basophil_Count",
        "Basophil_Percentage", "Lymphocyte_Count", "Monocyte_Count",
        "Neutrophil_Count", "Eosinophil_Count", "Mean_Corpuscular_Hemoglobin",
        "Mean_Corpuscular_Hemoglobin_Concentration", "Mean_Corpuscular_Volume",
        "Plateletcrit", "Platelet_Distribution_Width", "Reticulocyte_Count",
        "Reticulocyte_Percentage",
    ],
    "chemistry": [
        "Albumin", "Globulin", "Albumin_Globulin_Ratio",
        "Alanine_Aminotransferase", "Aspartate_Aminotransferase",
        "Alkaline_Phosphatase", "Creatinine", "Blood_Urea_Nitrogen",
        "Total_Bilirubin", "Direct_Bilirubin", "Indirect_Bilirubin",
        "Gamma_Glutamyl_Transferase", "Total_Protein", "Uric_Acid", "Sodium",
        "Potassium", "Chloride", "Glucose", "Lactate_Dehydrogenase",
        "Total_Calcium", "Magnesium", "Phosphorus",
        "Carbon_Dioxide_Combining_Power", "Total_Carbon_Dioxide", "Total_Bile_Acid",
    ],
    "tumor_markers": [
        "Alpha_Fetoprotein", "Cancer_Antigen_125", "Cancer_Antigen_153",
        "Cancer_Antigen_199", "Cancer_Antigen_724", "Carcinoembryonic_Antigen",
        "Cytokeratin_19_Fragment", "Neuron_Specific_Enolase",
        "Human_Epididymis_Protein_4", "Protein_Induced_by_Vitamin_K_Absence",
        "Pro_Gastrin_Releasing_Peptide", "Free_Prostate_Specific_Antigen",
        "Total_Prostate_Specific_Antigen", "Free_to_Total_PSA_Ratio",
        "Beta_Human_Chorionic_Gonadotropin", "Beta_2_Microglobulin",
        "Risk_of_Ovarian_Malignancy_Algorithm",
    ],
    "coagulation": [
        "Activated_Partial_Thromboplastin_Time", "Fibrinogen",
        "International_Normalized_Ratio", "Prothrombin_Time",
        "Prothrombin_Activity", "Prothrombin_Ratio", "Thrombin_Time", "D_Dimer",
    ],
    "lipids": [
        "Apolipoprotein_A1", "Apolipoprotein_B", "HDL_Cholesterol",
        "LDL_Cholesterol", "Lipoprotein_A", "Total_Cholesterol", "Triglycerides",
    ],
    "other_assays": [
        "Anti_Thyroglobulin_Antibody", "Anti_Thyroid_Microsomal_Antibody",
        "Cortisol", "Free_Triiodothyronine", "Free_Thyroxine",
        "Thyroid_Stimulating_Hormone", "Vitamin_B12", "Folate", "Ferritin",
        "Serum_Iron", "Prealbumin", "Glycated_Hemoglobin", "High_Sensitivity_CRP",
        "Cystatin_C", "Alpha_Amylase", "Cholinesterase", "Creatine_Kinase",
        "Cardiac_Troponin_I", "Alpha_Hydroxybutyrate_Dehydrogenase",
    ],
}

MIN_OBSERVED = {
    "cbc": 3,
    "chemistry": 3,
    "tumor_markers": 2,
    "coagulation": 2,
    "lipids": 2,
    "other_assays": 2,
}


OUTCOME_COMPONENTS = {
    "Outcome_NutriMetab": [195, 171],
    "Outcome_TumorBurden": [210, 212, 213, 215, 211],
    "Outcome_TreatComp": [161, 190, 186],
}


def fail_unless_release_valid() -> dict:
    validation = json.loads((RELEASE / "00_release_validation.json").read_text(encoding="utf-8"))
    if not validation.get("gate_passed"):
        raise RuntimeError("Formal database release validation did not pass")
    expected = {
        "patient_core_rows": 75248,
        "encounter_general_rows": 117157,
        "lab_observation_rows": 6996895,
    }
    for key, value in expected.items():
        if validation.get(key) != value:
            raise RuntimeError(f"Release count drift for {key}: {validation.get(key)} != {value}")
    return validation


def numeric_frame(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.apply(pd.to_numeric, errors="coerce").round(8)


def unique_signature_map(current: pd.DataFrame, formal_wide: pd.DataFrame, columns: list[str]) -> pd.Series:
    current_num = numeric_frame(current[columns])
    formal_num = formal_wide[[LAB_MAP[c] for c in columns]].round(8)
    sentinel = -9.87654321e99
    h_current = pd.util.hash_pandas_object(current_num.fillna(sentinel), index=False)
    h_formal = pd.util.hash_pandas_object(formal_num.fillna(sentinel), index=False)
    lookup = pd.DataFrame({"signature": h_formal.to_numpy(), "patient_key": formal_wide.index})
    lookup = lookup.loc[~lookup["signature"].duplicated(keep=False)].set_index("signature")["patient_key"]
    mapped = h_current.map(lookup)
    enough = current_num.notna().sum(axis=1) >= MIN_OBSERVED[next(k for k, v in PANELS.items() if v == columns)]
    return mapped.where(enough)


def make_outcome(core: pd.DataFrame, orders: list[int]) -> pd.Series:
    cols = [f"flag_legacy_final__{order}" for order in orders]
    values = core[cols].apply(pd.to_numeric, errors="coerce")
    complete = values.notna().all(axis=1)
    positive = values.eq(1).any(axis=1).astype(float)
    return positive.where(complete)


def patient_hash_split(patient_key: str) -> str:
    bucket = int(hashlib.sha256(patient_key.encode("utf-8")).hexdigest()[:8], 16) % 100
    if bucket < 70:
        return "train"
    if bucket < 85:
        return "validation"
    return "test"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    release_validation = fail_unless_release_valid()

    current = pd.read_csv(CURRENT_BASELINE, low_memory=False)
    current.insert(0, "current_row_id", np.arange(1, len(current) + 1))

    required_lab_ids = sorted(set(LAB_MAP.values()))
    summary_cols = [
        "patient_key", "final_variable_id", "first_numeric_value",
        "first_numeric_datetime", "first_record_key", "last_numeric_datetime",
    ]
    lab_long = pd.read_parquet(
        RELEASE / "162_patient_lab_summary_long_v1.parquet",
        filters=[("final_variable_id", "in", required_lab_ids)],
        columns=summary_cols,
    )
    formal_wide = lab_long.pivot(
        index="patient_key", columns="final_variable_id", values="first_numeric_value"
    )

    method_matches = {}
    for panel_name, panel_columns in PANELS.items():
        method_matches[panel_name] = unique_signature_map(current, formal_wide, panel_columns)
    method_df = pd.DataFrame(method_matches)
    candidate_count = method_df.notna().sum(axis=1)
    distinct_candidate_count = method_df.apply(lambda row: row.dropna().nunique(), axis=1)
    consensus = method_df.bfill(axis=1).iloc[:, 0].where(distinct_candidate_count == 1)

    linkage_status = np.select(
        [
            distinct_candidate_count > 1,
            (distinct_candidate_count == 1) & (candidate_count >= 2),
            (distinct_candidate_count == 1) & (candidate_count == 1),
        ],
        ["CONFLICT", "LINKED_MULTI_PANEL", "LINKED_SINGLE_PANEL"],
        default="UNLINKED",
    )

    core_min = pd.read_csv(
        RELEASE / "142_patient_nonlab_core_01na_v1.csv",
        usecols=[
            "patient_key", "first_observed_test_date", "last_observed_test_date",
            "observation_span_days", "encounter_count", "test_record_count",
            "age_at_first_observed_check", "sex_status",
        ],
        low_memory=False,
    ).set_index("patient_key")
    linked_core = core_min.reindex(consensus.to_numpy()).reset_index(drop=True)
    linked_age = pd.to_numeric(linked_core["age_at_first_observed_check"], errors="coerce")
    current_age = pd.to_numeric(current["Age"], errors="coerce")
    formal_sex = linked_core["sex_status"].map({"FEMALE": 0, "MALE": 1})
    current_sex = pd.to_numeric(current["Sex"], errors="coerce")

    linkage = pd.DataFrame({
        "current_row_id": current["current_row_id"],
        "patient_key": consensus,
        "linkage_status": linkage_status,
        "supporting_panel_count": candidate_count,
        "supporting_panels": method_df.apply(
            lambda row: ";".join(row.index[row.notna()].tolist()), axis=1
        ),
        "age_at_first_observed_check_matches": (
            (current_age == linked_age) | (current_age.isna() & linked_age.isna())
        ).where(consensus.notna()),
        "sex_matches": (
            (current_sex == formal_sex) | (current_sex.isna() & formal_sex.isna())
        ).where(consensus.notna()),
        "first_observed_test_date": linked_core["first_observed_test_date"],
        "last_observed_test_date": linked_core["last_observed_test_date"],
        "observation_span_days": linked_core["observation_span_days"],
        "encounter_count": linked_core["encounter_count"],
        "test_record_count": linked_core["test_record_count"],
    })
    linkage.to_csv(OUT / "current_table_to_formal_patient_linkage.csv", index=False, encoding="utf-8-sig")

    # Build a clean one-row-per-patient modeling table directly from the formal release.
    outcome_orders = sorted({order for orders in OUTCOME_COMPONENTS.values() for order in orders})
    core_cols = [
        "patient_key", "first_observed_test_date", "last_observed_test_date",
        "observation_span_days", "encounter_count", "test_record_count",
        *[f"flag_legacy_final__{order}" for order in outcome_orders],
    ]
    core = pd.read_csv(RELEASE / "142_patient_nonlab_core_01na_v1.csv", usecols=core_cols, low_memory=False)
    model_lab_wide = formal_wide[[LAB_MAP[feature] for feature in MODEL_FEATURES]].rename(
        columns={LAB_MAP[feature]: feature for feature in MODEL_FEATURES}
    )
    model = core.merge(model_lab_wide.reset_index(), on="patient_key", how="left", validate="one_to_one")
    for outcome, orders in OUTCOME_COMPONENTS.items():
        model[outcome] = make_outcome(model, orders)

    model["first_observed_test_date"] = pd.to_datetime(model["first_observed_test_date"], errors="coerce")
    model["last_observed_test_date"] = pd.to_datetime(model["last_observed_test_date"], errors="coerce")
    model["split_patient_hash"] = model["patient_key"].map(patient_hash_split)

    dated = model.loc[model["first_observed_test_date"].notna(), ["patient_key", "first_observed_test_date"]]
    dated = dated.sort_values(["first_observed_test_date", "patient_key"]).reset_index()
    cut70 = dated["first_observed_test_date"].quantile(0.70)
    cut85 = dated["first_observed_test_date"].quantile(0.85)
    dated["split_calendar_entry"] = np.select(
        [
            dated["first_observed_test_date"] <= cut70,
            dated["first_observed_test_date"] <= cut85,
        ],
        ["train", "validation"],
        default="test",
    )
    model["split_calendar_entry"] = "missing_date_excluded"
    model.loc[dated["index"], "split_calendar_entry"] = dated["split_calendar_entry"].to_numpy()

    keep = [
        "patient_key", "first_observed_test_date", "last_observed_test_date",
        "observation_span_days", "encounter_count", "test_record_count",
        "split_patient_hash", "split_calendar_entry", *MODEL_FEATURES,
        *OUTCOME_COMPONENTS.keys(),
    ]
    model[keep].to_csv(OUT / "hospital_patient_model_ready.csv", index=False, encoding="utf-8-sig")

    model_long = lab_long.loc[
        lab_long["final_variable_id"].isin([LAB_MAP[f] for f in MODEL_FEATURES]),
        [
            "patient_key", "final_variable_id", "first_numeric_value",
            "first_numeric_datetime", "last_numeric_datetime", "first_record_key",
        ],
    ].copy()
    reverse_map = {LAB_MAP[feature]: feature for feature in MODEL_FEATURES}
    model_long.insert(2, "model_feature", model_long["final_variable_id"].map(reverse_map))
    model_long.to_parquet(OUT / "hospital_predictor_measurement_times.parquet", index=False)

    encounters = pd.read_csv(
        RELEASE / "122_encounter_general_master_v1.csv",
        usecols=[
            "encounter_key", "patient_key", "first_test_date", "last_test_date",
            "test_record_count", "patient_source", "care_setting_from_department",
            "department_canonical",
        ],
        low_memory=False,
    )
    encounters.to_csv(OUT / "encounter_observation_windows.csv", index=False, encoding="utf-8-sig")

    plan = pd.read_csv(
        RELEASE / "supporting_materials" / "02_dictionaries_and_definitions"
        / "123_legacy60_status_and_candidate_variable_plan_v3.csv",
        low_memory=False,
    )
    provenance_rows = []
    for outcome, orders in OUTCOME_COMPONENTS.items():
        for order in orders:
            row = plan.loc[plan["original_order"] == order].iloc[0]
            provenance_rows.append({
                "outcome": outcome,
                "original_order": order,
                "original_variable": row["original_variable"],
                "candidate_variable": row["formal_or_candidate_new_variable"],
                "candidate_readiness": row["candidate_readiness"],
                "formal_research_use": row["old_flag_formal_research_use"],
                "validated_event_date_available": False,
                "time_interpretation": "Patient-ever legacy flag; no validated event date",
            })
    provenance = pd.DataFrame(provenance_rows)
    provenance.to_csv(OUT / "outcome_provenance_gate.csv", index=False, encoding="utf-8-sig")

    linkage_counts = pd.Series(linkage_status).value_counts().to_dict()
    linked = linkage.loc[linkage["patient_key"].notna()]
    duplicate_patient_rows = int(linked.duplicated("patient_key", keep=False).sum())
    duplicate_patient_keys = int(linked.loc[linked.duplicated("patient_key", keep=False), "patient_key"].nunique())
    audit = {
        "release_validation": release_validation,
        "current_table_rows": int(len(current)),
        "linkage_status_counts": {str(k): int(v) for k, v in linkage_counts.items()},
        "linked_distinct_patient_keys": int(linked["patient_key"].nunique()),
        "linked_rows_in_duplicate_patient_groups": duplicate_patient_rows,
        "duplicate_formal_patient_keys": duplicate_patient_keys,
        "age_agreement_among_linked": float(linkage.loc[linkage["patient_key"].notna(), "age_at_first_observed_check_matches"].mean()),
        "sex_agreement_among_linked": float(linkage.loc[linkage["patient_key"].notna(), "sex_matches"].mean()),
        "formal_patient_rows": int(len(model)),
        "calendar_split_counts": {str(k): int(v) for k, v in model["split_calendar_entry"].value_counts().items()},
        "patient_hash_split_counts": {str(k): int(v) for k, v in model["split_patient_hash"].value_counts().items()},
        "outcomes": {
            outcome: {
                "nonmissing": int(model[outcome].notna().sum()),
                "positive": int(model[outcome].eq(1).sum()),
                "prevalence": float(model[outcome].mean()),
            }
            for outcome in OUTCOME_COMPONENTS
        },
        "availability_gate": {
            "pseudonymous_patient_key": True,
            "encounter_key": True,
            "lab_collection_datetime": True,
            "encounter_first_last_test_date": True,
            "true_admission_datetime": False,
            "true_discharge_datetime": False,
            "validated_outcome_event_datetime": False,
            "calendar_entry_split_is_patient_disjoint": True,
            "calendar_entry_split_is_prospective_outcome_validation": False,
        },
        "interpretation": (
            "Patient identity and laboratory timing are recoverable. Encounter windows are test-observation "
            "windows, not admission/discharge. Current outcomes are patient-ever legacy flags without "
            "validated event dates and are prohibited as formal research outcomes by the locked release "
            "dictionary. Retraining is therefore an audit/sensitivity analysis, not a prospective model."
        ),
    }
    (OUT / "hospital_data_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
