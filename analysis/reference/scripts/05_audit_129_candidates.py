from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


RELEASE = Path("hospital_source")
PACKAGE = Path("reference")
OUT = PACKAGE / "outputs"

ROLE_FILE = RELEASE / "analysis" / "coupling_atlas_nc_v1" / "config" / "lab_variable_roles_v1.csv"
SUMMARY_FILE = RELEASE / "162_patient_lab_summary_long_v1.parquet"
BASE_MODEL_FILE = OUT / "hospital_patient_model_ready.csv"

EXPECTED_PATIENTS = 75_248
EXPECTED_LABS = 129
SPLITS = ("train", "validation", "test")
OUTCOMES = ("Outcome_NutriMetab", "Outcome_TumorBurden", "Outcome_TreatComp")


def missing_tier(missing_pct: float, is_numeric: bool) -> str:
    if not is_numeric or pd.isna(missing_pct):
        return "C"
    if missing_pct <= 30.0:
        return "A"
    if missing_pct <= 50.0:
        return "B"
    return "C"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    roles = pd.read_csv(ROLE_FILE, low_memory=False)
    if len(roles) != EXPECTED_LABS or roles["final_variable_id"].nunique() != EXPECTED_LABS:
        raise RuntimeError("The locked laboratory registry is not 129 unique variables")

    base = pd.read_csv(BASE_MODEL_FILE, low_memory=False)
    if len(base) != EXPECTED_PATIENTS or base["patient_key"].nunique() != EXPECTED_PATIENTS:
        raise RuntimeError("The patient-level base table is not 75,248 unique patients")
    if base["patient_key"].isna().any():
        raise RuntimeError("Missing patient_key in patient-level base table")

    lab = pd.read_parquet(
        SUMMARY_FILE,
        columns=[
            "patient_key",
            "final_variable_id",
            "first_numeric_value",
            "first_numeric_datetime",
            "first_value_status",
        ],
    )
    if lab.duplicated(["patient_key", "final_variable_id"]).any():
        raise RuntimeError("Duplicate patient-variable pairs in laboratory summary")
    unknown_ids = sorted(set(lab["final_variable_id"]) - set(roles["final_variable_id"]))
    if unknown_ids:
        raise RuntimeError(f"Laboratory summary contains unregistered IDs: {unknown_ids[:5]}")

    wide = lab.pivot(
        index="patient_key", columns="final_variable_id", values="first_numeric_value"
    ).reindex(columns=roles["final_variable_id"].tolist())
    wide.columns.name = None
    expanded = base[
        [
            "patient_key",
            "first_observed_test_date",
            "last_observed_test_date",
            "observation_span_days",
            "encounter_count",
            "test_record_count",
            "split_patient_hash",
            "split_calendar_entry",
            *OUTCOMES,
        ]
    ].merge(wide.reset_index(), on="patient_key", how="left", validate="one_to_one")

    expanded.to_csv(OUT / "hospital_patient_model_ready_129.csv", index=False, encoding="utf-8-sig")

    split_rows: list[pd.DataFrame] = []
    for split_name in ("overall", *SPLITS):
        part = expanded if split_name == "overall" else expanded.loc[
            expanded["split_calendar_entry"] == split_name
        ]
        nonmissing = part[roles["final_variable_id"]].notna().sum(axis=0)
        split_rows.append(
            pd.DataFrame(
                {
                    "split": split_name,
                    "final_variable_id": roles["final_variable_id"],
                    "n_patients": len(part),
                    "n_nonmissing_first_numeric": nonmissing.reindex(
                        roles["final_variable_id"]
                    ).to_numpy(),
                }
            )
        )
    split_long = pd.concat(split_rows, ignore_index=True)
    split_long["missing_pct"] = (
        100.0
        * (split_long["n_patients"] - split_long["n_nonmissing_first_numeric"])
        / split_long["n_patients"].replace(0, np.nan)
    )
    split_long.to_csv(
        OUT / "hospital_129_missingness_by_split_long.csv", index=False, encoding="utf-8-sig"
    )

    split_wide = split_long.pivot(
        index="final_variable_id", columns="split", values="missing_pct"
    ).rename(columns=lambda x: f"{x}_missing_pct")
    split_n = split_long.pivot(
        index="final_variable_id", columns="split", values="n_nonmissing_first_numeric"
    ).rename(columns=lambda x: f"{x}_nonmissing_n")
    audit = roles.merge(split_wide.reset_index(), on="final_variable_id", how="left")
    audit = audit.merge(split_n.reset_index(), on="final_variable_id", how="left")
    audit["value_type_is_numeric"] = audit["value_type"].eq("NUMERIC")
    audit["hospital_missingness_tier"] = [
        missing_tier(m, bool(n))
        for m, n in zip(audit["train_missing_pct"], audit["value_type_is_numeric"])
    ]
    audit["calendar_missingness_range_pp"] = audit[
        ["train_missing_pct", "validation_missing_pct", "test_missing_pct"]
    ].max(axis=1) - audit[
        ["train_missing_pct", "validation_missing_pct", "test_missing_pct"]
    ].min(axis=1)
    audit["calendar_missingness_drift_gt20pp"] = (
        audit["calendar_missingness_range_pp"] > 20.0
    )
    audit.to_csv(OUT / "lab_129_hospital_audit.csv", index=False, encoding="utf-8-sig")

    outcome_rows: list[pd.DataFrame] = []
    train = expanded.loc[expanded["split_calendar_entry"] == "train"]
    for outcome in OUTCOMES:
        eligible = train.loc[train[outcome].notna()]
        nonmissing = eligible[roles["final_variable_id"]].notna().sum(axis=0)
        frame = pd.DataFrame(
            {
                "outcome": outcome,
                "final_variable_id": roles["final_variable_id"],
                "n_outcome_evaluable_train": len(eligible),
                "n_nonmissing_first_numeric": nonmissing.reindex(
                    roles["final_variable_id"]
                ).to_numpy(),
            }
        )
        frame["missing_pct"] = 100.0 * (
            frame["n_outcome_evaluable_train"] - frame["n_nonmissing_first_numeric"]
        ) / frame["n_outcome_evaluable_train"]
        outcome_rows.append(frame)
    pd.concat(outcome_rows, ignore_index=True).to_csv(
        OUT / "hospital_129_outcome_train_missingness.csv",
        index=False,
        encoding="utf-8-sig",
    )

    qa = {
        "patient_rows": int(len(expanded)),
        "unique_patient_keys": int(expanded["patient_key"].nunique()),
        "registered_laboratory_variables": int(len(roles)),
        "laboratory_summary_rows": int(len(lab)),
        "patient_variable_pair_duplicates": int(
            lab.duplicated(["patient_key", "final_variable_id"]).sum()
        ),
        "split_counts": {
            str(k): int(v)
            for k, v in expanded["split_calendar_entry"].value_counts().items()
        },
        "hospital_missingness_tier_counts": {
            str(k): int(v)
            for k, v in audit["hospital_missingness_tier"].value_counts().items()
        },
        "numeric_variable_count": int(audit["value_type_is_numeric"].sum()),
        "calendar_missingness_drift_gt20pp_count": int(
            audit["calendar_missingness_drift_gt20pp"].sum()
        ),
        "eligibility_basis": (
            "Patient-level first numeric value; A <=30% missing in training, "
            "B >30%-50%, C >50% or nonnumeric"
        ),
    }
    (OUT / "lab_129_hospital_audit_qa.json").write_text(
        json.dumps(qa, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(qa, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
