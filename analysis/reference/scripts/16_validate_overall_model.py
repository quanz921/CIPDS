from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


PACKAGE_DIR = Path(
    os.environ.get(
        "CIPDS_PACKAGE_DIR_PY",
        "reference",
    )
)
OUT = PACKAGE_DIR / "outputs"


def require(condition: bool, name: str, detail: str, checks: list[dict]) -> None:
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


def close(a: float, b: float, tol: float = 1e-10) -> bool:
    return math.isfinite(a) and math.isfinite(b) and abs(a - b) <= tol


def ranked_positive_average_precision(y, score):
    return float(average_precision_score(y, score))


required_files = [
    "overall_oof_reproduction_checks.csv",
    "overall_frozen_component_reproduction_checks.csv",
    "overall_meta_training_oof.csv",
    "overall_hospital_outcome_state_counts.csv",
    "overall_calendar_performance.csv",
    "overall_calendar_predictions.csv",
    "overall_calendar_paired_comparison.csv",
    "overall_calendar_calibration.csv",
    "overall_stacked_model_coefficients.csv",
    "overall_stacked_ordinal_model.rds",
    "overall_stacked_run_manifest.json",
    "overall_nhanes_scores.csv",
    "overall_nhanes_coverage.csv",
    "overall_nhanes_allcause.csv",
    "overall_nhanes_cancer_history_interaction.csv",
    "overall_nhanes_phenoage_incremental_tests.csv",
    "overall_nhanes_phenoage_incremental_effects.csv",
    "overall_nhanes_component_incremental_tests.csv",
    "overall_nhanes_component_incremental_effects.csv",
    "overall_nhanes_cause_specific.csv",
    "overall_nhanes_design_audit.csv",
]

checks: list[dict] = []
missing = [name for name in required_files if not (OUT / name).exists()]
require(not missing, "required_outputs_present", f"missing={missing}", checks)
if missing:
    raise SystemExit(f"Missing required outputs: {missing}")

manifest = json.loads((OUT / "overall_stacked_run_manifest.json").read_text(encoding="utf-8-sig"))
require(
    manifest.get("base_models_retrained") is False,
    "base_models_remained_frozen",
    f"base_models_retrained={manifest.get('base_models_retrained')}",
    checks,
)
require(
    manifest.get("phenoage_used_as_input") is False
    and manifest.get("nhanes_used_for_training") is False,
    "no_phenoage_or_nhanes_training_input",
    f"pheno={manifest.get('phenoage_used_as_input')}; nhanes={manifest.get('nhanes_used_for_training')}",
    checks,
)
require(
    manifest.get("no_auc_weighted_composite") is True,
    "no_auc_weighted_composite",
    f"value={manifest.get('no_auc_weighted_composite')}",
    checks,
)
require(
    int(manifest.get("thread_count", -1)) == 23,
    "thread_count_23",
    f"thread_count={manifest.get('thread_count')}",
    checks,
)

oof_check = pd.read_csv(OUT / "overall_oof_reproduction_checks.csv")
require(
    len(oof_check) == 15 and oof_check["reproduced_within_tolerance"].astype(bool).all(),
    "all_15_outer_predictions_reproduced",
    f"rows={len(oof_check)}; max_auc_diff={oof_check['absolute_auc_difference'].max():.3g}",
    checks,
)

frozen_check = pd.read_csv(OUT / "overall_frozen_component_reproduction_checks.csv")
require(
    len(frozen_check) == 3 and frozen_check["within_tolerance"].astype(bool).all(),
    "frozen_test_component_scores_reproduced",
    f"rows={len(frozen_check)}; max_diff={frozen_check['max_absolute_probability_difference'].max():.3g}",
    checks,
)

meta = pd.read_csv(OUT / "overall_meta_training_oof.csv")
expected_meta_counts = {0: 39248, 1: 8126, 2: 1547, 3: 161}
observed_meta_counts = meta["burden"].value_counts().sort_index().to_dict()
require(
    len(meta) == 49082 and meta["patient_key"].nunique() == 49082,
    "meta_training_patient_grain",
    f"n={len(meta)}; unique={meta['patient_key'].nunique()}",
    checks,
)
require(
    observed_meta_counts == expected_meta_counts,
    "meta_training_burden_counts",
    f"observed={observed_meta_counts}",
    checks,
)
require(
    meta[["p_NM", "p_TB", "p_TC"]].notna().all().all()
    and ((meta[["p_NM", "p_TB", "p_TC"]] >= 0) & (meta[["p_NM", "p_TB", "p_TC"]] <= 1)).all().all(),
    "meta_component_probability_validity",
    "all OOF component probabilities are finite and within [0,1]",
    checks,
)

pred = pd.read_csv(OUT / "overall_calendar_predictions.csv")
test = pred.loc[pred["split"] == "calendar_test"].copy()
validation = pred.loc[pred["split"] == "calendar_validation"].copy()
test_counts = test["observed_burden"].value_counts().sort_index().to_dict()
require(
    len(test) == 11125 and test["patient_key"].nunique() == 11125,
    "calendar_test_patient_grain",
    f"n={len(test)}; unique={test['patient_key'].nunique()}",
    checks,
)
require(
    test_counts == {0: 9128, 1: 1717, 2: 259, 3: 21},
    "calendar_test_burden_counts",
    f"observed={test_counts}",
    checks,
)
prob_cols = [f"probability_burden_{i}" for i in range(4)]
prob_sum_error = np.abs(test[prob_cols].sum(axis=1).to_numpy() - 1).max()
expected_recomputed = sum(i * test[f"probability_burden_{i}"] for i in range(4))
expected_error = np.abs(expected_recomputed - test["overall_expected_burden"]).max()
require(
    prob_sum_error < 1e-10 and expected_error < 1e-10,
    "ordinal_probability_internal_consistency",
    f"max_probability_sum_error={prob_sum_error:.3g}; max_expected_burden_error={expected_error:.3g}",
    checks,
)

any_y = (test["observed_burden"] >= 1).astype(int)
multi_y = (test["observed_burden"] >= 2).astype(int)
any_auc = roc_auc_score(any_y, test["probability_any_domain"])
any_ap = ranked_positive_average_precision(any_y, test["probability_any_domain"])
any_ap_sklearn = average_precision_score(any_y, test["probability_any_domain"])
multi_auc = roc_auc_score(multi_y, test["probability_multidomain"])
multi_ap = ranked_positive_average_precision(multi_y, test["probability_multidomain"])
multi_ap_sklearn = average_precision_score(multi_y, test["probability_multidomain"])
perf = pd.read_csv(OUT / "overall_calendar_performance.csv")
reported_test = perf.loc[perf["split"] == "calendar_test"].iloc[0]
require(
    close(any_auc, float(reported_test["any_auc"]))
    and close(any_ap, float(reported_test["any_average_precision"]))
    and close(multi_auc, float(reported_test["multidomain_auc"]))
    and close(multi_ap, float(reported_test["multidomain_average_precision"])),
    "calendar_test_metrics_independently_recomputed",
    f"any_auc={any_auc:.6f}; any_ap={any_ap:.6f}; multi_auc={multi_auc:.6f}; multi_ap={multi_ap:.6f}",
    checks,
)

paired = pd.read_csv(OUT / "overall_calendar_paired_comparison.csv")
require(
    (paired["auc_difference_learned_minus_equal_sum"] < 0).all(),
    "no_false_superiority_over_equal_sum",
    "learned stack did not outperform equal-sum comparator; superiority language is prohibited",
    checks,
)

nhanes = pd.read_csv(OUT / "overall_nhanes_scores.csv")
require(
    len(nhanes) == 44772 and nhanes["SEQN"].nunique() == 44772,
    "nhanes_score_patient_grain",
    f"n={len(nhanes)}; unique={nhanes['SEQN'].nunique()}",
    checks,
)
valid_overall = nhanes["Overall_expected_burden"].dropna()
require(
    len(valid_overall) > 0 and valid_overall.between(0, 3).all(),
    "nhanes_overall_score_range",
    f"nonmissing={len(valid_overall)}; min={valid_overall.min():.6g}; max={valid_overall.max():.6g}",
    checks,
)

coverage = pd.read_csv(OUT / "overall_nhanes_coverage.csv")
age60_cov = coverage.loc[coverage["domain"] == "age60"].iloc[0]
age65_cov = coverage.loc[coverage["domain"] == "age65"].iloc[0]
require(
    int(age60_cov["source_n"]) == 15048
    and int(age60_cov["overall_score_nonmissing_n"]) == 14715
    and int(age60_cov["overall_score_missing_n"]) == 333,
    "age60_coverage_gate",
    age60_cov.to_json(),
    checks,
)
require(
    int(age65_cov["source_n"]) == 11013
    and int(age65_cov["overall_score_nonmissing_n"]) == 10763
    and int(age65_cov["overall_score_missing_n"]) == 250,
    "age65_coverage_gate",
    age65_cov.to_json(),
    checks,
)

allcause = pd.read_csv(OUT / "overall_nhanes_allcause.csv")
primary_overall = allcause.loc[
    (allcause["domain"] == "age60")
    & (allcause["analysis_role"] == "PRIMARY")
    & (allcause["population"] == "Overall")
].iloc[0]
require(
    int(primary_overall["n"]) == 9753 and int(primary_overall["events"]) == 3252,
    "primary_age60_analytic_counts",
    f"n={primary_overall['n']}; events={primary_overall['events']}",
    checks,
)
require(
    float(primary_overall["hazard_ratio_per_domain_weighted_sd"]) > 1
    and float(primary_overall["ci_lower"]) > 1
    and float(primary_overall["p_value"]) < 0.05,
    "primary_age60_overall_association",
    f"HR={primary_overall['hazard_ratio_per_domain_weighted_sd']}; CI={primary_overall['ci_lower']}-{primary_overall['ci_upper']}; p={primary_overall['p_value']}",
    checks,
)

interaction = pd.read_csv(OUT / "overall_nhanes_cancer_history_interaction.csv")
age60_int = interaction.loc[interaction["domain"] == "age60"].iloc[0]
require(
    int(age60_int["n"]) == 9753
    and math.isfinite(float(age60_int["interaction_p"])),
    "cancer_history_interaction_estimable",
    f"interaction_p={age60_int['interaction_p']}",
    checks,
)

pheno_tests = pd.read_csv(OUT / "overall_nhanes_phenoage_incremental_tests.csv")
pheno_primary = pheno_tests.loc[
    (pheno_tests["domain"] == "age60") & (pheno_tests["population"] == "Overall")
]
require(
    len(pheno_primary) == 4
    and (pheno_primary["n"] == 6987).all()
    and (pheno_primary["events"] == 2650).all(),
    "phenoage_common_complete_case_counts",
    f"rows={len(pheno_primary)}; n={sorted(pheno_primary['n'].unique())}; events={sorted(pheno_primary['events'].unique())}",
    checks,
)
p_overall_beyond_pheno = float(
    pheno_primary.loc[pheno_primary["comparison"] == "M1_to_M3", "design_adjusted_block_wald_p"].iloc[0]
)
p_pheno_beyond_overall = float(
    pheno_primary.loc[pheno_primary["comparison"] == "M2_to_M3", "design_adjusted_block_wald_p"].iloc[0]
)
require(
    p_overall_beyond_pheno < 0.05 and p_pheno_beyond_overall < 0.05,
    "overall_and_phenoage_bidirectional_incremental_value",
    f"overall_beyond_pheno_p={p_overall_beyond_pheno:.6g}; pheno_beyond_overall_p={p_pheno_beyond_overall:.6g}",
    checks,
)

component_tests = pd.read_csv(OUT / "overall_nhanes_component_incremental_tests.csv")
component_primary = component_tests.loc[
    (component_tests["domain"] == "age60") & (component_tests["population"] == "Overall")
]
p_overall_beyond_components = float(
    component_primary.loc[
        component_primary["comparison"] == "COMPONENTS_to_BOTH",
        "design_adjusted_block_wald_p",
    ].iloc[0]
)
p_components_beyond_overall = float(
    component_primary.loc[
        component_primary["comparison"] == "OVERALL_to_BOTH",
        "design_adjusted_block_wald_p",
    ].iloc[0]
)
require(
    p_components_beyond_overall < 0.05,
    "three_components_retain_information_beyond_overall",
    f"components_beyond_overall_p={p_components_beyond_overall:.6g}",
    checks,
)
require(
    p_overall_beyond_components >= 0.05,
    "overall_does_not_replace_components",
    f"overall_beyond_components_p={p_overall_beyond_components:.6g}; overall must be reported with components",
    checks,
)

causes = pd.read_csv(OUT / "overall_nhanes_cause_specific.csv")
key_causes = causes.loc[
    (causes["domain"] == "age60")
    & (causes["population"] == "Overall")
    & causes["mortality_outcome"].isin(["Death_Cancer", "Death_CVD"])
]
require(
    len(key_causes) == 2
    and np.isfinite(key_causes[["hazard_ratio_per_domain_weighted_sd", "ci_lower", "ci_upper", "p_value"]].to_numpy()).all()
    and (key_causes["hazard_ratio_per_domain_weighted_sd"] > 1).all()
    and (key_causes["p_value"] < 0.05).all(),
    "key_cause_specific_results_valid",
    key_causes[["mortality_outcome", "n", "events", "hazard_ratio_per_domain_weighted_sd", "p_value"]].to_json(orient="records"),
    checks,
)

design = pd.read_csv(OUT / "overall_nhanes_design_audit.csv")
design_map = dict(zip(design["item"], design["value"]))
require(
    design_map.get("full_design_before_age_domain") == "YES"
    and design_map.get("pooled_weight_formula") == "WTMEC2YR/9",
    "complex_survey_domain_order",
    f"full_design={design_map.get('full_design_before_age_domain')}; weight={design_map.get('pooled_weight_formula')}",
    checks,
)
require(
    design_map.get("phenoage_used_as_overall_input") == "FALSE"
    and design_map.get("nhanes_outcomes_used_for_training") == "FALSE",
    "overall_model_provenance_lock",
    "PhenoAge and NHANES outcomes excluded from overall-model training",
    checks,
)

model_sha256 = hashlib.sha256((OUT / "overall_stacked_ordinal_model.rds").read_bytes()).hexdigest()
all_passed = all(row["passed"] for row in checks)

with (OUT / "overall_qa_checks.csv").open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=["check", "passed", "detail"])
    writer.writeheader()
    writer.writerows(checks)

summary = {
    "all_checks_passed": all_passed,
    "checks_total": len(checks),
    "checks_passed": sum(row["passed"] for row in checks),
    "checks_failed": sum(not row["passed"] for row in checks),
    "model_sha256": model_sha256,
    "hospital_calendar_test": {
        "n": len(test),
        "burden_counts": {str(k): int(v) for k, v in test_counts.items()},
        "any_domain_auc": any_auc,
        "any_domain_average_precision": any_ap,
        "any_domain_average_precision_sklearn_tie_definition": any_ap_sklearn,
        "multidomain_auc": multi_auc,
        "multidomain_average_precision": multi_ap,
        "multidomain_average_precision_sklearn_tie_definition": multi_ap_sklearn,
    },
    "primary_age60_overall": primary_overall.to_dict(),
    "primary_age60_cancer_history_interaction": age60_int.to_dict(),
    "primary_age60_phenoage_incremental": {
        "overall_beyond_phenoage_p": p_overall_beyond_pheno,
        "phenoage_beyond_overall_p": p_pheno_beyond_overall,
    },
    "primary_age60_component_reduction": {
        "overall_beyond_three_components_p": p_overall_beyond_components,
        "three_components_beyond_overall_p": p_components_beyond_overall,
    },
    "interpretation_locks": [
        "Overall score is trained only from hospital cross-fitted NM, TB and TC probabilities.",
        "PhenoAge and NHANES mortality are not overall-model inputs.",
        "The learned stack did not outperform the equal-sum comparator on hospital calendar test discrimination.",
        "The overall scalar summarizes multidomain burden but does not replace separate NM, TB and TC reporting.",
        "No AUC-weighted composite was generated.",
        "Hospital legacy outcomes remain undated full-observation-window flags; no future-event prediction claim is allowed.",
    ],
}
(OUT / "overall_qa_summary.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2, default=str),
    encoding="utf-8",
)

print(json.dumps(summary, ensure_ascii=False, indent=2, default=str))
if not all_passed:
    failed = [row for row in checks if not row["passed"]]
    raise SystemExit(f"Overall-model QA failed: {failed}")
