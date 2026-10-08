from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE_DIR = Path(
    os.environ.get(
        "CIPDS_PACKAGE_DIR_PY",
        "reference",
    )
)
OUT = PACKAGE_DIR / "outputs"

checks: list[dict[str, object]] = []


def require(condition: bool, name: str, detail: str) -> None:
    checks.append({"check": name, "passed": bool(condition), "detail": detail})


required_files = [
    "nhanes_paired_discrimination_primary_cohort.csv",
    "nhanes_paired_discrimination_estimates.csv",
    "nhanes_paired_discrimination_comparisons.csv",
    "nhanes_paired_discrimination_cohort_audit.csv",
    "nhanes_paired_discrimination_timeROC_crosscheck.csv",
    "nhanes_paired_discrimination_replicates.rds",
    "nhanes_paired_discrimination_manifest.json",
]
missing = [name for name in required_files if not (OUT / name).exists()]
require(not missing, "required_outputs_present", f"missing={missing}")
if missing:
    raise SystemExit(f"Missing paired-discrimination outputs: {missing}")

manifest = json.loads(
    (OUT / "nhanes_paired_discrimination_manifest.json").read_text(encoding="utf-8-sig")
)
estimates = pd.read_csv(OUT / "nhanes_paired_discrimination_estimates.csv")
comparisons = pd.read_csv(OUT / "nhanes_paired_discrimination_comparisons.csv")
cohorts = pd.read_csv(OUT / "nhanes_paired_discrimination_cohort_audit.csv")
crosscheck = pd.read_csv(OUT / "nhanes_paired_discrimination_timeROC_crosscheck.csv")
primary = pd.read_csv(OUT / "nhanes_paired_discrimination_primary_cohort.csv")

require(
    manifest["bootstrap_replicates"] == 1000 and 1 <= manifest["threads"] <= 23,
    "formal_1000_replicate_run_within_23_worker_budget",
    f"replicates={manifest['bootstrap_replicates']}; threads={manifest['threads']}",
)
require(
    manifest["full_survey_design_created_before_domain_restriction"] is True,
    "full_design_precedes_domain_restriction",
    str(manifest["full_survey_design_created_before_domain_restriction"]),
)
require(
    manifest["score_models_retrained"] is False,
    "frozen_scores_not_retrained",
    str(manifest["score_models_retrained"]),
)

require(
    len(cohorts) == 6
    and set(cohorts["domain"]) == {"age60", "age65"}
    and set(cohorts["population"]) == {"Overall", "Cancer", "Noncancer"},
    "all_domains_and_subgroups_present",
    f"rows={len(cohorts)}",
)
primary_audit = cohorts.loc[
    (cohorts["domain"] == "age60") & (cohorts["population"] == "Overall")
].iloc[0]
require(
    int(primary_audit["n"]) == 11410 and int(primary_audit["events"]) == 5272,
    "primary_common_cohort_counts",
    f"n={primary_audit['n']}; events={primary_audit['events']}",
)
require(
    len(primary) == 11410
    and primary["SEQN"].nunique() == 11410
    and int(primary["Death_AllCause"].sum()) == 5272,
    "primary_participant_grain",
    f"rows={len(primary)}; unique={primary['SEQN'].nunique()}; events={primary['Death_AllCause'].sum()}",
)
score_columns = [
    "Overall_expected_burden",
    "PhenoAge_acceleration",
    "NM_component",
    "TB_component",
    "TC_component",
]
require(
    primary[score_columns].notna().all().all()
    and np.isfinite(primary[score_columns].to_numpy()).all(),
    "all_primary_scores_finite_on_common_cohort",
    "all five scores are complete and finite",
)

require(
    len(estimates) == 150
    and estimates[["domain", "population", "metric", "model"]].drop_duplicates().shape[0]
    == 150,
    "estimate_grid_complete",
    f"rows={len(estimates)}",
)
require(
    estimates["estimate"].between(0, 1).all()
    and estimates["ci_lower"].between(0, 1).all()
    and estimates["ci_upper"].between(0, 1).all()
    and (estimates["ci_lower"] <= estimates["estimate"]).all()
    and (estimates["estimate"] <= estimates["ci_upper"]).all(),
    "estimate_ranges_and_ci_containment",
    f"estimate_range={estimates['estimate'].min():.6f}-{estimates['estimate'].max():.6f}",
)
require(
    len(comparisons) == 210
    and comparisons[[
        "domain", "population", "metric", "comparison_family", "model_a", "model_b"
    ]].drop_duplicates().shape[0]
    == 210,
    "paired_comparison_grid_complete",
    f"rows={len(comparisons)}",
)
max_difference_error = np.abs(
    comparisons["paired_difference"]
    - (comparisons["estimate_a"] - comparisons["estimate_b"])
).max()
require(
    max_difference_error < 1e-12,
    "paired_difference_orientation_recomputed",
    f"max_error={max_difference_error:.3g}",
)
require(
    comparisons["paired_p_value"].between(0, 1).all()
    and (comparisons["difference_ci_lower"] <= comparisons["paired_difference"]).all()
    and (comparisons["paired_difference"] <= comparisons["difference_ci_upper"]).all(),
    "paired_p_values_and_ci_valid",
    "all P values and confidence intervals are valid",
)
require(
    float(crosscheck["absolute_difference"].max()) < 1e-4,
    "time_dependent_auc_matches_timeROC",
    f"max_absolute_difference={crosscheck['absolute_difference'].max():.8f}",
)


def weighted_censoring_survival(
    time: np.ndarray, event: np.ndarray, weight: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(time, kind="stable")
    time_sorted = time[order]
    event_sorted = event[order]
    weight_sorted = weight[order]
    unique_time, group = np.unique(time_sorted, return_inverse=True)
    total = np.bincount(group, weights=weight_sorted)
    censor = np.bincount(group, weights=weight_sorted * (event_sorted == 0))
    risk = np.cumsum(total[::-1])[::-1]
    hazard = np.divide(censor, risk, out=np.zeros_like(censor), where=risk > 0)
    after = np.cumprod(1 - np.clip(hazard, 0, 1))
    before = np.r_[1.0, after[:-1]]
    return unique_time, before, after


def weighted_auc(
    time: np.ndarray,
    event: np.ndarray,
    score: np.ndarray,
    weight: np.ndarray,
    horizon: float,
) -> float:
    cases = (event == 1) & (time < horizon)
    controls = time > horizon
    unique_time, before, after = weighted_censoring_survival(time, event, weight)
    case_index = np.searchsorted(unique_time, time[cases])
    case_weight = weight[cases] / np.maximum(before[case_index], 1e-8)
    horizon_index = np.searchsorted(unique_time, horizon, side="right") - 1
    g_horizon = 1.0 if horizon_index < 0 else after[horizon_index]
    control_weight = weight[controls] / max(g_horizon, 1e-8)
    control_score = score[controls]
    order = np.argsort(control_score, kind="stable")
    control_score = control_score[order]
    control_weight = control_weight[order]
    cumulative = np.cumsum(control_weight)
    case_score = score[cases]
    left = np.searchsorted(control_score, case_score, side="left")
    right = np.searchsorted(control_score, case_score, side="right")
    less = np.where(left > 0, cumulative[np.maximum(left - 1, 0)], 0.0)
    less_or_equal = np.where(right > 0, cumulative[np.maximum(right - 1, 0)], 0.0)
    equal = less_or_equal - less
    wins = less + 0.5 * equal
    return float(
        np.sum(case_weight * wins) / (np.sum(case_weight) * np.sum(control_weight))
    )


def weighted_concordance(
    time: np.ndarray,
    event: np.ndarray,
    score: np.ndarray,
    weight: np.ndarray,
    *,
    ymax: float | None,
    uno: bool,
) -> float:
    unique_time, before, _ = weighted_censoring_survival(time, event, weight)
    event_index = np.flatnonzero(event == 1)
    if ymax is not None:
        event_index = event_index[time[event_index] <= ymax]
    numerator = 0.0
    denominator = 0.0
    for index in event_index:
        comparable = time > time[index]
        if not np.any(comparable):
            continue
        pair_weight = weight[index] * weight[comparable]
        if uno:
            time_index = np.searchsorted(unique_time, time[index])
            pair_weight = pair_weight / max(before[time_index] ** 2, 1e-16)
        comparison = score[index] - score[comparable]
        numerator += float(
            np.sum(pair_weight * ((comparison > 0) + 0.5 * (comparison == 0)))
        )
        denominator += float(np.sum(pair_weight))
    return numerator / denominator


time = primary["Follow_Up_Years"].to_numpy(float)
event = primary["Death_AllCause"].to_numpy(int)
weight = primary["pooled_mec_weight"].to_numpy(float)
independent_rows: list[dict[str, object]] = []
for model, column in {
    "OVERALL": "Overall_expected_burden",
    "PHENO": "PhenoAge_acceleration",
}.items():
    score = primary[column].to_numpy(float)
    for horizon in (5, 10, 15):
        value = weighted_auc(time, event, score, weight, horizon)
        reported = estimates.loc[
            (estimates["domain"] == "age60")
            & (estimates["population"] == "Overall")
            & (estimates["metric"] == f"auc{horizon}")
            & (estimates["model"] == model),
            "estimate",
        ].iloc[0]
        independent_rows.append(
            {
                "model": model,
                "metric": f"auc{horizon}",
                "independent_value": value,
                "reported_value": reported,
                "absolute_difference": abs(value - reported),
            }
        )
    for metric, ymax, uno in (
        ("uno_c15", 15.0, True),
        ("harrell_c", None, False),
    ):
        value = weighted_concordance(
            time, event, score, weight, ymax=ymax, uno=uno
        )
        reported = estimates.loc[
            (estimates["domain"] == "age60")
            & (estimates["population"] == "Overall")
            & (estimates["metric"] == metric)
            & (estimates["model"] == model),
            "estimate",
        ].iloc[0]
        independent_rows.append(
            {
                "model": model,
                "metric": metric,
                "independent_value": value,
                "reported_value": reported,
                "absolute_difference": abs(value - reported),
            }
        )
independent = pd.DataFrame(independent_rows)
require(
    independent.loc[independent["metric"].str.startswith("auc"), "absolute_difference"].max()
    < 1e-10,
    "primary_weighted_auc_independently_recomputed",
    "AUC calculations agree to machine precision",
)
require(
    independent.loc[
        independent["metric"].isin(["uno_c15", "harrell_c"]), "absolute_difference"
    ].max()
    < 1e-3,
    "primary_weighted_c_index_pairwise_crosscheck",
    (
        "Independent simple comparable-pair calculation agrees within 0.001; "
        "small residual reflects survival::concordance handling of tied event/censor times"
    ),
)
independent.to_csv(
    OUT / "nhanes_paired_discrimination_independent_recomputation.csv",
    index=False,
    encoding="utf-8-sig",
)

primary_pair = comparisons.loc[
    (comparisons["domain"] == "age60")
    & (comparisons["population"] == "Overall")
    & (comparisons["comparison_family"] == "PRIMARY_OVERALL_VS_PHENO")
].set_index("metric")
require(
    set(primary_pair.index) == {"uno_c15", "harrell_c", "auc5", "auc10", "auc15"},
    "primary_metric_set_complete",
    f"metrics={sorted(primary_pair.index.tolist())}",
)
require(
    primary_pair.loc["uno_c15", "difference_ci_lower"] <= 0
    <= primary_pair.loc["uno_c15", "difference_ci_upper"],
    "primary_uno_c_difference_not_resolved",
    primary_pair.loc["uno_c15"].to_json(force_ascii=False),
)
require(
    all(primary_pair.loc[m, "difference_ci_upper"] < 0 for m in ("auc5", "auc10", "auc15"))
    and all(primary_pair.loc[m, "paired_p_holm"] < 0.05 for m in ("auc5", "auc10", "auc15")),
    "phenoage_higher_at_all_three_time_auc_horizons",
    primary_pair.loc[["auc5", "auc10", "auc15"]].to_json(force_ascii=False),
)

summary = {
    "status": "PASS" if all(bool(item["passed"]) for item in checks) else "FAIL",
    "checks_passed": sum(bool(item["passed"]) for item in checks),
    "checks_total": len(checks),
    "primary_n": int(primary_audit["n"]),
    "primary_events": int(primary_audit["events"]),
    "bootstrap_replicates": manifest["bootstrap_replicates"],
    "primary_comparison": "OVERALL minus PHENO",
    "primary_results": primary_pair.reset_index().to_dict(orient="records"),
    "interpretation": (
        "The primary 15-year Uno C-index does not resolve a difference. "
        "PhenoAge Acceleration has higher 5-, 10-, and 15-year time-dependent AUCs "
        "after Holm correction. Therefore the evidence does not support universal Overall-score "
        "superiority; it supports horizon-specific PhenoAge discrimination advantages."
    ),
    "checks": checks,
}
(OUT / "nhanes_paired_discrimination_qa_summary.json").write_text(
    json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
)
pd.DataFrame(checks).to_csv(
    OUT / "nhanes_paired_discrimination_qa_checks.csv", index=False, encoding="utf-8-sig"
)

if summary["status"] != "PASS":
    failed = [item for item in checks if not bool(item["passed"])]
    raise SystemExit(f"Paired-discrimination QA failed: {failed}")

print(
    f"Paired-discrimination QA PASS: {summary['checks_passed']}/{summary['checks_total']}; "
    f"primary n/events={summary['primary_n']}/{summary['primary_events']}"
)
