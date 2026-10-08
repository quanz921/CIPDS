from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "supplement_v2"
REPORT = ROOT / "reports" / "supplement_v2"
REPORT.mkdir(parents=True, exist_ok=True)

checks: list[dict[str, str]] = []


def check(name: str, condition: bool, detail: str) -> None:
    checks.append({"check": name, "status": "PASS" if condition else "FAIL", "detail": detail})
    if not condition:
        raise AssertionError(f"{name}: {detail}")


def read_csv(name: str) -> pd.DataFrame:
    path = OUT / name
    check(f"file:{name}", path.exists() and path.stat().st_size > 0, str(path))
    return pd.read_csv(path, encoding="utf-8-sig")


def finite_frame(df: pd.DataFrame, columns: list[str], label: str) -> None:
    values = df[columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    check(f"finite:{label}", np.isfinite(values).all(), f"rows={len(df)}, columns={columns}")


def probability_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    values = df[columns].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    check(f"probability:{label}", ((values >= 0) & (values <= 1)).all(), str(columns))


def interval_contains(df: pd.DataFrame, estimate: str, lower: str, upper: str, label: str) -> None:
    e = pd.to_numeric(df[estimate], errors="coerce")
    lo = pd.to_numeric(df[lower], errors="coerce")
    hi = pd.to_numeric(df[upper], errors="coerce")
    ok = np.isfinite(e).all() and np.isfinite(lo).all() and np.isfinite(hi).all()
    ok = ok and (lo <= e + 1e-12).all() and (e <= hi + 1e-12).all()
    check(f"interval:{label}", bool(ok), f"rows={len(df)}")


required = [
    "selection_participant_flow.csv",
    "selection_weighted_baseline_smd.csv",
    "selection_ipw_audit.csv",
    "selection_ipw_association_sensitivity.csv",
    "ph_assumption_tests.csv",
    "restricted_cubic_spline_tests.csv",
    "restricted_cubic_spline_curves.csv",
    "incremental_prediction_estimates.csv",
    "incremental_prediction_paired_comparisons.csv",
    "cycle_specific_score_effects.csv",
    "cycle_stability_interactions.csv",
    "cycle_specific_demographic_adjusted_effects.csv",
    "cycle_stability_demographic_adjusted_interactions.csv",
    "fair_target_cohort_audit.csv",
    "fair_target_temporal_validation_estimates.csv",
    "fair_target_paired_comparisons.csv",
    "temporal_calibration_cohort_audit.csv",
    "temporal_calibration_metrics.csv",
    "temporal_calibration_groups.csv",
    "competing_risk_finegray_effects.csv",
    "competing_risk_cif_q4_q1_contrasts.csv",
    "competing_risk_temporal_absolute_calibration.csv",
]
for filename in required:
    path = OUT / filename
    check(f"required:{filename}", path.exists() and path.stat().st_size > 0, str(path))

flow = read_csv("selection_participant_flow.csv")
expected_flow = {
    "source_age60": (15048, 6191),
    "overall_available": (14715, 6041),
    "component_fully_adjusted": (9753, 3252),
    "paired_score_only": (11410, 5272),
    "incremental_common_complete": (6987, 2650),
}
observed_flow = {r.node: (int(r.n), int(r.deaths)) for r in flow.itertuples()}
check("selection:locked_counts", observed_flow == expected_flow, repr(observed_flow))

smd = read_csv("selection_weighted_baseline_smd.csv")
finite_frame(smd, ["absolute_smd"], "selection SMD")
check("selection:smd_nonnegative", (smd["absolute_smd"] >= 0).all(), f"max={smd['absolute_smd'].max():.3f}")

ipw_audit = read_csv("selection_ipw_audit.csv")
finite_frame(
    ipw_audit,
    [
        "survey_weighted_selection_prevalence",
        "selection_probability_min",
        "selection_probability_p01",
        "selection_probability_median",
        "selection_probability_p99",
        "selection_probability_max",
        "ipw_trim_lower",
        "ipw_trim_upper",
        "ipw_effective_sample_size",
    ],
    "IPW audit",
)
check(
    "selection:finite_probabilities_in_supported_cycles",
    float(ipw_audit.loc[0, "selection_probability_min"]) > 0
    and float(ipw_audit.loc[0, "selection_probability_max"]) < 1,
    f"range={ipw_audit.loc[0, 'selection_probability_min']:.3f}-{ipw_audit.loc[0, 'selection_probability_max']:.3f}",
)

ipw_assoc = read_csv("selection_ipw_association_sensitivity.csv")
ipw_terms = ipw_assoc[ipw_assoc["term"] != "__MODEL__"].copy()
finite_frame(ipw_terms, ["hazard_ratio", "ci_lower", "ci_upper", "p_value"], "IPW associations")
interval_contains(ipw_terms, "hazard_ratio", "ci_lower", "ci_upper", "IPW associations")

ph = read_csv("ph_assumption_tests.csv")
check("PH:five_scores", len(ph) == 5 and set(ph["score"]) == {"NM", "TB", "TC", "OVERALL", "PHENO"}, repr(ph["score"].tolist()))
finite_frame(ph, ["hazard_ratio_per_weighted_sd", "ci_lower", "ci_upper"], "PH effects")
probability_columns(
    ph,
    [
        "linear_association_p",
        "schoenfeld_score_p",
        "schoenfeld_global_p",
        "score_log_time_interaction_p",
        "schoenfeld_score_p_holm",
        "score_log_time_interaction_p_holm",
    ],
    "PH tests",
)

rcs = read_csv("restricted_cubic_spline_tests.csv")
curves = read_csv("restricted_cubic_spline_curves.csv")
check("RCS:five_scores", len(rcs) == 5 and curves["score"].nunique() == 5, f"tests={len(rcs)}, curves={len(curves)}")
probability_columns(rcs, ["spline_overall_p", "spline_nonlinearity_p", "spline_nonlinearity_p_holm"], "RCS tests")
finite_frame(curves, ["score_value", "hazard_ratio", "ci_lower", "ci_upper"], "RCS curves")
interval_contains(curves, "hazard_ratio", "ci_lower", "ci_upper", "RCS curves")

inc = read_csv("incremental_prediction_estimates.csv")
inc_cmp = read_csv("incremental_prediction_paired_comparisons.csv")
check("incremental:shape", len(inc) == 48 and len(inc_cmp) == 88, f"estimates={len(inc)}, comparisons={len(inc_cmp)}")
check("incremental:common_cohort", set(inc["n"]) == {6987} and set(inc["events"]) == {2650}, "N=6987; deaths=2650")
finite_frame(inc, ["estimate", "standard_error", "ci_lower", "ci_upper", "finite_replicate_fraction"], "incremental estimates")
interval_contains(inc, "estimate", "ci_lower", "ci_upper", "incremental estimates")
check("incremental:replicate_success", (inc["finite_replicate_fraction"] >= 0.98).all(), f"minimum={inc['finite_replicate_fraction'].min():.3f}")
delta_error = np.abs(
    inc_cmp["paired_difference_a_minus_b"] - (inc_cmp["estimate_a"] - inc_cmp["estimate_b"])
).max()
check("incremental:paired_arithmetic", delta_error < 1e-12, f"maximum_error={delta_error:.3g}")
probability_columns(inc_cmp, ["paired_p_value", "paired_p_holm"], "incremental paired tests")

cycle = read_csv("cycle_stability_interactions.csv")
cycle_sens = read_csv("cycle_stability_demographic_adjusted_interactions.csv")
check("cycle:five_scores", len(cycle) == 5 and len(cycle_sens) == 5, f"primary={len(cycle)}, sensitivity={len(cycle_sens)}")
probability_columns(
    cycle,
    ["cycle_interaction_p", "early_late_interaction_p", "cycle_interaction_p_holm", "early_late_interaction_p_holm"],
    "cycle primary",
)
probability_columns(cycle_sens, ["cycle_interaction_p", "cycle_interaction_p_holm"], "cycle sensitivity")

fair_audit = read_csv("fair_target_cohort_audit.csv")
fair = read_csv("fair_target_temporal_validation_estimates.csv")
fair_cmp = read_csv("fair_target_paired_comparisons.csv")
check(
    "fair:cohort",
    list(fair_audit["n"]) == [2548, 3045] and list(fair_audit["deaths"]) == [1400, 1085],
    fair_audit[["role", "n", "deaths"]].to_dict("records").__repr__(),
)
check("fair:shape", len(fair) == 24 and len(fair_cmp) == 36, f"estimates={len(fair)}, comparisons={len(fair_cmp)}")
check("fair:no_15y_extrapolation", not fair["metric"].astype(str).str.contains("15").any(), repr(sorted(fair["metric"].unique())))
interval_contains(fair, "estimate", "ci_lower", "ci_upper", "fair estimates")
check("fair:replicate_success", (fair["finite_replicate_fraction"] >= 0.98).all(), f"minimum={fair['finite_replicate_fraction'].min():.3f}")
fair_delta_error = np.abs(
    fair_cmp["paired_difference_a_minus_b"] - (fair_cmp["estimate_a"] - fair_cmp["estimate_b"])
).max()
check("fair:paired_arithmetic", fair_delta_error < 1e-12, f"maximum_error={fair_delta_error:.3g}")

cal_audit = read_csv("temporal_calibration_cohort_audit.csv")
cal = read_csv("temporal_calibration_metrics.csv")
cal_groups = read_csv("temporal_calibration_groups.csv")
check(
    "calibration:cohort",
    list(cal_audit["n"]) == [2579, 3048] and list(cal_audit["deaths"]) == [1417, 1085],
    cal_audit[["role", "n", "deaths"]].to_dict("records").__repr__(),
)
check("calibration:shape", len(cal) == 36 and len(cal_groups) == 120, f"metrics={len(cal)}, grouped={len(cal_groups)}")
check("calibration:no_15y_extrapolation", not cal["metric"].astype(str).str.contains("15").any(), repr(sorted(cal["metric"].unique())))
check("calibration:replicate_success", (cal["finite_replicate_fraction"] >= 0.98).all(), f"minimum={cal['finite_replicate_fraction'].min():.3f}")
check(
    "calibration:no_stale_error_file",
    not (OUT / "temporal_calibration_bootstrap_errors.csv").exists(),
    "stale failed-run audit removed",
)

fg = read_csv("competing_risk_finegray_effects.csv")
cif_delta = read_csv("competing_risk_cif_q4_q1_contrasts.csv")
cif_cal = read_csv("competing_risk_temporal_absolute_calibration.csv")
check("competing:FineGray_shape", len(fg) == 10, f"rows={len(fg)}")
interval_contains(fg, "subdistribution_hazard_ratio", "ci_lower", "ci_upper", "Fine-Gray effects")
probability_columns(fg, ["p_value", "p_holm"], "Fine-Gray tests")
check("competing:CIF_delta_shape", len(cif_delta) == 30, f"rows={len(cif_delta)}")
interval_contains(
    cif_delta,
    "q4_minus_q1_absolute_risk_difference",
    "ci_lower",
    "ci_upper",
    "CIF risk differences",
)
probability_columns(cif_delta, ["p_value", "p_holm"], "CIF contrasts")
check(
    "competing:temporal_calibration_shape",
    len(cif_cal) == 40 and set(cif_cal["horizon_years"]) == {5, 10},
    f"rows={len(cif_cal)}, horizons={sorted(cif_cal['horizon_years'].unique())}",
)

for manifest_name in [
    "selection_manifest.json",
    "ph_nonlinearity_manifest.json",
    "incremental_prediction_manifest.json",
    "cycle_stability_manifest.json",
    "fair_target_manifest.json",
    "temporal_calibration_manifest.json",
    "competing_risk_manifest.json",
]:
    manifest_path = OUT / manifest_name
    check(f"manifest:{manifest_name}", manifest_path.exists(), str(manifest_path))
    with manifest_path.open("r", encoding="utf-8-sig") as handle:
        json.load(handle)

hash_rows = []
for path in sorted(OUT.iterdir(), key=lambda p: p.name.lower()):
    if not path.is_file() or path.name == "supplement_v2_sha256.csv":
        continue
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    hash_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()})

pd.DataFrame(hash_rows).to_csv(OUT / "supplement_v2_sha256.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(checks).to_csv(REPORT / "supplement_v2_validation_checks.csv", index=False, encoding="utf-8-sig")
qa = {
    "status": "PASS",
    "checks_passed": len(checks),
    "checks_failed": 0,
    "scope": "Seven prespecified NHANES supplementary analysis modules",
    "locked_counts": expected_flow,
    "warnings": [
        "Fine-Gray regression is a weighted PSU-robust sensitivity analysis, not a full stratified multistage variance estimator.",
        "The matched-target experiment and temporal calibration are internal temporal validations within NHANES, not independent external cohorts.",
        "Fifteen-year temporal validation/calibration was not estimated because 2007-2010 follow-up ends at 13.17 years.",
    ],
}
with (REPORT / "supplement_v2_validation_summary.json").open("w", encoding="utf-8") as handle:
    json.dump(qa, handle, ensure_ascii=False, indent=2)

print(f"PASS: {len(checks)} checks")
print(f"Validation table: {REPORT / 'supplement_v2_validation_checks.csv'}")
print(f"SHA256 registry: {OUT / 'supplement_v2_sha256.csv'}")

