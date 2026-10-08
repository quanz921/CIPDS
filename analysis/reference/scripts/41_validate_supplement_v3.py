from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "supplement_v3"
OUT2 = ROOT / "outputs" / "supplement_v2"
REPORT = ROOT / "reports" / "supplement_v3"
REPORT.mkdir(parents=True, exist_ok=True)

checks: list[dict[str, str]] = []


def check(name: str, condition: bool, evidence: str) -> None:
    if not condition:
        raise AssertionError(f"{name}: {evidence}")
    checks.append({"check": name, "status": "PASS", "evidence": evidence})


def read(name: str) -> pd.DataFrame:
    path = OUT / name
    check(f"file:{name}", path.exists() and path.stat().st_size > 0, str(path))
    return pd.read_csv(path, encoding="utf-8-sig")


def finite(frame: pd.DataFrame, columns: list[str], label: str) -> None:
    arr = frame[columns].to_numpy(dtype=float)
    check(f"finite:{label}", np.isfinite(arr).all(), f"shape={arr.shape}")


def intervals(frame: pd.DataFrame, est: str, lo: str, hi: str, label: str) -> None:
    ok = (frame[lo] <= frame[est]).all() and (frame[est] <= frame[hi]).all()
    check(f"interval:{label}", bool(ok), f"rows={len(frame)}")


landmark = read("early_death_landmark_sensitivity.csv")
check("landmark:shape", len(landmark) == 15, f"rows={len(landmark)}")
check("landmark:scores", set(landmark.score) == {"NM", "TB", "TC", "OVERALL", "PHENO"}, str(sorted(landmark.score.unique())))
check("landmark:cohorts", landmark.groupby("landmark_years").n.first().to_dict() == {0: 6987, 1: 6821, 2: 6626}, str(landmark.groupby("landmark_years").n.first().to_dict()))
check("landmark:events", landmark.groupby("landmark_years").events_after_landmark.first().to_dict() == {0: 2650, 1: 2484, 2: 2289}, str(landmark.groupby("landmark_years").events_after_landmark.first().to_dict()))
finite(landmark, ["hazard_ratio", "ci_lower", "ci_upper", "p_value", "p_holm_within_landmark"], "landmark")
intervals(landmark, "hazard_ratio", "ci_lower", "ci_upper", "landmark")
check("landmark:overall_persistent", (landmark.loc[landmark.score == "OVERALL", "ci_lower"] > 1).all(), landmark.loc[landmark.score == "OVERALL", ["landmark_years", "hazard_ratio", "ci_lower", "ci_upper"]].to_dict("records").__repr__())

age = read("age_time_scale_sensitivity.csv")
check("age:shape", len(age) == 10 and set(age.time_scale) == {"FOLLOW_UP_TIME", "ATTAINED_AGE"}, f"rows={len(age)}")
check("age:locked_cohort", set(age.n) == {6987} and set(age.events) == {2650}, "N=6987; deaths=2650")
finite(age, ["hazard_ratio", "ci_lower", "ci_upper", "p_value", "p_holm_within_time_scale"], "age time scale")
intervals(age, "hazard_ratio", "ci_lower", "ci_upper", "age time scale")
check("age:overall_persistent", (age.loc[age.score == "OVERALL", "ci_lower"] > 1).all(), age.loc[age.score == "OVERALL", ["time_scale", "hazard_ratio", "ci_lower", "ci_upper"]].to_dict("records").__repr__())

piece = read("time_varying_piecewise_hr.csv")
hetero = read("time_varying_interval_heterogeneity.csv")
curve = read("time_varying_continuous_hr_curves.csv")
trend = read("time_varying_logtime_tests.csv")
check("time:piecewise_shape", len(piece) == 15 and piece.interval.nunique() == 3, f"rows={len(piece)}")
check("time:event_partition", piece.groupby("score").events_in_interval.sum().eq(2650).all(), piece.groupby("score").events_in_interval.sum().to_dict().__repr__())
finite(piece, ["hazard_ratio", "ci_lower", "ci_upper", "p_value", "p_holm_within_interval"], "piecewise HR")
intervals(piece, "hazard_ratio", "ci_lower", "ci_upper", "piecewise HR")
check("time:heterogeneity_shape", len(hetero) == 5, f"rows={len(hetero)}")
finite(hetero, ["wald_chisq", "interval_heterogeneity_p", "interval_heterogeneity_p_holm"], "time heterogeneity")
check("time:curve_shape", len(curve) == 295 and curve.groupby("score").size().eq(59).all(), f"rows={len(curve)}")
finite(curve, ["follow_up_years", "hazard_ratio_per_weighted_sd", "ci_lower", "ci_upper"], "continuous HR")
intervals(curve, "hazard_ratio_per_weighted_sd", "ci_lower", "ci_upper", "continuous HR")
check("time:trend_shape", len(trend) == 5, f"rows={len(trend)}")

mapping = read("residual_mapping_quality.csv")
identity = read("residual_fidelity_identity_audit.csv")
est = read("residual_fidelity_validation_estimates.csv")
cmp = read("residual_fidelity_paired_comparisons.csv")
corr = read("residual_dimension_lab_correlations.csv")
check("residual:mapping_roles", set(mapping.role) == {"TRAIN_OOF", "TRAIN_APPARENT", "TEMPORAL_VALIDATION"}, str(mapping.role.tolist()))
finite(mapping, ["weighted_r_squared", "weighted_rmse", "weighted_correlation", "observed_weighted_sd", "residual_weighted_sd"], "mapping quality")
check("residual:temporal_r2_valid", mapping.loc[mapping.role == "TEMPORAL_VALIDATION", "weighted_r_squared"].between(0, 1).all(), mapping.to_dict("records").__repr__())
check("residual:identity", len(identity) == 4 and identity.passed.astype(bool).all(), identity.to_dict("records").__repr__())
check("residual:identity_precision", identity.maximum_absolute_error.max() < 1e-10, f"max={identity.maximum_absolute_error.max():.3g}")
check("residual:estimate_shape", len(est) == 30 and est.model.nunique() == 5, f"rows={len(est)}")
check("residual:paired_cohort", set(est.n_validation) == {3045} and set(est.validation_deaths) == {1085}, "N=3045; deaths=1085")
finite(est, ["estimate", "standard_error", "ci_lower", "ci_upper", "finite_replicate_fraction"], "residual estimates")
intervals(est, "estimate", "ci_lower", "ci_upper", "residual estimates")
check("residual:bootstrap_success", (est.finite_replicate_fraction >= 0.98).all(), f"minimum={est.finite_replicate_fraction.min():.3f}")
check("residual:comparison_shape", len(cmp) == 36, f"rows={len(cmp)}")
finite(cmp, ["estimate_a", "estimate_b", "paired_difference_a_minus_b", "difference_ci_lower", "difference_ci_upper", "paired_p_value", "paired_p_holm"], "residual comparisons")
paired_error = np.abs(cmp.paired_difference_a_minus_b - (cmp.estimate_a - cmp.estimate_b)).max()
check("residual:paired_arithmetic", paired_error < 1e-12, f"maximum_error={paired_error:.3g}")

pivot = est.pivot(index="metric", columns="model", values="estimate")
check("residual:hybrid_exact_metrics", np.abs(pivot.HYBRID4 - pivot.LAB26).max() < 1e-12, f"max={np.abs(pivot.HYBRID4 - pivot.LAB26).max():.3g}")
for metric in ["uno_c10", "auc5", "auc10"]:
    check(f"residual:recovered_{metric}", pivot.loc[metric, "LAB26"] > pivot.loc[metric, "PROJECTED3"], f"LAB26={pivot.loc[metric, 'LAB26']:.6f}; projected={pivot.loc[metric, 'PROJECTED3']:.6f}")
for metric in ["brier5", "brier10", "ibs10"]:
    check(f"residual:recovered_{metric}", pivot.loc[metric, "LAB26"] < pivot.loc[metric, "PROJECTED3"], f"LAB26={pivot.loc[metric, 'LAB26']:.6f}; projected={pivot.loc[metric, 'PROJECTED3']:.6f}")
check("residual:lab_correlations", len(corr) == 26 and corr.descriptive_rank.nunique() == 26, f"rows={len(corr)}")
finite(corr, ["glmnet_coefficient_in_lab26_model", "weighted_correlation_with_residual", "absolute_weighted_correlation"], "residual lab correlations")

# Reconcile the reconstructed legacy models against the already validated v2
# matched-target outputs.
old = pd.read_csv(OUT2 / "fair_target_temporal_validation_estimates.csv", encoding="utf-8-sig")
merged = est[est.model.isin(["LAB26", "COMPONENT3", "PHENO9"])].merge(
    old[["metric", "model", "estimate"]], on=["metric", "model"], suffixes=("_v3", "_v2"), how="inner"
)
reconcile_error = np.abs(merged.estimate_v3 - merged.estimate_v2).max()
check("reconcile:v2_fair_target", len(merged) == 18 and reconcile_error < 1e-12, f"rows={len(merged)}; max_error={reconcile_error:.3g}")

for manifest_name in ["time_reverse_age_manifest.json", "residual_fidelity_manifest.json"]:
    path = OUT / manifest_name
    check(f"manifest:{manifest_name}", path.exists(), str(path))
    with path.open("r", encoding="utf-8-sig") as handle:
        json.load(handle)

hash_rows = []
for path in sorted(OUT.iterdir(), key=lambda p: p.name.lower()):
    if not path.is_file() or path.name == "supplement_v3_sha256.csv":
        continue
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    hash_rows.append({"file": path.name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()})

pd.DataFrame(hash_rows).to_csv(OUT / "supplement_v3_sha256.csv", index=False, encoding="utf-8-sig")
pd.DataFrame(checks).to_csv(REPORT / "supplement_v3_validation_checks.csv", index=False, encoding="utf-8-sig")
summary = {
    "status": "PASS",
    "checks_passed": len(checks),
    "checks_failed": 0,
    "scope": [
        "time-varying effects",
        "one- and two-year landmark analyses",
        "attained-age time scale",
        "target-specific residual-fidelity decomposition",
    ],
    "warnings": [
        "The attained-age analysis inherits NHANES public-release age top-coding.",
        "Piecewise and smooth time-varying HRs use sampling-weighted PSU-robust Cox models; landmark and age-time-scale analyses use the full NHANES complex survey design.",
        "The residual dimension is target-specific and is not lossless compression of the 26 raw laboratory values.",
        "The residual dimension is not a fourth biological domain or aging archetype.",
        "The residual-fidelity analysis is internal temporal validation within NHANES and does not modify the hospital-trained frozen scores.",
    ],
}
with (REPORT / "supplement_v3_validation_summary.json").open("w", encoding="utf-8") as handle:
    json.dump(summary, handle, ensure_ascii=False, indent=2)

print(f"PASS: {len(checks)} checks")
print(REPORT / "supplement_v3_validation_checks.csv")

