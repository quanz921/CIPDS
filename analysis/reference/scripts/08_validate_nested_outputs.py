from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[1]
OUT = PACKAGE / "outputs"
SCENARIO = os.environ.get("CIPDS_QA_SCENARIO", "a_primary").lower()
if SCENARIO not in {"a_primary", "ab_sensitivity"}:
    raise RuntimeError("CIPDS_QA_SCENARIO must be a_primary or ab_sensitivity")
PREFIX = f"nested_{SCENARIO}_"
OUTCOMES = ("Outcome_NutriMetab", "Outcome_TumorBurden", "Outcome_TreatComp")


def auc_rank(y: np.ndarray, p: np.ndarray) -> float:
    y = np.asarray(y, dtype=int)
    p = np.asarray(p, dtype=float)
    n1 = int((y == 1).sum())
    n0 = int((y == 0).sum())
    ranks = pd.Series(p).rank(method="average").to_numpy()
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def average_precision(y: np.ndarray, p: np.ndarray) -> float:
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(y,p))


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, evidence: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "evidence": str(evidence)})

    audit = pd.read_csv(OUT / "lab_129_crosscohort_eligibility.csv", low_memory=False)
    base = pd.read_csv(
        OUT / "hospital_patient_model_ready_129.csv",
        usecols=["patient_key", "split_calendar_entry", *OUTCOMES],
        low_memory=False,
    )
    outer = pd.read_csv(OUT / f"{PREFIX}outer_performance.csv")
    inner = pd.read_csv(OUT / f"{PREFIX}inner_tuning.csv")
    final_tuning = pd.read_csv(OUT / f"{PREFIX}final_tuning.csv")
    votes = pd.read_csv(OUT / f"{PREFIX}feature_votes.csv")
    final_features = pd.read_csv(OUT / f"{PREFIX}final_features.csv")
    metrics = pd.read_csv(OUT / f"{PREFIX}calendar_metrics.csv")
    predictions = pd.read_csv(OUT / f"{PREFIX}calendar_test_predictions.csv")
    calibration = pd.read_csv(OUT / f"{PREFIX}calibration_and_model_registry.csv")
    manifest = json.loads((OUT / f"{PREFIX}run_manifest.json").read_text(encoding="utf-8"))

    check("patient grain", len(base) == 75248 and base["patient_key"].nunique() == 75248,
          f"rows={len(base)}, unique={base['patient_key'].nunique()}")
    split_sets = {
        s: set(base.loc[base["split_calendar_entry"] == s, "patient_key"])
        for s in ("train", "validation", "test")
    }
    overlap = (
        len(split_sets["train"] & split_sets["validation"])
        + len(split_sets["train"] & split_sets["test"])
        + len(split_sets["validation"] & split_sets["test"])
    )
    check("calendar split patient disjointness", overlap == 0, f"overlap={overlap}")

    check("129-variable registry", len(audit) == 129 and audit["final_variable_id"].nunique() == 129,
          f"rows={len(audit)}")
    if SCENARIO == "a_primary":
        allowed = set(audit.loc[audit["primary_crosscohort_eligible"] == True, "final_variable_id"])
        expected_allowed = 32
        allowed_label = "A-primary"
    else:
        allowed = set(audit.loc[audit["full_period_sensitivity_eligible"] == True, "final_variable_id"])
        expected_allowed = 36
        allowed_label = "A+B sensitivity"
    check(f"{allowed_label} candidate count", len(allowed) == expected_allowed,
          f"count={len(allowed)}")

    outer_counts = outer.groupby("outcome")["outer_fold"].nunique().to_dict()
    check("five outer folds per outcome", all(outer_counts.get(o) == 5 for o in OUTCOMES), outer_counts)
    check("outer performance rows", len(outer) == 15, f"rows={len(outer)}")
    check("outer AUC range", outer["outer_auc"].between(0, 1).all(),
          f"min={outer['outer_auc'].min()}, max={outer['outer_auc'].max()}")

    expected_inner = len(OUTCOMES) * 5 * 12 * 5
    expected_final_inner = len(OUTCOMES) * 12 * 5
    check("outer-loop inner tuning completeness", len(inner) == expected_inner,
          f"rows={len(inner)}, expected={expected_inner}")
    check("full-training inner tuning completeness", len(final_tuning) == expected_final_inner,
          f"rows={len(final_tuning)}, expected={expected_final_inner}")

    selected_vote_ok = final_features["votes"].ge(3).all()
    check("final Boruta stability rule", selected_vote_ok,
          final_features.groupby("outcome").size().to_dict())
    selected_outside_allowed = set(final_features["feature"]) - allowed
    check(f"final features restricted to {allowed_label}", not selected_outside_allowed,
          sorted(selected_outside_allowed))

    hard_excluded = {
        "Outcome_NutriMetab": {
            "lab_alb", "lab_a_ratio_g", "lab_glob", "lab_tp", "lab_hgb", "lab_hct",
            "lab_rbc", "lab_rdw", "lab_mcv", "lab_mch", "lab_mchc",
        },
        "Outcome_TumorBurden": set(),
        "Outcome_TreatComp": {
            "lab_baso", "lab_baso_pct", "lab_eos", "lab_eos_pct", "lab_hct", "lab_hgb",
            "lab_lym", "lab_lym_pct", "lab_mch", "lab_mchc", "lab_mcv", "lab_mono",
            "lab_mono_pct", "lab_neu", "lab_neu_pct", "lab_plt", "lab_rbc", "lab_rdw",
            "lab_wbc", "lab_mpv",
        },
    }
    violations = {}
    for outcome, excluded in hard_excluded.items():
        selected = set(final_features.loc[final_features["outcome"] == outcome, "feature"])
        bad = sorted(selected & excluded)
        if bad:
            violations[outcome] = bad
    check("outcome-overlap exclusions", not violations, violations)

    check("manifest 5x5", manifest["nested_cv"] == {"outer_folds": 5, "inner_folds": 5},
          manifest["nested_cv"])
    check("manifest 23 threads", manifest["catboost"]["thread_count"] == 23,
          manifest["catboost"]["thread_count"])
    check("analysis status caveat", manifest["analysis_status"] == "AUDIT_ONLY_LEGACY_UNDATED_OUTCOMES",
          manifest["analysis_status"])

    test_metrics = metrics.loc[metrics["split"] == "calendar_test"].set_index("outcome")
    prediction_counts = predictions.groupby("outcome").size().to_dict()
    prediction_unique = predictions.groupby("outcome")["patient_key"].nunique().to_dict()
    check("test prediction row counts", all(prediction_counts.get(o) == 11125 for o in OUTCOMES),
          prediction_counts)
    check("test patient uniqueness", all(prediction_unique.get(o) == 11125 for o in OUTCOMES),
          prediction_unique)
    check("prediction probability validity",
          predictions["calibrated_probability"].notna().all()
          and predictions["calibrated_probability"].between(0, 1).all(),
          f"range={predictions['calibrated_probability'].min()}..{predictions['calibrated_probability'].max()}")

    recomputed_rows = []
    all_metric_match = True
    for outcome in OUTCOMES:
        d = predictions.loc[predictions["outcome"] == outcome]
        y = d["observed"].to_numpy(dtype=int)
        p = d["calibrated_probability"].to_numpy(dtype=float)
        auc = auc_rank(y, p)
        ap = average_precision(y, p)
        brier = float(np.mean((p - y) ** 2))
        reported = test_metrics.loc[outcome]
        matched = (
            abs(auc - reported["auc"]) < 1e-10
            and abs(ap - reported["average_precision"]) < 1e-10
            and abs(brier - reported["brier"]) < 1e-10
        )
        all_metric_match &= matched
        recomputed_rows.append(
            {
                "outcome": outcome,
                "auc_recomputed": auc,
                "auc_reported": reported["auc"],
                "average_precision_recomputed": ap,
                "average_precision_reported": reported["average_precision"],
                "brier_recomputed": brier,
                "brier_reported": reported["brier"],
                "matched": matched,
            }
        )
    check("independent metric recomputation", all_metric_match, recomputed_rows)

    model_status = {}
    for _, row in calibration.iterrows():
        p = OUT / row["model_file"]
        model_status[row["outcome"]] = {"exists": p.exists(), "bytes": p.stat().st_size if p.exists() else 0}
    check("three nonempty CatBoost model files",
          len(model_status) == 3 and all(v["exists"] and v["bytes"] > 0 for v in model_status.values()),
          model_status)

    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(OUT / f"{PREFIX}qa_checks.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(recomputed_rows).to_csv(
        OUT / f"{PREFIX}metric_recomputation.csv", index=False, encoding="utf-8-sig"
    )
    outer_summary = outer.groupby("outcome").agg(
        outer_auc_mean=("outer_auc", "mean"),
        outer_auc_sd=("outer_auc", "std"),
        outer_auc_min=("outer_auc", "min"),
        outer_auc_max=("outer_auc", "max"),
        selected_features_mean=("n_selected", "mean"),
    ).reset_index()
    outer_summary.to_csv(
        OUT / f"{PREFIX}outer_performance_summary.csv", index=False, encoding="utf-8-sig"
    )

    summary = {
        "all_checks_passed": bool(checks_df["passed"].all()),
        "checks_total": int(len(checks_df)),
        "checks_passed": int(checks_df["passed"].sum()),
        "checks_failed": int((~checks_df["passed"]).sum()),
        "outer_performance": outer_summary.to_dict(orient="records"),
        "calendar_test_metrics": test_metrics.reset_index()[
            ["outcome", "n", "positives", "auc", "auc_ci_lower", "auc_ci_upper",
             "average_precision", "brier", "calibration_intercept", "calibration_slope"]
        ].to_dict(orient="records"),
        "selected_features": {
            o: final_features.loc[final_features["outcome"] == o, "feature"].tolist()
            for o in OUTCOMES
        },
        "required_caveat": (
            "Hospital diagnoses preceded the sampled laboratory measurements, as confirmed by the data owner. "
            "Hospital models classify established deterioration states; NHANES evaluates subsequent mortality."
        ),
    }
    (OUT / f"{PREFIX}qa_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not summary["all_checks_passed"]:
        raise SystemExit("Nested output QA failed")


if __name__ == "__main__":
    main()
