from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(
    os.environ.get(
        "CIPDS_PACKAGE_DIR_PY",
        "reference",
    )
)
OUT = PACKAGE / "outputs"
SCENARIOS = ("a_primary", "ab_sensitivity")
EXPECTED_SCENARIO_LABELS = {
    "a_primary": "A_PRIMARY",
    "ab_sensitivity": "AB_SENSITIVITY",
}
COMPONENT_TERMS = {"NM_z", "TB_z", "TC_z"}


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, evidence: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "evidence": str(evidence)})

    scenario_summaries: dict[str, object] = {}
    score_tables: dict[str, pd.DataFrame] = {}

    for slug in SCENARIOS:
        prefix = f"nested_{slug}_nhanes_"
        label = EXPECTED_SCENARIO_LABELS[slug]
        scores = pd.read_csv(OUT / f"{prefix}component_scores.csv")
        registry = pd.read_csv(OUT / f"{prefix}score_registry.csv")
        scaling = pd.read_csv(OUT / f"{prefix}survey_scaling.csv")
        separate = pd.read_csv(OUT / f"{prefix}survey_component_separate.csv")
        joint = pd.read_csv(OUT / f"{prefix}survey_components_joint.csv")
        wald = pd.read_csv(OUT / f"{prefix}survey_joint_wald.csv")
        design = pd.read_csv(OUT / f"{prefix}survey_design_audit.csv")
        final_features = pd.read_csv(OUT / f"nested_{slug}_final_features.csv")
        nested_qa = json.loads(
            (OUT / f"nested_{slug}_qa_summary.json").read_text(encoding="utf-8")
        )
        score_tables[slug] = scores

        check(
            f"{label}: NHANES participant grain",
            len(scores) == 44772 and scores["SEQN"].nunique() == 44772,
            f"rows={len(scores)}, unique_SEQN={scores['SEQN'].nunique()}",
        )
        check(
            f"{label}: nine cycles represented",
            scores["CYCLE"].nunique() == 9,
            sorted(scores["CYCLE"].astype(str).unique()),
        )
        probability_columns = ["NM_component", "TB_component", "TC_component"]
        probability_ok = (
            scores[probability_columns].notna().all().all()
            and scores[probability_columns].apply(lambda x: x.between(0, 1).all()).all()
        )
        check(
            f"{label}: component score validity",
            probability_ok,
            {
                c: [float(scores[c].min()), float(scores[c].max()), int(scores[c].isna().sum())]
                for c in probability_columns
            },
        )

        feature_counts = final_features.groupby("outcome").size().to_dict()
        registry_counts = registry.set_index("outcome")["n_features"].astype(int).to_dict()
        registry_ok = (
            len(registry) == 3
            and registry["outcome"].nunique() == 3
            and registry["score_missing_n"].eq(0).all()
            and feature_counts == registry_counts
        )
        check(
            f"{label}: score registry and selected features agree",
            registry_ok,
            {"selected": feature_counts, "registry": registry_counts},
        )
        check(
            f"{label}: survey result table completeness",
            len(scaling) == 4 and len(separate) == 24 and len(joint) == 18 and len(wald) == 6,
            {
                "scaling": len(scaling),
                "separate": len(separate),
                "joint": len(joint),
                "wald": len(wald),
            },
        )

        numeric_columns = [
            "hazard_ratio_per_survey_weighted_sd",
            "ci_lower",
            "ci_upper",
            "p_value",
        ]
        result_numbers_ok = True
        for table in (separate, joint):
            vals = table[numeric_columns].to_numpy(dtype=float)
            result_numbers_ok &= bool(np.isfinite(vals).all())
            result_numbers_ok &= bool((table["hazard_ratio_per_survey_weighted_sd"] > 0).all())
            result_numbers_ok &= bool((table["ci_lower"] > 0).all())
            result_numbers_ok &= bool((table["ci_upper"] >= table["ci_lower"]).all())
            result_numbers_ok &= bool(table["p_value"].between(0, 1).all())
        check(f"{label}: finite coherent Cox estimates", result_numbers_ok, "separate and joint tables")

        primary_sep = separate[
            (separate["adjustment"] == "fully_adjusted")
            & (separate["term"].isin(COMPONENT_TERMS))
        ]
        primary_joint = joint[joint["adjustment"] == "fully_adjusted"]
        n_event_pairs_sep = {
            pop: sorted(set(zip(d["n"].astype(int), d["events"].astype(int))))
            for pop, d in primary_sep.groupby("population")
        }
        n_event_pairs_joint = {
            pop: sorted(set(zip(d["n"].astype(int), d["events"].astype(int))))
            for pop, d in primary_joint.groupby("population")
        }
        check(
            f"{label}: common complete-case cohort for separate and joint components",
            n_event_pairs_sep == n_event_pairs_joint
            and all(len(v) == 1 for v in n_event_pairs_sep.values()),
            {"separate": n_event_pairs_sep, "joint": n_event_pairs_joint},
        )

        audit = dict(zip(design["item"].astype(str), design["value"].astype(str)))
        design_ok = (
            audit.get("analysis_rows") == "44772"
            and audit.get("cycles") == "9"
            and audit.get("weight_variable") == "WTMEC2YR"
            and audit.get("pooled_weight_formula") == "WTMEC2YR/9"
            and audit.get("lonely_psu_handling") == "adjust"
            and audit.get("weighted_composite_generated") == "FALSE"
        )
        check(f"{label}: complex-survey design lock", design_ok, audit)
        check(
            f"{label}: upstream 5x5 nested QA passed",
            nested_qa.get("all_checks_passed") is True
            and nested_qa.get("checks_passed") == nested_qa.get("checks_total") == 20,
            {
                "passed": nested_qa.get("all_checks_passed"),
                "checks": f"{nested_qa.get('checks_passed')}/{nested_qa.get('checks_total')}",
            },
        )

        fully_adjusted = primary_sep[
            [
                "population",
                "term",
                "n",
                "events",
                "hazard_ratio_per_survey_weighted_sd",
                "ci_lower",
                "ci_upper",
                "p_value",
            ]
        ].to_dict(orient="records")
        scenario_summaries[label] = {
            "score_feature_counts": registry_counts,
            "fully_adjusted_separate_component_results": fully_adjusted,
            "joint_wald_fully_adjusted": wald[
                wald["adjustment"] == "fully_adjusted"
            ][["population", "n", "events", "joint_wald_p"]].to_dict(orient="records"),
        }

    same_seqn_order = score_tables["a_primary"]["SEQN"].equals(
        score_tables["ab_sensitivity"]["SEQN"]
    )
    max_abs_differences = {
        c: float(
            np.max(
                np.abs(
                    score_tables["a_primary"][c].to_numpy()
                    - score_tables["ab_sensitivity"][c].to_numpy()
                )
            )
        )
        for c in ("NM_component", "TB_component", "TC_component")
    }
    check(
        "A and A+B scored on identical NHANES participants",
        same_seqn_order,
        f"same_order={same_seqn_order}",
    )
    check(
        "A+B sensitivity scores are not accidental copies of A-primary",
        all(v > 0 for v in max_abs_differences.values()),
        max_abs_differences,
    )

    comparison = pd.read_csv(OUT / "nested_a_vs_ab_paired_test_comparison.csv")
    missingness = pd.read_csv(OUT / "nested_b_tier_missingness_only_diagnostic.csv")
    check(
        "paired calendar-test comparison completeness",
        len(comparison) == 3
        and comparison["n"].eq(11125).all()
        and comparison["auc_difference_ci_orientation"].eq(
            "AB_SENSITIVITY minus A_PRIMARY"
        ).all(),
        comparison[["outcome", "auc_difference_ab_minus_a", "paired_delong_p"]].to_dict(
            orient="records"
        ),
    )
    test_missingness = missingness[missingness["split"] == "calendar_test"]
    check(
        "B-tier missingness-only diagnostic completeness",
        len(test_missingness) == 3
        and test_missingness["n"].eq(11125).all()
        and test_missingness["missingness_only_auc"].between(0.5, 1).all(),
        test_missingness[["outcome", "missingness_only_auc"]].to_dict(orient="records"),
    )

    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(
        OUT / "nested_nhanes_qa_checks.csv", index=False, encoding="utf-8-sig"
    )
    summary = {
        "all_checks_passed": bool(checks_df["passed"].all()),
        "checks_total": int(len(checks_df)),
        "checks_passed": int(checks_df["passed"].sum()),
        "checks_failed": int((~checks_df["passed"]).sum()),
        "scenario_summaries": scenario_summaries,
        "max_abs_score_difference_ab_vs_a": max_abs_differences,
        "paired_calendar_test_comparison": comparison.to_dict(orient="records"),
        "b_tier_missingness_only_calendar_test": test_missingness.to_dict(orient="records"),
        "interpretation_lock": (
            "A_PRIMARY is primary. AB_SENSITIVITY is sensitivity only because B-tier "
            "missingness indicators alone discriminate all three outcomes; do not attribute "
            "A+B gains solely to additional physiology."
        ),
        "composite_lock": "No AUC-weighted composite was generated.",
    }
    (OUT / "nested_nhanes_qa_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not summary["all_checks_passed"]:
        failed = checks_df.loc[~checks_df["passed"], ["check", "evidence"]]
        raise SystemExit("NHANES nested QA failed:\n" + failed.to_string(index=False))


if __name__ == "__main__":
    main()
