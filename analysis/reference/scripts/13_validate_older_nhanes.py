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
SOURCE_CSV = Path("external_inputs/NHANES_cleaned.csv")
SCRIPT_R = PACKAGE / "scripts" / "12_nhanes_older_adult_primary.R"


def main() -> None:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, evidence: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "evidence": str(evidence)})

    source_all = pd.read_csv(
        SOURCE_CSV,
        usecols=[
            "SEQN",
            "Age",
            "Cancer_Diagnosed",
            "Dead",
            "UCOD_LEADING",
            "CYCLE",
            "Follow_Up_Years",
        ],
        low_memory=False,
    )
    expected_cycles = {
        "1999-2000",
        "2001-2002",
        "2003-2004",
        "2005-2006",
        "2007-2008",
        "2009-2010",
        "2011-2012",
        "2013-2014",
        "2015-2016",
    }
    source = source_all[source_all["CYCLE"].isin(expected_cycles)].copy()
    check(
        "independent source contains one later non-analysis cycle",
        len(source_all) == 49774
        and source_all["CYCLE"].nunique() == 10
        and set(source_all["CYCLE"].unique()) - expected_cycles == {"2017-2018"},
        {
            "source_rows": len(source_all),
            "source_cycles": sorted(source_all["CYCLE"].unique()),
            "locked_analysis_cycles": sorted(expected_cycles),
        },
    )
    later_cycle = source_all[source_all["CYCLE"] == "2017-2018"]
    check(
        "2017-2018 excluded because survival time is unavailable",
        len(later_cycle) == 5002 and later_cycle["Follow_Up_Years"].notna().sum() == 0,
        {
            "rows": len(later_cycle),
            "nonmissing_follow_up_years": int(later_cycle["Follow_Up_Years"].notna().sum()),
        },
    )
    check(
        "independent source grain",
        len(source) == 44772 and source["SEQN"].nunique() == 44772,
        f"rows={len(source)}, unique_SEQN={source['SEQN'].nunique()}",
    )
    check("nine NHANES cycles", source["CYCLE"].nunique() == 9, source["CYCLE"].nunique())

    expected_counts: dict[tuple[str, str], dict[str, int]] = {}
    for domain, threshold in (("age60", 60), ("age65", 65)):
        domain_mask = source["Age"].notna() & source["Age"].ge(threshold)
        populations = {
            "Overall": domain_mask,
            "Cancer": domain_mask & source["Cancer_Diagnosed"].eq(1),
            "Noncancer": domain_mask & source["Cancer_Diagnosed"].eq(0),
        }
        for population, mask in populations.items():
            d = source.loc[mask]
            expected_counts[(domain, population)] = {
                "n": len(d),
                "Death_AllCause": int(d["Dead"].eq(1).sum()),
                "Death_Cancer": int((d["Dead"].eq(1) & d["UCOD_LEADING"].eq(2)).sum()),
                "Death_CVD": int(
                    (d["Dead"].eq(1) & d["UCOD_LEADING"].isin([1, 5]) & d["CYCLE"].ne("2015-2016")).sum()
                ),
                "Death_CLRD": int((d["Dead"].eq(1) & d["UCOD_LEADING"].eq(3)).sum()),
                "Death_Alzheimer": int(
                    (d["Dead"].eq(1) & d["UCOD_LEADING"].eq(6)).sum()
                ),
                "Death_Diabetes_UCOD": int(
                    (d["Dead"].eq(1) & d["UCOD_LEADING"].eq(7)).sum()
                ),
                "Death_Kidney": int(
                    (d["Dead"].eq(1) & d["UCOD_LEADING"].eq(9)).sum()
                ),
            }

    audit = pd.read_csv(OUT / "older_nhanes_cohort_event_audit.csv")
    audit_match = True
    mismatches: list[dict[str, object]] = []
    for (domain, population), expected in expected_counts.items():
        row = audit[(audit["domain"] == domain) & (audit["population"] == population)]
        if len(row) != 1:
            audit_match = False
            mismatches.append({"domain": domain, "population": population, "rows": len(row)})
            continue
        actual = row.iloc[0]
        for field, value in expected.items():
            if int(actual[field]) != value:
                audit_match = False
                mismatches.append(
                    {
                        "domain": domain,
                        "population": population,
                        "field": field,
                        "expected": value,
                        "actual": int(actual[field]),
                    }
                )
    check("independent cohort and cause-count reconciliation", audit_match, mismatches or expected_counts)
    check(
        "age >=60 headline cohort",
        expected_counts[("age60", "Overall")]
        == {
            "n": 15048,
            "Death_AllCause": 6191,
            "Death_Cancer": 1253,
            "Death_CVD": 2033,
            "Death_CLRD": 360,
            "Death_Alzheimer": 282,
            "Death_Diabetes_UCOD": 214,
            "Death_Kidney": 148,
        },
        expected_counts[("age60", "Overall")],
    )

    script_text = SCRIPT_R.read_text(encoding="utf-8")
    check(
        "correct cardiovascular mortality code lock",
        'UCOD_LEADING %in% c(1, 5)' in script_text
        and 'UCOD_LEADING %in% c(1, 3)' not in script_text,
        "CVD=UCOD 1+5; UCOD 3 is CLRD",
    )

    design = pd.read_csv(OUT / "older_nhanes_design_audit.csv")
    design_map = dict(zip(design["item"].astype(str), design["value"].astype(str)))
    design_ok = (
        design_map.get("full_design_rows_before_domain_subset") == "44772"
        and design_map.get("full_design_cycles") == "9"
        and design_map.get("analytic_cycle_scope")
        == "1999-2016; 9 cycles; 44,772 participants with valid follow-up time"
        and design_map.get("weight_variable") == "WTMEC2YR"
        and design_map.get("pooled_weight_formula") == "WTMEC2YR/9"
        and design_map.get("primary_domain") == "Age >= 60 years"
        and design_map.get("age65_sensitivity_domain") == "Age >= 65 years"
        and design_map.get("primary_outcome") == "All-cause mortality"
        and design_map.get("phenoage_acceleration_definition")
        == "survey-weighted residual from PhenoAge regressed on chronological age within each domain"
        and design_map.get("weighted_composite_generated") == "FALSE"
        and design_map.get("hospital_model_retrained") == "FALSE"
        and design_map.get("all_age_nhanes_role") == "SENSITIVITY"
    )
    check("survey-domain design and analysis locks", design_ok, design_map)

    scaling = pd.read_csv(OUT / "older_nhanes_domain_scaling.csv")
    check(
        "domain-specific survey scaling completeness",
        len(scaling) == 8
        and set(scaling["domain"]) == {"age60", "age65"}
        and scaling["survey_weighted_sd"].gt(0).all()
        and np.isfinite(scaling[["survey_weighted_mean", "survey_weighted_sd"]]).all().all(),
        scaling.groupby("domain").size().to_dict(),
    )
    pheno_scaling = scaling[scaling["standardized_variable"].str.startswith("PHENO_z_")]
    check(
        "PhenoAge acceleration uses weighted age-regression residuals",
        len(pheno_scaling) == 2
        and pheno_scaling["age_regression_intercept"].notna().all()
        and pheno_scaling["age_regression_slope"].notna().all()
        and pheno_scaling["survey_weighted_mean"].abs().lt(1e-8).all()
        and pheno_scaling["standardization_definition"].str.contains(
            "residual from PhenoAge ~ chronological age", regex=False
        ).all(),
        pheno_scaling[
            [
                "domain",
                "survey_weighted_mean",
                "survey_weighted_sd",
                "age_regression_intercept",
                "age_regression_slope",
            ]
        ].to_dict(orient="records"),
    )

    coverage = pd.read_csv(OUT / "older_nhanes_final_feature_coverage.csv")
    component_coverage = pd.read_csv(OUT / "older_nhanes_component_input_coverage.csv")
    all_cycle_coverage = coverage[coverage["cycle"] == "ALL"]
    check(
        "frozen feature coverage audit completeness",
        len(coverage) == 520
        and len(all_cycle_coverage) == 52
        and coverage["missing_pct"].between(0, 100).all(),
        {"rows": len(coverage), "all_cycle_rows": len(all_cycle_coverage)},
    )
    zero_map = {
        (r.domain, r.outcome): int(r.zero_observed_n)
        for r in component_coverage.itertuples(index=False)
    }
    check(
        "individual component observation gate audited",
        zero_map[("age60", "Outcome_NutriMetab")] == 0
        and zero_map[("age60", "Outcome_TumorBurden")] == 0
        and zero_map[("age60", "Outcome_TreatComp")] == 333
        and zero_map[("age65", "Outcome_TreatComp")] == 250,
        zero_map,
    )

    separate = pd.read_csv(OUT / "older_nhanes_allcause_component_separate.csv")
    joint = pd.read_csv(OUT / "older_nhanes_allcause_components_joint.csv")
    joint_wald = pd.read_csv(OUT / "older_nhanes_allcause_joint_wald.csv")
    check(
        "all-cause result table completeness",
        len(separate) == 48 and len(joint) == 36 and len(joint_wald) == 12,
        {"separate": len(separate), "joint": len(joint), "joint_wald": len(joint_wald)},
    )

    def result_values_valid(frame: pd.DataFrame) -> bool:
        numeric = frame[
            [
                "hazard_ratio_per_domain_weighted_sd",
                "ci_lower",
                "ci_upper",
                "p_value",
            ]
        ].astype(float)
        return bool(
            np.isfinite(numeric.to_numpy()).all()
            and numeric["hazard_ratio_per_domain_weighted_sd"].gt(0).all()
            and numeric["ci_lower"].gt(0).all()
            and numeric["ci_upper"].ge(numeric["ci_lower"]).all()
            and numeric["p_value"].between(0, 1).all()
        )

    check(
        "finite coherent all-cause estimates",
        result_values_valid(separate) and result_values_valid(joint),
        "separate and joint",
    )
    primary_separate = separate[
        (separate["domain"] == "age60") & (separate["analysis_role"] == "PRIMARY")
    ]
    check(
        "primary age60 stratified all-cause results",
        len(primary_separate) == 12
        and set(primary_separate["population"]) == {"Overall", "Cancer", "Noncancer"}
        and set(primary_separate["component"]) == {"NM", "TB", "TC", "PHENO"},
        primary_separate.groupby("population").size().to_dict(),
    )
    tc_n = int(
        primary_separate[
            (primary_separate["population"] == "Overall")
            & (primary_separate["component"] == "TC")
        ]["n"].iloc[0]
    )
    nm_n = int(
        primary_separate[
            (primary_separate["population"] == "Overall")
            & (primary_separate["component"] == "NM")
        ]["n"].iloc[0]
    )
    check(
        "zero-input TC scores excluded from modeling",
        tc_n < nm_n,
        {"TC_n": tc_n, "NM_n": nm_n},
    )

    interactions = pd.read_csv(OUT / "older_nhanes_cancer_history_interactions.csv")
    age60_interactions = interactions[interactions["domain"] == "age60"]
    check(
        "cancer-history interaction completeness",
        len(interactions) == 6
        and len(age60_interactions) == 3
        and interactions["interaction_p"].between(0, 1).all()
        and interactions["interaction_fdr_bh"].between(0, 1).all(),
        age60_interactions[
            ["component", "interaction_ratio_of_hrs", "interaction_p", "interaction_fdr_bh"]
        ].to_dict(orient="records"),
    )

    incremental = pd.read_csv(OUT / "older_nhanes_phenoage_incremental_tests.csv")
    expected_comparisons = {"M0_to_M1", "M0_to_M2", "M1_to_M3", "M2_to_M3"}
    same_cohort_ok = True
    for _, d in incremental.groupby(["domain", "population"]):
        same_cohort_ok &= (
            set(d["comparison"]) == expected_comparisons
            and d["n"].nunique() == 1
            and d["events"].nunique() == 1
            and d["common_complete_case_cohort"].astype(str).str.lower().eq("true").all()
        )
    check(
        "PhenoAge incremental comparisons use paired common cohorts",
        len(incremental) == 24
        and same_cohort_ok
        and incremental["design_adjusted_block_wald_p"].between(0, 1).all(),
        incremental.groupby(["domain", "population"])["n"].first().to_dict(),
    )

    cause_separate = pd.read_csv(OUT / "older_nhanes_cause_specific_component_separate.csv")
    cause_joint = pd.read_csv(OUT / "older_nhanes_cause_specific_components_joint.csv")
    cause_wald = pd.read_csv(OUT / "older_nhanes_cause_specific_joint_wald.csv")
    cause_ok = (
        len(cause_separate) == 48
        and len(cause_joint) == 48
        and len(cause_wald) == 16
        and result_values_valid(cause_separate)
        and result_values_valid(cause_joint)
        and cause_separate["p_fdr_bh"].between(0, 1).all()
        and cause_joint["p_fdr_bh"].between(0, 1).all()
    )
    check(
        "cause-specific result completeness and FDR",
        cause_ok,
        {
            "separate": len(cause_separate),
            "joint": len(cause_joint),
            "wald": len(cause_wald),
        },
    )
    check(
        "key secondary causes correctly labeled",
        set(
            cause_separate.loc[
                cause_separate["analysis_role"] == "KEY_SECONDARY", "mortality_outcome"
            ]
        )
        == {"Death_Cancer", "Death_CVD"},
        sorted(cause_separate["mortality_outcome"].unique()),
    )

    all_age_files_exist = all(
        (OUT / name).exists()
        for name in (
            "nested_a_primary_nhanes_survey_component_separate.csv",
            "nested_a_primary_nhanes_survey_components_joint.csv",
            "nested_a_primary_nhanes_survey_joint_wald.csv",
        )
    )
    check("prior all-age NHANES sensitivity retained", all_age_files_exist, all_age_files_exist)

    checks_df = pd.DataFrame(checks)
    checks_df.to_csv(OUT / "older_nhanes_qa_checks.csv", index=False, encoding="utf-8-sig")

    primary_joint = joint[
        (joint["domain"] == "age60") & (joint["analysis_role"] == "PRIMARY")
    ][
        [
            "population",
            "component",
            "n",
            "events",
            "hazard_ratio_per_domain_weighted_sd",
            "ci_lower",
            "ci_upper",
            "p_value",
        ]
    ]
    primary_incremental = incremental[incremental["domain"] == "age60"]
    primary_key_causes = cause_separate[
        (cause_separate["domain"] == "age60")
        & (cause_separate["population"] == "Overall")
        & (cause_separate["analysis_role"] == "KEY_SECONDARY")
    ][
        [
            "mortality_outcome",
            "component",
            "n",
            "events",
            "hazard_ratio_per_domain_weighted_sd",
            "ci_lower",
            "ci_upper",
            "p_value",
            "p_fdr_bh",
        ]
    ]
    summary = {
        "all_checks_passed": bool(checks_df["passed"].all()),
        "checks_total": int(len(checks_df)),
        "checks_passed": int(checks_df["passed"].sum()),
        "checks_failed": int((~checks_df["passed"]).sum()),
        "age60_source_counts": expected_counts[("age60", "Overall")],
        "age65_source_counts": expected_counts[("age65", "Overall")],
        "primary_age60_joint_allcause": primary_joint.to_dict(orient="records"),
        "primary_age60_interactions": age60_interactions.to_dict(orient="records"),
        "primary_age60_phenoage_incremental_tests": primary_incremental.to_dict(
            orient="records"
        ),
        "primary_age60_key_cause_specific": primary_key_causes.to_dict(orient="records"),
        "interpretation_locks": [
            "Frozen all-age pan-cancer models were transported to older NHANES adults; models were not retrained in NHANES.",
            "Scores identify mortality risk or laboratory vulnerability, not cancer, metastasis, or treatment complications.",
            "All-cause mortality is the sole primary outcome; cancer and cardiovascular mortality are key secondary outcomes.",
            "CVD mortality is UCOD 1 plus 5; UCOD 3 is chronic lower respiratory disease.",
            "Incremental value is supported by design-adjusted block Wald tests on common complete-case cohorts; no unpaired C-index superiority claim is made.",
            "PhenoAge Acceleration is a survey-weighted age-regression residual within each older-adult domain, not the simple PhenoAge minus age difference.",
            "No AUC-weighted composite was generated.",
        ],
    }
    (OUT / "older_nhanes_qa_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if not summary["all_checks_passed"]:
        failed = checks_df.loc[~checks_df["passed"], ["check", "evidence"]]
        raise SystemExit("Older NHANES QA failed:\n" + failed.to_string(index=False))


if __name__ == "__main__":
    main()
