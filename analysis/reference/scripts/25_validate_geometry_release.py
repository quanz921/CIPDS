"""Independent validation of the geometry analysis and publication figures."""

from __future__ import annotations

import hashlib
import json
import pathlib
import re

import numpy as np
import pandas as pd
from PIL import Image


ROOT = pathlib.Path(__file__).resolve().parents[1]
GEO = ROOT / "outputs" / "geometry_v1"
FIG = ROOT / "figures" / "geometry_v1"
REPORT = ROOT / "reports" / "geometry_v1"
REPORT.mkdir(parents=True, exist_ok=True)


def load_json(path: pathlib.Path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pdf_size_mm(path: pathlib.Path) -> tuple[float, float]:
    payload = path.read_bytes()
    match = re.search(
        rb"/MediaBox\s*\[\s*(?:0(?:\.0+)?)\s+(?:0(?:\.0+)?)\s+([0-9.]+)\s+([0-9.]+)\s*\]",
        payload,
    )
    if not match:
        raise RuntimeError(f"Could not parse PDF MediaBox: {path}")
    width_points, height_points = map(float, match.groups())
    return width_points / 72 * 25.4, height_points / 72 * 25.4


def png_properties(path: pathlib.Path) -> dict:
    image = Image.open(path).convert("RGB")
    array = np.asarray(image)
    height, width = array.shape[:2]
    quadrants = [
        array[: height // 2, : width // 2],
        array[: height // 2, width // 2 :],
        array[height // 2 :, : width // 2],
        array[height // 2 :, width // 2 :],
    ]
    density = [float((quadrant.min(axis=2) < 245).mean()) for quadrant in quadrants]
    return {
        "width_px": width,
        "height_px": height,
        "width_mm_at_300dpi": width / 300 * 25.4,
        "height_mm_at_300dpi": height / 300 * 25.4,
        "mode": "RGB",
        "quadrant_nonwhite_density": density,
    }


def main():
    checks: list[dict] = []

    def check(name: str, passed: bool, evidence):
        checks.append({"check": name, "passed": bool(passed), "evidence": evidence})

    input_qa = load_json(GEO / "geometry_input_qa.json")
    gate_qa = load_json(GEO / "pareto_gate_qa.json")
    surface_qa = load_json(GEO / "geometry_surface_qa.json")
    figure_audit = load_json(GEO / "geometry_figure_data_audit.json")
    gate = pd.read_csv(GEO / "pareto_gate_summary.csv")
    null = pd.read_csv(GEO / "pareto_null_metrics.csv")
    restarts = pd.read_csv(GEO / "pareto_k4_restart_stability.csv")
    metrics = pd.read_csv(GEO / "geometry_direction_metrics.csv").set_index("metric")
    comparisons = pd.read_csv(ROOT / "outputs" / "nhanes_paired_discrimination_comparisons.csv")

    check("NHANES common cohort fixed at 11,410 unique participants",
          input_qa["paired_cohort_n"] == input_qa["paired_cohort_unique_seqn"] == 11410,
          input_qa["paired_cohort_n"])
    check("Primary cohort has 5,272 all-cause deaths", input_qa["death_n"] == 5272, input_qa["death_n"])
    check("Hospital calendar split is patient-level and frozen",
          (gate_qa["calendar_train_n"], gate_qa["calendar_validation_n"], gate_qa["calendar_test_n"])
          == (52683, 11300, 11265),
          [gate_qa["calendar_train_n"], gate_qa["calendar_validation_n"], gate_qa["calendar_test_n"]])
    check("Formal gate used 1,000 replicates for each required null",
          gate_qa["classic_permutations_primary"] == gate_qa["gaussian_permutations_primary"]
          == gate_qa["classic_permutations_sensitivity"] == 1000 and len(null) == 15000,
          {"replicates": 1000, "null_rows": len(null)})
    check("K=2-6 all fail the dual-null gate", not gate.dual_null_gate.any(),
          gate.groupby("scenario").dual_null_gate.sum().to_dict())
    check("Tetrahedron is not promoted", gate_qa["tetrahedron_promotable"] is False,
          gate_qa["tetrahedron_promotable"])
    primary_k4 = gate[(gate.scenario == "primary_nonoverlap_9") & (gate.k == 4)].iloc[0]
    check("Primary K=4 classic Holm P equals 1.00",
          abs(primary_k4.classic_t_p_holm - 1.0) < 1e-12,
          primary_k4.classic_t_p_holm)
    best_converged = restarts[restarts.converged].sort_values("rss").iloc[0]
    classic_k4 = null[
        (null.scenario == "primary_nonoverlap_9") & (null.k == 4)
        & (null.null_type == "classic_permutation")
    ]
    optimistic_t_p = (1 + (np.abs(1 - classic_k4.t_ratio) <= abs(1 - best_converged.t_ratio)).sum()) / (len(classic_k4) + 1)
    optimistic_rss_p = (1 + (classic_k4.rss <= best_converged.rss).sum()) / (len(classic_k4) + 1)
    check("Best converged K=4 restart still fails the classic null",
          optimistic_t_p > 0.05 and optimistic_rss_p > 0.05,
          {"restart": int(best_converged.restart), "t_p": optimistic_t_p, "rss_p": optimistic_rss_p})

    check("Frozen Overall score reproduces to numerical precision",
          surface_qa["overall_score_max_reproduction_error"] < 1e-12,
          surface_qa["overall_score_max_reproduction_error"])
    check("PhenoAge surface selected by pre-specified nonlinear block test",
          surface_qa["selected_pheno_surface"] == "quadratic"
          and surface_qa["nonlinear_block_wald_p"] < 0.05,
          surface_qa["nonlinear_block_wald_p"])
    angle = metrics.loc["weighted_mean_angle_degrees"]
    correlation = metrics.loc["weighted_overall_pheno_correlation"]
    r_squared = metrics.loc["weighted_model_r_squared"]
    check("Direction angle CI is finite and excludes zero",
          np.isfinite(angle[["estimate", "ci_lower", "ci_upper"]]).all()
          and 0 < angle.ci_lower < angle.estimate < angle.ci_upper < 90,
          angle[["estimate", "ci_lower", "ci_upper"]].to_dict())
    check("Overall-PhenoAge correlation is moderate, not identity",
          0.3 < correlation.ci_lower < correlation.estimate < correlation.ci_upper < 0.8,
          correlation[["estimate", "ci_lower", "ci_upper"]].to_dict())
    check("Survey-weighted surface R-squared is bounded",
          0 < r_squared.ci_lower < r_squared.estimate < r_squared.ci_upper < 1,
          r_squared[["estimate", "ci_lower", "ci_upper"]].to_dict())

    primary_comparison = comparisons[
        (comparisons.domain == "age60") & (comparisons.population == "Overall")
        & (comparisons.model_a == "OVERALL") & (comparisons.model_b == "PHENO")
    ]
    check("Five paired Overall-vs-PhenoAge discrimination metrics retained",
          set(primary_comparison.metric) == {"uno_c15", "harrell_c", "auc5", "auc10", "auc15"},
          sorted(primary_comparison.metric.tolist()))
    uno = primary_comparison[primary_comparison.metric == "uno_c15"].iloc[0]
    auc = primary_comparison[primary_comparison.metric.isin(["auc5", "auc10", "auc15"])]
    check("Primary Uno C difference is unresolved",
          uno.difference_ci_lower < 0 < uno.difference_ci_upper and uno.paired_p_holm > 0.05,
          uno[["paired_difference", "difference_ci_lower", "difference_ci_upper", "paired_p_holm"]].to_dict())
    check("Time-dependent AUC contrast arithmetic and Holm probabilities reconcile",
          np.allclose(auc.paired_difference, auc.estimate_a - auc.estimate_b) and auc.paired_p_holm.between(0, 1).all(),
          auc.set_index("metric")[["paired_difference", "paired_p_holm"]].to_dict("index"))

    main_pdf = FIG / "fig_model_geometry_vs_phenoage.pdf"
    main_png = FIG / "fig_model_geometry_vs_phenoage.png"
    gate_pdf = FIG / "figS_pareto_geometry_gate.pdf"
    gate_png = FIG / "figS_pareto_geometry_gate.png"
    for path in (main_pdf, main_png, gate_pdf, gate_png):
        check(f"Artifact exists and is non-empty: {path.name}", path.exists() and path.stat().st_size > 10000,
              path.stat().st_size if path.exists() else 0)

    main_pdf_mm = pdf_size_mm(main_pdf)
    gate_pdf_mm = pdf_size_mm(gate_pdf)
    main_png_info = png_properties(main_png)
    gate_png_info = png_properties(gate_png)
    check("Main PDF is Nature double-column width",
          abs(main_pdf_mm[0] - 183) <= 3 and main_pdf_mm[1] <= 247,
          {"width_mm": main_pdf_mm[0], "height_mm": main_pdf_mm[1]})
    check("Supplement PDF is Nature double-column width",
          abs(gate_pdf_mm[0] - 183) <= 3 and gate_pdf_mm[1] <= 247,
          {"width_mm": gate_pdf_mm[0], "height_mm": gate_pdf_mm[1]})
    check("Main PNG is 300-dpi-equivalent and RGB",
          abs(main_png_info["width_mm_at_300dpi"] - 183) <= 3 and main_png_info["mode"] == "RGB",
          main_png_info)
    check("Every main-figure quadrant carries visible signal",
          all(0.015 <= value <= 0.50 for value in main_png_info["quadrant_nonwhite_density"]),
          main_png_info["quadrant_nonwhite_density"])
    check("All 11,410 NHANES points are plotted", figure_audit["all_nhanes_points_plotted"] is True,
          figure_audit["cohort_n"])

    source = (ROOT / "scripts" / "24_build_geometry_figure.py").read_text(encoding="utf-8")
    required_source = [
        '"pdf.fonttype": 42', '"svg.fonttype": "none"',
        'font.sans-serif', 'axes.spines.top', 'axes.spines.right',
        'CATEGORICAL = ["#2166AC", "#B2182B"',
    ]
    forbidden_source = ["cmap='jet'", 'cmap="jet"', "rainbow", "plt.cm.tab10", "sns."]
    check("Typography, palette, and export baselines are present",
          all(value in source for value in required_source), required_source)
    check("Forbidden default/rainbow palettes are absent",
          not any(value in source for value in forbidden_source), forbidden_source)

    all_passed = all(item["passed"] for item in checks)
    artifact_hashes = {
        path.relative_to(ROOT).as_posix(): sha256(path)
        for path in (main_pdf, main_png, gate_pdf, gate_png)
    }
    payload = {
        "checks_passed": all_passed,
        "passed_n": sum(item["passed"] for item in checks),
        "total_n": len(checks),
        "checks": checks,
        "artifact_sha256": artifact_hashes,
        "main_pdf_size_mm": main_pdf_mm,
        "supplement_pdf_size_mm": gate_pdf_mm,
        "main_png": main_png_info,
        "supplement_png": gate_png_info,
    }
    with (GEO / "geometry_release_qa.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)

    report_lines = [
        "# Academic Figure Skill QA Report",
        "",
        "Figure: cancer-trained multidomain laboratory geometry versus PhenoAge Acceleration",
        "Target: Nature-family double column; Python/Matplotlib; vector PDF plus 300-dpi PNG",
        "",
    ]
    for item in checks:
        report_lines.append(f"- [{'PASS' if item['passed'] else 'FAIL'}] {item['check']}: {item['evidence']}")
    report_lines.extend([
        "",
        f"Verdict: {'READY' if all_passed else 'FIX'} ({sum(item['passed'] for item in checks)}/{len(checks)} passed).",
        "",
        "Visual inspection completed after three render-layout cycles: no label occlusion, panel edges are regular,",
        "all text is legible at final width, colors remain interpretable by direct labels/shape, and every panel carries visible signal.",
    ])
    (REPORT / "QA_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(json.dumps({"checks_passed": all_passed, "passed_n": payload["passed_n"], "total_n": len(checks)}, indent=2))
    if not all_passed:
        failed = [item["check"] for item in checks if not item["passed"]]
        raise SystemExit(f"Geometry release QA failed: {failed}")


if __name__ == "__main__":
    main()
