# Asset Confirmation Table
# Figure | Plot type | Production assets inspected | Adaptation decision
# S1 | participant flow + SMD + forest | SankeyDiagram, BarComparison | use compact branch flow and lollipop/forest; no decorative Sankey links
# S2 | dose-response small multiples | LineTrend, multipanel | reuse consistent line/ribbon grammar with shared axes and direct P-value labels
# S3/S5 | paired performance differences | BarComparison, AUROC | use CI forest plots because paired uncertainty, not raw bars or ROC curves, supports the claim
# S4 | cycle-specific forest + heterogeneity | LineTrend, multipanel | use connected cycle effects plus explicit interaction P values
# S6 | temporal calibration | LineTrend, multipanel | observed-versus-predicted quintiles with identity line and uncertainty
# S7 | competing-risk forest | BarComparison, multipanel | use two aligned CI forests for relative and absolute effects

from __future__ import annotations

import csv
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Academic Figure Skill Typography Baseline — COPY VERBATIM, place at TOP of script
import matplotlib as mpl
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 8,
    "figure.titlesize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "legend.frameon": False,
})

# Academic Figure Skill Nature/Cell/Science Color Palette -- COPY VERBATIM
CATEGORICAL = ["#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666"]
CATEGORICAL_EXTENDED = [
    "#2166AC", "#B2182B", "#1B7837", "#F1A340", "#762A83", "#666666",
    "#4393C3", "#D6604D", "#5AAE61", "#B35806", "#9970AB", "#999999",
]
DIVERGING   = ["#2166AC", "#F7F7F7", "#B2182B"]
SEQUENTIAL  = ["#F7FBFF", "#6BAED6", "#08306B"]
ACCENT_RED  = "#B2182B"
GREY        = "#999999"
BLACK       = "#222222"

# Academic Figure Skill Export Baseline — COPY VERBATIM
mpl.rcParams.update({
    "pdf.fonttype": 42,         # TrueType font embedding
    "svg.fonttype": "none",     # editable text in SVG
    "savefig.bbox": "tight",    # trim whitespace
    "savefig.dpi": 300,
})

def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs" / "supplement_v2"
FIG = ROOT / "figures" / "supplement_v2"
FIG.mkdir(parents=True, exist_ok=True)

COLORS = {
    "clinical": "#666666",
    "overall": "#B2182B",
    "pheno": "#2166AC",
    "component": "#F1A340",
    "lab26": "#1B7837",
    "nm": "#762A83",
    "tb": "#B35806",
    "tc": "#5AAE61",
    "cancer": "#B2182B",
    "cvd": "#2166AC",
    "light": "#E9EDF2",
    "ink": "#202124",
}
mpl.rcParams.update({"savefig.facecolor": "white", "figure.facecolor": "white"})


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(OUT / name, encoding="utf-8-sig")


def clean_axis(ax: plt.Axes, grid_axis: str | None = None) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color="#D9DDE2", linewidth=0.55, alpha=0.75)
        ax.set_axisbelow(True)


def panel_label(ax: plt.Axes, letter: str) -> None:
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")


manifest_rows: list[dict[str, str | int | float]] = []


def save(fig: plt.Figure, stem: str, claim: str) -> None:
    png = FIG / f"{stem}.png"
    pdf = FIG / f"{stem}.pdf"
    save_cns_figure(fig, str(FIG / stem))
    fig.savefig(FIG / f"{stem}_800dpi.png", dpi=800, bbox_inches="tight", facecolor="white")
    width, height = fig.get_size_inches()
    manifest_rows.append(
        {
            "figure": stem,
            "png": png.name,
            "pdf": pdf.name,
            "width_inches": width,
            "height_inches": height,
            "png_dpi": 300,
            "claim": claim,
        }
    )
    plt.close(fig)


def figure_s1() -> None:
    flow = read("selection_participant_flow.csv")
    smd = read("selection_weighted_baseline_smd.csv")
    assoc = read("selection_ipw_association_sensitivity.csv")
    fig = plt.figure(figsize=(7.2, 7.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.22])

    ax = fig.add_subplot(gs[0, :])
    ax.set_axis_off()
    source = flow.loc[flow.node == "source_age60"].iloc[0]
    source_box = (0.36, 0.70, 0.28, 0.20)
    box = FancyBboxPatch(
        (source_box[0], source_box[1]), source_box[2], source_box[3],
        boxstyle="round,pad=0.012,rounding_size=0.012", linewidth=1.1,
        edgecolor=COLORS["ink"], facecolor="#F5F6F7", transform=ax.transAxes,
    )
    ax.add_patch(box)
    ax.text(0.50, 0.80, f"NHANES age ≥60 years\nN={int(source.n):,}; deaths={int(source.deaths):,}",
            ha="center", va="center", transform=ax.transAxes, fontsize=9, fontweight="bold")
    branches = [
        ("overall_available", "Overall\navailable", 0.02, COLORS["overall"]),
        ("component_fully_adjusted", "Three-component\nassociation", 0.265, COLORS["component"]),
        ("paired_score_only", "Direct paired\ndiscrimination", 0.51, COLORS["pheno"]),
        ("incremental_common_complete", "Common incremental\ncohort", 0.755, COLORS["lab26"]),
    ]
    for node, label, x0, color in branches:
        row = flow.loc[flow.node == node].iloc[0]
        ax.add_patch(FancyArrowPatch((0.50, 0.70), (x0 + 0.105, 0.48), arrowstyle="-|>",
                                     mutation_scale=8, linewidth=0.8, color="#7A7D81",
                                     transform=ax.transAxes, connectionstyle="arc3,rad=0"))
        ax.add_patch(FancyBboxPatch((x0, 0.20), 0.21, 0.25, boxstyle="round,pad=0.010,rounding_size=0.01",
                                    linewidth=1.0, edgecolor=color, facecolor=mpl.colors.to_rgba(color, 0.09),
                                    transform=ax.transAxes))
        ax.text(x0 + 0.105, 0.385, label, ha="center", va="center", transform=ax.transAxes,
                fontsize=7.4, fontweight="bold", color=color)
        ax.text(x0 + 0.105, 0.270,
                f"N={int(row.n):,}\nDeaths={int(row.deaths):,}\nRetained={float(row.retained_percent):.1f}%",
                ha="center", va="center", transform=ax.transAxes, fontsize=7)
    ax.text(0.01, 0.03, "Branches are analysis-specific, not sequential exclusions.", transform=ax.transAxes,
            fontsize=7, color="#555555")
    panel_label(ax, "a")

    ax = fig.add_subplot(gs[1, 0])
    smd2 = smd.sort_values("absolute_smd", ascending=False).head(14).copy().iloc[::-1]
    labels = []
    for r in smd2.itertuples():
        level = "" if pd.isna(r.level) else f": {r.level}"
        labels.append(f"{r.variable_label}{level}")
    colors = [COLORS["overall"] if v >= 0.1 else COLORS["clinical"] for v in smd2.absolute_smd]
    ax.hlines(np.arange(len(smd2)), 0, smd2.absolute_smd, color=colors, linewidth=1.2)
    ax.scatter(smd2.absolute_smd, np.arange(len(smd2)), s=20, color=colors, zorder=3)
    ax.axvline(0.10, color="#444444", linestyle="--", linewidth=0.8)
    ax.set_yticks(np.arange(len(smd2)), labels)
    ax.set_xlabel("Absolute standardized mean difference")
    ax.set_title("Largest included–excluded imbalances")
    clean_axis(ax, "x")
    panel_label(ax, "b")

    ax = fig.add_subplot(gs[1, 1])
    keep = assoc[
        ((assoc.model == "M3_CLINICAL_OVERALL_PHENO") & assoc.term.isin(["OVERALL_z_age60", "PHENO_z_age60"]))
        | ((assoc.model == "M5_CLINICAL_COMPONENTS_PHENO") & assoc.term.isin(["NM_z_age60", "TB_z_age60", "TC_z_age60"]))
    ].copy()
    order = ["OVERALL_z_age60", "PHENO_z_age60", "NM_z_age60", "TB_z_age60", "TC_z_age60"]
    label_map = {
        "OVERALL_z_age60": "Overall",
        "PHENO_z_age60": "PhenoAge acceleration",
        "NM_z_age60": "NM",
        "TB_z_age60": "TB",
        "TC_z_age60": "TC",
    }
    positions = {term: len(order) - 1 - i for i, term in enumerate(order)}
    for j, (weighting, marker, color) in enumerate(
        [("ORIGINAL_SURVEY", "o", COLORS["clinical"]), ("COMPLETE_CASE_IPW", "D", COLORS["overall"])]
    ):
        d = keep[keep.weighting == weighting]
        y = np.array([positions[t] for t in d.term]) + (-0.10 if j == 0 else 0.10)
        ax.errorbar(d.hazard_ratio, y,
                    xerr=np.vstack([d.hazard_ratio - d.ci_lower, d.ci_upper - d.hazard_ratio]),
                    fmt=marker, color=color, markersize=4.2, capsize=2, linewidth=0.9,
                    label="Original survey weights" if j == 0 else "Complete-case IPW")
    ax.axvline(1, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_yticks([positions[t] for t in order], [label_map[t] for t in order])
    ax.set_xlabel("Adjusted hazard ratio per SD")
    ax.set_xscale("log")
    ax.set_xlim(0.9, 1.55)
    ax.set_xticks([0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
    ax.set_xticklabels(["0.9", "1.0", "1.1", "1.2", "1.3", "1.4", "1.5"])
    ax.set_title("Selection-weight sensitivity")
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=COLORS["clinical"], linewidth=0.9,
                   markersize=4, label="Original survey weights"),
            Line2D([0], [0], marker="D", color=COLORS["overall"], linewidth=0.9,
                   markersize=4, label="Complete-case IPW"),
        ],
        frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2,
    )
    clean_axis(ax, "x")
    panel_label(ax, "c")
    save(fig, "Figure_S1_selection_bias", "Complete-case selection is substantial and cycle-structured, but IPW does not materially change score associations.")


def figure_s2() -> None:
    tests = read("restricted_cubic_spline_tests.csv").set_index("score")
    ph = read("ph_assumption_tests.csv").set_index("score")
    curves = read("restricted_cubic_spline_curves.csv")
    order = ["NM", "TB", "TC", "OVERALL", "PHENO"]
    labels = {
        "NM": "NM",
        "TB": "TB",
        "TC": "TC",
        "OVERALL": "Overall",
        "PHENO": "PhenoAge acceleration",
    }
    score_colors = {"NM": COLORS["nm"], "TB": COLORS["tb"], "TC": COLORS["tc"],
                    "OVERALL": COLORS["overall"], "PHENO": COLORS["pheno"]}
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0), constrained_layout=True)
    for i, score in enumerate(order):
        ax = axes.flat[i]
        d = curves[curves.score == score].sort_values("score_value")
        color = score_colors[score]
        ax.fill_between(d.score_value, d.ci_lower, d.ci_upper, color=color, alpha=0.16, linewidth=0)
        ax.plot(d.score_value, d.hazard_ratio, color=color, linewidth=1.6)
        ax.axhline(1, color="#666666", linestyle="--", linewidth=0.7)
        ax.axvline(float(d.reference_value.iloc[0]), color="#AAAAAA", linestyle=":", linewidth=0.6)
        ax.set_yscale("log")
        ax.set_ylim(0.50, 5.4)
        ax.set_yticks([0.6, 1, 2, 3, 4])
        ax.set_yticklabels(["0.6", "1", "2", "3", "4"])
        ax.set_title(labels[score], color=color, fontweight="bold")
        ax.text(0.03, 0.96,
                f"P(nonlin), Holm={tests.loc[score, 'spline_nonlinearity_p_holm']:.3g}\n"
                f"P(PH), Holm={ph.loc[score, 'schoenfeld_score_p_holm']:.3g}",
                transform=ax.transAxes, va="top", fontsize=6.7)
        ax.set_xlabel("Score (survey-weighted SD)")
        if i % 3 == 0:
            ax.set_ylabel("Hazard ratio (95% CI)")
        clean_axis(ax, "y")
        panel_label(ax, chr(ord("a") + i))
    axes.flat[5].set_axis_off()
    axes.flat[5].text(0.05, 0.86, "Interpretation", fontsize=9, fontweight="bold", transform=axes.flat[5].transAxes)
    axes.flat[5].text(
        0.05, 0.72,
        "NM and PhenoAge show nonlinear\nassociations after Holm correction.\n\n"
        "TB, TC, Overall and PhenoAge\nshow time-varying hazard effects;\nNM does not.",
        fontsize=7.5, linespacing=1.45, va="top", transform=axes.flat[5].transAxes,
    )
    save(fig, "Figure_S2_ph_nonlinearity", "NM and PhenoAge show nonlinear associations; TB, TC, Overall and PhenoAge show time-varying associations.")


def _paired_forest(ax: plt.Axes, data: pd.DataFrame, metrics: list[str], title: str) -> None:
    comparison_order = [
        "Overall beyond clinical base",
        "PhenoAge beyond clinical base",
        "Overall beyond PhenoAge",
        "PhenoAge beyond Overall",
        "Three components versus Overall when both include PhenoAge",
        "Three components beyond PhenoAge",
    ]
    comp_colors = {
        "Overall beyond clinical base": COLORS["overall"],
        "PhenoAge beyond clinical base": COLORS["pheno"],
        "Overall beyond PhenoAge": "#C77C79",
        "PhenoAge beyond Overall": "#6E9DBD",
        "Three components versus Overall when both include PhenoAge": COLORS["component"],
        "Three components beyond PhenoAge": "#4B7C6B",
    }
    sub = data[data.metric.isin(metrics) & data.comparison.isin(comparison_order)].copy()
    ylabels, y, point, low, high, color, significant = [], [], [], [], [], [], []
    k = 0
    for metric in metrics:
        for comp in comparison_order:
            row = sub[(sub.metric == metric) & (sub.comparison == comp)].iloc[0]
            direction = 1 if bool(row.higher_is_better) else -1
            metric_label = {"uno_c15": "Uno C (15 y)", "auc5": "AUC (5 y)",
                            "auc10": "AUC (10 y)", "auc15": "AUC (15 y)"}[metric]
            comp_label = {
                "Overall beyond clinical base": "Overall vs clinical",
                "PhenoAge beyond clinical base": "PhenoAge vs clinical",
                "Overall beyond PhenoAge": "Overall added to PhenoAge",
                "PhenoAge beyond Overall": "PhenoAge added to Overall",
                "Three components versus Overall when both include PhenoAge": "Components vs Overall; both +PhenoAge",
                "Three components beyond PhenoAge": "Components added to PhenoAge",
            }[comp]
            ylabels.append(f"{metric_label} · {comp_label}")
            y.append(k)
            point.append(direction * row.paired_difference_a_minus_b * 100)
            a = direction * row.difference_ci_lower * 100
            b = direction * row.difference_ci_upper * 100
            low.append(min(a, b)); high.append(max(a, b))
            color.append(comp_colors[comp]); significant.append(row.paired_p_holm < 0.05)
            k += 1
        k += 0.55
    y = np.asarray(y)
    for yi, p, lo, hi, c, sig in zip(y, point, low, high, color, significant):
        ax.plot([lo, hi], [yi, yi], color=c, linewidth=1.0)
        ax.scatter([p], [yi], s=25, marker="s" if sig else "o", facecolor=c if sig else "white",
                   edgecolor=c, linewidth=0.9, zorder=3)
    ax.axvline(0, color="#3C4043", linewidth=0.8)
    ax.set_yticks(y, ylabels)
    ax.invert_yaxis()
    ax.set_xlabel("Benefit-oriented paired difference (percentage points)")
    ax.set_title(title)
    clean_axis(ax, "x")


def figure_s3() -> None:
    inc = read("incremental_prediction_estimates.csv")
    cmp = read("incremental_prediction_paired_comparisons.csv")
    fig = plt.figure(figsize=(7.2, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[0.92, 1.8])
    ax = fig.add_subplot(gs[0])
    metric_order = ["uno_c15", "auc5", "auc10", "auc15"]
    model_order = ["M0_CLINICAL", "M1_CLINICAL_OVERALL", "M2_CLINICAL_PHENO",
                   "M3_CLINICAL_OVERALL_PHENO", "M4_CLINICAL_COMPONENTS", "M5_CLINICAL_COMPONENTS_PHENO"]
    model_short = {
        "M0_CLINICAL": "Clinical",
        "M1_CLINICAL_OVERALL": "+Overall",
        "M2_CLINICAL_PHENO": "+PhenoAge",
        "M3_CLINICAL_OVERALL_PHENO": "+Overall+PhenoAge",
        "M4_CLINICAL_COMPONENTS": "+NM+TB+TC",
        "M5_CLINICAL_COMPONENTS_PHENO": "+NM+TB+TC+PhenoAge",
    }
    model_color = {
        "M0_CLINICAL": COLORS["clinical"], "M1_CLINICAL_OVERALL": COLORS["overall"],
        "M2_CLINICAL_PHENO": COLORS["pheno"], "M3_CLINICAL_OVERALL_PHENO": "#6D4C7D",
        "M4_CLINICAL_COMPONENTS": COLORS["component"], "M5_CLINICAL_COMPONENTS_PHENO": "#4B7C6B",
    }
    xbase = np.arange(len(metric_order))
    offsets = np.linspace(-0.25, 0.25, len(model_order))
    for off, model in zip(offsets, model_order):
        d = inc[inc.model == model].set_index("metric").loc[metric_order]
        ax.errorbar(xbase + off, d.estimate,
                    yerr=np.vstack([d.estimate - d.ci_lower, d.ci_upper - d.estimate]),
                    fmt="o", markersize=3.6, capsize=1.7, linewidth=0.8,
                    color=model_color[model], label=model_short[model])
    ax.set_xticks(xbase, ["Uno C (15 y)", "AUC (5 y)", "AUC (10 y)", "AUC (15 y)"])
    ax.set_ylabel("Discrimination (95% CI)")
    ax.set_ylim(0.72, 0.85)
    ax.set_title("Apparent performance on the common cohort (N=9,319; deaths=4,177)")
    ax.legend(frameon=False, ncol=3, loc="lower right")
    clean_axis(ax, "y")
    panel_label(ax, "a")

    ax = fig.add_subplot(gs[1])
    _paired_forest(ax, cmp, ["uno_c15", "auc5", "auc10", "auc15"],
                   "Paired discrimination increments; filled squares indicate Holm P<0.05")
    panel_label(ax, "b")
    save(fig, "Figure_S3_incremental_prediction", "Overall and the three components provide incremental discrimination beyond the clinical model and PhenoAge.")


def figure_s4() -> None:
    effects = read("cycle_specific_score_effects.csv")
    tests = read("cycle_stability_interactions.csv").set_index("score")
    sens = read("cycle_stability_demographic_adjusted_interactions.csv").set_index("score")
    d = effects[(effects.standardization == "global age>=60 survey-weighted SD") & (effects.estimable == True)].copy()
    cycle_order = sorted(d.cycle.unique())
    score_order = ["NM", "TB", "TC", "OVERALL", "PHENO"]
    score_colors = {"NM": COLORS["nm"], "TB": COLORS["tb"], "TC": COLORS["tc"],
                    "OVERALL": COLORS["overall"], "PHENO": COLORS["pheno"]}
    fig = plt.figure(figsize=(7.2, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(5, 2, width_ratios=[3.1, 1.25], wspace=0.12)
    for i, score in enumerate(score_order):
        ax = fig.add_subplot(gs[i, 0])
        z = d[d.score == score].set_index("cycle").reindex(cycle_order).dropna(subset=["hazard_ratio"])
        x = [cycle_order.index(v) for v in z.index]
        ax.plot(x, z.hazard_ratio, color=score_colors[score], linewidth=1.0, alpha=0.7)
        ax.errorbar(x, z.hazard_ratio,
                    yerr=np.vstack([z.hazard_ratio - z.ci_lower, z.ci_upper - z.hazard_ratio]),
                    fmt="o", color=score_colors[score], markersize=3.8, capsize=1.8, linewidth=0.9)
        ax.axhline(1, color="#666666", linestyle="--", linewidth=0.65)
        ax.set_xlim(-0.4, len(cycle_order)-0.6)
        ax.set_ylabel(score, rotation=0, ha="right", va="center", fontweight="bold", color=score_colors[score])
        if i == len(score_order) - 1:
            ax.set_xticks(range(len(cycle_order)), [c.replace("-", "–") for c in cycle_order], rotation=40, ha="right")
            ax.set_xlabel("NHANES cycle")
        else:
            ax.set_xticks(range(len(cycle_order)), [])
        clean_axis(ax, "y")
        if i == 0:
            ax.set_title("Cycle-specific adjusted HR per global SD")
            panel_label(ax, "a")

    ax = fig.add_subplot(gs[:, 1])
    y = np.arange(len(score_order))[::-1]
    primary = -np.log10(np.clip([tests.loc[s, "cycle_interaction_p_holm"] for s in score_order], 1e-12, 1))
    lighter = -np.log10(np.clip([sens.loc[s, "cycle_interaction_p_holm"] for s in score_order], 1e-12, 1))
    ax.scatter(primary, y + 0.10, color=COLORS["ink"], marker="o", s=24, label="Fully adjusted")
    ax.scatter(lighter, y - 0.10, color=COLORS["overall"], marker="D", s=22, label="Demographic-adjusted sensitivity")
    ax.axvline(-np.log10(0.05), color="#555555", linestyle="--", linewidth=0.8)
    ax.set_yticks(y, score_order)
    ax.set_xlabel("−log10(Holm P interaction)")
    ax.set_title("Temporal heterogeneity")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=1)
    clean_axis(ax, "x")
    panel_label(ax, "b")
    save(fig, "Figure_S4_cycle_stability", "Cancer-trained score interactions were nonsignificant across nine cycles after Holm correction; PhenoAge Acceleration showed heterogeneity across its seven available cycles.")


def figure_s5() -> None:
    est = read("fair_target_temporal_validation_estimates.csv")
    cmp = read("fair_target_paired_comparisons.csv")
    model_order = ["CLINICAL", "COMPONENT3", "PHENO9", "LAB26"]
    model_label = {"CLINICAL": "Clinical", "COMPONENT3": "Clinical + 3 components",
                   "PHENO9": "Clinical + 9 PhenoAge labs", "LAB26": "Clinical + 26 labs"}
    model_color = {"CLINICAL": COLORS["clinical"], "COMPONENT3": COLORS["component"],
                   "PHENO9": COLORS["pheno"], "LAB26": COLORS["lab26"]}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.3), constrained_layout=True, gridspec_kw={"width_ratios": [1.0, 1.1]})
    ax = axes[0]
    metrics = ["uno_c10", "auc5", "auc10"]
    xbase = np.arange(len(metrics))
    offsets = [-0.24, -0.08, 0.08, 0.24]
    for off, model in zip(offsets, model_order):
        d = est[est.model == model].set_index("metric").loc[metrics]
        ax.errorbar(xbase + off, d.estimate,
                    yerr=np.vstack([d.estimate - d.ci_lower, d.ci_upper - d.estimate]),
                    fmt="o", markersize=4, capsize=2, linewidth=0.9,
                    color=model_color[model], label=model_label[model])
    ax.set_xticks(xbase, ["Uno C (10 y)", "AUC (5 y)", "AUC (10 y)"])
    ax.set_ylabel("Temporal validation estimate")
    ax.set_ylim(0.73, 0.86)
    ax.set_title("Same target, algorithm and time split")
    ax.legend(frameon=False, loc="lower right")
    clean_axis(ax, "y")
    panel_label(ax, "a")

    ax = axes[1]
    pair_order = ["26 labs versus nine PhenoAge labs", "Three components versus nine PhenoAge labs"]
    pair_color = {pair_order[0]: COLORS["lab26"], pair_order[1]: COLORS["component"]}
    rows = []
    for metric in metrics:
        for pair in pair_order:
            r = cmp[(cmp.metric == metric) & (cmp.comparison == pair)].iloc[0]
            rows.append((metric, pair, r.paired_difference_a_minus_b * 100,
                         r.difference_ci_lower * 100, r.difference_ci_upper * 100, r.paired_p_holm))
    for i, (metric, pair, point, low, high, p) in enumerate(rows):
        ax.plot([low, high], [i, i], color=pair_color[pair], linewidth=1.0)
        ax.scatter(point, i, s=27, marker="s" if p < 0.05 else "o",
                   facecolor=pair_color[pair] if p < 0.05 else "white",
                   edgecolor=pair_color[pair], linewidth=0.9)
    ax.axvline(0, color="#444444", linewidth=0.8)
    metric_label = {"uno_c10": "Uno C (10 y)", "auc5": "AUC (5 y)", "auc10": "AUC (10 y)"}
    ax.set_yticks(range(len(rows)), [f"{metric_label[m]} · {'26 vs 9 labs' if p == pair_order[0] else '3 components vs 9 labs'}" for m, p, *_ in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Paired difference (percentage points)")
    ax.set_title("Positive values favor first model")
    clean_axis(ax, "x")
    panel_label(ax, "b")
    save(fig, "Figure_S5_fair_target", "Under a matched mortality target, 26 laboratories modestly outperform nine PhenoAge laboratories, while the three-component compression performs worse.")


def figure_s6() -> None:
    groups = read("temporal_calibration_groups.csv")
    metrics = read("temporal_calibration_metrics.csv")
    models = ["M1_CLINICAL_OVERALL", "M2_CLINICAL_PHENO", "M3_CLINICAL_OVERALL_PHENO", "M5_CLINICAL_COMPONENTS_PHENO"]
    labels = {"M1_CLINICAL_OVERALL": "Overall", "M2_CLINICAL_PHENO": "PhenoAge",
              "M3_CLINICAL_OVERALL_PHENO": "Overall + PhenoAge", "M5_CLINICAL_COMPONENTS_PHENO": "Components + PhenoAge"}
    colors = {"M1_CLINICAL_OVERALL": COLORS["overall"], "M2_CLINICAL_PHENO": COLORS["pheno"],
              "M3_CLINICAL_OVERALL_PHENO": "#6D4C7D", "M5_CLINICAL_COMPONENTS_PHENO": "#4B7C6B"}
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.0), constrained_layout=True, gridspec_kw={"width_ratios": [1, 1, 1.05]})
    for j, horizon in enumerate([5, 10]):
        ax = axes[j]
        maxv = 0
        for model in models:
            d = groups[(groups.model == model) & (groups.horizon_years == horizon)]
            obs = d[d.value_type == "observed"].sort_values("group")
            exp = d[d.value_type == "expected"].sort_values("group")
            maxv = max(maxv, obs.ci_upper.max(), exp.estimate.max())
            ax.errorbar(exp.estimate, obs.estimate,
                        yerr=np.vstack([obs.estimate - obs.ci_lower, obs.ci_upper - obs.estimate]),
                        color=colors[model], marker="o", markersize=3.2, linewidth=0.8, capsize=1.5,
                        label=labels[model])
        lim = min(0.75, maxv * 1.08)
        ax.plot([0, lim], [0, lim], color="#555555", linestyle="--", linewidth=0.8)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel("Predicted mortality risk")
        ax.set_ylabel("Observed weighted risk")
        ax.set_title(f"{horizon}-year temporal calibration")
        clean_axis(ax, "both")
        panel_label(ax, chr(ord("a") + j))
        if j == 1:
            ax.legend(frameon=False, loc="upper left", fontsize=6.3)

    ax = axes[2]
    all_models = list(metrics.model.unique())
    short = {"M0_CLINICAL": "Clinical", "M1_CLINICAL_OVERALL": "Overall", "M2_CLINICAL_PHENO": "PhenoAge",
             "M3_CLINICAL_OVERALL_PHENO": "Overall+PhenoAge", "M4_CLINICAL_COMPONENTS": "Components",
             "M5_CLINICAL_COMPONENTS_PHENO": "Components+PhenoAge"}
    for i, model in enumerate(all_models):
        for h, marker, off in [(5, "o", -0.11), (10, "D", 0.11)]:
            r = metrics[(metrics.model == model) & (metrics.metric == f"slope{h}")].iloc[0]
            y = len(all_models) - 1 - i + off
            color = COLORS["pheno"] if "PHENO" in model else (COLORS["overall"] if "OVERALL" in model else COLORS["clinical"])
            ax.errorbar(r.estimate, y, xerr=[[r.estimate-r.ci_lower], [r.ci_upper-r.estimate]],
                        fmt=marker, color=color, markersize=3.5, capsize=1.5, linewidth=0.8)
    ax.axvline(1, color="#444444", linestyle="--", linewidth=0.8)
    ax.set_yticks(range(len(all_models))[::-1], [short[m] for m in all_models])
    ax.set_xlabel("Calibration slope (95% CI)")
    ax.set_title("Temporal calibration slopes\n(circle: 5 y; diamond: 10 y)")
    clean_axis(ax, "x")
    panel_label(ax, "c")
    save(fig, "Figure_S6_temporal_calibration", "Early-cycle models show slopes near one in later cycles, with modest risk-level miscalibration and lower Brier scores for PhenoAge-containing models.")


def figure_s7() -> None:
    fg = read("competing_risk_finegray_effects.csv")
    delta = read("competing_risk_cif_q4_q1_contrasts.csv")
    score_order = ["PHENO", "OVERALL", "TB", "TC", "NM"]
    score_label = {"PHENO": "PhenoAge", "OVERALL": "Overall", "TB": "TB", "TC": "TC", "NM": "NM"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.4), constrained_layout=True, gridspec_kw={"width_ratios": [0.9, 1.3]})
    ax = axes[0]
    y0 = np.arange(len(score_order))[::-1]
    for cause, off, marker, color in [("CANCER", 0.11, "o", COLORS["cancer"]), ("CVD", -0.11, "D", COLORS["cvd"])]:
        d = fg[fg.cause == cause].set_index("score").loc[score_order]
        ax.errorbar(d.subdistribution_hazard_ratio, y0 + off,
                    xerr=np.vstack([d.subdistribution_hazard_ratio-d.ci_lower, d.ci_upper-d.subdistribution_hazard_ratio]),
                    fmt=marker, color=color, markersize=4, capsize=2, linewidth=0.9, label=cause)
    ax.axvline(1, color="#555555", linestyle="--", linewidth=0.8)
    ax.set_yticks(y0, [score_label[s] for s in score_order])
    ax.set_xlabel("Subdistribution HR per SD")
    ax.set_title("Fine–Gray subdistribution associations")
    ax.legend(frameon=False, loc="lower right")
    clean_axis(ax, "x")
    panel_label(ax, "a")

    ax = axes[1]
    sub = delta[delta.score.isin(score_order)].copy()
    rows = []
    for score in score_order:
        for cause in ["CANCER", "CVD"]:
            for h in [5, 10, 15]:
                r = sub[(sub.score == score) & (sub.cause == cause) & (sub.horizon_years == h)].iloc[0]
                rows.append((score, cause, h, r.q4_minus_q1_absolute_risk_difference * 100,
                             r.ci_lower * 100, r.ci_upper * 100, r.p_holm))
    y = np.arange(len(rows))
    for yi, (score, cause, h, p, lo, hi, pval) in zip(y, rows):
        color = COLORS["cancer"] if cause == "CANCER" else COLORS["cvd"]
        ax.plot([lo, hi], [yi, yi], color=color, linewidth=0.85)
        ax.scatter(p, yi, s=20, marker="s" if pval < 0.05 else "o",
                   facecolor=color if pval < 0.05 else "white", edgecolor=color, linewidth=0.8)
    ax.axvline(0, color="#555555", linewidth=0.8)
    ax.set_yticks(y, [f"{score_label[s]} · {c if c == 'CVD' else 'Cancer'} · {h} y" for s, c, h, *_ in rows])
    ax.invert_yaxis()
    ax.set_xlabel("Q4−Q1 absolute cumulative-incidence difference (%)")
    ax.set_title("Weighted risk contrasts; filled squares: Holm P<0.05")
    clean_axis(ax, "x")
    panel_label(ax, "b")
    save(fig, "Figure_S7_competing_risks", "Overall, TB and TC generalize to cancer and cardiovascular cumulative incidence, whereas NM is not independently informative.")


if __name__ == "__main__":
    figure_s1()
    figure_s2()
    figure_s3()
    figure_s4()
    figure_s5()
    figure_s6()
    figure_s7()
    pd.DataFrame(manifest_rows).to_csv(FIG / "figure_manifest.csv", index=False, encoding="utf-8-sig")
    print(f"Created {len(manifest_rows)} publication figures in {FIG}")
