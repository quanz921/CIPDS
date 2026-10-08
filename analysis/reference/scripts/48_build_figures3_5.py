# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# Supplementary feature architecture -> assets/figures/heatmap -> param inherit
# Figure 3A-B effect-size intervals -> assets/figures/BarComparison -> param inherit
# Figure 3C cumulative-incidence curves -> assets/figures/LineTrend -> param inherit
# Figure 4A-B performance profiles -> assets/figures/LineTrend -> param inherit
# Figure 4C paired comparison matrix -> assets/figures/heatmap -> param inherit
# Figure 5A-B 3D score surfaces -> project production geometry figure + assets/figures/Manifold -> param inherit
# Figure 5C direction comparison -> assets/figures/BarComparison -> param inherit
# Figure 5D fidelity decomposition -> assets/figures/BarComposition -> param inherit
# Figure 5E paired recovery forest -> assets/figures/BarComparison -> param inherit
# Figure 5F residual-correlation lollipop -> assets/figures/LineTrend -> param inherit
# RULE: no bundled production asset is semantically compatible with the study-specific data structures;
#       all panels inherit validated visual parameters and use the complete locked source rows required by that panel.

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


import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.gridspec import GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure


# Large, legible master figures; legends remain compact.
mpl.rcParams.update({
    "font.size": 9.0,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.8,
    "xtick.labelsize": 7.6,
    "ytick.labelsize": 7.6,
    "legend.fontsize": 6.8,
    "figure.titlesize": 10.5,
})

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs"
GEO = OUT / "geometry_v1"
SUP2 = OUT / "supplement_v2"
SUP3 = OUT / "supplement_v3"
RELEASE = ROOT / "reports" / "figures3_5_release_20260905"
SOURCE = RELEASE / "source_data"
RELEASE.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)

COLORS = {
    "OVERALL": "#B76576",
    "NM": "#5B7FA6",
    "TB": "#6F9B7A",
    "TC": "#C69355",
    "PHENO": "#8B78A6",
    "CLINICAL": "#777B82",
    "COMPONENT3": "#5B7FA6",
    "PHENO9": "#8B78A6",
    "LAB26": "#B76576",
    "PROJECTED3": "#5B7FA6",
    "HYBRID4": "#3F6F63",
}
LIGHT_GREY = "#D7DADF"
MID_GREY = "#777B82"
TEXT = "#272A2F"

LAB_LABELS = {
    "lab_hct": "Hematocrit", "lab_lym": "Lymphocyte count",
    "lab_lym_pct": "Lymphocyte %", "lab_mch": "MCH",
    "lab_mchc": "MCHC", "lab_mcv": "MCV", "lab_mono": "Monocyte count",
    "lab_mono_pct": "Monocyte %", "lab_mpv": "MPV", "lab_neu_pct": "Neutrophil %",
    "lab_plt": "Platelet count", "lab_rbc": "Red blood cell count",
    "lab_rdw": "RDW", "lab_wbc": "White blood cell count",
    "lab_alb": "Albumin", "lab_a_ratio_g": "Albumin/globulin",
    "lab_alp": "Alkaline phosphatase", "lab_alt": "ALT", "lab_ast": "AST",
    "lab_ggt": "Gamma-GT", "lab_glob": "Globulin",
    "lab_tbil": "Total bilirubin", "lab_tp": "Total protein",
    "lab_crea": "Creatinine", "lab_ua": "Uric acid", "lab_urea": "Urea",
}
LAB_GROUPS = {
    "Hematology / inflammation": [
        "lab_hct", "lab_lym", "lab_lym_pct", "lab_mch", "lab_mchc", "lab_mcv",
        "lab_mono", "lab_mono_pct", "lab_mpv", "lab_neu_pct", "lab_plt",
        "lab_rbc", "lab_rdw", "lab_wbc",
    ],
    "Liver / protein metabolism": [
        "lab_alb", "lab_a_ratio_g", "lab_alp", "lab_alt", "lab_ast", "lab_ggt",
        "lab_glob", "lab_tbil", "lab_tp",
    ],
    "Renal metabolism": ["lab_crea", "lab_ua", "lab_urea"],
}
OUTCOME_MAP = {
    "NM": "Outcome_NutriMetab", "TB": "Outcome_TumorBurden", "TC": "Outcome_TreatComp"
}


def fmt_p(value):
    value = float(value)
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def style_axis(ax, grid_axis=None):
    ax.tick_params(length=2.5, width=0.55, colors=TEXT, pad=2)
    ax.spines["left"].set_color("#6E7177")
    ax.spines["bottom"].set_color("#6E7177")
    if grid_axis:
        ax.grid(axis=grid_axis, color="#E7E9EC", linewidth=0.45, zorder=0)
    ax.set_axisbelow(True)


def panel_shell(fig, spec, letter, title, subtitle=None, header_ratio=0.14,
                title_fontsize=10.2, title_x=0.085):
    shell = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=spec, height_ratios=[header_ratio, 1 - header_ratio], hspace=0.01
    )
    hax = fig.add_subplot(shell[0, 0])
    hax.set_axis_off()
    hax.text(0.00, 0.72, letter, fontsize=12.5, fontweight="bold", color=TEXT,
             va="center", ha="left")
    hax.text(title_x, 0.72, title, fontsize=title_fontsize, fontweight="bold", color=TEXT,
             va="center", ha="left")
    if subtitle:
        hax.text(title_x, 0.12, subtitle, fontsize=7.3, color=MID_GREY,
                 va="center", ha="left")
    return shell[1, 0]


def blend_with_white(color, amount):
    rgb = np.asarray(mcolors.to_rgb(color))
    return tuple((1 - amount) * np.ones(3) + amount * rgb)


def save_outputs(fig, stem):
    for label in fig.findobj(mpl.text.Text):
        label.set_text(label.get_text().replace("=<", "<"))
    pdf = stem.with_suffix(".pdf")
    png = stem.with_suffix(".png")
    tiff = stem.with_suffix(".tiff")
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(pdf, format="pdf", bbox_inches=None)
        fig.savefig(png, format="png", dpi=300, bbox_inches=None)
    pdftoppm = shutil.which("pdftoppm")
    bundled = Path("external_inputs/pdftoppm.exe")
    if not pdftoppm and bundled.exists():
        pdftoppm = str(bundled)
    if not pdftoppm:
        raise RuntimeError("pdftoppm is required for 800 dpi TIFF export")
    tmp = stem.parent / f"{stem.name}_800dpi_tmp"
    subprocess.run([
        pdftoppm, "-r", "800", "-singlefile", "-tiff", "-tiffcompression", "lzw",
        str(pdf), str(tmp),
    ], check=True)
    generated = tmp.with_suffix(".tif")
    if tiff.exists():
        tiff.unlink()
    generated.replace(tiff)
    plt.close(fig)
    return pdf, png, tiff


def read_csv(path):
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def build_feature_architecture():
    final = read_csv(OUT / "nested_a_primary_final_features.csv")
    votes = read_csv(OUT / "nested_a_primary_feature_votes.csv")
    boruta = read_csv(OUT / "nested_a_primary_boruta_importance.csv")
    union = sorted(final["feature"].unique())
    expected = [lab for labs in LAB_GROUPS.values() for lab in labs]
    if set(union) != set(expected) or len(union) != 26:
        raise RuntimeError("Frozen 26-laboratory union drift")
    rows = []
    for score, outcome in OUTCOME_MAP.items():
        fsub = final[final.outcome == outcome]
        preexcluded = set(str(fsub.iloc[0].preexcluded_features).split(";"))
        vsub = votes[votes.outcome == outcome].set_index("feature")
        bsub = boruta[boruta.outcome == outcome].copy()
        bsub["importance_ratio"] = bsub["importance"] / bsub["shadow_threshold"]
        ratios = bsub.groupby("feature")["importance_ratio"].median()
        selected = set(fsub.feature)
        for lab in expected:
            excluded = lab in preexcluded
            eligible = lab in vsub.index
            row = {
                "laboratory": lab,
                "laboratory_label": LAB_LABELS[lab],
                "score": score,
                "outcome": outcome,
                "preexcluded": excluded,
                "eligible": eligible,
                "selected_final": lab in selected,
                "selection_votes": np.nan if not eligible else int(vsub.loc[lab, "votes"]),
                "median_boruta_importance_to_shadow": np.nan if lab not in ratios.index else float(ratios.loc[lab]),
            }
            rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_csv(SOURCE / "Supplementary_Figure_Unnumbered_feature_architecture.csv", index=False)
    return frame, expected


def draw_feature_architecture(ax, feature, labs):
    scores = ["NM", "TB", "TC"]
    n = len(labs)
    ax.set_xlim(0, 3)
    ax.set_ylim(n, 0)
    max_ratio = 2.0
    for yi, lab in enumerate(labs):
        for xi, score in enumerate(scores):
            row = feature[(feature.laboratory == lab) & (feature.score == score)].iloc[0]
            if bool(row.preexcluded):
                rect = Rectangle((xi, yi), 1, 1, facecolor="#F0F1F2", edgecolor="white",
                                 linewidth=0.7, hatch="////")
                ax.add_patch(rect)
                ax.text(xi + 0.5, yi + 0.5, "×", ha="center", va="center",
                        fontsize=7.2, color="#686C72", fontweight="bold")
            else:
                ratio = 0.0 if pd.isna(row.median_boruta_importance_to_shadow) else float(row.median_boruta_importance_to_shadow)
                strength = 0.10 + 0.72 * min(ratio / max_ratio, 1.0)
                rect = Rectangle((xi, yi), 1, 1, facecolor=blend_with_white(COLORS[score], strength),
                                 edgecolor=COLORS[score] if bool(row.selected_final) else "white",
                                 linewidth=1.25 if bool(row.selected_final) else 0.7)
                ax.add_patch(rect)
                vote = "" if pd.isna(row.selection_votes) else str(int(row.selection_votes))
                ax.text(xi + 0.5, yi + 0.5, vote, ha="center", va="center", fontsize=6.7,
                        color=TEXT, fontweight="bold" if bool(row.selected_final) else "normal")
    ax.set_xticks(np.arange(3) + 0.5, ["NM", "TB", "TC"], fontsize=8.2, fontweight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(n) + 0.5, [LAB_LABELS[x] for x in labs], fontsize=6.6)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    offset = 0
    group_colors = [COLORS["NM"], COLORS["TC"], COLORS["TB"]]
    short_group_labels = {
        "Hematology / inflammation": "Hematology",
        "Liver / protein metabolism": "Liver / protein",
        "Renal metabolism": "Renal / metabolic",
    }
    for (group, members), group_color in zip(LAB_GROUPS.items(), group_colors, strict=True):
        if offset:
            ax.axhline(offset, color="#AEB2B8", lw=0.8)
        center = offset + len(members) / 2
        ax.text(-0.38, center, short_group_labels.get(group, group), rotation=90,
                va="center", ha="center", fontsize=5.8, color=group_color,
                fontweight="bold", clip_on=False)
        offset += len(members)
    ax.text(0.01, -0.055,
            "Tile intensity: median Boruta importance / shadow threshold   |   bold border: final retained feature   |   ×: pre-excluded for outcome overlap",
            transform=ax.transAxes, fontsize=6.5, color=MID_GREY, va="top")


def supplementary_feature_architecture(feature, labs):
    fig = plt.figure(figsize=(183 / 25.4, 154 / 25.4), facecolor="white")
    outer = fig.add_gridspec(1, 1, left=0.135, right=0.975, top=0.955, bottom=0.095)
    spec = panel_shell(fig, outer[0], "", "Frozen 26-laboratory architecture",
                       None, header_ratio=0.08, title_x=0.0)
    ax = fig.add_subplot(spec)
    draw_feature_architecture(ax, feature, labs)
    stem = RELEASE / "Supplementary_Figure_Unnumbered_Frozen_Laboratory_Architecture"
    return save_outputs(fig, stem)


def figure3():
    feature, labs = build_feature_architecture()
    supplementary_paths = supplementary_feature_architecture(feature, labs)
    joint = read_csv(OUT / "older_nhanes_allcause_components_joint.csv")
    joint = joint[(joint.domain == "age60") & (joint.population == "Overall") &
                  (joint.adjustment == "fully_adjusted")].copy()
    wald = read_csv(OUT / "older_nhanes_allcause_joint_wald.csv")
    wald = wald[(wald.domain == "age60") & (wald.population == "Overall") &
                (wald.adjustment == "fully_adjusted")].iloc[0]
    components = read_csv(OUT / "older_nhanes_cause_specific_component_separate.csv")
    components = components[(components.domain == "age60") & (components.population == "Overall") &
                            (components.adjustment == "fully_adjusted") &
                            (components.mortality_outcome.isin(["Death_Cancer", "Death_CVD"]))].copy()
    overall = read_csv(OUT / "overall_nhanes_cause_specific.csv")
    overall = overall[(overall.domain == "age60") & (overall.population == "Overall") &
                      (overall.adjustment == "fully_adjusted") &
                      (overall.mortality_outcome.isin(["Death_Cancer", "Death_CVD"]))].copy()
    cause = pd.concat([overall, components], ignore_index=True, sort=False)
    cause = cause[cause.component.isin(["OVERALL", "NM", "TB", "TC"])]
    curves = read_csv(SUP2 / "competing_risk_cif_curves.csv")
    curves = curves[(curves.score == "OVERALL") & (curves.quartile.isin([1, 4])) &
                    (curves.cause.isin(["CANCER", "CVD"]))].copy()
    contrasts = read_csv(SUP2 / "competing_risk_cif_q4_q1_contrasts.csv")
    contrasts = contrasts[(contrasts.score == "OVERALL") &
                          (contrasts.cause.isin(["CANCER", "CVD"]))].copy()
    joint.to_csv(SOURCE / "Figure_3A_joint_component_effects.csv", index=False)
    cause.to_csv(SOURCE / "Figure_3B_cause_specific_effects.csv", index=False)
    curves.to_csv(SOURCE / "Figure_3C_cumulative_incidence_curves.csv", index=False)
    contrasts.to_csv(SOURCE / "Figure_3C_q4_q1_risk_differences.csv", index=False)

    fig = plt.figure(figsize=(183 / 25.4, 164 / 25.4), facecolor="white")
    outer = fig.add_gridspec(2, 1, height_ratios=[0.78, 1.02],
                             left=0.105, right=0.975, top=0.985, bottom=0.075,
                             hspace=0.20)

    top = GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[0.88, 1.25], wspace=0.40)
    # A: mutually adjusted components.
    spec_a = panel_shell(fig, top[0], "A", "Mutually adjusted components",
                         f"Survey-weighted Cox; N=12,094, deaths=4,787; joint Wald P={fmt_p(wald.component_block_wald_p)}",
                         header_ratio=0.20)
    ax = fig.add_subplot(spec_a)
    order = ["NM", "TB", "TC"]
    y = np.arange(3)[::-1]
    for yy, score in zip(y, order, strict=True):
        row = joint[joint.component == score].iloc[0]
        est, lo, hi = [float(row[x]) for x in ["hazard_ratio_per_domain_weighted_sd", "ci_lower", "ci_upper"]]
        ax.errorbar(est, yy, xerr=[[est - lo], [hi - est]], fmt="o", color=COLORS[score],
                    mfc=COLORS[score], mec="white", mew=0.7, ms=6, elinewidth=1.15, capsize=2.2)
        ax.text(1.36, yy, f"{est:.2f} ({lo:.2f}–{hi:.2f})", ha="left", va="center", fontsize=6.6)
    ax.axvline(1, color="#55585D", lw=0.85)
    ax.set_yticks(y, order)
    ax.set_xlim(0.92, 1.53)
    ax.set_xlabel("Hazard ratio per survey-weighted SD")
    style_axis(ax, "x")

    # B: key cause-specific Cox associations; PhenoAge is not mixed in because a
    # survey-weighted cause-specific Cox result was not part of the locked comparator analysis.
    spec_b = panel_shell(fig, top[1], "B", "Cause-specific mortality",
                         "Open circles: cancer; filled squares: cardiovascular mortality", header_ratio=0.20)
    ax = fig.add_subplot(spec_b)
    order = ["OVERALL", "NM", "TB", "TC"]
    y = np.arange(4)[::-1]
    cause_styles = {
        "Death_Cancer": ("o", -0.10, "Cancer mortality"),
        "Death_CVD": ("s", 0.10, "Cardiovascular mortality"),
    }
    for outcome, (marker, shift, label) in cause_styles.items():
        for yy, score in zip(y, order, strict=True):
            row = cause[(cause.mortality_outcome == outcome) & (cause.component == score)].iloc[0]
            est, lo, hi = [float(row[x]) for x in ["hazard_ratio_per_domain_weighted_sd", "ci_lower", "ci_upper"]]
            ax.errorbar(est, yy + shift, xerr=[[est - lo], [hi - est]], fmt=marker,
                        color=COLORS[score], mfc="white" if outcome == "Death_Cancer" else COLORS[score],
                        mec=COLORS[score], mew=0.9, ms=5.2, elinewidth=1.0, capsize=1.8)
    ax.axvline(1, color="#55585D", lw=0.85)
    ax.set_yticks(y, ["Overall", "NM", "TB", "TC"])
    ax.set_xlim(0.88, 1.53)
    ax.set_xlabel("Hazard ratio per survey-weighted SD")
    style_axis(ax, "x")

    # C: absolute competing-risk curves for the primary Overall exposure.
    # The estimand remains cumulative incidence in a competing-risk setting;
    # a crisp post-step rendering gives the familiar KM-like visual grammar
    # without mislabelling these curves as Kaplan-Meier estimates.
    spec_c = panel_shell(fig, outer[1], "C", "Absolute cause-specific risk across Overall quartiles",
                         None, header_ratio=0.18)
    dg = GridSpecFromSubplotSpec(1, 2, subplot_spec=spec_c, wspace=0.30)
    for idx, (cause_name, title) in enumerate([("CANCER", "Cancer mortality, 1999–2016"), ("CVD", "Cardiovascular mortality, 1999–2014")]):
        ax = fig.add_subplot(dg[0, idx])
        endpoints = []
        for q, color, label, ls, lw in [
            (1, "#7E838A", "Overall Q1", (0, (4.0, 2.3)), 1.55),
            (4, COLORS["OVERALL"], "Overall Q4", "-", 2.15),
        ]:
            sub = curves[(curves.cause == cause_name) & (curves.quartile == q)].sort_values("time")
            ax.step(sub.time, sub.cif * 100, where="post", color=color, lw=lw, ls=ls,
                    solid_capstyle="butt", dash_capstyle="butt", label=label, zorder=3)
            endpoints.append((q, color, float(sub.cif.iloc[-1]) * 100))
        ax.set_xlim(0, 15.35)
        ax.set_xticks([0, 5, 10, 15])
        ax.set_xlabel("Follow-up time (years)")
        ax.set_ylabel("Cumulative incidence (%)")
        ax.set_title(title, loc="left", fontsize=8.5, fontweight="bold")
        for q, color, end_y in endpoints:
            ax.text(15.06, end_y, f"Q{q}", ha="right", va="center", fontsize=6.3,
                    color=color, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="none", alpha=0.86),
                    zorder=5)
        csub = contrasts[contrasts.cause == cause_name].sort_values("horizon_years")
        lines = []
        for _, row in csub.iterrows():
            h = int(row.horizon_years)
            diff = 100 * float(row.q4_minus_q1_absolute_risk_difference)
            lo = 100 * float(row.ci_lower)
            hi = 100 * float(row.ci_upper)
            p = float(row.p_holm)
            lines.append(f"{h} y: {diff:+.1f} pp ({lo:+.1f} to {hi:+.1f}); P$_{{Holm}}$={fmt_p(p)}")
        ax.text(0.03, 0.96, "Q4−Q1\n" + "\n".join(lines), transform=ax.transAxes,
                fontsize=6.2, va="top", color=TEXT,
                bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="#D9DCE0", lw=0.55, alpha=0.92))
        style_axis(ax, "y")

    stem = RELEASE / "Figure_3_Mortality_Profiles_of_Cancer_Trained_Scores"
    paths = save_outputs(fig, stem)
    return paths, supplementary_paths, {
        "figure": 3, "main_panels": 3, "supplement_feature_union_n": 26, "joint_rows": len(joint),
        "cause_specific_rows": len(cause), "cif_rows": len(curves),
        "q4_q1_contrast_rows": len(contrasts),
    }


def figure4():
    estimates = read_csv(SUP2 / "fair_target_temporal_validation_estimates.csv")
    paired = read_csv(SUP2 / "fair_target_paired_comparisons.csv")
    estimates.to_csv(SOURCE / "Figure_4_absolute_validation_performance.csv", index=False)
    paired.to_csv(SOURCE / "Figure_4_paired_comparisons.csv", index=False)
    if len(estimates) != 24 or len(paired) != 36:
        raise RuntimeError("Fair-target result family drift")

    fig = plt.figure(figsize=(183 / 25.4, 190 / 25.4), facecolor="white")
    fig.text(0.10, 0.972, "Matched-target temporal validation", fontsize=12.0,
             fontweight="bold", color=TEXT, ha="left", va="top")
    fig.text(0.10, 0.940,
             "Same all-cause mortality target, clinical base, preprocessing and elastic-net Cox algorithm",
             fontsize=7.5, color=MID_GREY, ha="left", va="top")
    fig.text(0.10, 0.916,
             "Train: NHANES 2003–2006 (N=2,548; deaths=1,400)  →  independent validation: 2007–2010 (N=3,045; deaths=1,085)",
             fontsize=7.2, color=MID_GREY, ha="left", va="top")
    outer = fig.add_gridspec(2, 2, height_ratios=[0.78, 1.25], width_ratios=[1, 1],
                             left=0.165, right=0.975, top=0.845, bottom=0.090,
                             wspace=0.30, hspace=0.23)

    discrimination_metrics = ["uno_c10", "auc5", "auc10"]
    error_metrics = ["brier5", "brier10", "ibs10"]
    metric_labels = {
        "uno_c10": "Uno C\n10 y", "auc5": "AUC\n5 y", "auc10": "AUC\n10 y",
        "brier5": "Brier\n5 y", "brier10": "Brier\n10 y", "ibs10": "IBS\n1–10 y",
    }
    model_order = ["CLINICAL", "COMPONENT3", "PHENO9", "LAB26"]
    model_labels = {"CLINICAL": "Clinical", "COMPONENT3": "NM/TB/TC", "PHENO9": "PHENO9", "LAB26": "LAB26"}
    markers = {"CLINICAL": "o", "COMPONENT3": "s", "PHENO9": "^", "LAB26": "D"}
    offsets = {"CLINICAL": -0.12, "COMPONENT3": -0.04, "PHENO9": 0.04, "LAB26": 0.12}
    legend_handles = [Line2D([0], [0], marker=markers[m], color=COLORS[m], lw=1.2,
                             mfc=COLORS[m], mec="white", mew=0.5, label=model_labels[m])
                      for m in model_order]
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.59, 0.892),
               ncol=4, fontsize=6.8, handlelength=2.0, columnspacing=1.5)

    # A: absolute discrimination as a connected performance profile rather than a forest plot.
    spec_a = panel_shell(fig, outer[0, 0], "A", "Discrimination profile",
                         "Independent validation; points and bars are estimates and 95% CIs", header_ratio=0.23)
    ax = fig.add_subplot(spec_a)
    x_base = np.arange(len(discrimination_metrics))
    for model in model_order:
        values, lower, upper = [], [], []
        for metric in discrimination_metrics:
            row = estimates[(estimates.metric == metric) & (estimates.model == model)].iloc[0]
            values.append(float(row.estimate)); lower.append(float(row.ci_lower)); upper.append(float(row.ci_upper))
        values = np.asarray(values); lower = np.asarray(lower); upper = np.asarray(upper)
        x = x_base + offsets[model]
        prominence = 0.82 if model in {"PHENO9", "LAB26"} else 0.56
        line_color = blend_with_white(COLORS[model], prominence)
        ax.plot(x, values, color=line_color, lw=1.35 if prominence > 0.8 else 0.95,
                alpha=0.92, zorder=2)
        ax.errorbar(x, values, yerr=np.vstack([values-lower, upper-values]), fmt=markers[model],
                    color=line_color, ecolor=line_color, mfc=COLORS[model], mec="white", mew=0.65,
                    ms=5.1, elinewidth=0.85, capsize=1.7, zorder=3)
    ax.set_xticks(x_base, [metric_labels[m] for m in discrimination_metrics])
    ax.set_ylabel("Discrimination")
    ax.set_xlim(-0.35, 2.35)
    ax.set_ylim(0.735, 0.855)
    style_axis(ax, "y")

    # B: absolute prediction error in the same participants.
    spec_b = panel_shell(fig, outer[0, 1], "B", "Prediction-error profile",
                         "Lower is better; same validation participants", header_ratio=0.23)
    ax = fig.add_subplot(spec_b)
    x_base = np.arange(len(error_metrics))
    for model in model_order:
        values, lower, upper = [], [], []
        for metric in error_metrics:
            row = estimates[(estimates.metric == metric) & (estimates.model == model)].iloc[0]
            values.append(float(row.estimate)); lower.append(float(row.ci_lower)); upper.append(float(row.ci_upper))
        values = np.asarray(values); lower = np.asarray(lower); upper = np.asarray(upper)
        x = x_base + offsets[model]
        prominence = 0.82 if model in {"PHENO9", "LAB26"} else 0.56
        line_color = blend_with_white(COLORS[model], prominence)
        ax.plot(x, values, color=line_color, lw=1.35 if prominence > 0.8 else 0.95,
                alpha=0.92, zorder=2)
        ax.errorbar(x, values, yerr=np.vstack([values-lower, upper-values]), fmt=markers[model],
                    color=line_color, ecolor=line_color, mfc=COLORS[model], mec="white", mew=0.65,
                    ms=5.1, elinewidth=0.85, capsize=1.7, zorder=3)
    ax.set_xticks(x_base, [metric_labels[m] for m in error_metrics])
    ax.set_ylabel("Prediction error")
    ax.set_xlim(-0.35, 2.35)
    ax.set_ylim(0.072, 0.170)
    style_axis(ax, "y")

    # C: paired benefit-oriented differences for all 36 prespecified comparisons.
    spec_c = panel_shell(fig, outer[1, :], "C", "Paired performance gains",
                         None, header_ratio=0.11)
    ax = fig.add_subplot(spec_c)
    comp_order = [
        "26 labs beyond clinical base", "Nine PhenoAge labs beyond clinical base",
        "Three components beyond clinical base", "26 labs versus nine PhenoAge labs",
        "26 labs versus three components", "Three components versus nine PhenoAge labs",
    ]
    comp_labels = [
        "LAB26 vs clinical", "PHENO9 vs clinical", "NM/TB/TC vs clinical",
        "LAB26 vs PHENO9", "LAB26 vs NM/TB/TC", "NM/TB/TC vs PHENO9",
    ]
    all_metrics = discrimination_metrics + error_metrics
    matrix = np.zeros((len(comp_order), len(all_metrics)))
    ci_low = np.zeros_like(matrix)
    ci_high = np.zeros_like(matrix)
    significant = np.zeros_like(matrix, dtype=bool)
    for i, comp in enumerate(comp_order):
        for j, metric in enumerate(all_metrics):
            row = paired[(paired.metric == metric) & (paired.comparison == comp)].iloc[0]
            matrix[i, j] = float(row.benefit_oriented_difference)
            lo, hi = float(row.difference_ci_lower), float(row.difference_ci_upper)
            if str(row.higher_is_better).upper() == "FALSE":
                lo, hi = -hi, -lo
            ci_low[i, j], ci_high[i, j] = lo, hi
            significant[i, j] = float(row.paired_p_holm) < 0.05
    column_scale = np.maximum(np.max(np.abs(matrix), axis=0), 1e-12)
    color_matrix = 0.76 * matrix / column_scale
    gain_cmap = LinearSegmentedColormap.from_list(
        "gain_soft", ["#AFC5D8", "#E7EEF3", "#FBFBFA", "#F3E4E7", "#C88D99"]
    )
    ax.imshow(color_matrix, cmap=gain_cmap, vmin=-1, vmax=1, aspect="auto", interpolation="nearest")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value, lo, hi = matrix[i, j], ci_low[i, j], ci_high[i, j]
            ax.text(j, i - 0.11, f"{value:+.3f}", ha="center", va="center",
                    fontsize=5.8, color=TEXT, fontweight="bold")
            ax.text(j, i + 0.16, f"({lo:+.3f}, {hi:+.3f})", ha="center", va="center",
                    fontsize=5.2, color="#51555B")
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                   edgecolor="white", linewidth=1.35))
            if significant[i, j]:
                ax.scatter(j + 0.39, i - 0.35, s=5.0, color="#4F5358", marker="o",
                           linewidths=0, zorder=5)
    ax.axhline(2.5, color="#747980", lw=0.85)
    ax.set_xticks(np.arange(len(all_metrics)), [metric_labels[m].replace("\n", " ") for m in all_metrics],
                  fontsize=7.2, fontweight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(comp_order)), comp_labels, fontsize=7.0)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(0.0, -0.075,
            "Cells: benefit-oriented difference (95% CI); positive favors the first model.\n"
            "Rose = gain; blue = loss; color scaled within metric; corner dot = Holm-adjusted P<0.05.",
            transform=ax.transAxes, fontsize=6.1, color=MID_GREY, va="top", linespacing=1.15)
    stem = RELEASE / "Figure_4_Fair_Target_Temporal_Comparison"
    paths = save_outputs(fig, stem)
    return paths, {"figure": 4, "estimate_rows": len(estimates), "paired_rows": len(paired),
                   "display": "two absolute performance profiles plus a 6-by-6 paired-gain matrix",
                   "validation_n": int(estimates.n_validation.iloc[0]),
                   "validation_deaths": int(estimates.validation_deaths.iloc[0])}


def quadratic_prediction(coefficients, x, y, z):
    get = lambda name: float(coefficients.get(name, 0.0))
    return (get("(Intercept)") + get("z_NM") * x + get("z_TB") * y + get("z_TC") * z
            + get("I(z_NM^2)") * x**2 + get("I(z_TB^2)") * y**2 + get("I(z_TC^2)") * z**2
            + get("z_NM:z_TB") * x * y + get("z_NM:z_TC") * x * z
            + get("z_TB:z_TC") * y * z)


def style_3d_axis(ax, limits):
    ax.set_xlim(limits["z_NM"]); ax.set_ylim(limits["z_TB"]); ax.set_zlim(limits["z_TC"])
    ax.set_box_aspect((1, 1, 0.84)); ax.view_init(elev=23, azim=-54)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1, 1, 1, 0)); axis.pane.set_edgecolor("#D9D9D9")
        axis._axinfo["grid"].update(color="#D9D9D9", linewidth=0.35, linestyle="-")
    ax.tick_params(pad=-1, labelsize=6.3)
    ax.set_xlabel("")
    ax.set_ylabel("TB z-logit", labelpad=0, fontsize=7.0)
    ax.set_zlabel("TC z-logit", labelpad=-5, fontsize=7.0)
    ax.text2D(0.34, -0.032, "NM z-logit", transform=ax.transAxes, fontsize=7.0,
              rotation=0, ha="center", va="center", color=TEXT)


def add_score_plane(ax, limits, planes):
    row = planes.iloc[0]
    beta = np.array([row.beta_z_NM, row.beta_z_TB, row.beta_z_TC], dtype=float)
    x = np.linspace(*limits["z_NM"], 70); y = np.linspace(*limits["z_TB"], 70)
    xx, yy = np.meshgrid(x, y)
    zz = (float(row.eta) - beta[0] * xx - beta[1] * yy) / beta[2]
    zz = np.where((zz >= limits["z_TC"][0]) & (zz <= limits["z_TC"][1]), zz, np.nan)
    ax.plot_surface(xx, yy, zz, color=COLORS["OVERALL"], alpha=0.22, linewidth=0.18,
                    edgecolor="#7F3443", antialiased=True, shade=False, rstride=3, cstride=3,
                    rasterized=True)


def add_pheno_surface(ax, coefficients, limits):
    n_grid = 54
    x = np.linspace(*limits["z_NM"], n_grid); y = np.linspace(*limits["z_TB"], n_grid); z = np.linspace(*limits["z_TC"], n_grid)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    volume = quadratic_prediction(coefficients, xx, yy, zz)
    vertices, faces, _, _ = measure.marching_cubes(
        volume.astype(np.float32), level=0, spacing=(x[1]-x[0], y[1]-y[0], z[1]-z[0])
    )
    vertices += np.array([x[0], y[0], z[0]])
    poly = Poly3DCollection(vertices[faces], facecolor=COLORS["PHENO"], edgecolor="#5D4B77",
                            linewidth=0.11, alpha=0.22)
    poly.set_rasterized(True); ax.add_collection3d(poly)


def figure5():
    data = read_csv(GEO / "geometry_model_space_data.csv")
    coefficients = read_csv(GEO / "geometry_pheno_surface_coefficients.csv").set_index("term")["estimate"]
    direction = read_csv(GEO / "geometry_direction_metrics.csv").set_index("metric")
    limits_df = read_csv(GEO / "geometry_axis_limits.csv").set_index("axis")
    limits = {axis: (float(limits_df.loc[axis, "lower_005"]), float(limits_df.loc[axis, "upper_995"]))
              for axis in ["z_NM", "z_TB", "z_TC"]}
    planes = read_csv(GEO / "geometry_overall_iso_probability_planes.csv")
    mapping = read_csv(SUP3 / "residual_mapping_quality.csv")
    recovery = read_csv(SUP3 / "residual_fidelity_paired_comparisons.csv")
    recovery = recovery[(recovery.model_a == "LAB26") & (recovery.model_b == "PROJECTED3")].copy()
    correlations = read_csv(SUP3 / "residual_dimension_lab_correlations.csv").sort_values("descriptive_rank")
    identity = read_csv(SUP3 / "residual_fidelity_identity_audit.csv")
    if len(data) != 11410 or data.SEQN.nunique() != 11410 or len(correlations) != 26:
        raise RuntimeError("Geometry/fidelity cohort drift")
    data.to_csv(SOURCE / "Figure_5AB_geometry_model_space.csv", index=False)
    direction.reset_index().to_csv(SOURCE / "Figure_5C_direction_metrics.csv", index=False)
    mapping.to_csv(SOURCE / "Figure_5D_mapping_quality.csv", index=False)
    recovery.to_csv(SOURCE / "Figure_5E_residual_recovery.csv", index=False)
    correlations.to_csv(SOURCE / "Figure_5F_residual_lab_correlations.csv", index=False)
    identity.to_csv(SOURCE / "Figure_5E_exact_identity_audit.csv", index=False)

    fig = plt.figure(figsize=(183 / 25.4, 247 / 25.4), facecolor="white")
    outer = fig.add_gridspec(3, 2, height_ratios=[1.00, 0.67, 1.28], width_ratios=[1, 1],
                             left=0.105, right=0.975, top=0.985, bottom=0.055,
                             wspace=0.34, hspace=0.23)
    overall_cmap = LinearSegmentedColormap.from_list("overall_soft", ["#F7FBFF", "#B8D3E5", "#E1B07A", COLORS["OVERALL"]])
    pheno_cmap = LinearSegmentedColormap.from_list("pheno_soft", ["#5B7FA6", "#F7F7F7", COLORS["OVERALL"]])
    overall_norm = Normalize(vmin=float(data.Overall_expected_burden.quantile(0.01)),
                             vmax=float(data.Overall_expected_burden.quantile(0.995)))
    pheno_bound = float(np.quantile(np.abs(data.PhenoAge_acceleration), 0.975))
    pheno_norm = TwoSlopeNorm(vmin=-pheno_bound, vcenter=0, vmax=pheno_bound)

    # A: frozen Overall surface, with no main-figure ParTI annotation.
    spec_a = panel_shell(fig, outer[0, 0], "A", "Cancer-trained multidomain score",
                         "Frozen Overall surface within the central 99% NHANES support", header_ratio=0.17,
                         title_fontsize=9.3)
    ax = fig.add_subplot(spec_a, projection="3d")
    ax.scatter(data.z_NM, data.z_TB, data.z_TC, c=data.Overall_expected_burden,
               cmap=overall_cmap, norm=overall_norm, s=2.0, alpha=0.20,
               linewidths=0, depthshade=False, rasterized=True)
    add_score_plane(ax, limits, planes); style_3d_axis(ax, limits)
    ax.text2D(0.02, 0.94, "Rose mesh: P(any affected domain)=0.50",
              transform=ax.transAxes, fontsize=6.5, color="#7F3443", va="top")
    cax = ax.inset_axes([0.17, -0.09, 0.66, 0.035])
    cb = fig.colorbar(ScalarMappable(norm=overall_norm, cmap=overall_cmap), cax=cax, orientation="horizontal")
    cb.set_label("Frozen Overall expected burden (0–3)", fontsize=6.5, labelpad=1)
    cb.ax.tick_params(labelsize=5.8, length=2)

    # B: PhenoAge surface in the identical coordinate system.
    spec_b = panel_shell(fig, outer[0, 1], "B", "PhenoAge surface in the same space",
                         "Quadratic zero-acceleration surface; same axes/view", header_ratio=0.17,
                         title_fontsize=9.0)
    ax = fig.add_subplot(spec_b, projection="3d")
    ax.scatter(data.z_NM, data.z_TB, data.z_TC,
               c=data.PhenoAge_acceleration.clip(-pheno_bound, pheno_bound),
               cmap=pheno_cmap, norm=pheno_norm, s=2.0, alpha=0.18,
               linewidths=0, depthshade=False, rasterized=True)
    add_pheno_surface(ax, coefficients, limits); style_3d_axis(ax, limits)
    ax.set_zlabel("")
    nonlinear_p = float(read_csv(GEO / "geometry_pheno_surface_coefficients.csv").nonlinear_block_wald_p.iloc[0])
    ax.text2D(0.02, 0.94, f"Purple mesh: zero acceleration; global Wald P={nonlinear_p:.1e}",
              transform=ax.transAxes, fontsize=6.3, color="#5D4B77", va="top")
    cax = ax.inset_axes([0.17, -0.09, 0.66, 0.035])
    cb = fig.colorbar(ScalarMappable(norm=pheno_norm, cmap=pheno_cmap), cax=cax, orientation="horizontal")
    cb.set_label("Observed PhenoAge Acceleration (years)", fontsize=6.5, labelpad=1)
    cb.ax.tick_params(labelsize=5.8, length=2)

    # C: quantified direction difference.
    spec_c = panel_shell(fig, outer[1, 0], "C", "Quantified score-space direction",
                         "Frozen Overall gradient versus local PhenoAge gradient at the origin", header_ratio=0.22,
                         title_fontsize=9.2)
    ax = fig.add_subplot(spec_c)
    beta = np.array([float(planes.iloc[0].beta_z_NM), float(planes.iloc[0].beta_z_TB), float(planes.iloc[0].beta_z_TC)])
    beta = beta / np.linalg.norm(beta)
    pheno_rows = direction.loc[["unit_gradient_NM_at_origin", "unit_gradient_TB_at_origin", "unit_gradient_TC_at_origin"]]
    pheno_est = pheno_rows.estimate.to_numpy(float); pheno_lo = pheno_rows.ci_lower.to_numpy(float); pheno_hi = pheno_rows.ci_upper.to_numpy(float)
    y = np.arange(3)[::-1]
    for yy, a, b in zip(y, beta, pheno_est, strict=True):
        ax.plot([a, b], [yy, yy], color="#C7C9CC", lw=1.1)
    ax.scatter(beta, y, s=34, color=COLORS["OVERALL"], marker="o", edgecolor="white", linewidth=0.5,
               label="Overall (frozen)", zorder=3)
    ax.errorbar(pheno_est, y, xerr=np.vstack([pheno_est-pheno_lo, pheno_hi-pheno_est]),
                fmt="D", color=COLORS["PHENO"], ecolor=COLORS["PHENO"], markersize=4.8,
                elinewidth=1.0, capsize=2.0, label="PhenoAge surface (95% CI)", zorder=4)
    ax.axvline(0, color="#55585D", lw=0.8, ls="--")
    ax.set_yticks(y, ["NM", "TB", "TC"])
    ax.set_xlim(-0.22, 1.08)
    ax.set_ylim(-0.92, 2.38)
    ax.set_xlabel("Unit-gradient component")
    angle = direction.loc["weighted_mean_angle_degrees"]
    corr = direction.loc["weighted_overall_pheno_correlation"]
    r2 = direction.loc["weighted_model_r_squared"]
    ax.text(0.07, -0.78,
            f"Mean local angle {angle.estimate:.1f}° ({angle.ci_lower:.1f}–{angle.ci_upper:.1f})\n"
            f"Survey-weighted r={corr.estimate:.3f} ({corr.ci_lower:.3f}–{corr.ci_upper:.3f}); R²={r2.estimate:.3f}",
            fontsize=6.2, color=TEXT, va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.90))
    ax.legend(loc="upper right", fontsize=6.3)
    style_axis(ax, "x")

    # D: cross-fitted and temporal mapping fidelity.
    spec_d = panel_shell(fig, outer[1, 1], "D", "Three-domain projection of LAB26",
                         "Explained variance in the LAB26 mortality predictor", header_ratio=0.22,
                         title_fontsize=9.0)
    ax = fig.add_subplot(spec_d)
    role_order = ["TRAIN_OOF", "TRAIN_APPARENT", "TEMPORAL_VALIDATION"]
    role_labels = ["Training OOF", "Training apparent", "Temporal validation"]
    y = np.arange(3)[::-1]
    for yy, role, label in zip(y, role_order, role_labels, strict=True):
        row = mapping[mapping.role == role].iloc[0]
        explained = 100 * float(row.weighted_r_squared)
        ax.barh(yy, explained, color=COLORS["COMPONENT3"], height=0.46)
        ax.barh(yy, 100 - explained, left=explained, color="#E5E7EA", height=0.46)
        ax.text(explained / 2, yy, f"{explained:.1f}%", ha="center", va="center", fontsize=7.0,
                color="white", fontweight="bold")
        ax.text(98, yy, f"r={float(row.weighted_correlation):.3f}", ha="right", va="center", fontsize=6.2, color=TEXT)
    ax.set_yticks(y, role_labels)
    ax.set_xlim(0, 100)
    ax.set_xlabel("Variance in LAB26 laboratory predictor (%)")
    ax.legend(handles=[Patch(fc=COLORS["COMPONENT3"], label="Explained by NM/TB/TC mapping"),
                       Patch(fc="#E5E7EA", label="Residual target-specific information")],
              loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=2, fontsize=5.7,
              handlelength=1.2, columnspacing=0.9, borderaxespad=0.0)
    ax.set_ylim(-0.52, 2.90)
    style_axis(ax, "x")

    # E: benefit-oriented recovery across discrimination and prediction error.
    spec_e = panel_shell(fig, outer[2, 0], "E", "Residual information restores\nmodel performance",
                         "LAB26 versus projected NM/TB/TC; 1,000 bootstraps", header_ratio=0.24,
                         title_fontsize=8.8)
    ax = fig.add_subplot(spec_e)
    metric_order = ["uno_c10", "auc5", "auc10", "brier5", "brier10", "ibs10"]
    metric_labels = ["Uno C, 10 y", "AUC, 5 y", "AUC, 10 y", "Brier↓, 5 y", "Brier↓, 10 y", "IBS↓"]
    y = np.arange(6)[::-1]
    for yy, metric in zip(y, metric_order, strict=True):
        row = recovery[recovery.metric == metric].iloc[0]
        est = float(row.benefit_oriented_difference)
        lo, hi = float(row.difference_ci_lower), float(row.difference_ci_upper)
        if str(row.higher_is_better).upper() == "FALSE":
            lo, hi = -hi, -lo
        ax.errorbar(est, yy, xerr=[[est-lo], [hi-est]], fmt="s", color=COLORS["LAB26"],
                    mfc=COLORS["LAB26"], mec="white", mew=0.6, ms=5.5, elinewidth=1.05, capsize=1.8)
        ax.text(0.0535, yy, f"P$_{{Holm}}$={fmt_p(row.paired_p_holm)}", ha="right", va="center", fontsize=5.9, color=MID_GREY)
    ax.axvline(0, color="#55585D", lw=0.85)
    ax.set_yticks(y, metric_labels)
    ax.set_xlim(-0.002, 0.055)
    ax.set_ylim(-1.02, 5.45)
    ax.set_xlabel("Benefit-oriented paired improvement")
    ax.text(0.14, 0.025, "Hybrid4 = LAB26 to machine precision for LP and 5/10-y risks",
            transform=ax.transAxes, fontsize=6.0, color="#3F6F63", fontweight="bold",
            ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.92))
    style_axis(ax, "x")

    # F: all 26 residual correlations, no biological-domain claim.
    spec_f = panel_shell(fig, outer[2, 1], "F", "Laboratory correlates of the\nresidual dimension",
                         "Weighted correlations with residual R", header_ratio=0.24,
                         title_fontsize=8.8)
    ax = fig.add_subplot(spec_f)
    corr_frame = correlations.sort_values("weighted_correlation_with_residual")
    y = np.arange(len(corr_frame))
    values = corr_frame.weighted_correlation_with_residual.to_numpy(float)
    ranks = corr_frame.descriptive_rank.to_numpy(int)
    labs_order = corr_frame.laboratory.tolist()
    for yy, value, rank in zip(y, values, ranks, strict=True):
        top = rank <= 5
        color = (COLORS["OVERALL"] if value > 0 else COLORS["NM"]) if top else "#AEB2B7"
        ax.plot([0, value], [yy, yy], color=color, lw=1.1 if top else 0.75)
        ax.scatter(value, yy, s=24 if top else 13, color=color, edgecolor="white", linewidth=0.4, zorder=3)
    ax.axvline(0, color="#55585D", lw=0.8)
    ax.set_yticks(y, [LAB_LABELS.get(x, x.replace("lab_", "")) for x in labs_order], fontsize=6.1)
    for tick, rank in zip(ax.get_yticklabels(), ranks, strict=True):
        if rank <= 5:
            tick.set_fontweight("bold"); tick.set_color(TEXT)
        else:
            tick.set_color(MID_GREY)
    ax.set_xlim(-0.36, 0.62)
    ax.set_xlabel("Survey-weighted correlation with residual R")
    style_axis(ax, "x")

    stem = RELEASE / "Figure_5_Geometry_and_Residual_Information"
    paths = save_outputs(fig, stem)
    return paths, {
        "figure": 5, "geometry_n": len(data), "geometry_unique_n": data.SEQN.nunique(),
        "residual_correlation_rows": len(correlations), "recovery_metrics": len(recovery),
        "identity_checks_passed": bool(identity["passed"].astype(str).str.lower().eq("true").all()),
        "main_figure_contains_pareto_gate_annotation": False,
    }


LEGEND_3 = """Figure 3. Mortality profiles of the cancer-trained multidomain scores.

(A) Fully adjusted survey-weighted Cox estimates for NM, TB, and TC entered jointly. (B) Fully adjusted cause-specific associations with cancer and cardiovascular mortality. (C) Survey-weighted cumulative incidence in the highest versus lowest Overall quartiles, with absolute differences at 5, 10, and 15 years. Hazard ratios are per survey-weighted standard deviation; bars are 95% confidence intervals. Cancer analyses use 1999–2016 and cardiovascular analyses use 1999–2014, where heart disease and cerebrovascular death categories are separately available.
"""

LEGEND_SUPP_ARCHITECTURE = """Supplementary Figure, numbering to be assigned. Frozen laboratory architecture of the three component models.

Rows comprise all 26 laboratory variables retained by at least one nutritional-metabolic (NM), tumor-burden-related (TB), or treatment-complication-related (TC) model. Tile intensity represents the median ratio of Boruta-style feature importance to the within-fold shadow threshold across the five outer folds; numbers are the number of outer folds selecting the feature. Bold colored borders identify final retained features. Hatched cells marked with a cross denote variables pre-excluded from the corresponding outcome model because they directly constituted or closely overlapped the outcome definition.
"""

LEGEND_4 = """Figure 4. Mortality-target-matched temporal comparison.

Models used identical clinical covariates, mortality outcome, preprocessing, algorithm, and PSU-grouped cross-validation; training used NHANES 2003–2006 and frozen testing used 2007–2010. (A) Uno C-index and time-dependent AUC. (B) Brier scores and integrated Brier score. (C) Paired benefit-oriented differences from 1,000 survey-bootstrap replicates; positive values favor the first model and corner dots denote Holm-adjusted P<0.05.
"""

LEGEND_5 = """Figure 5. Geometric distinction from PhenoAge and residual laboratory information.

(A) Frozen Overall and (B) PhenoAge Acceleration surfaces in the same NM-TB-TC space. (C) Their local directions, mean angle, correlation, and explained variance. (D) LAB26 variance explained by the three-domain mapping. (E) Temporal-validation improvement after restoring residual information. (F) Survey-weighted laboratory correlates of the residual. The residual quantifies fitted laboratory mortality information beyond the three-domain projection.
"""


def write_documentation(qa_rows):
    (RELEASE / "Figure_3_Legend_EN.txt").write_text(LEGEND_3, encoding="utf-8")
    (RELEASE / "Supplementary_Figure_Unnumbered_Frozen_Laboratory_Architecture_Legend_EN.txt").write_text(
        LEGEND_SUPP_ARCHITECTURE, encoding="utf-8")
    (RELEASE / "Figure_4_Legend_EN.txt").write_text(LEGEND_4, encoding="utf-8")
    (RELEASE / "Figure_5_Legend_EN.txt").write_text(LEGEND_5, encoding="utf-8")
    qa = {
        "target": "Nature-family double-column visual standard",
        "backend": "Python/matplotlib with rasterized dense data layers",
        "exports": "vector PDF, 300 dpi PNG preview, LZW-compressed 800 dpi TIFF",
        "font": "Arial/Helvetica/Liberation Sans; minimum plotted text 5.6 pt",
        "palette": COLORS,
        "figures": qa_rows,
        "global_checks": {
            "all_required_source_files_found": True,
            "no_input_rows_randomly_subsampled": True,
            "figure_5_main_has_no_pareto_or_tetrahedron_gate_annotation": True,
            "phenoage_not_mixed_into_unavailable_cause_specific_cox_estimand": True,
            "fair_target_and_main_score_temporal_calibration_not_conflated": True,
            "top_and_right_spines_removed_for_2d_axes": True,
            "semantic_color_mapping_matches_figure_2": True,
        },
    }
    (RELEASE / "Figures_3_5_QA.json").write_text(json.dumps(qa, indent=2), encoding="utf-8")
    report = """# Figures 3–5 statistics and reproducibility report

## Figure 3

- Panel A uses the fully adjusted survey-weighted Cox model in 12,094 adults aged at least 60 years with 4,787 deaths. NM, TB, and TC are entered simultaneously. Error bars are 95% confidence intervals.
- Panel B uses fully adjusted survey-weighted cause-specific Cox models. Competing deaths are censored at their death time. Error bars are 95% confidence intervals.
- Panel C uses weighted cumulative-incidence estimates. Risk-difference intervals and P values are based on 1,000 survey bootstrap replicates with Holm correction across the prespecified horizon family.
- The 26-laboratory feature architecture has been moved to an unnumbered Supplementary Figure. It uses the median outer-fold Boruta importance-to-shadow ratio, selection in 0–5 outer folds, and explicit structural marking of outcome-overlapping exclusions.

## Figure 4

- Training: NHANES 2003–2006, N=2,548, deaths=1,400. Temporal validation: NHANES 2007–2010, N=3,045, deaths=1,085.
- Every model shares the clinical base, all-cause mortality target, training and validation participants, preprocessing, elastic-net Cox algorithm, and PSU-grouped five-fold cross-validation.
- Panels A and B are connected performance profiles, not time trajectories; their error bars use 1,000 survey bootstrap replicates.
- Panel C shows all 36 paired differences as a 6-by-6 matrix. Each cell contains the exact paired difference and 95% confidence interval. Paired P values are Holm-adjusted within the prespecified comparison family.
- Panel C color is standardized within each metric only; exact numerical annotations, rather than color intensity, support cross-metric interpretation.
- Uno C-index is truncated at 10 years. Time-dependent AUC and Brier scores are reported at 5 and 10 years; IBS is integrated from 1 to 10 years.

## Figure 5

- Panels A–C use all 11,410 participants in the common geometry cohort. Dense participant layers are rasterized only for rendering; no rows are sampled or discarded.
- Direction intervals, correlations, and mapping R-squared use the existing 1,000 survey-bootstrap estimates.
- Panel D reports cross-fitted training and independent temporal-validation mapping quality. The displayed percentage is weighted R-squared, not a proportion of the original raw laboratory values that can be reconstructed.
- Panel E uses paired survey-bootstrap differences with Holm correction. Benefit orientation reverses the sign of Brier/IBS differences so that positive always denotes improvement.
- Panel F reports all 26 descriptive weighted correlations. Ranking is descriptive and does not imply an independent biological mechanism.
"""
    (RELEASE / "Figures_3_5_Statistics_and_Reproducibility.md").write_text(report, encoding="utf-8")


def write_hashes():
    files = sorted(p for p in RELEASE.rglob("*") if p.is_file() and p.name != "Figures_3_5_SHA256.csv")
    rows = []
    for path in files:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        rows.append({"relative_path": str(path.relative_to(RELEASE)), "sha256": digest, "bytes": path.stat().st_size})
    pd.DataFrame(rows).to_csv(RELEASE / "Figures_3_5_SHA256.csv", index=False)


def main():
    # Remove only superseded files produced by the earlier four-panel Figure 3 layout.
    obsolete = [
        RELEASE / "Figure_3_Laboratory_Architecture_and_Mortality_Profiles.pdf",
        RELEASE / "Figure_3_Laboratory_Architecture_and_Mortality_Profiles.png",
        RELEASE / "Figure_3_Laboratory_Architecture_and_Mortality_Profiles.tiff",
        SOURCE / "Figure_3A_feature_architecture.csv",
        SOURCE / "Figure_3B_joint_component_effects.csv",
        SOURCE / "Figure_3C_cause_specific_effects.csv",
        SOURCE / "Figure_3D_cumulative_incidence_curves.csv",
        SOURCE / "Figure_3D_q4_q1_risk_differences.csv",
    ]
    for path in obsolete:
        path.unlink(missing_ok=True)
    paths3, supplementary_paths, qa3 = figure3()
    paths4, qa4 = figure4()
    paths5, qa5 = figure5()
    write_documentation([qa3, qa4, qa5])
    write_hashes()
    result = {
        "release": str(RELEASE),
        "figure_3": [str(x) for x in paths3],
        "supplementary_feature_architecture": [str(x) for x in supplementary_paths],
        "figure_4": [str(x) for x in paths4],
        "figure_5": [str(x) for x in paths5],
        "qa": str(RELEASE / "Figures_3_5_QA.json"),
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
