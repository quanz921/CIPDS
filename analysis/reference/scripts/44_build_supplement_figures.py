# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (a-d) Line/forest sensitivity panels -> LineTrend + BarComparison parameter inherit -> param inherit
# (a-d) Residual-fidelity diagnostic panels -> BarComparison parameter inherit -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.

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

from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "supplementary_materials_release_20260907" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

SCORE_ORDER = ["OVERALL", "NM", "TB", "TC", "PHENO"]
SCORE_LABEL = {
    "OVERALL": "Overall",
    "NM": "NM",
    "TB": "TB",
    "TC": "TC",
    "PHENO": "PhenoAge acceleration",
}
SCORE_COLOR = {
    "OVERALL": ACCENT_RED,
    "NM": CATEGORICAL[0],
    "TB": CATEGORICAL[2],
    "TC": CATEGORICAL[3],
    "PHENO": CATEGORICAL[4],
}


def panel_label(ax, label):
    ax.text(-0.13, 1.08, label, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left", color=BLACK)


def clean_axis(ax, grid_axis="x"):
    ax.spines["left"].set_color("#444444")
    ax.spines["bottom"].set_color("#444444")
    ax.grid(axis=grid_axis, color="#D9DEE5", linewidth=0.45, alpha=0.75)
    ax.set_axisbelow(True)


def build_time_reverse_age_figure():
    src = ROOT / "outputs" / "supplement_v3"
    piece = pd.read_csv(src / "time_varying_piecewise_hr.csv")
    curve = pd.read_csv(src / "time_varying_continuous_hr_curves.csv")
    landmark = pd.read_csv(src / "early_death_landmark_sensitivity.csv")
    ages = pd.read_csv(src / "age_time_scale_sensitivity.csv")
    logtests = pd.read_csv(src / "time_varying_logtime_tests.csv")

    fig = plt.figure(figsize=(183 / 25.4, 205 / 25.4), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.985, bottom=0.08, top=0.965,
                          wspace=0.32, hspace=0.40)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # A: interval-specific effects
    interval_order = ["0-5 years", "5-10 years", ">10 years"]
    offsets = {"0-5 years": 0.22, "5-10 years": 0.0, ">10 years": -0.22}
    markers = {"0-5 years": "o", "5-10 years": "s", ">10 years": "D"}
    ybase = {score: len(SCORE_ORDER) - 1 - i for i, score in enumerate(SCORE_ORDER)}
    for score in SCORE_ORDER:
        sub = piece[piece["score"] == score]
        for interval in interval_order:
            row = sub[sub["interval"] == interval].iloc[0]
            y = ybase[score] + offsets[interval]
            ax_a.errorbar(row["hazard_ratio"], y,
                          xerr=[[row["hazard_ratio"] - row["ci_lower"]],
                                [row["ci_upper"] - row["hazard_ratio"]]],
                          fmt=markers[interval], color=SCORE_COLOR[score],
                          markerfacecolor="white" if interval == "5-10 years" else SCORE_COLOR[score],
                          markersize=4.2, linewidth=1.0, capsize=1.8, zorder=3)
    ax_a.axvline(1, color="#666666", linestyle="--", linewidth=0.75)
    ax_a.set_yticks([ybase[s] for s in SCORE_ORDER], [SCORE_LABEL[s] for s in SCORE_ORDER])
    ax_a.set_xlabel("Hazard ratio per survey-weighted SD")
    ax_a.set_title("Interval-specific mortality associations", loc="left", fontweight="bold")
    ax_a.set_xlim(0.94, max(1.90, piece["ci_upper"].max() * 1.03))
    clean_axis(ax_a)
    handles = [Line2D([0], [0], marker=markers[i], color="#555555", lw=0,
                      markerfacecolor="white" if i == "5-10 years" else "#555555",
                      markersize=4.2, label=i) for i in interval_order]
    ax_a.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=3,
                handletextpad=0.4, columnspacing=0.8)
    panel_label(ax_a, "A")

    # B: continuous time-varying HR curves
    for score in SCORE_ORDER:
        sub = curve[curve["score"] == score].sort_values("follow_up_years")
        ax_b.plot(sub["follow_up_years"], sub["hazard_ratio_per_weighted_sd"],
                  color=SCORE_COLOR[score], linewidth=1.45, label=SCORE_LABEL[score])
        ax_b.fill_between(sub["follow_up_years"].to_numpy(), sub["ci_lower"].to_numpy(),
                          sub["ci_upper"].to_numpy(), color=SCORE_COLOR[score], alpha=0.09,
                          linewidth=0)
    ax_b.axhline(1, color="#666666", linestyle="--", linewidth=0.75)
    ax_b.set_xlim(0.5, 15)
    ax_b.set_xlabel("Follow-up time (years)")
    ax_b.set_ylabel("Time-varying hazard ratio per SD")
    ax_b.set_title("Continuous score-by-log-time models", loc="left", fontweight="bold")
    clean_axis(ax_b, grid_axis="both")
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=3,
                handlelength=1.8, columnspacing=0.8)
    for score in SCORE_ORDER:
        p = logtests.loc[logtests["score"] == score, "score_log_time_p_holm"].iloc[0]
        if p < 0.05:
            pass
    panel_label(ax_b, "B")

    # C: reverse-causation landmarks
    lm_order = [0, 1, 2]
    lm_label = {0: "Baseline", 1: "1-year landmark", 2: "2-year landmark"}
    lm_marker = {0: "o", 1: "s", 2: "D"}
    lm_offset = {0: 0.22, 1: 0.0, 2: -0.22}
    for score in SCORE_ORDER:
        sub = landmark[landmark["score"] == score]
        for lm in lm_order:
            row = sub[sub["landmark_years"] == lm].iloc[0]
            y = ybase[score] + lm_offset[lm]
            ax_c.errorbar(row["hazard_ratio"], y,
                          xerr=[[row["hazard_ratio"] - row["ci_lower"]],
                                [row["ci_upper"] - row["hazard_ratio"]]],
                          fmt=lm_marker[lm], color=SCORE_COLOR[score],
                          markerfacecolor="white" if lm == 1 else SCORE_COLOR[score],
                          markersize=4.2, linewidth=1.0, capsize=1.8)
    ax_c.axvline(1, color="#666666", linestyle="--", linewidth=0.75)
    ax_c.set_yticks([ybase[s] for s in SCORE_ORDER], [SCORE_LABEL[s] for s in SCORE_ORDER])
    ax_c.set_xlabel("Hazard ratio per survey-weighted SD")
    ax_c.set_title("Landmark analyses excluding early deaths", loc="left", fontweight="bold")
    ax_c.set_xlim(0.94, max(1.82, landmark["ci_upper"].max() * 1.03))
    clean_axis(ax_c)
    handles = [Line2D([0], [0], marker=lm_marker[lm], color="#555555", lw=0,
                      markerfacecolor="white" if lm == 1 else "#555555", markersize=4.2,
                      label=lm_label[lm]) for lm in lm_order]
    ax_c.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=3,
                handletextpad=0.4, columnspacing=0.8)
    panel_label(ax_c, "C")

    # D: attained-age time scale
    ts_order = ["FOLLOW_UP_TIME", "ATTAINED_AGE"]
    ts_label = {"FOLLOW_UP_TIME": "Follow-up time + baseline age", "ATTAINED_AGE": "Attained age"}
    ts_marker = {"FOLLOW_UP_TIME": "o", "ATTAINED_AGE": "D"}
    ts_offset = {"FOLLOW_UP_TIME": 0.13, "ATTAINED_AGE": -0.13}
    for score in SCORE_ORDER:
        sub = ages[ages["score"] == score]
        rows = []
        for ts in ts_order:
            row = sub[sub["time_scale"] == ts].iloc[0]
            rows.append(row)
            ax_d.errorbar(row["hazard_ratio"], ybase[score] + ts_offset[ts],
                          xerr=[[row["hazard_ratio"] - row["ci_lower"]],
                                [row["ci_upper"] - row["hazard_ratio"]]],
                          fmt=ts_marker[ts], color=SCORE_COLOR[score],
                          markerfacecolor="white" if ts == "ATTAINED_AGE" else SCORE_COLOR[score],
                          markersize=4.2, linewidth=1.0, capsize=1.8)
        ax_d.plot([rows[0]["hazard_ratio"], rows[1]["hazard_ratio"]],
                  [ybase[score] + ts_offset[ts_order[0]], ybase[score] + ts_offset[ts_order[1]]],
                  color="#BFC3C8", linewidth=0.65, zorder=0)
    ax_d.axvline(1, color="#666666", linestyle="--", linewidth=0.75)
    ax_d.set_yticks([ybase[s] for s in SCORE_ORDER], [SCORE_LABEL[s] for s in SCORE_ORDER])
    ax_d.set_xlabel("Hazard ratio per survey-weighted SD")
    ax_d.set_title("Sensitivity using attained age as the time scale", loc="left", fontweight="bold")
    ax_d.set_xlim(0.94, max(1.82, ages["ci_upper"].max() * 1.03))
    clean_axis(ax_d)
    handles = [Line2D([0], [0], marker=ts_marker[ts], color="#555555", lw=0,
                      markerfacecolor="white" if ts == "ATTAINED_AGE" else "#555555",
                      markersize=4.2, label=ts_label[ts]) for ts in ts_order]
    ax_d.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.16), ncol=2,
                handletextpad=0.4, columnspacing=0.8)
    panel_label(ax_d, "D")

    fig.text(0.10, 0.012,
             "Common cohort: N=9,319; 4,177 deaths. Error bars and shaded bands are 95% confidence intervals.",
             fontsize=7, color="#555555")
    stem = OUT / "Figure_S8_time_reverse_age_sensitivity"
    save_cns_figure(fig, stem)
    fig.savefig(f"{stem}_800dpi.png", bbox_inches="tight", dpi=800, facecolor="white")
    plt.close(fig)


def metric_name(code):
    return {
        "uno_c10": "Uno C, 10 y",
        "auc5": "AUC, 5 y",
        "auc10": "AUC, 10 y",
        "brier5": "Brier, 5 y",
        "brier10": "Brier, 10 y",
        "ibs10": "IBS, 1-10 y",
    }.get(code, code)


def build_residual_fidelity_figure():
    src = ROOT / "outputs" / "supplement_v3"
    mapping = pd.read_csv(src / "residual_mapping_quality.csv")
    paired = pd.read_csv(src / "residual_fidelity_paired_comparisons.csv")
    cor = pd.read_csv(src / "residual_dimension_lab_correlations.csv")
    ident = pd.read_csv(src / "residual_fidelity_identity_audit.csv")

    fig = plt.figure(figsize=(183 / 25.4, 195 / 25.4), facecolor="white")
    gs = fig.add_gridspec(2, 2, left=0.10, right=0.985, bottom=0.09, top=0.965,
                          wspace=0.36, hspace=0.40)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # A: variance accounted for by three-domain mapping
    role_label = {"TRAIN_OOF": "Training OOF", "TRAIN_APPARENT": "Training apparent",
                  "TEMPORAL_VALIDATION": "Temporal validation"}
    m = mapping.iloc[::-1].copy()
    y = np.arange(len(m))
    explained = 100 * m["weighted_r_squared"].to_numpy()
    ax_a.barh(y, explained, color=CATEGORICAL[0], alpha=0.78, height=0.48,
              label="Explained by NM/TB/TC mapping")
    ax_a.barh(y, 100 - explained, left=explained, color="#E1E4E8", height=0.48,
              label="Residual target-specific information")
    for yi, val, r in zip(y, explained, m["weighted_correlation"]):
        ax_a.text(val / 2, yi, f"{val:.1f}%", ha="center", va="center",
                  fontsize=7, color="white", fontweight="bold")
        ax_a.text(98, yi, f"r={r:.3f}", ha="right", va="center", fontsize=6.6, color=BLACK)
    ax_a.set_yticks(y, [role_label[x] for x in m["role"]])
    ax_a.set_xlim(0, 100)
    ax_a.set_xlabel("Variance in LAB26 laboratory predictor (%)")
    ax_a.set_title("Three-domain projection of LAB26", loc="left", fontweight="bold")
    clean_axis(ax_a)
    ax_a.legend(loc="upper left", bbox_to_anchor=(0, -0.17), ncol=1)
    panel_label(ax_a, "A")

    # B: benefit-oriented performance restored by residual information
    p = paired[(paired["model_a"] == "LAB26") & (paired["model_b"] == "PROJECTED3")].copy()
    order = ["uno_c10", "auc5", "auc10", "brier5", "brier10", "ibs10"]
    p["ord"] = p["metric"].map({v: i for i, v in enumerate(order)})
    p = p.sort_values("ord", ascending=False)
    yy = np.arange(len(p))
    scale = 100.0
    est = p["benefit_oriented_difference"].to_numpy() * scale
    raw_lo = p["difference_ci_lower"].to_numpy() * scale
    raw_hi = p["difference_ci_upper"].to_numpy() * scale
    higher = p["higher_is_better"].astype(bool).to_numpy()
    lo = np.where(higher, raw_lo, -raw_hi)
    hi = np.where(higher, raw_hi, -raw_lo)
    ax_b.errorbar(est, yy, xerr=[est - lo, hi - est], fmt="s", color=ACCENT_RED,
                  markersize=4.3, linewidth=1.05, capsize=1.8)
    ax_b.axvline(0, color="#555555", linewidth=0.75)
    ax_b.set_yticks(yy, [metric_name(x) for x in p["metric"]])
    ax_b.set_xlabel("Benefit-oriented paired improvement (percentage points)")
    ax_b.set_title("Residual information restores performance", loc="left", fontweight="bold")
    clean_axis(ax_b)
    panel_label(ax_b, "B")

    # C: coefficient versus residual correlation
    ax_c.axhline(0, color="#BFC3C8", linewidth=0.6)
    ax_c.axvline(0, color="#BFC3C8", linewidth=0.6)
    sizes = 18 + 55 * cor["absolute_weighted_correlation"].to_numpy()
    colors = np.where(cor["weighted_correlation_with_residual"] >= 0, CATEGORICAL[0], ACCENT_RED)
    ax_c.scatter(cor["glmnet_coefficient_in_lab26_model"],
                 cor["weighted_correlation_with_residual"], s=sizes, c=colors,
                 alpha=0.82, edgecolor="white", linewidth=0.35)
    top = cor.nlargest(7, "absolute_weighted_correlation")
    for _, row in top.iterrows():
        ax_c.annotate(row["laboratory"].replace("lab_", ""),
                      (row["glmnet_coefficient_in_lab26_model"], row["weighted_correlation_with_residual"]),
                      xytext=(4, 4), textcoords="offset points", fontsize=6.2, color=BLACK)
    ax_c.set_xlabel("Elastic-net coefficient in LAB26 model")
    ax_c.set_ylabel("Survey-weighted correlation with residual R")
    ax_c.set_title("Laboratory correlates of residual information", loc="left", fontweight="bold")
    clean_axis(ax_c, grid_axis="both")
    panel_label(ax_c, "C")

    # D: numerical identity relative to prespecified tolerance
    ident = ident.copy()
    ident["ratio"] = ident["maximum_absolute_error"] / ident["tolerance"]
    ident["log_ratio"] = np.log10(ident["ratio"].replace(0, np.nan))
    short = ["LAB26 LP decomposition", "Hybrid4 LP", "Hybrid4 5-y risk", "Hybrid4 10-y risk"]
    y = np.arange(len(ident))[::-1]
    vals = ident["log_ratio"].fillna(-16).to_numpy()
    ax_d.barh(y, vals, color=CATEGORICAL[2], alpha=0.78, height=0.46)
    ax_d.axvline(0, color="#555555", linestyle="--", linewidth=0.75)
    ax_d.set_yticks(y, short)
    ax_d.set_xlabel("log10(maximum absolute error / tolerance)")
    ax_d.set_title("Machine-precision reconstruction audit", loc="left", fontweight="bold")
    ax_d.set_xlim(min(-8.5, np.nanmin(vals) - 0.5), 0.5)
    clean_axis(ax_d)
    for yi, err in zip(y, ident["maximum_absolute_error"]):
        ax_d.text(-0.15, yi, f"{err:.2e}", ha="right", va="center", fontsize=6.2, color=BLACK)
    panel_label(ax_d, "D")

    fig.text(0.10, 0.018,
             "Training: NHANES 2003-2006; temporal validation: 2007-2010 (N=3,045; 1,085 deaths). "
             "Intervals use 1,000 survey-bootstrap replicates.", fontsize=7, color="#555555")
    stem = OUT / "Figure_S9_residual_fidelity_diagnostics"
    save_cns_figure(fig, stem)
    fig.savefig(f"{stem}_800dpi.png", bbox_inches="tight", dpi=800, facecolor="white")
    plt.close(fig)


def write_statistics_manifest():
    manifest = {
        "Figure_S8": {
            "scientific_claim": "Cancer-trained scores remain associated with mortality after excluding early deaths and using attained age, while several effects decline over follow-up.",
            "archetype": "quantitative_grid",
            "cohort": "NHANES age >=60 common complete cases, N=9319, deaths=4177",
            "intervals": "95% confidence intervals from sampling-weighted Cox models with PSU-cluster robust variance or the full survey design as specified",
            "multiplicity": "Holm correction across five scores within each analysis family",
            "sources": ["time_varying_piecewise_hr.csv", "time_varying_continuous_hr_curves.csv", "early_death_landmark_sensitivity.csv", "age_time_scale_sensitivity.csv"],
        },
        "Figure_S9": {
            "scientific_claim": "The three interpretable domains retain only part of the LAB26 mortality signal; a target-specific residual restores that fitted prediction exactly.",
            "archetype": "quantitative_grid",
            "training": "NHANES 2003-2006, N=2548, deaths=1400",
            "validation": "NHANES 2007-2010, N=3045, deaths=1085",
            "intervals": "95% confidence intervals from 1000 survey-bootstrap replicates",
            "nonclaim": "The residual is not lossless compression of raw laboratory measurements and is not a fourth biological domain.",
            "sources": ["residual_mapping_quality.csv", "residual_fidelity_paired_comparisons.csv", "residual_dimension_lab_correlations.csv", "residual_fidelity_identity_audit.csv"],
        },
    }
    (OUT.parent / "supplementary_figure_statistics_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    build_time_reverse_age_figure()
    build_residual_fidelity_figure()
    write_statistics_manifest()
    print(f"Wrote supplementary figures to {OUT}")
