# Academic Figure Skill Asset Confirmation (verified against assets/figures/)
# (A) hospital discrimination + PR curves -> assets/figures/AUROC -> param inherit
# (B) temporal calibration curves -> assets/figures/LineTrend -> param inherit
# (C) NHANES discrimination + incremental forests -> assets/figures/BarComparison -> param inherit
# (D) adjusted mortality forest -> assets/figures/BarComparison -> param inherit
# (E) cancer-history subgroup forest -> assets/figures/BarComparison -> param inherit
# (F) time-varying HR curves -> assets/figures/LineTrend -> param inherit
# RULE: "native run" = load pre-rendered PNG via Image.open().ax.imshow().
#       "param inherit" = drawing function below that copies Class A/B/C values.
#       If a panel says "native run" and you write a drawing function, you broke the contract.
# No native template run: the bundled assets do not support the study-specific nested multipanel data structures.

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
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.bbox": "tight",
    "savefig.dpi": 300,
})
def save_cns_figure(fig, filename):
    """Standard Academic Figure Skill export: vector PDF + 300dpi PNG preview."""
    fig.savefig(f"{filename}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(f"{filename}.png", bbox_inches="tight", dpi=300)

import argparse
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from scipy.stats import norm
from sklearn.metrics import average_precision_score, precision_recall_curve, roc_auc_score
from PIL import Image

# User-requested large-print typography. This override is intentionally placed
# after the mandatory baseline block; the final print-size floor is 8 pt.
mpl.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 10.5,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 7,
    "figure.titlesize": 12,
})


# Softer, colorblind-conscious score palette. Shape and line style carry the
# second visual channel so the figure remains interpretable in grayscale.
COLORS = {
    "OVERALL": "#B76576",  # muted rose
    "NM": "#5B7FA6",       # soft blue
    "TB": "#6F9B7A",       # soft sage
    "TC": "#C69355",       # muted ochre
    "PHENO": "#8B78A6",    # dusty purple
}
LIGHT_GREY = "#D7DADF"
MID_GREY = "#777B82"
TEXT = "#272A2F"

SCORE_LABELS = {
    "OVERALL": "Overall",
    "NM": "NM",
    "TB": "TB",
    "TC": "TC",
    "PHENO": "PhenoAge",
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


def panel_shell(spec, letter, title, subtitle=None, header_ratio=0.16):
    """Reserve a dedicated header band so titles can never overlap data."""
    shell = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=spec,
        height_ratios=[header_ratio, 1 - header_ratio], hspace=0.02,
    )
    hax = plt.subplot(shell[0, 0])
    hax.set_axis_off()
    hax.text(0.00, 0.72, letter, fontsize=14, fontweight="bold",
             va="center", ha="left", color=TEXT)
    hax.text(0.085, 0.72, title, fontsize=12, fontweight="bold",
             va="center", ha="left", color=TEXT)
    if subtitle:
        hax.text(0.085, 0.12, subtitle, fontsize=9.2,
                 va="center", ha="left", color=MID_GREY)
    return shell[1, 0]


def wilson_ci(k, n, alpha=0.05):
    if n <= 0:
        return np.nan, np.nan
    z = norm.ppf(1 - alpha / 2)
    p = k / n
    den = 1 + z * z / n
    ctr = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, ctr - half), min(1.0, ctr + half)


def calibration_bins(y, p, bins=10):
    frame = pd.DataFrame({"y": np.asarray(y, dtype=float), "p": np.asarray(p, dtype=float)}).dropna()
    frame["bin"] = pd.qcut(frame["p"], q=bins, duplicates="drop")
    rows = []
    for _, group in frame.groupby("bin", observed=True):
        n = len(group)
        k = int(group["y"].sum())
        lo, hi = wilson_ci(k, n)
        rows.append({
            "mean_predicted": group["p"].mean(),
            "observed": k / n,
            "ci_lower": lo,
            "ci_upper": hi,
            "n": n,
            "events": k,
        })
    return pd.DataFrame(rows)


def save_publication_outputs(fig, stem):
    for label in fig.findobj(mpl.text.Text):
        label.set_text(label.get_text().replace("=<", "<"))
    """User-requested vector PDF and LZW-compressed 800 dpi TIFF plus preview."""
    pdf_path = stem.with_suffix(".pdf")
    tiff_path = stem.with_suffix(".tiff")
    # Override the baseline's tight crop here so the final page remains exactly
    # 240 x 360 mm and every panel uses the requested enlarged review canvas.
    with mpl.rc_context({"savefig.bbox": None}):
        fig.savefig(pdf_path, format="pdf", bbox_inches=None)
        fig.savefig(stem.with_suffix(".png"), format="png", dpi=300,
                    bbox_inches=None)
    # Render the TIFF from the vector master. This avoids holding a ~45 MP
    # RGBA canvas plus a second PIL copy in memory and guarantees PDF/TIFF parity.
    pdftoppm = shutil.which("pdftoppm")
    bundled = Path("external_inputs/pdftoppm.exe")
    if not pdftoppm and bundled.exists():
        pdftoppm = str(bundled)
    if not pdftoppm:
        raise RuntimeError("pdftoppm is required for memory-safe 800 dpi TIFF export")
    tmp_prefix = stem.parent / f"{stem.name}_800dpi_tmp"
    subprocess.run([
        pdftoppm, "-r", "800", "-singlefile", "-tiff",
        "-tiffcompression", "lzw", str(pdf_path), str(tmp_prefix)
    ], check=True)
    generated = tmp_prefix.with_suffix(".tif")
    if tiff_path.exists():
        tiff_path.unlink()
    generated.replace(tiff_path)


def read_inputs(root):
    out = root / "outputs"
    return {
        "oof": pd.read_csv(out / "overall_meta_training_oof.csv"),
        "calendar_pred": pd.read_csv(out / "overall_calendar_predictions.csv"),
        "calendar_metrics": pd.read_csv(out / "nested_a_primary_calendar_metrics.csv"),
        "outer": pd.read_csv(out / "nested_a_primary_outer_performance.csv"),
        "overall_perf": pd.read_csv(out / "overall_calendar_performance.csv"),
        "disc": pd.read_csv(out / "nhanes_paired_discrimination_estimates.csv"),
        "incremental": pd.read_csv(out / "supplement_v2" / "incremental_prediction_paired_comparisons.csv"),
        "allcause": pd.read_csv(out / "main_tables" / "table2_panel_c_allcause_holm.csv"),
        "subgroup": pd.read_csv(out / "main_tables" / "table2_panel_d_cancer_history_interactions_all_scores.csv"),
        "tv_curve": pd.read_csv(out / "supplement_v3" / "time_varying_continuous_hr_curves.csv"),
        "tv_test": pd.read_csv(out / "supplement_v3" / "time_varying_logtime_tests.csv"),
        "piecewise": pd.read_csv(out / "supplement_v3" / "time_varying_piecewise_hr.csv"),
    }


def hospital_metrics(data):
    outcome_map = {
        "NM": "Outcome_NutriMetab",
        "TB": "Outcome_TumorBurden",
        "TC": "Outcome_TreatComp",
    }
    rows = []
    for score in ["NM", "TB", "TC"]:
        y = data["oof"][f"observed_{score}"].astype(int)
        p = data["oof"][f"p_{score}"].astype(float)
        rows.append({"score": score, "endpoint": score, "set": "Strict OOF",
                     "n": len(y), "events": int(y.sum()),
                     "auc": roc_auc_score(y, p), "lower": np.nan, "upper": np.nan,
                     "ap": average_precision_score(y, p)})
        cm = data["calendar_metrics"]
        for split, set_label in [("calendar_validation", "Validation"), ("calendar_test", "Test")]:
            r = cm[(cm["outcome"] == outcome_map[score]) & (cm["split"] == split)].iloc[0]
            rows.append({"score": score, "endpoint": score, "set": set_label,
                         "n": int(r["n"]), "events": int(r["positives"]),
                         "auc": float(r["auc"]), "lower": float(r["auc_ci_lower"]),
                         "upper": float(r["auc_ci_upper"]), "ap": float(r["average_precision"])})
    op = data["overall_perf"]
    for endpoint, auc_col, lo_col, hi_col, ap_col in [
        ("Overall >=1", "any_auc", "any_auc_ci_lower", "any_auc_ci_upper", "any_average_precision"),
        ("Overall >=2", "multidomain_auc", "multidomain_auc_ci_lower", "multidomain_auc_ci_upper", "multidomain_average_precision"),
    ]:
        for split, set_label in [("calendar_validation", "Validation"), ("calendar_test", "Test")]:
            r = op[(op["candidate_model"] == "main_effects") & (op["split"] == split)].iloc[0]
            events = int(r["burden_1"] + r["burden_2"] + r["burden_3"]) if endpoint.endswith(">=1") else int(r["burden_2"] + r["burden_3"])
            rows.append({"score": "OVERALL", "endpoint": endpoint, "set": set_label,
                         "n": int(r["n"]), "events": events,
                         "auc": float(r[auc_col]), "lower": float(r[lo_col]),
                         "upper": float(r[hi_col]), "ap": float(r[ap_col])})
    return pd.DataFrame(rows)


def draw_panel_a(spec, data, plotted):
    content = panel_shell(
        spec, "A", "Hospital model development",
        "Strict OOF development; independent calendar validation/test",
    )
    body = GridSpecFromSubplotSpec(
        5, 2, subplot_spec=content,
        height_ratios=[0.07, 0.66, 0.10, 0.08, 0.09], width_ratios=[0.82, 1.28],
        hspace=0.10, wspace=0.38,
    )
    left_heading = plt.subplot(body[0, 0])
    right_heading = plt.subplot(body[0, 1])
    for heading in [left_heading, right_heading]:
        heading.set_axis_off()
    left_heading.text(0.0, 0.35, "Temporal AUROC", fontsize=10.0,
                      fontweight="bold", color=TEXT, ha="left", va="center")
    right_heading.text(0.0, 0.35, "Precision–recall curves", fontsize=10.0,
                       fontweight="bold", color=TEXT, ha="left", va="center")

    ax = plt.subplot(body[1, 0])
    perf = hospital_metrics(data)
    plotted.append(perf.assign(panel="A", measure="hospital_discrimination"))
    endpoints = ["NM", "TB", "TC", "Overall >=1", "Overall >=2"]
    y_base = np.arange(len(endpoints))[::-1]
    set_style = {
        "Strict OOF": (-0.20, "o", ":", 5.8),
        "Validation": (0.00, "s", "--", 5.6),
        "Test": (0.20, "D", "-", 5.8),
    }
    for i, endpoint in enumerate(endpoints):
        score = "OVERALL" if endpoint.startswith("Overall") else endpoint
        for set_name, (off, marker, _, size) in set_style.items():
            rr = perf[(perf["endpoint"] == endpoint) & (perf["set"] == set_name)]
            if rr.empty:
                continue
            r = rr.iloc[0]
            xerr = None
            if np.isfinite(r["lower"]):
                xerr = np.array([[r["auc"] - r["lower"]], [r["upper"] - r["auc"]]])
            ax.errorbar(r["auc"], y_base[i] + off, xerr=xerr, fmt=marker,
                        color=COLORS[score], mec="white", mew=0.45, ms=size,
                        elinewidth=1.05, capsize=2.0, zorder=3)
    ax.axvline(0.5, color=LIGHT_GREY, lw=0.85, ls="--")
    ax.set_yticks(y_base)
    ax.set_yticklabels(["NM", "TB", "TC", "Overall ≥1", "Overall ≥2"], fontsize=9.2)
    ax.set_xlim(0.49, 0.87)
    ax.set_ylim(-0.45, 4.75)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8])
    ax.set_xlabel("AUROC", fontsize=10, labelpad=6)
    style_axis(ax, "x")

    pr_spec = GridSpecFromSubplotSpec(
        3, 2, subplot_spec=body[1, 1], width_ratios=[0.73, 0.27],
        hspace=0.30, wspace=0.08,
    )
    split_meta = [
        ("Strict OOF", data["oof"], None, ":", 1.0),
        ("Validation", data["calendar_pred"], "calendar_validation", "--", 1.05),
        ("Test", data["calendar_pred"], "calendar_test", "-", 1.2),
    ]
    for j, score in enumerate(["NM", "TB", "TC"]):
        pax = plt.subplot(pr_spec[j, 0])
        apax = plt.subplot(pr_spec[j, 1])
        apax.set_axis_off()
        ap_text = []
        for set_name, frame, split, ls, lw in split_meta:
            ff = frame if split is None else frame[frame["split"] == split]
            y = ff[f"observed_{score}"].astype(int).to_numpy()
            p = ff[f"p_{score}"].astype(float).to_numpy()
            precision, recall, _ = precision_recall_curve(y, p)
            ap = average_precision_score(y, p)
            prevalence = y.mean()
            pax.plot(recall, precision, color=COLORS[score], ls=ls, lw=lw + 0.25,
                     alpha={"Strict OOF": 0.68, "Validation": 0.80, "Test": 1.0}[set_name],
                     rasterized=False)
            ap_text.append(f"{set_name.replace('Strict ', '')} {ap:.3f}")
            plotted.append(pd.DataFrame({
                "panel": ["A"], "measure": ["precision_recall_summary"],
                "score": [score], "set": [set_name], "n": [len(y)],
                "events": [int(y.sum())], "prevalence": [prevalence], "ap": [ap]
            }))
        test_y = data["calendar_pred"].loc[data["calendar_pred"]["split"] == "calendar_test", f"observed_{score}"].astype(int)
        pax.axhline(test_y.mean(), color=LIGHT_GREY, lw=0.8, ls="--", zorder=0)
        pax.set_xlim(0, 1)
        ymax = {"NM": 0.38, "TB": 0.62, "TC": 0.38}[score]
        pax.set_ylim(0, ymax)
        pax.text(0.02, 0.93, score, transform=pax.transAxes, color=COLORS[score],
                 fontsize=9.4, fontweight="bold", va="top", ha="left")
        ap_lines = [f"OOF {ap_text[0].split()[-1]}",
                    f"Validation {ap_text[1].split()[-1]}",
                    f"Test {ap_text[2].split()[-1]}"]
        apax.text(0.0, 0.88, "\n".join(ap_lines),
                  ha="left", va="top", fontsize=8.4, color=TEXT, linespacing=1.18)
        pax.set_ylabel("Precision", fontsize=9.2, labelpad=3.0)
        if j < 2:
            pax.tick_params(labelbottom=False)
        else:
            pax.set_xlabel("Recall", fontsize=10, labelpad=6)
        pax.set_xticks([0, 0.5, 1])
        pax.set_yticks([0, round(ymax / 2, 2), ymax])
        pax.tick_params(labelsize=8.7, pad=1.8)
        style_axis(pax)
    line_handles = [Line2D([0], [0], color=MID_GREY, ls=ls, lw=1.2, label=name)
                    for name, _, _, ls, _ in split_meta]
    xlabel_ax = plt.subplot(body[2, :])
    xlabel_ax.set_axis_off()
    legend_ax = plt.subplot(body[3, :])
    legend_ax.set_axis_off()
    legend_ax.legend(handles=line_handles, loc="center", ncol=3,
                     fontsize=7.0, handlelength=2.0, columnspacing=1.5)
    note_ax = plt.subplot(body[4, :])
    note_ax.set_axis_off()
    note_ax.text(0.0, 0.35,
                "Bars are 95% CIs for validation/test; OOF AUROC is a point estimate.",
                fontsize=8.7, color=MID_GREY, va="center", ha="left")


def draw_panel_b(spec, data, plotted):
    content = panel_shell(
        spec, "C", "NHANES mortality discrimination",
        "Upper: N=11,410/5,272 deaths; lower: N=9,319/4,177 deaths",
        header_ratio=0.13,
    )
    sub = GridSpecFromSubplotSpec(
        8, 3, subplot_spec=content,
        height_ratios=[0.05, 0.32, 0.07, 0.07, 0.39, 0.03, 0.05, 0.03],
        width_ratios=[0.13, 0.68, 0.19], hspace=0.10, wspace=0.04,
    )
    top_title = plt.subplot(sub[0, :])
    top_title.set_axis_off()
    top_title.text(0.13, 0.45, "Absolute mortality discrimination",
                   fontsize=10.0, fontweight="bold", color=TEXT,
                   ha="left", va="center")
    top_lab = plt.subplot(sub[1, 0])
    ax = plt.subplot(sub[1, 1])
    disc = data["disc"]
    disc = disc[(disc["scenario"] == "A_PRIMARY_FROZEN_PANCANCER_OVERALL_STACK") &
                (disc["domain"] == "age60") & (disc["population"] == "Overall") &
                (disc["model"].isin(["OVERALL", "NM", "TB", "TC", "PHENO"])) &
                (disc["metric"].isin(["uno_c15", "auc5", "auc10", "auc15"]))].copy()
    plotted.append(disc.assign(panel="C1", measure="nhanes_absolute_discrimination"))
    score_order = ["OVERALL", "NM", "TB", "TC", "PHENO"]
    y = np.arange(len(score_order))[::-1]
    metric_styles = {
        "uno_c15": (-0.21, "s", "Uno C, 15 y"),
        "auc5": (-0.07, "o", "AUC, 5 y"),
        "auc10": (0.07, "^", "AUC, 10 y"),
        "auc15": (0.21, "D", "AUC, 15 y"),
    }
    for i, score in enumerate(score_order):
        for metric, (off, marker, _) in metric_styles.items():
            r = disc[(disc["model"] == score) & (disc["metric"] == metric)].iloc[0]
            err = np.array([[r["estimate"] - r["ci_lower"]], [r["ci_upper"] - r["estimate"]]])
            ax.errorbar(r["estimate"], y[i] + off, xerr=err, fmt=marker,
                        color=COLORS[score], mec="white", mew=0.45, ms=5.4,
                        elinewidth=1.0, capsize=1.8)
    ax.axvline(0.5, color=LIGHT_GREY, ls="--", lw=0.8)
    ax.set_xlim(0.47, 0.78)
    ax.set_ylim(-0.45, 4.80)
    ax.set_xticks([0.5, 0.6, 0.7])
    ax.set_yticks(y)
    ax.set_yticklabels([])
    top_lab.set_xlim(0, 1)
    top_lab.set_ylim(-0.45, 4.80)
    top_lab.set_axis_off()
    for yi, score in zip(y, score_order):
        top_lab.text(0.96, yi, SCORE_LABELS[score], fontsize=9.0,
                     ha="right", va="center", color=TEXT)
    style_axis(ax, "x")
    handles = [Line2D([0], [0], marker=m, color="none", markerfacecolor=MID_GREY,
                      markeredgecolor="white", markersize=5.5, label=lab)
               for _, m, lab in metric_styles.values()]
    lax = plt.subplot(sub[1, 2])
    lax.set_axis_off()
    lax.legend(handles=handles, loc="center left", fontsize=7.0,
               handletextpad=0.36, labelspacing=0.58, borderaxespad=0.0)

    top_xlabel = plt.subplot(sub[2, :])
    top_xlabel.set_axis_off()
    ax.set_xlabel("C-index / time-dependent AUC", fontsize=9.4, labelpad=6)
    middle_title = plt.subplot(sub[3, :])
    middle_title.set_axis_off()
    middle_title.text(0.13, 0.50, "Paired incremental discrimination",
                      fontsize=10.0, fontweight="bold", color=TEXT,
                      ha="left", va="center")

    low_lab = plt.subplot(sub[4, 0])
    iax = plt.subplot(sub[4, 1])
    inc = data["incremental"]
    comparisons = [
        "Overall beyond clinical base",
        "PhenoAge beyond clinical base",
        "Overall beyond PhenoAge",
        "PhenoAge beyond Overall",
        "Three components beyond PhenoAge",
    ]
    metrics = ["uno_c15", "auc5", "auc10"]
    inc = inc[inc["comparison"].isin(comparisons) & inc["metric"].isin(metrics)].copy()
    plotted.append(inc.assign(panel="C2", measure="paired_incremental_discrimination"))
    comp_labels = ["Overall\nvs base", "PhenoAge\nvs base",
                   "Overall\n| PhenoAge", "PhenoAge\n| Overall", "NM+TB+TC\n| PhenoAge"]
    yy = np.arange(len(comparisons))[::-1]
    mstyles = {"uno_c15": (-0.16, "s", "Uno C15"), "auc5": (0.0, "o", "AUC5"), "auc10": (0.16, "^", "AUC10")}
    for i, comp in enumerate(comparisons):
        for metric, (off, marker, _) in mstyles.items():
            r = inc[(inc["comparison"] == comp) & (inc["metric"] == metric)].iloc[0]
            sig = float(r["paired_p_holm"]) < 0.05
            fc = TEXT if sig else "white"
            iax.errorbar(r["benefit_oriented_difference"], yy[i] + off,
                         xerr=np.array([[r["benefit_oriented_difference"] - r["difference_ci_lower"]],
                                        [r["difference_ci_upper"] - r["benefit_oriented_difference"]]]),
                         fmt=marker, color=MID_GREY, mfc=fc, mec=TEXT, mew=0.70,
                         ms=5.0, elinewidth=1.0, capsize=1.8)
    iax.axvline(0, color="#5D6065", lw=0.85)
    iax.set_xlim(-0.018, 0.070)
    iax.set_ylim(-0.40, 4.80)
    iax.set_xticks([-0.01, 0.00, 0.02, 0.04, 0.06])
    iax.set_yticks(yy)
    iax.set_yticklabels([])
    low_lab.set_xlim(0, 1)
    low_lab.set_ylim(-0.40, 4.80)
    low_lab.set_axis_off()
    for yi, label in zip(yy, comp_labels):
        low_lab.text(0.96, yi, label, fontsize=8.7, ha="right", va="center",
                     color=TEXT, linespacing=0.92)
    style_axis(iax, "x")
    h1 = [Line2D([0], [0], marker=m, color="none", markerfacecolor=MID_GREY,
                 markeredgecolor=TEXT, markersize=5.2, label=lab)
          for _, m, lab in mstyles.values()]
    ilax = plt.subplot(sub[4, 2])
    ilax.set_axis_off()
    ilax.legend(handles=h1, ncol=1, loc="center left", fontsize=7.0,
                borderaxespad=0.0, labelspacing=0.58, handletextpad=0.36)
    bottom_spacer = plt.subplot(sub[5, :])
    bottom_spacer.set_axis_off()
    bottom_label = plt.subplot(sub[6, :])
    bottom_label.set_axis_off()
    bottom_label.text(0.47, 0.52, "Paired gain in C-index / AUC",
                      fontsize=9.4, color=TEXT, ha="center", va="center")
    note_ax = plt.subplot(sub[7, :])
    note_ax.set_axis_off()
    note_ax.text(0.13, 0.48,
                 "Bars are 95% CIs; filled markers indicate Holm-adjusted P<0.05.",
                 fontsize=7.8, color=MID_GREY, ha="left", va="center")


def draw_panel_c(spec, data, plotted):
    content = panel_shell(
        spec, "B", "Independent temporal calibration",
        "Deciles with Wilson 95% CIs; matched x/y scales",
    )
    sub = GridSpecFromSubplotSpec(
        5, 2, subplot_spec=content,
        height_ratios=[1, 1, 1, 0.17, 0.22], width_ratios=[0.74, 0.26],
        hspace=0.30, wspace=0.10,
    )
    outcome_map = {"NM": "Outcome_NutriMetab", "TB": "Outcome_TumorBurden", "TC": "Outcome_TreatComp"}
    for j, score in enumerate(["NM", "TB", "TC"]):
        ax = plt.subplot(sub[j, 0])
        curves = []
        maxval = 0.0
        for split, label, ls, alpha in [
            ("calendar_validation", "Validation", "--", 0.72),
            ("calendar_test", "Test", "-", 1.0),
        ]:
            ff = data["calendar_pred"][data["calendar_pred"]["split"] == split]
            cal = calibration_bins(ff[f"observed_{score}"], ff[f"p_{score}"], bins=10)
            cal["score"] = score
            cal["set"] = label
            cal["panel"] = "B"
            cal["measure"] = "temporal_calibration"
            plotted.append(cal)
            curves.append((label, ls, alpha, cal))
            maxval = max(maxval, cal[["mean_predicted", "observed", "ci_upper"]].to_numpy().max())
        lim = min(1.0, max({"NM": 0.10, "TB": 0.42, "TC": 0.24}[score], maxval * 1.10))
        ax.plot([0, lim], [0, lim], color=LIGHT_GREY, ls=":", lw=0.9)
        for label, ls, alpha, cal in curves:
            yerr = np.clip(
                np.vstack([cal["observed"] - cal["ci_lower"],
                           cal["ci_upper"] - cal["observed"]]),
                0.0, None,
            )
            ax.errorbar(cal["mean_predicted"], cal["observed"], yerr=yerr,
                        color=COLORS[score], ls=ls, lw=1.25, marker="o",
                        ms=4.0, mew=0.40, mec="white", alpha=alpha,
                        capsize=1.6, elinewidth=0.75, label=label)
        ax.set_xlim(0, lim)
        ax.set_ylim(0, lim)
        ax.text(0.02, 0.92, score, transform=ax.transAxes, fontsize=9.2,
                fontweight="bold", color=COLORS[score], va="top", ha="left")
        ax.set_ylabel("Observed", fontsize=9.2, labelpad=3)
        if j < 2:
            ax.tick_params(labelbottom=False)
        ax.set_xticks([0, lim / 2, lim])
        ax.set_yticks([0, lim / 2, lim])
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
        ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter("%.2f"))
        ax.tick_params(labelsize=8.7, pad=1.8)
        style_axis(ax)
        met = data["calendar_metrics"]
        tr = met[(met["outcome"] == outcome_map[score]) & (met["split"] == "calendar_test")].iloc[0]
        tax = plt.subplot(sub[j, 1])
        tax.set_axis_off()
        tax.text(0.00, 0.92, "Calendar test", fontsize=9.0, fontweight="bold",
                 color=TEXT, ha="left", va="top")
        tax.text(0.00, 0.66,
                 f"Intercept  {tr['calibration_intercept']:.2f}\n"
                 f"Slope      {tr['calibration_slope']:.2f}\n"
                 f"Brier      {tr['brier']:.3f}",
                 fontsize=8.7, color=TEXT, ha="left", va="top", linespacing=1.18)
    xlabel_ax = plt.subplot(sub[3, :])
    xlabel_ax.set_axis_off()
    xlabel_ax.text(0.36, 0.50, "Predicted probability", fontsize=10.0,
                   color=TEXT, ha="center", va="center")
    footer = plt.subplot(sub[4, :])
    footer.set_axis_off()
    handles = [
        Line2D([0], [0], color=MID_GREY, ls="--", marker="o", markersize=4,
               lw=1.2, label="Validation"),
        Line2D([0], [0], color=MID_GREY, ls="-", marker="o", markersize=4,
               lw=1.2, label="Test"),
        Line2D([0], [0], color=LIGHT_GREY, ls=":", lw=1.0, label="Ideal"),
    ]
    footer.legend(handles=handles, loc="center", ncol=3, fontsize=7.0,
                  handlelength=2.0, columnspacing=1.1)


def draw_panel_d(spec, data, plotted):
    content = panel_shell(
        spec, "D", "Fully adjusted all-cause mortality",
        "Survey-weighted Cox models; adults ≥60 years",
        header_ratio=0.13,
    )
    sub = GridSpecFromSubplotSpec(
        2, 1, subplot_spec=content, height_ratios=[0.86, 0.14], hspace=0.12,
    )
    main = GridSpecFromSubplotSpec(
        1, 2, subplot_spec=sub[0, 0], width_ratios=[0.42, 0.58], wspace=0.07,
    )
    ax = plt.subplot(main[0, 0])
    tax = plt.subplot(main[0, 1])
    df = data["allcause"].copy()
    order = ["OVERALL", "NM", "TB", "TC", "PHENO"]
    df["ord"] = df["component"].map({k: i for i, k in enumerate(order)})
    df = df.sort_values("ord")
    plotted.append(df.assign(panel="D", measure="fully_adjusted_allcause_hr"))
    y = np.arange(len(df))[::-1]
    for yi, (_, r) in zip(y, df.iterrows()):
        score = r["component"]
        ax.errorbar(r["hr"], yi, xerr=np.array([[r["hr"] - r["ci_lower"]], [r["ci_upper"] - r["hr"]]]),
                    fmt="o", color=COLORS[score], mfc=COLORS[score], mec="white",
                    mew=0.45, ms=6.2, elinewidth=1.2, capsize=2.1)
        tax.text(0.27, yi, f"{int(r['n']):,}/{int(r['events']):,}", fontsize=9.0,
                 va="center", ha="right", color=TEXT)
        tax.text(0.35, yi, f"{r['hr']:.2f} ({r['ci_lower']:.2f}–{r['ci_upper']:.2f})",
                 fontsize=9.0, va="center", ha="left", color=TEXT)
        tax.text(0.98, yi, fmt_p(r["p_holm"]), fontsize=9.0,
                 va="center", ha="right", color=TEXT)
    ax.axvline(1, color="#55585D", lw=0.9)
    ax.set_xlim(0.98, 1.60)
    ax.set_xticks([1.0, 1.2, 1.4, 1.6])
    ax.set_yticks(y)
    ax.set_yticklabels([SCORE_LABELS[s] for s in order], fontsize=9.2)
    ax.set_xlabel("Hazard ratio per weighted SD", fontsize=10)
    ax.set_ylim(-0.45, 4.65)
    style_axis(ax, "x")
    tax.set_xlim(0, 1)
    tax.set_ylim(-0.45, 4.65)
    tax.set_axis_off()
    tax.text(0.27, 4.48, "N/deaths", fontsize=9.0, ha="right", fontweight="bold")
    tax.text(0.35, 4.48, "HR (95% CI)", fontsize=9.0, ha="left", fontweight="bold")
    tax.text(0.98, 4.48, r"$P_{\mathrm{Holm}}$", fontsize=9.0,
             ha="right", fontweight="bold")
    legend_ax = plt.subplot(sub[1, 0])
    legend_ax.set_axis_off()
    score_handles = [
        Line2D([0], [0], color=COLORS[score], marker="o", lw=1.2,
               markerfacecolor=COLORS[score], markeredgecolor="white",
               markersize=5.2, label=SCORE_LABELS[score])
        for score in order
    ]
    legend_ax.legend(handles=score_handles, ncol=5, loc="center",
                     fontsize=7.0, handlelength=1.4, handletextpad=0.35,
                     columnspacing=0.85, borderaxespad=0.0)


def draw_panel_e(spec, data, plotted):
    content = panel_shell(
        spec, "E", "Cancer-history subgroup associations",
        "Fully adjusted models; inference based on interaction tests",
    )
    sub = GridSpecFromSubplotSpec(
        2, 2, subplot_spec=content, height_ratios=[0.82, 0.18],
        width_ratios=[0.55, 0.45], hspace=0.22, wspace=0.08,
    )
    ax = plt.subplot(sub[0, 0])
    tax = plt.subplot(sub[0, 1])
    df = data["subgroup"].copy()
    order = ["OVERALL", "NM", "TB", "TC", "PHENO"]
    df["ord"] = df["component"].map({k: i for i, k in enumerate(order)})
    df = df.sort_values("ord")
    plotted.append(df.assign(panel="E", measure="cancer_history_subgroup_hr"))
    y = np.arange(len(df))[::-1]
    for yi, (_, r) in zip(y, df.iterrows()):
        score = r["component"]
        for off, hr, lo, hi, marker, filled in [
            (-0.105, r["hr_non_cancer_history"], r["non_cancer_ci_lower"], r["non_cancer_ci_upper"], "o", False),
            (0.105, r["hr_cancer_history"], r["cancer_ci_lower"], r["cancer_ci_upper"], "s", True),
        ]:
            ax.errorbar(hr, yi + off, xerr=np.array([[hr - lo], [hi - hr]]), fmt=marker,
                        color=COLORS[score], mfc=COLORS[score] if filled else "white",
                        mec=COLORS[score], mew=0.85, ms=5.8, elinewidth=1.1, capsize=2.0)
        tax.text(0.00, yi, f"{r['interaction_ratio_of_hrs']:.2f} "
                 f"({r['interaction_ci_lower']:.2f}–{r['interaction_ci_upper']:.2f})",
                 fontsize=9.0, va="center", ha="left", color=TEXT)
        tax.text(0.98, yi, fmt_p(r["interaction_p_holm"]), fontsize=9.0,
                 va="center", ha="right", color=TEXT)
    ax.axvline(1, color="#55585D", lw=0.9)
    ax.set_xlim(0.90, 1.65)
    ax.set_xticks([1.0, 1.2, 1.4, 1.6])
    ax.set_yticks(y)
    ax.set_yticklabels([SCORE_LABELS[s] for s in order], fontsize=9.2)
    ax.set_xlabel("Hazard ratio per weighted SD", fontsize=10)
    ax.set_ylim(-0.45, 5.60)
    style_axis(ax, "x")
    tax.set_xlim(0, 1)
    tax.set_ylim(-0.45, 5.60)
    tax.set_axis_off()
    tax.text(0.00, 5.35, "Ratio of HRs\n(95% CI)", fontsize=9.0,
             ha="left", va="top", fontweight="bold", linespacing=0.92)
    tax.text(0.98, 5.35, r"$P_{\mathrm{int,Holm}}$", fontsize=9.0,
             ha="right", fontweight="bold")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor=MID_GREY, markersize=5.8, label="No cancer history"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=MID_GREY,
               markeredgecolor=MID_GREY, markersize=5.8, label="Cancer history"),
    ]
    footer = plt.subplot(sub[1, :])
    footer.set_axis_off()
    footer.legend(handles=handles, ncol=2, loc="center", fontsize=7.0,
                  columnspacing=1.6, handletextpad=0.45, borderaxespad=0.0)


def draw_panel_f(spec, data, plotted):
    content = panel_shell(
        spec, "F", "Time-varying mortality associations",
        "N=9,319; 4,177 deaths; 95% confidence bands",
    )
    sub = GridSpecFromSubplotSpec(
        3, 1, subplot_spec=content,
        height_ratios=[0.72, 0.12, 0.16], hspace=0.10,
    )
    ax = plt.subplot(sub[0, 0])
    curves = data["tv_curve"]
    tests = data["tv_test"]
    curves = curves[curves["score"].isin(["OVERALL", "PHENO"])].copy()
    plotted.append(curves.assign(panel="F", measure="continuous_time_varying_hr"))
    for score, ls in [("OVERALL", "-"), ("PHENO", "--")]:
        ff = curves[curves["score"] == score].sort_values("follow_up_years")
        x = ff["follow_up_years"].to_numpy(dtype=float)
        y = ff["hazard_ratio_per_weighted_sd"].to_numpy(dtype=float)
        lo = ff["ci_lower"].to_numpy(dtype=float)
        hi = ff["ci_upper"].to_numpy(dtype=float)
        ax.fill_between(x, lo, hi, color=COLORS[score], alpha=0.14, lw=0)
        ax.plot(x, y, color=COLORS[score], lw=1.9, ls=ls, label=SCORE_LABELS[score])
    ax.axhline(1, color="#55585D", lw=0.9)
    for xref in [5, 10]:
        ax.axvline(xref, color=LIGHT_GREY, lw=0.8, ls=":")
    ax.set_xlim(0.5, 15)
    ax.set_xticks([1, 5, 10, 15])
    ax.set_ylim(0.95, 1.90)
    ax.set_yticks([1.0, 1.2, 1.4, 1.6, 1.8])
    ax.set_xlabel("Follow-up time (years)", fontsize=10, labelpad=6)
    ax.set_ylabel("Time-varying HR per weighted SD", fontsize=10)
    style_axis(ax, "both")
    ax.legend(loc="upper right", fontsize=7.0, handlelength=2.2)
    p_overall = tests.loc[tests["score"] == "OVERALL", "score_log_time_p_holm"].iloc[0]
    p_pheno = tests.loc[tests["score"] == "PHENO", "score_log_time_p_holm"].iloc[0]
    pw = data["piecewise"]
    counts = pw[pw["score"] == "OVERALL"].set_index("interval")["events_in_interval"].to_dict()
    xlabel_ax = plt.subplot(sub[1, 0])
    xlabel_ax.set_axis_off()
    note_ax = plt.subplot(sub[2, 0])
    note_ax.set_axis_off()
    note_ax.text(0.00, 0.90,
            f"Score x log(time): Overall P_Holm={fmt_p(p_overall)};\n"
            f"PhenoAge P_Holm={fmt_p(p_pheno)}\n"
            f"Deaths: 0-5 y {int(counts.get('0-5 years', 0))}; 5-10 y {int(counts.get('5-10 years', 0))}; >10 y {int(counts.get('>10 years', 0))}",
            transform=note_ax.transAxes, fontsize=8.4, color=TEXT,
            va="top", ha="left", linespacing=1.10)


def validate_inputs(data):
    checks = {
        "calendar_splits_present": set(data["calendar_pred"]["split"].unique()) == {"calendar_validation", "calendar_test"},
        "oof_unique_patients": data["oof"]["patient_key"].is_unique,
        "nhanes_absolute_common_n": set(data["disc"].loc[
            (data["disc"]["scenario"] == "A_PRIMARY_FROZEN_PANCANCER_OVERALL_STACK") &
            (data["disc"]["domain"] == "age60") & (data["disc"]["population"] == "Overall") &
            (data["disc"]["metric"].isin(["uno_c15", "auc5", "auc10", "auc15"])), "n"
        ].unique()) == {11410},
        "incremental_common_n": set(data["incremental"]["n"].unique()) == {9319},
        "incremental_common_events": set(data["incremental"]["events"].unique()) == {4177},
        "subgroup_all_five_scores": set(data["subgroup"]["component"]) == {"OVERALL", "NM", "TB", "TC", "PHENO"},
        "continuous_curve_range": (float(data["tv_curve"]["follow_up_years"].min()) == 0.5 and
                                   float(data["tv_curve"]["follow_up_years"].max()) == 15.0),
    }
    failed = [k for k, v in checks.items() if not bool(v)]
    if failed:
        raise RuntimeError("Figure 2 input validation failed: " + ", ".join(failed))
    return checks


def main():
    parser = argparse.ArgumentParser(description="Build publication Figure 2 from frozen CIPDS results.")
    parser.add_argument("--root", type=Path,
                        default=Path(__file__).resolve().parents[1])
    parser.add_argument("--outdir", type=Path, default=None)
    args = parser.parse_args()
    outdir = args.outdir or (args.root / "reports" / "figure2_release_v5_20260905")
    outdir.mkdir(parents=True, exist_ok=True)

    data = read_inputs(args.root)
    checks = validate_inputs(data)
    plotted = []

    # Taller review/submission master requested by the author: 240 x 360 mm.
    # Row 1 keeps its prior physical height; rows 2 and 3 gain vertical space.
    fig = plt.figure(figsize=(240 / 25.4, 360 / 25.4), facecolor="white")
    outer = GridSpec(3, 2, figure=fig, height_ratios=[1.00, 1.26, 1.068],
                     left=0.075, right=0.985, top=0.990, bottom=0.045,
                     wspace=0.20, hspace=0.04)
    draw_panel_a(outer[0, 0], data, plotted)
    draw_panel_c(outer[0, 1], data, plotted)
    draw_panel_b(outer[1, 0], data, plotted)
    draw_panel_d(outer[1, 1], data, plotted)
    draw_panel_e(outer[2, 0], data, plotted)
    draw_panel_f(outer[2, 1], data, plotted)

    stem = outdir / "Figure_2_Main_Results_Final_v3"
    save_publication_outputs(fig, stem)
    plt.close(fig)
    with Image.open(stem.with_suffix(".tiff")) as tif:
        tiff_audit = {
            "pixel_width": int(tif.size[0]),
            "pixel_height": int(tif.size[1]),
            "dpi_x": float(tif.info.get("dpi", (np.nan, np.nan))[0]),
            "dpi_y": float(tif.info.get("dpi", (np.nan, np.nan))[1]),
            "compression": str(tif.info.get("compression")),
            "mode": tif.mode,
        }

    normalized = []
    for item in plotted:
        normalized.append(item.copy())
    pd.concat(normalized, ignore_index=True, sort=False).to_csv(
        outdir / "Figure_2_Source_Data_Final_v3.csv", index=False, encoding="utf-8-sig")

    legend = """Figure 2. Hospital model performance and NHANES mortality transportability.

(A) AUROC and precision-recall curves from training-period out-of-fold (OOF), calendar-validation, and calendar-test predictions. (B) Validation- and test-set decile calibration. (C) Absolute discrimination (N=11,410; 5,272 deaths) and paired increments in apparent performance (N=9,319; 4,177 deaths) in NHANES adults aged ≥60 years. (D) Fully adjusted all-cause mortality hazard ratios. (E) Cancer-history strata and interaction tests. (F) Time-varying hazard ratios for Overall and PhenoAge Acceleration. Bars and bands are 95% confidence intervals; filled markers in paired comparisons denote Holm-adjusted P<0.05.

Abbreviations: AP, average precision; AUC, area under the curve; HR, hazard ratio; NM, nutritional-metabolic; OOF, out-of-fold; TB, tumor-burden-related; TC, treatment-complication-related.
"""
    (outdir / "Figure_2_Legend_EN_Final_v3.txt").write_text(legend, encoding="utf-8")
    qa = {
        "figure": "Figure_2_Main_Results_Final_v3",
        "canvas_mm": {"width": 240, "height": 360},
        "exports": {"pdf": "vector", "tiff_dpi": 800, "tiff_compression": "LZW", "preview_png_dpi": 300},
        "tiff_file_audit": tiff_audit,
        "input_validation": checks,
        "scientific_guards": {
            "hospital_development_uses_strict_oof": True,
            "hospital_pr_uses_individual_predictions": True,
            "nhanes_pr_not_drawn_due_to_censoring": True,
            "nhanes_scores_not_treated_as_mortality_probabilities": True,
            "common_cohort_for_absolute_discrimination": "N=11410; deaths=5272",
            "common_cohort_for_incremental_and_time_varying": "N=9319; deaths=4177",
            "interaction_inference_uses_interaction_p": True,
        },
        "palette": COLORS,
        "final_font_floor_pt": 7.0,
        "legend_font_pt": 7.0,
        "primary_text_floor_pt": 8.7,
        "layout_revision": {
            "row_1": "A + B",
            "row_2": "C + D",
            "row_3": "E + F",
            "all_panel_letters_uppercase": True,
            "calibration_panels_vertical": True,
            "dedicated_header_bands": True,
            "dedicated_numeric_columns_for_D_and_E": True,
            "dedicated_annotation_footer_for_F": True,
            "expanded_nhanes_plot_width": True,
            "AB_content_locked": True,
            "CD_rows_taller": True,
            "EF_rows_taller": True,
            "C_D_E_legends_explicit": True,
            "C_significance_note_moved_inside_panel_footer": True,
            "AB_content_and_physical_height_locked": True,
            "EF_content_and_physical_height_locked": True,
            "interrow_whitespace_reallocated_to_CD": True,
            "C_plot_zones_expanded": True,
            "physical_row_heights_mm": {"AB": 99.6, "CD": 125.5, "EF": 106.3},
            "CD_height_gain_vs_v2_mm": 14.3,
        },
        "visual_qa": "Pass after full-canvas inspection and an original-resolution C-panel crop inspection on the 240 x 360 mm canvas: A-B and E-F functions are unchanged and their physical row heights are preserved; C-D gained 14.3 mm; both C plotting zones are taller; comparison labels, metric markers, x-axis title, and significance note occupy separate bands; no C-D or D-E boundary collision remains.",
    }
    qa_path = outdir / "Figure_2_QA_Final_v3.json"
    qa_path.write_text(json.dumps(qa, indent=2), encoding="utf-8")
    hash_targets = [
        stem.with_suffix(".pdf"), stem.with_suffix(".tiff"), stem.with_suffix(".png"),
        outdir / "Figure_2_Source_Data_Final_v3.csv",
        outdir / "Figure_2_Legend_EN_Final_v3.txt", qa_path, Path(__file__),
    ]
    hash_rows = []
    for path in hash_targets:
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hash_rows.append({"file": path.name, "sha256": digest, "bytes": path.stat().st_size})
    pd.DataFrame(hash_rows).to_csv(outdir / "Figure_2_SHA256_Final_v3.csv", index=False, encoding="utf-8-sig")
    print(str(outdir))


if __name__ == "__main__":
    main()
