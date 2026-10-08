# Asset Confirmation Table (Academic Figure Skill)
# Panel (a) cancer-trained 3D score space -> assets/figures/Manifold/plot_hole_manifold.py -> parameter inheritance: transparent panes, fine surface mesh, fixed camera, restrained depth
# Panel (b) PhenoAge 3D response surface -> assets/figures/Manifold/plot_hole_manifold.py -> parameter inheritance: translucent surface, identical camera and box aspect
# Panel (c) directional comparison -> assets/figures/BarComparison/plot_comparison_GeneRegulatory.py -> parameter inheritance: compact paired marks, direct labels, zero reference
# Panel (d) paired discrimination forest -> assets/figures/BarComparison/plot_comparison_Trajectory.py -> parameter inheritance: ordered estimates with confidence intervals
# Asset decision: no native-run template was semantically valid; only the listed visual parameters were inherited.

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


import json
import math
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure


ROOT = pathlib.Path(__file__).resolve().parents[1]
GEO_DIR = ROOT / "outputs" / "geometry_v1"
FIG_DIR = ROOT / "figures" / "geometry_v1"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MODEL_SPACE = GEO_DIR / "geometry_model_space_data.csv"
SURFACE_COEFFICIENTS = GEO_DIR / "geometry_pheno_surface_coefficients.csv"
DIRECTION_METRICS = GEO_DIR / "geometry_direction_metrics.csv"
AXIS_LIMITS = GEO_DIR / "geometry_axis_limits.csv"
PLANES = GEO_DIR / "geometry_overall_iso_probability_planes.csv"
SURFACE_QA = GEO_DIR / "geometry_surface_qa.json"
GATE_SUMMARY = GEO_DIR / "pareto_gate_summary.csv"
GATE_QA = GEO_DIR / "pareto_gate_qa.json"
NULL_METRICS = GEO_DIR / "pareto_null_metrics.csv"
PAIRED = ROOT / "outputs" / "nhanes_paired_discrimination_comparisons.csv"


def load_inputs():
    required = [
        MODEL_SPACE, SURFACE_COEFFICIENTS, DIRECTION_METRICS, AXIS_LIMITS,
        PLANES, SURFACE_QA, GATE_SUMMARY, GATE_QA, NULL_METRICS, PAIRED,
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing figure inputs: {missing}")
    data = pd.read_csv(MODEL_SPACE)
    coefficients = pd.read_csv(SURFACE_COEFFICIENTS).set_index("term")["estimate"]
    metrics = pd.read_csv(DIRECTION_METRICS).set_index("metric")
    limits_frame = pd.read_csv(AXIS_LIMITS).set_index("axis")
    limits = {
        axis: (float(limits_frame.loc[axis, "lower_005"]), float(limits_frame.loc[axis, "upper_995"]))
        for axis in ("z_NM", "z_TB", "z_TC")
    }
    planes = pd.read_csv(PLANES)
    with SURFACE_QA.open(encoding="utf-8") as handle:
        surface_qa = json.load(handle)
    with GATE_QA.open(encoding="utf-8") as handle:
        gate_qa = json.load(handle)
    gate = pd.read_csv(GATE_SUMMARY)
    null_metrics = pd.read_csv(NULL_METRICS)
    paired = pd.read_csv(PAIRED)
    return data, coefficients, metrics, limits, planes, surface_qa, gate_qa, gate, null_metrics, paired


def quadratic_prediction(coefficients: pd.Series, x, y, z):
    get = lambda name: float(coefficients.get(name, 0.0))
    return (
        get("(Intercept)") + get("z_NM") * x + get("z_TB") * y + get("z_TC") * z
        + get("I(z_NM^2)") * x**2 + get("I(z_TB^2)") * y**2 + get("I(z_TC^2)") * z**2
        + get("z_NM:z_TB") * x * y + get("z_NM:z_TC") * x * z
        + get("z_TB:z_TC") * y * z
    )


def style_3d_axis(ax, limits, show_ticklabels=True):
    ax.set_xlim(limits["z_NM"])
    ax.set_ylim(limits["z_TB"])
    ax.set_zlim(limits["z_TC"])
    ax.set_box_aspect((1, 1, 0.82))
    ax.view_init(elev=23, azim=-54)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1, 1, 1, 0))
        axis.pane.set_edgecolor("#D9D9D9")
        axis._axinfo["grid"].update(color="#D9D9D9", linewidth=0.35, linestyle="-")
        axis._axinfo["axisline"].update(color="#777777", linewidth=0.5)
    ax.tick_params(pad=-1, labelsize=6)
    if not show_ticklabels:
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
    ax.set_xlabel("NM z-logit", labelpad=0, fontsize=6.5)
    ax.set_ylabel("TB z-logit", labelpad=0, fontsize=6.5)
    ax.set_zlabel("TC z-logit", labelpad=-5, fontsize=6.5)


def plot_score_plane(ax, limits, planes):
    first = planes.iloc[0]
    beta = np.array([first.beta_z_NM, first.beta_z_TB, first.beta_z_TC], dtype=float)
    x = np.linspace(*limits["z_NM"], 70)
    y = np.linspace(*limits["z_TB"], 70)
    xx, yy = np.meshgrid(x, y)
    zz = (float(first.eta) - beta[0] * xx - beta[1] * yy) / beta[2]
    valid = (zz >= limits["z_TC"][0]) & (zz <= limits["z_TC"][1])
    zz = np.where(valid, zz, np.nan)
    ax.plot_surface(
        xx, yy, zz,
        color=ACCENT_RED, alpha=0.20, linewidth=0.20, edgecolor="#7F0000",
        antialiased=True, shade=False, rstride=3, cstride=3,
    )
    return int(valid.sum())


def plot_pheno_zero_surface(ax, coefficients, limits):
    n_grid = 56
    x = np.linspace(*limits["z_NM"], n_grid)
    y = np.linspace(*limits["z_TB"], n_grid)
    z = np.linspace(*limits["z_TC"], n_grid)
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    volume = quadratic_prediction(coefficients, xx, yy, zz)
    if not (float(np.nanmin(volume)) <= 0 <= float(np.nanmax(volume))):
        raise RuntimeError("PhenoAge zero surface does not cross the plotting support")
    spacing = (x[1] - x[0], y[1] - y[0], z[1] - z[0])
    vertices, faces, _, _ = measure.marching_cubes(volume.astype(np.float32), level=0, spacing=spacing)
    vertices += np.array([x[0], y[0], z[0]])
    collection = Poly3DCollection(
        vertices[faces], facecolor=CATEGORICAL[0], edgecolor="#164A7B",
        linewidth=0.12, alpha=0.20,
    )
    collection.set_rasterized(True)
    ax.add_collection3d(collection)
    return len(vertices), len(faces), float(np.nanmin(volume)), float(np.nanmax(volume))


def format_p(value: float) -> str:
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def add_panel_label(ax, label):
    if hasattr(ax, "text2D"):
        ax.text2D(-0.10, 1.02, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")
    else:
        ax.text(-0.10, 1.02, label, transform=ax.transAxes, fontsize=9, fontweight="bold", va="top")


def build_gate_figure(gate: pd.DataFrame, null_metrics: pd.DataFrame):
    primary_gate = gate[gate.scenario == "primary_nonoverlap_9"].sort_values("k")
    primary_null = null_metrics[null_metrics.scenario == "primary_nonoverlap_9"]
    summaries = (
        primary_null.groupby(["k", "null_type"])
        .agg(
            t_median=("t_ratio", "median"),
            t_low=("t_ratio", lambda values: values.quantile(0.025)),
            t_high=("t_ratio", lambda values: values.quantile(0.975)),
            rss_median=("rss", "median"),
        )
        .reset_index()
    )
    fig, axes = plt.subplots(
        1, 2, figsize=(179 / 25.4, 73 / 25.4),
        gridspec_kw={"left": 0.08, "right": 0.965, "bottom": 0.19, "top": 0.86, "wspace": 0.30},
    )
    style = {
        "classic_permutation": (CATEGORICAL[0], "Classic PC permutation null"),
        "gaussian_covariance": (CATEGORICAL[3], "Gaussian covariance null"),
    }
    observed_k = primary_gate.k.to_numpy()

    ax = axes[0]
    for null_type, (color, label) in style.items():
        subset = summaries[summaries.null_type == null_type].sort_values("k")
        ax.fill_between(subset.k, subset.t_low, subset.t_high, color=color, alpha=0.16, linewidth=0)
        ax.plot(subset.k, subset.t_median, marker="o", markersize=3.5, color=color,
                linewidth=1.0, label=f"{label} (median, 95% interval)")
    ax.plot(observed_k, primary_gate.observed_t_ratio, color=BLACK, marker="D", markersize=4,
            linewidth=1.1, label="Observed hospital data")
    ax.set_yscale("log")
    ax.set_xticks(observed_k)
    ax.set_xlabel("Candidate archetype count (K)")
    ax.set_ylabel("T-ratio (polytope volume / data-hull volume)")
    ax.set_title("Geometry across archetype counts", loc="left", fontweight="bold")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 0), fontsize=6.2)
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.35)
    add_panel_label(ax, "a")

    ax = axes[1]
    for null_type, (color, label) in style.items():
        subset = summaries[summaries.null_type == null_type].sort_values("k")
        ratios = primary_gate.observed_rss.to_numpy() / subset.rss_median.to_numpy()
        ax.plot(observed_k, ratios, marker="s", markersize=4, color=color,
                linewidth=1.1)
        ax.text(6.10, ratios[-1], label, color=color, fontsize=6.2, va="center", ha="left")
    ax.axhline(1, color=BLACK, linewidth=0.7, linestyle="--")
    ax.set_yscale("log")
    ax.set_xticks(observed_k)
    ax.set_xlim(1.8, 7.55)
    ax.set_xlabel("Candidate archetype count (K)")
    ax.set_ylabel("Observed RSS / null median RSS")
    ax.set_title("Residual fit under the two null models", loc="left", fontweight="bold")
    ax.text(
        0.02, 0.04,
        "Pass rule: both T-ratio and RSS significant after Holm correction\n"
        "Primary K=2–6: no dual-null pass; K=4 Holm P=1.00 (classic)",
        transform=ax.transAxes, fontsize=6.2, va="bottom", color="#444444",
    )
    ax.grid(axis="y", color="#E5E5E5", linewidth=0.35)
    add_panel_label(ax, "b")

    output_stem = FIG_DIR / "figS_pareto_geometry_gate"
    save_cns_figure(fig, str(output_stem))
    plt.close(fig)
    return str(output_stem.with_suffix(".pdf")), str(output_stem.with_suffix(".png"))


def main():
    data, coefficients, metrics, limits, planes, surface_qa, gate_qa, gate, null_metrics, paired = load_inputs()
    if len(data) != 11410 or data["SEQN"].nunique() != 11410:
        raise RuntimeError("Figure cohort grain drift")
    if surface_qa["selected_pheno_surface"] != "quadratic":
        raise RuntimeError("Figure is locked to the formally selected quadratic PhenoAge surface")
    if gate_qa["tetrahedron_promotable"]:
        raise RuntimeError("Figure route mismatch: tetrahedron unexpectedly passed the formal gate")
    k4 = gate[(gate.scenario == "primary_nonoverlap_9") & (gate.k == 4)].iloc[0]

    overall_cmap = LinearSegmentedColormap.from_list(
        "overall_vulnerability", ["#F7FBFF", "#9ECAE1", "#F1A340", "#B2182B"]
    )
    pheno_cmap = LinearSegmentedColormap.from_list("pheno_acceleration", DIVERGING)
    pheno_bound = float(np.nanquantile(np.abs(data["PhenoAge_acceleration"]), 0.975))

    fig = plt.figure(figsize=(183 / 25.4, 158 / 25.4), facecolor="white")
    grid = fig.add_gridspec(
        2, 2, height_ratios=[1.18, 0.82], width_ratios=[1, 1],
        left=0.07, right=0.955, bottom=0.09, top=0.965, wspace=0.28, hspace=0.50,
    )

    ax_a = fig.add_subplot(grid[0, 0], projection="3d")
    overall_norm = Normalize(
        vmin=float(np.nanquantile(data["Overall_expected_burden"], 0.01)),
        vmax=float(np.nanquantile(data["Overall_expected_burden"], 0.995)),
    )
    ax_a.scatter(
        data.z_NM, data.z_TB, data.z_TC,
        c=data.Overall_expected_burden, cmap=overall_cmap, norm=overall_norm,
        s=2.0, alpha=0.22, linewidths=0, depthshade=False, rasterized=True,
    )
    valid_plane_cells = plot_score_plane(ax_a, limits, planes)
    style_3d_axis(ax_a, limits)
    ax_a.set_title("Cancer-trained multidomain score", loc="left", fontweight="bold", pad=3)
    ax_a.text2D(
        0.02, 0.94,
        "Red mesh: P(any domain)=0.50\n"
        "P(≥2 domains)=0.50 lies beyond 99% NHANES support",
        transform=ax_a.transAxes, fontsize=6.2, color="#6B0F13", va="top",
    )
    ax_a.text2D(
        0.50, 0.86,
        f"ParTI K=4 gate: Holm P={k4.classic_t_p_holm:.2f}; no tetrahedron",
        transform=ax_a.transAxes, fontsize=6.2, color="#555555", va="top",
    )
    colorbar_axis_a = ax_a.inset_axes([0.16, -0.17, 0.68, 0.035])
    colorbar_a = fig.colorbar(ScalarMappable(norm=overall_norm, cmap=overall_cmap), cax=colorbar_axis_a,
                             orientation="horizontal")
    colorbar_a.set_label("Frozen Overall expected burden (0–3)", fontsize=6.5, labelpad=1)
    colorbar_a.ax.tick_params(labelsize=6, length=2)
    add_panel_label(ax_a, "a")

    ax_b = fig.add_subplot(grid[0, 1], projection="3d")
    pheno_norm = TwoSlopeNorm(vmin=-pheno_bound, vcenter=0, vmax=pheno_bound)
    ax_b.scatter(
        data.z_NM, data.z_TB, data.z_TC,
        c=data.PhenoAge_acceleration.clip(-pheno_bound, pheno_bound),
        cmap=pheno_cmap, norm=pheno_norm,
        s=2.0, alpha=0.20, linewidths=0, depthshade=False, rasterized=True,
    )
    surface_vertices, surface_faces, surface_min, surface_max = plot_pheno_zero_surface(
        ax_b, coefficients, limits
    )
    style_3d_axis(ax_b, limits)
    ax_b.set_title("PhenoAge Acceleration in the same space", loc="left", fontweight="bold", pad=3)
    nonlinear_p = float(surface_qa["nonlinear_block_wald_p"])
    ax_b.text2D(
        0.02, 0.94,
        "Blue mesh: predicted acceleration = 0 years\n"
        f"Quadratic surface selected; global Wald P={nonlinear_p:.1e}",
        transform=ax_b.transAxes, fontsize=6.2, color="#164A7B", va="top",
    )
    colorbar_axis_b = ax_b.inset_axes([0.16, -0.17, 0.68, 0.035])
    colorbar_b = fig.colorbar(ScalarMappable(norm=pheno_norm, cmap=pheno_cmap), cax=colorbar_axis_b,
                             orientation="horizontal")
    colorbar_b.set_label("Observed PhenoAge Acceleration (years)", fontsize=6.5, labelpad=1)
    colorbar_b.ax.tick_params(labelsize=6, length=2)
    add_panel_label(ax_b, "b")

    ax_c = fig.add_subplot(grid[1, 0])
    beta = np.array(list(surface_qa["overall_gradient"].values()), dtype=float)
    beta_unit = beta / np.linalg.norm(beta)
    labels = ["NM", "TB", "TC"]
    y = np.arange(3)[::-1]
    pheno_rows = metrics.loc[[f"unit_gradient_{label}_at_origin" for label in labels]]
    pheno_est = pheno_rows.estimate.to_numpy()
    pheno_low = pheno_rows.ci_lower.to_numpy()
    pheno_high = pheno_rows.ci_upper.to_numpy()
    for row_y, a_value, b_value in zip(y, beta_unit, pheno_est, strict=True):
        ax_c.plot([a_value, b_value], [row_y, row_y], color="#C7C7C7", linewidth=1.0, zorder=1)
    ax_c.scatter(beta_unit, y, s=30, marker="o", color=ACCENT_RED, edgecolor="white", linewidth=0.5,
                 label="Overall (frozen)", zorder=3)
    ax_c.errorbar(
        pheno_est, y,
        xerr=np.vstack([pheno_est - pheno_low, pheno_high - pheno_est]),
        fmt="D", markersize=4.3, color=CATEGORICAL[0], ecolor=CATEGORICAL[0],
        elinewidth=0.8, capsize=2, label="PhenoAge surface (95% CI)", zorder=4,
    )
    ax_c.axvline(0, color="#777777", linewidth=0.6, linestyle="--")
    ax_c.set_yticks(y, ["NM pattern", "TB pattern", "TC pattern"])
    ax_c.set_xlim(-0.22, 1.08)
    ax_c.set_xlabel("Unit-gradient component at the score-space origin")
    ax_c.set_title("Score-space directions differ", loc="left", fontweight="bold")
    angle = metrics.loc["weighted_mean_angle_degrees"]
    correlation = metrics.loc["weighted_overall_pheno_correlation"]
    ax_c.text(
        0.02, 0.38,
        f"Mean local angle {angle.estimate:.1f}° ({angle.ci_lower:.1f}–{angle.ci_upper:.1f})\n"
        f"Survey-weighted r={correlation.estimate:.3f} ({correlation.ci_lower:.3f}–{correlation.ci_upper:.3f})",
        transform=ax_c.transAxes, fontsize=6.5, va="center", color="#333333",
    )
    ax_c.legend(loc="upper right", bbox_to_anchor=(1.0, 0.98), fontsize=6.4, handletextpad=0.5)
    ax_c.tick_params(axis="y", length=0)
    ax_c.grid(axis="x", color="#E5E5E5", linewidth=0.35)
    add_panel_label(ax_c, "c")

    ax_d = fig.add_subplot(grid[1, 1])
    comparison = paired[
        (paired.domain == "age60") & (paired.population == "Overall")
        & (paired.model_a == "OVERALL") & (paired.model_b == "PHENO")
    ].copy()
    metric_order = ["uno_c15", "auc5", "auc10", "auc15", "harrell_c"]
    display_label = {
        "uno_c15": "Uno C (15 y; primary)",
        "auc5": "AUC (5 y)",
        "auc10": "AUC (10 y)",
        "auc15": "AUC (15 y)",
        "harrell_c": "Harrell C (sensitivity)",
    }
    comparison = comparison.set_index("metric").loc[metric_order].reset_index()
    y_d = np.arange(len(comparison))[::-1]
    significant = comparison.paired_p_holm < 0.05
    point_colors = np.where(significant, CATEGORICAL[4], GREY)
    for index, row in comparison.iterrows():
        yy = y_d[index]
        ax_d.plot([row.difference_ci_lower, row.difference_ci_upper], [yy, yy],
                  color=point_colors[index], linewidth=1.0, zorder=1)
        ax_d.plot([row.difference_ci_lower, row.difference_ci_lower], [yy - 0.08, yy + 0.08],
                  color=point_colors[index], linewidth=0.7)
        ax_d.plot([row.difference_ci_upper, row.difference_ci_upper], [yy - 0.08, yy + 0.08],
                  color=point_colors[index], linewidth=0.7)
        ax_d.scatter(row.paired_difference, yy, s=24, marker="s", color=point_colors[index],
                     edgecolor="white", linewidth=0.45, zorder=3)
        ax_d.text(
            0.019, yy,
            f"P$_{{Holm}}$={format_p(float(row.paired_p_holm))}",
            ha="right", va="center", fontsize=6.0, color="#444444",
        )
    ax_d.axvline(0, color=BLACK, linewidth=0.7)
    ax_d.set_yticks(y_d, [display_label[value] for value in metric_order])
    ax_d.set_xlim(-0.067, 0.021)
    ax_d.set_xlabel("Difference in C-index/AUC (Overall − PhenoAge)")
    ax_d.set_title(
        "Paired mortality discrimination\nOverall is not superior to PhenoAge",
        loc="left", fontweight="bold", linespacing=1.45,
    )
    ax_d.text(
        0.02, 0.18, "N=11,410; deaths=5,272; 1,000 survey bootstraps",
        transform=ax_d.transAxes, ha="left", va="center", fontsize=6.0, color="#555555",
    )
    ax_d.grid(axis="x", color="#E5E5E5", linewidth=0.35)
    ax_d.tick_params(axis="y", length=0)
    add_panel_label(ax_d, "d")

    output_stem = FIG_DIR / "fig_model_geometry_vs_phenoage"
    save_cns_figure(fig, str(output_stem))
    plt.close(fig)
    gate_pdf, gate_png = build_gate_figure(gate, null_metrics)

    audit = {
        "figure_width_mm": 183,
        "figure_height_mm": 158,
        "cohort_n": int(len(data)),
        "unique_seqn_n": int(data.SEQN.nunique()),
        "all_nhanes_points_plotted": True,
        "panel_a_any_domain_plane_valid_grid_cells": valid_plane_cells,
        "panel_b_zero_surface_vertices": int(surface_vertices),
        "panel_b_zero_surface_faces": int(surface_faces),
        "panel_b_predicted_range_on_grid": [surface_min, surface_max],
        "pareto_k4_classic_t_holm_p": float(k4.classic_t_p_holm),
        "tetrahedron_promotable": False,
        "selected_pheno_surface": "quadratic",
        "nonlinear_block_wald_p": nonlinear_p,
        "paired_comparisons_n": int(len(comparison)),
        "vector_master": str(output_stem.with_suffix(".pdf")),
        "png_preview": str(output_stem.with_suffix(".png")),
        "supplementary_gate_vector": gate_pdf,
        "supplementary_gate_png": gate_png,
        "checks_passed": True,
    }
    with (GEO_DIR / "geometry_figure_data_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
