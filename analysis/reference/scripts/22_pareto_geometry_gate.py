"""Run a ParTI-style, outcome-agnostic geometry gate for the CIPDS rebuild.

The gate uses the official ParTIpy 0.2.0 archetypal-analysis implementation
vendored from its signed PyPI wheel.  The code is loaded without importing the
single-cell-specific top-level dependencies that are irrelevant to this task.

Primary geometry excludes every CBC/anemia/marrow-suppression constituent and
albumin, so the polytope cannot be created by the concurrent outcome-defining
laboratories.  A 22-laboratory geometry is retained as a sensitivity analysis.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import pathlib
import sys
import time
from dataclasses import dataclass
from typing import Iterable

# Prevent 23 worker processes from each opening a full BLAS thread pool.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed, parallel_backend
from scipy.optimize import linear_sum_assignment
from scipy.spatial import ConvexHull, distance_matrix
from sklearn.decomposition import PCA


ROOT = pathlib.Path(__file__).resolve().parents[1]
PARTIPY_SOURCE = ROOT / ".vendor_partipy_src"
if str(PARTIPY_SOURCE) not in sys.path:
    sys.path.insert(0, str(PARTIPY_SOURCE))
os.environ["PYTHONPATH"] = str(PARTIPY_SOURCE) + os.pathsep + os.environ.get("PYTHONPATH", "")
OUTPUT_DIR = ROOT / "outputs" / "geometry_v1"
LOG_DIR = ROOT / "logs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(parents=True, exist_ok=True)


def _load_vendored_partipy():
    """Load the numerical AA modules without optional single-cell imports."""

    package_dir = PARTIPY_SOURCE / "partipy"
    if not package_dir.is_dir():
        raise FileNotFoundError(f"Vendored ParTIpy source not found: {package_dir}")

    from partipy.arch import AA
    from partipy.optim import _compute_A_projected_gradients

    return AA, _compute_A_projected_gradients


AA, _compute_A_projected_gradients = _load_vendored_partipy()


N_JOBS = int(os.getenv("CIPDS_THREADS", "23"))
N_PERM = int(os.getenv("CIPDS_GEOMETRY_PERMUTATIONS", "1000"))
N_GAUSSIAN = int(os.getenv("CIPDS_GEOMETRY_GAUSSIAN_PERMUTATIONS", str(N_PERM)))
N_SENSITIVITY = int(os.getenv("CIPDS_GEOMETRY_SENSITIVITY_PERMUTATIONS", str(N_PERM)))
N_BOOTSTRAP = int(os.getenv("CIPDS_GEOMETRY_BOOTSTRAPS", "1000"))
CORESET_SIZE = int(os.getenv("CIPDS_GEOMETRY_CORESET_SIZE", "5000"))
SEED = int(os.getenv("CIPDS_GEOMETRY_SEED", "20260831"))
K_VALUES = tuple(range(2, 7))
WINSOR_LOWER = 0.005
WINSOR_UPPER = 0.995

if N_JOBS < 1 or N_JOBS > 23:
    raise ValueError("CIPDS_THREADS must be between 1 and 23")
if min(N_PERM, N_GAUSSIAN, N_SENSITIVITY) < 50:
    raise ValueError("At least 50 null replicates are required")


LOGGER = logging.getLogger("pareto_geometry")
LOGGER.setLevel(logging.INFO)
LOGGER.handlers.clear()
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
file_handler = logging.FileHandler(LOG_DIR / "22_pareto_geometry_gate.log", mode="w", encoding="utf-8")
file_handler.setFormatter(formatter)
stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(formatter)
LOGGER.addHandler(file_handler)
LOGGER.addHandler(stream_handler)


OVERLAP_EXCLUDED = {
    "lab_baso_pct": "CBC index overlapping marrow-suppression/bleeding state",
    "lab_eos_pct": "CBC index overlapping marrow-suppression/bleeding state",
    "lab_hgb": "Direct anemia constituent",
    "lab_lym_pct": "CBC index overlapping marrow-suppression/bleeding state",
    "lab_mcv": "Erythrocyte index closely overlapping anemia/marrow state",
    "lab_mono_pct": "CBC index overlapping marrow-suppression/bleeding state",
    "lab_neu_pct": "CBC index overlapping marrow-suppression/bleeding state",
    "lab_plt": "CBC index overlapping marrow-suppression/bleeding state",
    "lab_rbc": "Erythrocyte index closely overlapping anemia",
    "lab_rdw": "Erythrocyte index closely overlapping anemia",
    "lab_wbc": "CBC index overlapping marrow-suppression state",
    "lab_mpv": "Platelet index overlapping marrow-suppression state",
    "lab_alb": "Direct low-albumin outcome constituent",
}


@dataclass
class Preprocessor:
    variables: list[str]
    lower: np.ndarray
    upper: np.ndarray
    median: np.ndarray
    mean: np.ndarray
    scale: np.ndarray

    @classmethod
    def fit(cls, frame: pd.DataFrame, variables: list[str]) -> "Preprocessor":
        raw = frame[variables].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        lower = np.nanquantile(raw, WINSOR_LOWER, axis=0)
        upper = np.nanquantile(raw, WINSOR_UPPER, axis=0)
        clipped = np.clip(raw, lower, upper)
        median = np.nanmedian(clipped, axis=0)
        imputed = np.where(np.isnan(clipped), median, clipped)
        mean = imputed.mean(axis=0)
        scale = imputed.std(axis=0, ddof=0)
        if np.any(~np.isfinite(scale) | (scale <= 0)):
            bad = [variables[i] for i in np.flatnonzero(~np.isfinite(scale) | (scale <= 0))]
            raise ValueError(f"Zero or invalid training scale: {bad}")
        return cls(variables, lower, upper, median, mean, scale)

    def transform(self, frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
        raw = frame[self.variables].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        missing = np.isnan(raw)
        clipped = np.clip(raw, self.lower, self.upper)
        imputed = np.where(missing, self.median, clipped)
        standardized = (imputed - self.mean) / self.scale
        return standardized.astype(np.float32), missing


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def holm_adjust(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted_ranked = np.maximum.accumulate((len(p) - np.arange(len(p))) * ranked)
    adjusted_ranked = np.minimum(adjusted_ranked, 1.0)
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = adjusted_ranked
    return adjusted


def t_ratio(X: np.ndarray, Z: np.ndarray) -> float:
    """ParTIpy definition: polytope volume divided by data convex-hull volume."""

    if X.ndim != 2 or Z.ndim != 2:
        raise ValueError("X and Z must be matrices")
    if X.shape[1] == 1:
        hull_volume = float(np.ptp(X[:, 0]))
        polytope_volume = float(np.ptp(Z[:, 0]))
    else:
        hull_volume = float(ConvexHull(X, qhull_options="QJ").volume)
        polytope_volume = float(ConvexHull(Z, qhull_options="QJ").volume)
    if hull_volume <= 0 or not np.isfinite(hull_volume):
        raise ValueError("Non-positive convex-hull volume")
    return polytope_volume / hull_volume


def fit_aa(X: np.ndarray, k: int, seed: int) -> dict:
    model = AA(
        n_archetypes=k,
        init="plus_plus",
        optim="projected_gradients",
        max_iter=250,
        rel_tol=1e-4,
        early_stopping=True,
        coreset_algorithm="standard",
        coreset_size=min(CORESET_SIZE, X.shape[0]),
        delta=0.0,
        centering=True,
        scaling=True,
        seed=int(seed),
        derivative_max_iter=80,
    ).fit(np.ascontiguousarray(X, dtype=np.float32))
    Z = np.asarray(model.Z, dtype=float)
    return {
        "Z": Z,
        "rss": float(model.RSS),
        "variance_explained": float(model.varexpl),
        "t_ratio": float(t_ratio(np.asarray(X, dtype=float), Z)),
        "converged": bool(model.fitting_info.get("conv")),
        "iterations": int(model.fitting_info.get("n_iter") or 0),
    }


def make_null(X: np.ndarray, null_type: str, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if null_type == "classic_permutation":
        return np.column_stack([rng.permutation(X[:, j]) for j in range(X.shape[1])]).astype(np.float32)
    if null_type == "gaussian_covariance":
        mean = X.mean(axis=0)
        covariance = np.cov(X, rowvar=False, ddof=1)
        if X.shape[1] == 1:
            scale = math.sqrt(float(np.asarray(covariance)))
            return rng.normal(float(mean[0]), scale, size=(X.shape[0], 1)).astype(np.float32)
        return rng.multivariate_normal(mean, covariance, size=X.shape[0]).astype(np.float32)
    raise ValueError(f"Unknown null type: {null_type}")


def run_null_batch(X: np.ndarray, k: int, null_type: str, seeds: list[int]) -> list[dict]:
    rows: list[dict] = []
    for seed in seeds:
        X_null = make_null(X, null_type, int(seed))
        fit = fit_aa(X_null, k=k, seed=int(seed) + 100_003)
        rows.append(
            {
                "replicate": int(seed),
                "t_ratio": fit["t_ratio"],
                "rss": fit["rss"],
                "variance_explained": fit["variance_explained"],
                "converged": fit["converged"],
                "iterations": fit["iterations"],
            }
        )
    return rows


def split_batches(values: np.ndarray, n_batches: int) -> list[list[int]]:
    return [list(map(int, batch)) for batch in np.array_split(values, min(n_batches, len(values))) if len(batch)]


def run_null_distribution(
    parallel: Parallel,
    X: np.ndarray,
    k: int,
    null_type: str,
    n_replicates: int,
    seed: int,
) -> pd.DataFrame:
    master = np.random.default_rng(seed)
    seeds = master.choice(np.arange(1, 2_000_000_000, dtype=np.int64), size=n_replicates, replace=False)
    batches = split_batches(seeds, N_JOBS)
    started = time.time()
    nested = parallel(delayed(run_null_batch)(X, k, null_type, batch) for batch in batches)
    rows = [row for batch in nested for row in batch]
    frame = pd.DataFrame(rows)
    frame["null_type"] = null_type
    frame["k"] = k
    LOGGER.info(
        "Completed %s null for k=%d: %d replicates in %.1f s",
        null_type,
        k,
        len(frame),
        time.time() - started,
    )
    return frame


def monte_carlo_p(observed: float, null: np.ndarray, metric: str) -> float:
    if metric == "t_ratio":
        more_extreme = np.abs(1.0 - null) <= abs(1.0 - observed)
    elif metric == "rss":
        more_extreme = null <= observed
    else:
        raise ValueError(metric)
    return float((1 + np.sum(more_extreme)) / (len(null) + 1))


def align_vertices(reference: np.ndarray, query: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    cost = distance_matrix(reference, query)
    row, col = linear_sum_assignment(cost)
    aligned = np.empty_like(query)
    aligned[row] = query[col]
    return aligned, cost[row, col]


def project_to_simplex(X: np.ndarray, Z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    A0 = np.full((X.shape[0], Z.shape[0]), 1.0 / Z.shape[0], dtype=np.float32)
    weights = _compute_A_projected_gradients(
        X=np.ascontiguousarray(X, dtype=np.float32),
        Z=np.ascontiguousarray(Z, dtype=np.float32),
        A=A0,
        derivative_max_iter=200,
        rel_tol_ls=1e-5,
        rel_tol_conv=1e-6,
    )
    reconstructed = weights @ Z
    error = np.sqrt(np.sum((X - reconstructed) ** 2, axis=1))
    return weights, error


def barycentric_coordinates(X: np.ndarray, Z: np.ndarray) -> np.ndarray:
    k = Z.shape[0]
    matrix = np.vstack([Z.T, np.ones(k)])
    rhs = np.column_stack([X, np.ones(X.shape[0])])
    return np.linalg.solve(matrix, rhs.T).T


def main() -> None:
    started = time.time()
    LOGGER.info(
        "Start geometry gate | workers=%d | classic=%d | gaussian=%d | sensitivity=%d | bootstrap=%d",
        N_JOBS,
        N_PERM,
        N_GAUSSIAN,
        N_SENSITIVITY,
        N_BOOTSTRAP,
    )

    wheel = ROOT / ".vendor_wheels" / "partipy-0.2.0-py3-none-any.whl"
    expected_wheel_hash = "33db20cb2e1917c818572ffcd2b3d21f5c281c8859234a3e9b333de64a2ef2bc"
    if not wheel.exists() or sha256(wheel) != expected_wheel_hash:
        raise RuntimeError("Vendored ParTIpy wheel missing or SHA256 mismatch")

    dictionary = pd.read_csv(OUTPUT_DIR / "geometry_lab_dictionary.csv")
    full_variables = dictionary["variable"].astype(str).tolist()
    if len(full_variables) != 22 or len(set(full_variables)) != 22:
        raise RuntimeError("Expected 22 unique cross-cohort geometry variables")
    primary_variables = [v for v in full_variables if v not in OVERLAP_EXCLUDED]
    if len(primary_variables) != 9:
        raise RuntimeError(f"Expected 9 non-overlapping primary variables, observed {len(primary_variables)}")

    hospital_path = ROOT / "outputs" / "hospital_patient_model_ready_129.csv"
    required = [
        "patient_key",
        "split_calendar_entry",
        "Outcome_NutriMetab",
        "Outcome_TumorBurden",
        "Outcome_TreatComp",
        *full_variables,
    ]
    hospital = pd.read_csv(hospital_path, usecols=required, low_memory=False)
    if len(hospital) != 75248 or hospital["patient_key"].nunique() != 75248:
        raise RuntimeError("Hospital patient grain drift")
    train = hospital.loc[hospital["split_calendar_entry"].eq("train")].copy()
    validation = hospital.loc[hospital["split_calendar_entry"].eq("validation")].copy()
    test = hospital.loc[hospital["split_calendar_entry"].eq("test")].copy()
    if (len(train), len(validation), len(test)) != (52683, 11300, 11265):
        raise RuntimeError("Hospital calendar split drift")

    scenarios = {
        "primary_nonoverlap_9": primary_variables,
        "sensitivity_all_A_22": full_variables,
    }
    preprocessor_map: dict[str, Preprocessor] = {}
    standardized_map: dict[str, dict[str, np.ndarray]] = {}
    pca_map: dict[tuple[str, int], PCA] = {}
    pc_map: dict[tuple[str, int], np.ndarray] = {}
    observed_map: dict[tuple[str, int], dict] = {}
    preprocessing_rows: list[dict] = []

    for scenario, variables in scenarios.items():
        preprocessor = Preprocessor.fit(train, variables)
        preprocessor_map[scenario] = preprocessor
        standardized_map[scenario] = {}
        for split_name, frame in (("train", train), ("validation", validation), ("test", test)):
            X, missing = preprocessor.transform(frame)
            standardized_map[scenario][split_name] = X
            if split_name == "train":
                for j, variable in enumerate(variables):
                    preprocessing_rows.append(
                        {
                            "scenario": scenario,
                            "variable": variable,
                            "winsor_lower": preprocessor.lower[j],
                            "winsor_upper": preprocessor.upper[j],
                            "training_median": preprocessor.median[j],
                            "training_mean_after_imputation": preprocessor.mean[j],
                            "training_sd_after_imputation": preprocessor.scale[j],
                            "training_missing_n": int(missing[:, j].sum()),
                            "training_missing_pct": float(100 * missing[:, j].mean()),
                            "excluded_from_primary": variable in OVERLAP_EXCLUDED,
                            "exclusion_reason": OVERLAP_EXCLUDED.get(variable, "Retained in non-overlap primary geometry"),
                        }
                    )

        X_train = standardized_map[scenario]["train"]
        for k in K_VALUES:
            pca = PCA(n_components=k - 1, svd_solver="full", random_state=SEED)
            pcs = pca.fit_transform(X_train).astype(np.float32)
            pca_map[(scenario, k)] = pca
            pc_map[(scenario, k)] = pcs
            fit = fit_aa(pcs, k=k, seed=SEED + 100 * k + (0 if scenario.startswith("primary") else 10_000))
            fit["pca_variance_cumulative"] = float(pca.explained_variance_ratio_.sum())
            observed_map[(scenario, k)] = fit
            LOGGER.info(
                "Observed %s k=%d | t=%.4f | RSS=%.1f | var=%.4f | PCvar=%.4f",
                scenario,
                k,
                fit["t_ratio"],
                fit["rss"],
                fit["variance_explained"],
                fit["pca_variance_cumulative"],
            )

    pd.DataFrame(preprocessing_rows).to_csv(OUTPUT_DIR / "pareto_preprocessing.csv", index=False)

    null_frames: list[pd.DataFrame] = []
    with parallel_backend("loky", inner_max_num_threads=1):
        with Parallel(n_jobs=N_JOBS, max_nbytes="20M", mmap_mode="r") as parallel:
            for scenario in scenarios:
                for k in K_VALUES:
                    X = pc_map[(scenario, k)]
                    n_classic = N_PERM if scenario.startswith("primary") else N_SENSITIVITY
                    frame = run_null_distribution(
                        parallel,
                        X,
                        k,
                        "classic_permutation",
                        n_classic,
                        SEED + 1_000_000 * (1 if scenario.startswith("primary") else 2) + 10_000 * k,
                    )
                    frame["scenario"] = scenario
                    null_frames.append(frame)
                    if scenario.startswith("primary"):
                        frame = run_null_distribution(
                            parallel,
                            X,
                            k,
                            "gaussian_covariance",
                            N_GAUSSIAN,
                            SEED + 3_000_000 + 10_000 * k,
                        )
                        frame["scenario"] = scenario
                        null_frames.append(frame)

    null_metrics = pd.concat(null_frames, ignore_index=True)
    null_metrics.to_csv(OUTPUT_DIR / "pareto_null_metrics.csv", index=False)

    summary_rows: list[dict] = []
    for scenario in scenarios:
        for k in K_VALUES:
            observed = observed_map[(scenario, k)]
            classic = null_metrics.loc[
                (null_metrics["scenario"] == scenario)
                & (null_metrics["k"] == k)
                & (null_metrics["null_type"] == "classic_permutation")
            ]
            gaussian = null_metrics.loc[
                (null_metrics["scenario"] == scenario)
                & (null_metrics["k"] == k)
                & (null_metrics["null_type"] == "gaussian_covariance")
            ]
            row = {
                "scenario": scenario,
                "k": k,
                "simplex_dimension": k - 1,
                "n_participants": len(train),
                "n_variables": len(scenarios[scenario]),
                "pca_variance_cumulative": observed["pca_variance_cumulative"],
                "observed_t_ratio": observed["t_ratio"],
                "observed_rss": observed["rss"],
                "observed_variance_explained": observed["variance_explained"],
                "observed_converged": observed["converged"],
                "classic_t_p": monte_carlo_p(observed["t_ratio"], classic["t_ratio"].to_numpy(), "t_ratio"),
                "classic_rss_p": monte_carlo_p(observed["rss"], classic["rss"].to_numpy(), "rss"),
                "classic_n": len(classic),
                "gaussian_t_p": np.nan,
                "gaussian_rss_p": np.nan,
                "gaussian_n": len(gaussian),
            }
            if len(gaussian):
                row["gaussian_t_p"] = monte_carlo_p(
                    observed["t_ratio"], gaussian["t_ratio"].to_numpy(), "t_ratio"
                )
                row["gaussian_rss_p"] = monte_carlo_p(observed["rss"], gaussian["rss"].to_numpy(), "rss")
            summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    for scenario in scenarios:
        mask = summary["scenario"].eq(scenario)
        for column in ("classic_t_p", "classic_rss_p", "gaussian_t_p", "gaussian_rss_p"):
            values = summary.loc[mask, column]
            if values.notna().all():
                summary.loc[mask, f"{column}_holm"] = holm_adjust(values)
            else:
                summary.loc[mask, f"{column}_holm"] = np.nan

    summary["classic_gate"] = (
        (summary["classic_t_p_holm"] < 0.05) & (summary["classic_rss_p_holm"] < 0.05)
    )
    summary["gaussian_gate"] = (
        (summary["gaussian_t_p_holm"] < 0.05) & (summary["gaussian_rss_p_holm"] < 0.05)
    )
    summary["dual_null_gate"] = summary["classic_gate"] & summary["gaussian_gate"]

    primary_mask = summary["scenario"].eq("primary_nonoverlap_9")
    supported_k = summary.loc[primary_mask & summary["dual_null_gate"], "k"].astype(int).tolist()
    minimal_supported_k = min(supported_k) if supported_k else None

    # Evaluate algorithmic stability for the candidate tetrahedron without using
    # these restarts to choose the primary observed fit.
    primary_k4_pca = pca_map[("primary_nonoverlap_9", 4)]
    primary_k4_X = pc_map[("primary_nonoverlap_9", 4)]
    reference_fit = observed_map[("primary_nonoverlap_9", 4)]
    reference_Z_pc = reference_fit["Z"]
    reference_Z_lab = primary_k4_pca.inverse_transform(reference_Z_pc)

    restart_rows: list[dict] = []
    restart_vertices: list[np.ndarray] = []
    for restart in range(20):
        restart_seed = SEED + 700_000 + restart
        fit = fit_aa(primary_k4_X, k=4, seed=restart_seed)
        query_lab = primary_k4_pca.inverse_transform(fit["Z"])
        aligned_lab, distances = align_vertices(reference_Z_lab, query_lab)
        aligned_pc = primary_k4_pca.transform(aligned_lab)
        restart_vertices.append(aligned_pc)
        separation = np.median(distance_matrix(reference_Z_lab, reference_Z_lab)[np.triu_indices(4, 1)])
        restart_rows.append(
            {
                "restart": restart + 1,
                "seed": restart_seed,
                "t_ratio": fit["t_ratio"],
                "rss": fit["rss"],
                "variance_explained": fit["variance_explained"],
                "converged": fit["converged"],
                "mean_matched_distance_lab_z": float(distances.mean()),
                "max_matched_distance_lab_z": float(distances.max()),
                "mean_distance_relative_to_vertex_separation": float(distances.mean() / separation),
                "max_distance_relative_to_vertex_separation": float(distances.max() / separation),
            }
        )
    restart_frame = pd.DataFrame(restart_rows)
    restart_frame.to_csv(OUTPUT_DIR / "pareto_k4_restart_stability.csv", index=False)
    restart_stable = bool(
        restart_frame["max_distance_relative_to_vertex_separation"].quantile(0.90) < 0.20
    )

    tetrahedron_promotable_prebootstrap = bool(
        4 in supported_k and not any(k < 4 for k in supported_k) and restart_stable
    )

    # Calendar validation/test transportability in the frozen primary K=4 space.
    split_rows: list[dict] = []
    for split_name, frame in (("train", train), ("validation", validation), ("test", test)):
        X_standardized = standardized_map["primary_nonoverlap_9"][split_name]
        X_pc = primary_k4_pca.transform(X_standardized).astype(np.float32)
        bary = barycentric_coordinates(X_pc, reference_Z_pc)
        weights, error = project_to_simplex(X_pc, reference_Z_pc)
        split_rows.append(
            {
                "split": split_name,
                "n": len(frame),
                "inside_simplex_pct_tolerance_0.01": float(100 * np.mean(np.min(bary, axis=1) >= -0.01)),
                "median_projection_distance": float(np.median(error)),
                "p90_projection_distance": float(np.quantile(error, 0.90)),
                "mean_projection_distance": float(np.mean(error)),
                "mean_max_archetype_weight": float(np.mean(np.max(weights, axis=1))),
            }
        )
    split_frame = pd.DataFrame(split_rows)
    split_frame.to_csv(OUTPUT_DIR / "pareto_k4_calendar_transport.csv", index=False)

    # Save observed archetype coordinates for every candidate geometry.
    archetype_rows: list[dict] = []
    for scenario, variables in scenarios.items():
        for k in K_VALUES:
            pca = pca_map[(scenario, k)]
            fit = observed_map[(scenario, k)]
            Z_lab = pca.inverse_transform(fit["Z"])
            for archetype in range(k):
                row = {
                    "scenario": scenario,
                    "k": k,
                    "archetype": archetype + 1,
                }
                for pc_index in range(k - 1):
                    row[f"PC{pc_index + 1}"] = fit["Z"][archetype, pc_index]
                for variable, value in zip(variables, Z_lab[archetype], strict=True):
                    row[f"z_{variable}"] = value
                archetype_rows.append(row)
    pd.DataFrame(archetype_rows).to_csv(OUTPUT_DIR / "pareto_archetypes.csv", index=False)

    bootstrap_frame = pd.DataFrame()
    bootstrap_stable = False
    if tetrahedron_promotable_prebootstrap:
        LOGGER.info("K=4 passed dual-null and restart gates; starting %d patient bootstraps", N_BOOTSTRAP)

        def bootstrap_batch(seeds: list[int]) -> list[dict]:
            rows: list[dict] = []
            X_full = standardized_map["primary_nonoverlap_9"]["train"]
            n = X_full.shape[0]
            for bootstrap_seed in seeds:
                rng = np.random.default_rng(bootstrap_seed)
                indices = rng.integers(0, n, size=n)
                X_boot = X_full[indices]
                pca_boot = PCA(n_components=3, svd_solver="full", random_state=bootstrap_seed)
                pc_boot = pca_boot.fit_transform(X_boot).astype(np.float32)
                fit = fit_aa(pc_boot, k=4, seed=bootstrap_seed + 500_000)
                query_lab = pca_boot.inverse_transform(fit["Z"])
                aligned_lab, distances = align_vertices(reference_Z_lab, query_lab)
                aligned_pc = primary_k4_pca.transform(aligned_lab)
                for archetype in range(4):
                    row = {
                        "bootstrap_seed": bootstrap_seed,
                        "archetype": archetype + 1,
                        "PC1": aligned_pc[archetype, 0],
                        "PC2": aligned_pc[archetype, 1],
                        "PC3": aligned_pc[archetype, 2],
                        "matched_distance_lab_z": distances[archetype],
                        "t_ratio": fit["t_ratio"],
                        "rss": fit["rss"],
                    }
                    for variable, value in zip(primary_variables, aligned_lab[archetype], strict=True):
                        row[f"z_{variable}"] = value
                    rows.append(row)
            return rows

        master = np.random.default_rng(SEED + 9_000_000)
        seeds = master.choice(np.arange(1, 2_000_000_000, dtype=np.int64), size=N_BOOTSTRAP, replace=False)
        batches = split_batches(seeds, N_JOBS)
        with parallel_backend("loky", inner_max_num_threads=1):
            with Parallel(n_jobs=N_JOBS, max_nbytes="20M", mmap_mode="r") as parallel:
                nested = parallel(delayed(bootstrap_batch)(batch) for batch in batches)
        bootstrap_frame = pd.DataFrame([row for batch in nested for row in batch])
        bootstrap_frame.to_csv(OUTPUT_DIR / "pareto_k4_bootstrap_vertices.csv", index=False)
        vertex_separation = np.median(
            distance_matrix(reference_Z_lab, reference_Z_lab)[np.triu_indices(4, 1)]
        )
        per_boot = bootstrap_frame.groupby("bootstrap_seed")["matched_distance_lab_z"].max()
        bootstrap_stable = bool(per_boot.quantile(0.90) / vertex_separation < 0.25)
        LOGGER.info(
            "Bootstrap vertex stability: p90 max displacement / separation = %.3f",
            per_boot.quantile(0.90) / vertex_separation,
        )

    tetrahedron_promotable = bool(tetrahedron_promotable_prebootstrap and bootstrap_stable)
    summary["minimal_supported_k_primary"] = minimal_supported_k
    summary["k4_restart_stable"] = restart_stable
    summary["k4_bootstrap_stable"] = bootstrap_stable
    summary["tetrahedron_promotable"] = tetrahedron_promotable
    summary.to_csv(OUTPUT_DIR / "pareto_gate_summary.csv", index=False)

    # Persist the frozen primary K=4 transform even if the tetrahedron fails;
    # it is useful for transparent diagnostic plots but cannot be named a true
    # tetrahedral state space unless tetrahedron_promotable is TRUE.
    np.savez_compressed(
        OUTPUT_DIR / "pareto_primary_k4_frozen_transform.npz",
        variables=np.asarray(primary_variables),
        winsor_lower=preprocessor_map["primary_nonoverlap_9"].lower,
        winsor_upper=preprocessor_map["primary_nonoverlap_9"].upper,
        median=preprocessor_map["primary_nonoverlap_9"].median,
        mean=preprocessor_map["primary_nonoverlap_9"].mean,
        scale=preprocessor_map["primary_nonoverlap_9"].scale,
        pca_mean=primary_k4_pca.mean_,
        pca_components=primary_k4_pca.components_,
        pca_explained_variance_ratio=primary_k4_pca.explained_variance_ratio_,
        archetypes_pc=reference_Z_pc,
        archetypes_lab_z=reference_Z_lab,
    )

    qa = {
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "hospital_total_n": int(len(hospital)),
        "hospital_unique_patient_n": int(hospital["patient_key"].nunique()),
        "calendar_train_n": int(len(train)),
        "calendar_validation_n": int(len(validation)),
        "calendar_test_n": int(len(test)),
        "primary_variable_n": len(primary_variables),
        "sensitivity_variable_n": len(full_variables),
        "candidate_k": list(K_VALUES),
        "classic_permutations_primary": N_PERM,
        "gaussian_permutations_primary": N_GAUSSIAN,
        "classic_permutations_sensitivity": N_SENSITIVITY,
        "bootstrap_replicates_requested": N_BOOTSTRAP,
        "bootstrap_replicates_completed": int(bootstrap_frame["bootstrap_seed"].nunique()) if len(bootstrap_frame) else 0,
        "coreset_size": CORESET_SIZE,
        "all_hospital_rows_used_for_preprocessing_pca_hull_and_rss": True,
        "archetype_optimizer": "ParTIpy 0.2.0 AA projected gradients with standard weighted coreset",
        "partipy_wheel_sha256": expected_wheel_hash,
        "primary_supported_k": supported_k,
        "minimal_supported_k": minimal_supported_k,
        "k4_restart_stable": restart_stable,
        "k4_bootstrap_stable": bootstrap_stable,
        "tetrahedron_promotable": tetrahedron_promotable,
        "checks_passed": True,
        "elapsed_seconds": time.time() - started,
    }
    with (OUTPUT_DIR / "pareto_gate_qa.json").open("w", encoding="utf-8") as handle:
        json.dump(qa, handle, indent=2, ensure_ascii=False)

    LOGGER.info(
        "Geometry gate complete in %.1f min | supported k=%s | tetrahedron promotable=%s",
        (time.time() - started) / 60,
        supported_k,
        tetrahedron_promotable,
    )


if __name__ == "__main__":
    main()
