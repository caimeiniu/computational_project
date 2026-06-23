#!/usr/bin/env python3
"""
Extended raw-PCA training-site selection experiment up to 1000 labels.

All outputs are written inside:
    project/point_selection/outputs/extended_budget/

This script does not modify the main-project results, figures, PPT, or CSV files.
It is resumable: existing method/budget/repeat rows in extended_budget_results.csv
are skipped unless --overwrite is passed.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from sklearn.metrics import pairwise_distances_argmin_min, r2_score


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PACKAGE_ROOT / "outputs" / "extended_budget"
DATA_PATH = PACKAGE_ROOT / "data" / "al_mg_pca_deltaE_verified.csv"
BASELINE_PATH = (
    PACKAGE_ROOT / "results" / "reference" / "all_data_10pc_linear_regression_baseline.csv"
)
RESULTS_PATH = OUTPUT_DIR / "extended_budget_results.csv"
SUMMARY_PATH = OUTPUT_DIR / "extended_budget_summary.csv"
REPORT_PATH = OUTPUT_DIR / "extended_budget_report.md"
SELECTED_DIR = OUTPUT_DIR / "selected_indices"
FIGURES_DIR = OUTPUT_DIR / "figures"

FEATURE_COLS = [f"pc{i}" for i in range(1, 11)]
TARGET_COL = "deltaE_kJmol"
BUDGETS = [20, 50, 70, 100, 150, 200, 500, 1000]
REPRESENTATIVE_BUDGETS = {100, 500, 1000}

KB_KJMOL_K = 0.008314462618
DEFAULT_T = 600.0
DEFAULT_X_TOT = 0.05
DEFAULT_F_GB = 0.10

METHOD_LABELS = {
    "random": "Random",
    "kmeans_raw_pca": "k-means",
    "active_bootstrap_raw_pca": "Bootstrap active",
    "tail_aware_active_raw_pca": "Tail-aware active",
    "xgb_aware_active_raw_pca": "XGB-aware active",
}
METHOD_ORDER = [
    "random",
    "kmeans_raw_pca",
    "active_bootstrap_raw_pca",
    "tail_aware_active_raw_pca",
    "xgb_aware_active_raw_pca",
]
METHOD_COLORS = {
    "random": "#6b7280",
    "kmeans_raw_pca": "#f59e0b",
    "active_bootstrap_raw_pca": "#2563eb",
    "tail_aware_active_raw_pca": "#16a34a",
    "xgb_aware_active_raw_pca": "#7c3aed",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overwrite", action="store_true", help="Rerun and overwrite extended experiment rows.")
    parser.add_argument("--budgets", nargs="*", type=int, default=BUDGETS, help="Optional budget subset.")
    return parser.parse_args()


def stable_sigmoid(z: np.ndarray) -> np.ndarray:
    out = np.empty_like(z, dtype=float)
    positive = z >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-z[positive]))
    exp_z = np.exp(z[~positive])
    out[~positive] = exp_z / (1.0 + exp_z)
    return out


def compute_occupation(deltaE: np.ndarray, X_c: float, T: float = DEFAULT_T) -> np.ndarray:
    deltaE = np.asarray(deltaE, dtype=float)
    X_c = float(np.clip(X_c, 1e-12, 1.0 - 1e-12))
    log_factor = np.log((1.0 - X_c) / X_c)
    z = -(log_factor + deltaE / (KB_KJMOL_K * T))
    return stable_sigmoid(z)


def solve_Xc(
    deltaE: np.ndarray,
    X_tot: float = DEFAULT_X_TOT,
    f_GB: float = DEFAULT_F_GB,
    T: float = DEFAULT_T,
) -> float:
    deltaE = np.asarray(deltaE, dtype=float)
    lower = 1e-12
    upper = float(X_tot) * (1.0 - 1e-12)

    def balance(X_c: float) -> float:
        p = compute_occupation(deltaE, X_c, T=T)
        return (1.0 - f_GB) * X_c + f_GB * float(np.mean(p)) - X_tot

    f_lower = balance(lower)
    f_upper = balance(upper)
    if f_lower > 0 or f_upper < 0:
        raise ValueError(
            "Mass-balance root is not bracketed in (1e-12, X_tot): "
            f"f(lower)={f_lower:.6g}, f(upper)={f_upper:.6g}"
        )
    try:
        from scipy.optimize import brentq

        return float(brentq(balance, lower, upper, xtol=1e-13, rtol=1e-12, maxiter=200))
    except ImportError:
        lo, hi = lower, upper
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if balance(mid) > 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-13:
                break
        return float(0.5 * (lo + hi))


def compute_xgb(
    deltaE: np.ndarray,
    X_tot: float = DEFAULT_X_TOT,
    f_GB: float = DEFAULT_F_GB,
    T: float = DEFAULT_T,
) -> tuple[float, float]:
    X_c = solve_Xc(deltaE, X_tot=X_tot, f_GB=f_GB, T=T)
    p = compute_occupation(deltaE, X_c, T=T)
    return float(np.mean(p)), X_c


def load_data() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    df = pd.read_csv(DATA_PATH)
    missing = [col for col in FEATURE_COLS + [TARGET_COL] if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if df[FEATURE_COLS + [TARGET_COL]].isna().any().any():
        raise ValueError("NaN values found in feature or target columns.")
    X = df[FEATURE_COLS].to_numpy(dtype=float)
    y = df[TARGET_COL].to_numpy(dtype=float)
    print(f"Loaded dataset: {DATA_PATH}")
    print(f"Dataset shape: {df.shape}")
    print(f"Features: {FEATURE_COLS}")
    print(f"Target: {TARGET_COL}")
    print("StandardScaler used: No")
    return df, X, y


def load_baseline() -> dict[str, float]:
    baseline = pd.read_csv(BASELINE_PATH).iloc[0]
    return {
        "mae": float(baseline["mae"]),
        "tail_mae": float(baseline["tail_mae"]),
        "xgb_abs_error": float(baseline["xgb_abs_error"]),
        "rmse": float(baseline["rmse"]),
        "r2": float(baseline["r2"]),
    }


def load_existing_results(overwrite: bool) -> pd.DataFrame:
    if overwrite or not RESULTS_PATH.exists():
        return pd.DataFrame()
    return pd.read_csv(RESULTS_PATH)


def save_results(results: pd.DataFrame) -> None:
    results = results.sort_values(["budget", "method", "repeat"], kind="mergesort")
    results.to_csv(RESULTS_PATH, index=False)


def completed_keys(results: pd.DataFrame) -> set[tuple[str, int, int]]:
    if results.empty:
        return set()
    return set(zip(results["method"].astype(str), results["budget"].astype(int), results["repeat"].astype(int)))


def selected_indices_frame(df: pd.DataFrame, selected: np.ndarray) -> pd.DataFrame:
    out = pd.DataFrame({"site_index": selected.astype(int)})
    if "row_id" in df.columns:
        out["row_id"] = df.iloc[selected]["row_id"].to_numpy()
    if "site_id_if_verified" in df.columns:
        out["site_id_if_verified"] = df.iloc[selected]["site_id_if_verified"].to_numpy()
    return out


def save_representative_indices(df: pd.DataFrame, selected: np.ndarray, method: str, budget: int, repeat: int) -> None:
    if budget not in REPRESENTATIVE_BUDGETS or repeat != 0:
        return
    if method == "kmeans_raw_pca":
        name = f"selected_indices_{method}_{budget}.csv"
    else:
        name = f"selected_indices_{method}_{budget}_repeat{repeat}.csv"
    selected_indices_frame(df, selected).to_csv(SELECTED_DIR / name, index=False)


def select_random(X: np.ndarray, n_labels: int, random_state: int) -> np.ndarray:
    rng = np.random.default_rng(random_state)
    return rng.choice(len(X), size=n_labels, replace=False).astype(int)


def select_kmeans(X: np.ndarray, n_labels: int, random_state: int = 42) -> np.ndarray:
    kmeans = KMeans(n_clusters=n_labels, random_state=random_state, n_init=10)
    kmeans.fit(X)
    nearest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, X)

    selected: list[int] = []
    seen = set()
    for idx in nearest.tolist():
        idx = int(idx)
        if idx not in seen:
            selected.append(idx)
            seen.add(idx)

    if len(selected) < n_labels:
        min_dist = kmeans.transform(X).min(axis=1)
        for idx in np.argsort(min_dist):
            idx = int(idx)
            if idx not in seen:
                selected.append(idx)
                seen.add(idx)
            if len(selected) == n_labels:
                break

    if len(selected) < n_labels:
        rng = np.random.default_rng(random_state)
        remaining = np.array([idx for idx in range(len(X)) if idx not in seen], dtype=int)
        selected.extend(rng.choice(remaining, size=n_labels - len(selected), replace=False).tolist())

    return np.array(selected[:n_labels], dtype=int)


def bootstrap_predictions(
    X: np.ndarray,
    y: np.ndarray,
    labeled: np.ndarray,
    query: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Bootstrap LinearRegression models trained only on labeled points."""
    n_features = X.shape[1]
    coefs = np.empty((n_bootstrap, n_features), dtype=float)
    intercepts = np.empty(n_bootstrap, dtype=float)

    for b in range(n_bootstrap):
        sample = rng.choice(labeled, size=len(labeled), replace=True)
        model = LinearRegression()
        model.fit(X[sample], y[sample])
        coefs[b] = model.coef_
        intercepts[b] = model.intercept_

    predictions = X[query] @ coefs.T + intercepts
    return predictions.mean(axis=1), predictions.std(axis=1)


def active_schedule(final_n: int) -> tuple[int, int, int]:
    if final_n <= 200:
        return 50, 20, 30
    return 50, 100, 10


def select_active_bootstrap(
    X: np.ndarray,
    y: np.ndarray,
    final_n: int,
    random_state: int,
) -> np.ndarray:
    initial_n, batch_size, n_bootstrap = active_schedule(final_n)
    rng = np.random.default_rng(random_state)
    if final_n <= 50:
        return rng.choice(len(X), size=final_n, replace=False).astype(int)

    labeled = rng.choice(len(X), size=initial_n, replace=False).astype(int)
    labeled_set = set(labeled.tolist())
    while len(labeled) < final_n:
        unlabeled = np.array([idx for idx in range(len(X)) if idx not in labeled_set], dtype=int)
        _, sigma = bootstrap_predictions(X, y, labeled, unlabeled, n_bootstrap, rng)
        acquire_n = min(batch_size, final_n - len(labeled))
        new_indices = unlabeled[np.argsort(sigma)[-acquire_n:]]
        labeled = np.concatenate([labeled, new_indices])
        labeled_set.update(new_indices.tolist())
        print(f"    active_bootstrap acquisition: labeled={len(labeled)}/{final_n}")
    return labeled.astype(int)


def select_tail_aware_active(
    X: np.ndarray,
    y: np.ndarray,
    final_n: int,
    random_state: int,
    lambda_uncertainty: float = 1.0,
) -> np.ndarray:
    initial_n, batch_size, n_bootstrap = active_schedule(final_n)
    rng = np.random.default_rng(random_state)
    if final_n <= 50:
        return rng.choice(len(X), size=final_n, replace=False).astype(int)

    labeled = rng.choice(len(X), size=initial_n, replace=False).astype(int)
    labeled_set = set(labeled.tolist())
    while len(labeled) < final_n:
        unlabeled = np.array([idx for idx in range(len(X)) if idx not in labeled_set], dtype=int)
        mu, sigma = bootstrap_predictions(X, y, labeled, unlabeled, n_bootstrap, rng)
        score = -mu + lambda_uncertainty * sigma
        acquire_n = min(batch_size, final_n - len(labeled))
        new_indices = unlabeled[np.argsort(score)[-acquire_n:]]
        labeled = np.concatenate([labeled, new_indices])
        labeled_set.update(new_indices.tolist())
        print(f"    tail_aware acquisition: labeled={len(labeled)}/{final_n}")
    return labeled.astype(int)


def select_xgb_aware_active(
    X: np.ndarray,
    y: np.ndarray,
    final_n: int,
    random_state: int,
) -> np.ndarray:
    initial_n, batch_size, n_bootstrap = active_schedule(final_n)
    rng = np.random.default_rng(random_state)
    if final_n <= 50:
        return rng.choice(len(X), size=final_n, replace=False).astype(int)

    all_indices = np.arange(len(X), dtype=int)
    labeled = rng.choice(len(X), size=initial_n, replace=False).astype(int)
    labeled_set = set(labeled.tolist())
    while len(labeled) < final_n:
        unlabeled = np.array([idx for idx in range(len(X)) if idx not in labeled_set], dtype=int)
        mu_all, sigma_all = bootstrap_predictions(X, y, labeled, all_indices, n_bootstrap, rng)
        xgb_pred, x_c = compute_xgb(mu_all)
        del xgb_pred
        occupation_unlabeled = compute_occupation(mu_all[unlabeled], x_c)
        sensitivity = occupation_unlabeled * (1.0 - occupation_unlabeled)
        score = sigma_all[unlabeled] * sensitivity
        acquire_n = min(batch_size, final_n - len(labeled))
        new_indices = unlabeled[np.argsort(score)[-acquire_n:]]
        labeled = np.concatenate([labeled, new_indices])
        labeled_set.update(new_indices.tolist())
        print(f"    xgb_aware acquisition: labeled={len(labeled)}/{final_n}")
    return labeled.astype(int)


def evaluate_selection(
    X: np.ndarray,
    y: np.ndarray,
    selected: np.ndarray,
    xgb_true: float,
) -> tuple[dict[str, float], np.ndarray]:
    model = LinearRegression()
    model.fit(X[selected], y[selected])
    y_pred = model.predict(X)
    residual = y_pred - y
    tail_cutoff = np.percentile(y, 10)
    tail_mask = y <= tail_cutoff
    xgb_pred, _ = compute_xgb(y_pred)
    metrics = {
        "n_selected": int(len(selected)),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "r2": float(r2_score(y, y_pred)),
        "tail_mae": float(np.mean(np.abs(residual[tail_mask]))),
        "xgb_true": float(xgb_true),
        "xgb_pred": float(xgb_pred),
        "xgb_abs_error": float(abs(xgb_pred - xgb_true)),
    }
    return metrics, y_pred


def random_repeats_for_budget(budget: int) -> int:
    return 20 if budget <= 200 else 10


def active_repeats_for_budget(budget: int) -> int:
    return 10 if budget <= 200 else 5


def method_runs_for_budget(budget: int) -> list[tuple[str, int, int]]:
    runs: list[tuple[str, int, int]] = []
    for repeat in range(random_repeats_for_budget(budget)):
        runs.append(("random", repeat, 100000 + budget * 100 + repeat))
    runs.append(("kmeans_raw_pca", 0, 42))
    for repeat in range(active_repeats_for_budget(budget)):
        runs.append(("active_bootstrap_raw_pca", repeat, 200000 + budget * 100 + repeat))
        runs.append(("tail_aware_active_raw_pca", repeat, 300000 + budget * 100 + repeat))
        runs.append(("xgb_aware_active_raw_pca", repeat, 400000 + budget * 100 + repeat))
    return runs


def select_for_method(method: str, X: np.ndarray, y: np.ndarray, budget: int, seed: int) -> np.ndarray:
    if method == "random":
        return select_random(X, budget, seed)
    if method == "kmeans_raw_pca":
        return select_kmeans(X, budget, seed)
    if method == "active_bootstrap_raw_pca":
        return select_active_bootstrap(X, y, budget, seed)
    if method == "tail_aware_active_raw_pca":
        return select_tail_aware_active(X, y, budget, seed, lambda_uncertainty=1.0)
    if method == "xgb_aware_active_raw_pca":
        return select_xgb_aware_active(X, y, budget, seed)
    raise ValueError(f"Unknown method: {method}")


def append_result(existing: pd.DataFrame, row: dict) -> pd.DataFrame:
    return pd.concat([existing, pd.DataFrame([row])], ignore_index=True)


def summarize_results(results: pd.DataFrame) -> pd.DataFrame:
    if results.empty:
        return pd.DataFrame()
    summary = (
        results.groupby(["method", "budget"], as_index=False)
        .agg(
            mae_mean=("mae", "mean"),
            mae_std=("mae", "std"),
            tail_mae_mean=("tail_mae", "mean"),
            tail_mae_std=("tail_mae", "std"),
            xgb_abs_error_mean=("xgb_abs_error", "mean"),
            xgb_abs_error_std=("xgb_abs_error", "std"),
            rmse_mean=("rmse", "mean"),
            rmse_std=("rmse", "std"),
            r2_mean=("r2", "mean"),
            r2_std=("r2", "std"),
            runtime_seconds_mean=("runtime_seconds", "mean"),
            runtime_seconds_std=("runtime_seconds", "std"),
            n_runs=("mae", "size"),
        )
        .sort_values(["budget", "method"], kind="mergesort")
    )
    return summary


def save_summary(results: pd.DataFrame) -> pd.DataFrame:
    summary = summarize_results(results)
    if not summary.empty:
        summary.to_csv(SUMMARY_PATH, index=False)
    return summary


def plot_metric(
    ax: plt.Axes,
    summary: pd.DataFrame,
    metric_mean: str,
    metric_std: str,
    baseline_value: float,
    baseline_label: str,
    ylabel: str,
    title: str,
) -> None:
    for method in METHOD_ORDER:
        data = summary[summary["method"] == method].sort_values("budget")
        if data.empty:
            continue
        x = data["budget"].to_numpy(dtype=float)
        y = data[metric_mean].to_numpy(dtype=float)
        color = METHOD_COLORS[method]
        ax.plot(x, y, marker="o", markersize=4.5, linewidth=2.0, color=color, label=METHOD_LABELS[method])
        if metric_std in data.columns:
            std = data[metric_std].fillna(0.0).to_numpy(dtype=float)
            if np.any(std > 0):
                ax.fill_between(x, y - std, y + std, color=color, alpha=0.14, linewidth=0)
    ax.axhline(baseline_value, color="#111827", linestyle="--", linewidth=1.7, label=baseline_label)
    ax.set_xscale("log")
    ax.set_xticks(BUDGETS)
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax.set_xlabel("Number of labeled sites")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, color="#d8dee9", linewidth=0.7, alpha=0.85)


def make_figures(summary: pd.DataFrame, baseline: dict[str, float]) -> list[Path]:
    if summary.empty:
        return []
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    specs = [
        (
            "extended_mae_vs_labels_to_1000.png",
            "mae_mean",
            "mae_std",
            baseline["mae"],
            f"All-data 10D LR baseline = {baseline['mae']:.3f} kJ/mol",
            "Global MAE in DeltaE prediction (kJ/mol)",
            "Global MAE vs label budget",
        ),
        (
            "extended_tail_mae_vs_labels_to_1000.png",
            "tail_mae_mean",
            "tail_mae_std",
            baseline["tail_mae"],
            f"All-data 10D LR baseline = {baseline['tail_mae']:.3f} kJ/mol",
            "Tail MAE, lowest-energy 10% (kJ/mol)",
            "Low-energy tail MAE vs label budget",
        ),
        (
            "extended_xgb_error_vs_labels_to_1000.png",
            "xgb_abs_error_mean",
            "xgb_abs_error_std",
            baseline["xgb_abs_error"],
            f"All-data 10D LR baseline = {baseline['xgb_abs_error']:.6f}",
            "XGB_abs_error",
            "XGB error vs label budget",
        ),
    ]
    out_paths: list[Path] = []
    for filename, mean_col, std_col, base_value, base_label, ylabel, title in specs:
        fig, ax = plt.subplots(figsize=(8.0, 5.4), dpi=300)
        plot_metric(ax, summary, mean_col, std_col, base_value, base_label, ylabel, title)
        ax.set_title(f"{title}\nExtended raw-PCA learning curves to 1000 labels")
        ax.legend(frameon=False, fontsize=9, loc="best")
        fig.tight_layout()
        out_path = FIGURES_DIR / filename
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        out_paths.append(out_path)

    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.2), dpi=300, constrained_layout=True)
    for ax, (_, mean_col, std_col, base_value, base_label, ylabel, title) in zip(axes, specs):
        plot_metric(ax, summary, mean_col, std_col, base_value, base_label, ylabel, title)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Extended raw-PCA learning curves to 1000 labels", fontsize=16, y=1.04)
    combined_path = FIGURES_DIR / "combined_extended_learning_curves_to_1000.png"
    fig.savefig(combined_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    out_paths.append(combined_path)
    return out_paths


def format_metric(summary: pd.DataFrame, method: str, budget: int, metric: str) -> str:
    row = summary[(summary["method"] == method) & (summary["budget"] == budget)]
    if row.empty:
        return "not completed"
    mean = float(row.iloc[0][f"{metric}_mean"])
    std_col = f"{metric}_std"
    std = float(row.iloc[0][std_col]) if std_col in row.columns and not pd.isna(row.iloc[0][std_col]) else 0.0
    return f"{mean:.6f} +/- {std:.6f}"


def write_report(summary: pd.DataFrame, baseline: dict[str, float], incomplete: list[str]) -> None:
    lines = [
        "# Extended Raw-PCA Label-Budget Experiment To 1000 Labels",
        "",
        "## A. What Was Done",
        "",
        "- Used raw `pc1` ... `pc10` from the verified Al-Mg dataset.",
        "- Did not use `StandardScaler`.",
        "- Did not fit a new PCA.",
        "- Used `LinearRegression()` as the final model for every method.",
        "- Extended label budgets to `[20, 50, 70, 100, 150, 200, 500, 1000]`.",
        "- Read the all-data 10D raw PCA `LinearRegression` baseline as a horizontal reference line.",
        "",
        "The all-data baseline is an in-sample reference for the same dataset and feature/model setup. Because `LinearRegression()` minimizes squared error rather than MAE or XGB_abs_error, it is not a guaranteed lower bound for every reported metric.",
        "",
        "## B. Results Summary",
        "",
        "| Method | Budget | Global MAE | Tail MAE | XGB_abs_error |",
        "|---|---:|---:|---:|---:|",
    ]
    for method in METHOD_ORDER:
        for budget in [100, 500, 1000]:
            row = summary[(summary["method"] == method) & (summary["budget"] == budget)]
            if row.empty:
                lines.append(f"| {METHOD_LABELS[method]} | {budget} | not completed | not completed | not completed |")
            else:
                r = row.iloc[0]
                lines.append(
                    f"| {METHOD_LABELS[method]} | {budget} | "
                    f"{r['mae_mean']:.6f} | {r['tail_mae_mean']:.6f} | {r['xgb_abs_error_mean']:.8f} |"
                )

    lines.extend(
        [
            "",
            "## C. Convergence To All-Data Baseline",
            "",
            f"All-data global MAE baseline: `{baseline['mae']:.6f}` kJ/mol.",
            f"All-data tail MAE baseline: `{baseline['tail_mae']:.6f}` kJ/mol.",
            f"All-data XGB_abs_error baseline: `{baseline['xgb_abs_error']:.8f}`.",
            "",
        ]
    )

    if not summary.empty and (summary["budget"] == 1000).any():
        at1000 = summary[summary["budget"] == 1000].copy()
        best_mae = at1000.loc[at1000["mae_mean"].idxmin()]
        best_tail = at1000.loc[at1000["tail_mae_mean"].idxmin()]
        best_xgb = at1000.loc[at1000["xgb_abs_error_mean"].idxmin()]
        lines.extend(
            [
                f"At 1000 labels, the best global MAE method is `{best_mae['method']}` with MAE `{best_mae['mae_mean']:.6f}`.",
                f"At 1000 labels, the best tail MAE method is `{best_tail['method']}` with tail MAE `{best_tail['tail_mae_mean']:.6f}`.",
                f"At 1000 labels, the best XGB_abs_error method is `{best_xgb['method']}` with XGB_abs_error `{best_xgb['xgb_abs_error_mean']:.8f}`.",
                "",
                "Distance to the all-data global MAE baseline at 1000 labels:",
                "",
                "| Method | MAE gap to all-data baseline | Tail MAE gap | XGB_abs_error gap |",
                "|---|---:|---:|---:|",
            ]
        )
        for _, r in at1000.sort_values("method").iterrows():
            lines.append(
                f"| {r['method']} | {r['mae_mean'] - baseline['mae']:.6f} | "
                f"{r['tail_mae_mean'] - baseline['tail_mae']:.6f} | "
                f"{r['xgb_abs_error_mean'] - baseline['xgb_abs_error']:.8f} |"
            )
    else:
        lines.append("No 1000-label results were completed yet.")

    lines.extend(
        [
            "",
            "## D. Bottleneck Interpretation",
            "",
            "If the extended curves flatten near the all-data global MAE line, then global-error improvement is likely limited by the 10D raw PCA + `LinearRegression` feature/model setup rather than selection alone. However, the all-data OLS fit minimizes squared error, not MAE, so a selected subset can occasionally have slightly lower MAE while still being the same linear model family. Tail MAE should be interpreted separately, because low-energy sites can remain difficult even when global MAE is close to the all-data fit.",
            "",
            "XGB_abs_error is nonlinear in the predicted energy spectrum and may not behave monotonically with global MAE. A method can have lower MAE but worse XGB_abs_error if the spectrum errors occur in occupation-sensitive regions.",
            "",
            "## E. Important Caution",
            "",
            "- The all-data baseline trains and evaluates on the same 82,645 sites.",
            "- It is an in-sample feature/model reference for the current dataset, not an independent generalization benchmark.",
            "- Because sklearn `LinearRegression()` optimizes squared error, the all-data fit is not guaranteed to be the best possible MAE or XGB_abs_error.",
            "- XGB_abs_error should be discussed as a downstream physical metric, not as a simple monotonic transformation of MAE.",
            "",
            "## F. Recommended Presentation Wording",
            "",
            "The extended-budget curves show how each raw-PCA selection strategy approaches the all-data 10D `LinearRegression` reference as the label budget increases. This all-data line is not an independent test result; it is an in-sample reference for the same feature/model setup. Since `LinearRegression()` optimizes squared error rather than MAE, a selected subset can occasionally fall slightly below the all-data MAE line, so the line should be interpreted as a practical reference rather than a strict MAE bound. Tail MAE and XGB error should still be interpreted separately because they emphasize low-energy sites and nonlinear segregation behavior.",
            "",
            "## G. Completion Notes",
            "",
        ]
    )
    if incomplete:
        lines.append("Incomplete or skipped runs:")
        lines.extend(f"- {item}" for item in incomplete)
    else:
        lines.append("All requested runs completed.")
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_experiment(args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    SELECTED_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df, X, y = load_data()
    xgb_true, _ = compute_xgb(y)
    results = load_existing_results(args.overwrite)
    done = completed_keys(results)
    incomplete: list[str] = []

    if args.overwrite and RESULTS_PATH.exists():
        RESULTS_PATH.unlink()

    for budget in args.budgets:
        print(f"\n=== Budget {budget} ===")
        for method, repeat, seed in method_runs_for_budget(budget):
            key = (method, int(budget), int(repeat))
            if key in done:
                print(f"Skipping completed: method={method}, budget={budget}, repeat={repeat}")
                continue
            print(f"Running method={method}, budget={budget}, repeat={repeat}, seed={seed}")
            start = time.perf_counter()
            try:
                selected = select_for_method(method, X, y, budget, seed)
                metrics, _ = evaluate_selection(X, y, selected, xgb_true)
                runtime = time.perf_counter() - start
                row = {
                    "method": method,
                    "budget": int(budget),
                    "repeat": int(repeat),
                    "seed": int(seed),
                    "n_selected": metrics["n_selected"],
                    "mae": metrics["mae"],
                    "rmse": metrics["rmse"],
                    "r2": metrics["r2"],
                    "tail_mae": metrics["tail_mae"],
                    "xgb_true": metrics["xgb_true"],
                    "xgb_pred": metrics["xgb_pred"],
                    "xgb_abs_error": metrics["xgb_abs_error"],
                    "runtime_seconds": float(runtime),
                }
                results = append_result(results, row)
                save_results(results)
                save_summary(results)
                save_representative_indices(df, selected, method, budget, repeat)
                done.add(key)
                print(
                    f"  done in {runtime:.1f}s | MAE={metrics['mae']:.4f}, "
                    f"tail={metrics['tail_mae']:.4f}, XGBerr={metrics['xgb_abs_error']:.6f}"
                )
            except KeyboardInterrupt:
                incomplete.append(f"Interrupted during {method}, budget={budget}, repeat={repeat}")
                save_results(results)
                raise
            except Exception as exc:
                msg = f"FAILED {method}, budget={budget}, repeat={repeat}: {exc}"
                print(msg)
                incomplete.append(msg)
                save_results(results)
                save_summary(results)
    return results, incomplete


def print_console_summary(summary: pd.DataFrame, baseline: dict[str, float], incomplete: list[str]) -> None:
    completed = summary.groupby("method")["budget"].apply(lambda s: sorted(set(map(int, s)))).to_dict() if not summary.empty else {}
    print("\nCompleted methods and budgets:")
    for method in METHOD_ORDER:
        print(f"  {method}: {completed.get(method, [])}")
    if incomplete:
        print("\nIncomplete/skipped runs:")
        for item in incomplete:
            print(f"  {item}")
    else:
        print("\nIncomplete/skipped runs: none")
    print("\nAll-data baseline values:")
    print(f"  global MAE = {baseline['mae']:.6f} kJ/mol")
    print(f"  tail MAE = {baseline['tail_mae']:.6f} kJ/mol")
    print(f"  XGB_abs_error = {baseline['xgb_abs_error']:.8f}")

    at1000 = summary[summary["budget"] == 1000] if not summary.empty else pd.DataFrame()
    if not at1000.empty:
        best_mae = at1000.loc[at1000["mae_mean"].idxmin()]
        best_tail = at1000.loc[at1000["tail_mae_mean"].idxmin()]
        best_xgb = at1000.loc[at1000["xgb_abs_error_mean"].idxmin()]
        print("\nBest methods at 1000 labels:")
        print(f"  global MAE: {best_mae['method']} ({best_mae['mae_mean']:.6f})")
        print(f"  tail MAE: {best_tail['method']} ({best_tail['tail_mae_mean']:.6f})")
        print(f"  XGB_abs_error: {best_xgb['method']} ({best_xgb['xgb_abs_error_mean']:.8f})")
    else:
        print("\nBest methods at 1000 labels: no 1000-label rows completed")
    print(f"\nReport path: {REPORT_PATH}")


def main() -> None:
    args = parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    baseline = load_baseline()
    incomplete: list[str] = []
    try:
        results, incomplete = run_experiment(args)
    except KeyboardInterrupt:
        print("\nInterrupted by user. Partial results have been saved.")
        results = load_existing_results(overwrite=False)
    summary = save_summary(results)
    figure_paths = make_figures(summary, baseline)
    write_report(summary, baseline, incomplete)
    print("\nGenerated figures:")
    for path in figure_paths:
        print(f"  {path}")
    print_console_summary(summary, baseline, incomplete)


if __name__ == "__main__":
    main()
