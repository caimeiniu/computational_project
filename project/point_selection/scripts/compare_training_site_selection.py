#!/usr/bin/env python3
"""
Compare raw-PCA training-site selection strategies for Al-Mg GB segregation energies.

Input:
    al_mg_pca_deltaE_verified.csv

Outputs:
    selection_results.csv
    selected_indices_random_100.csv
    selected_indices_kmeans_raw_pca_100.csv
    selected_indices_active_bootstrap_raw_pca_100.csv
    selected_indices_tail_aware_active_raw_pca_100.csv
    selected_indices_xgb_aware_active_raw_pca_100.csv
    predictions_kmeans_raw_pca_100.csv
    predictions_active_bootstrap_raw_pca_100.csv
    predictions_tail_aware_active_raw_pca_100.csv
    predictions_xgb_aware_active_raw_pca_100.csv
    figures/mae_vs_labels.png
    figures/tail_mae_vs_labels.png
    figures/xgb_error_vs_labels.png
    figures/pca_selected_points.png
    figures/pca_selected_points_all_methods.png
    figures/predicted_vs_true_100labels.png
    figures/predicted_vs_true_100labels_all_methods.png

The predictor is always LinearRegression. Selection methods may choose which
points are labeled, but final evaluation is always over all GB sites.
All methods use raw original-notebook-style PCA coordinates without StandardScaler.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.linear_model import LinearRegression
from sklearn.metrics import pairwise_distances_argmin_min, r2_score


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BUDGETS = [20, 50, 70, 90, 100, 150, 200]
KB_KJMOL_K = 0.008314462618
DEFAULT_T = 600.0
DEFAULT_X_TOT = 0.05
DEFAULT_F_GB = 0.10
TARGET_PATTERNS = [
    r"^deltae($|_)",
    r"^delta_e($|_)",
    r"^delta.*e",
    r"^segregation.*energy",
    r"^eseg($|_)",
    r"^e_seg($|_)",
    r"^target$",
]


@dataclass
class Evaluation:
    mae: float
    rmse: float
    r2: float
    tail_mae: float
    xgb_true: float
    xgb_pred: float
    xgb_abs_error: float


def normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")


def identify_target_column(df: pd.DataFrame) -> str:
    """Pick the segregation-energy target column or fail loudly."""
    normalized = {col: normalize_name(col) for col in df.columns}
    candidates = []
    for col, norm in normalized.items():
        if any(re.search(pattern, norm) for pattern in TARGET_PATTERNS):
            candidates.append(col)

    # The verified file uses deltaE_kJmol; prefer unit-bearing deltaE columns.
    if len(candidates) > 1:
        preferred = [col for col in candidates if "kj" in normalized[col] or "mol" in normalized[col]]
        if len(preferred) == 1:
            return preferred[0]
        print("Ambiguous target candidates:", candidates)
        raise ValueError("Could not uniquely identify the target column.")
    if not candidates:
        numeric = df.select_dtypes(include=[np.number]).columns.tolist()
        print("No target column matched known names.")
        print("Numeric candidate columns:", numeric)
        raise ValueError("Could not identify target column.")
    return candidates[0]


def identify_feature_columns(df: pd.DataFrame, target_col: str) -> list[str]:
    """Prefer pc1..pc10-like names; otherwise infer numeric non-target columns."""
    pc_matches: list[tuple[int, str]] = []
    for col in df.columns:
        match = re.match(r"^pc[_\s-]?(\d+)$", col, flags=re.IGNORECASE)
        if match:
            pc_matches.append((int(match.group(1)), col))

    if pc_matches:
        pc_matches.sort()
        return [col for _, col in pc_matches]

    id_like = {"row_id", "site_id", "site_id_if_verified", "atom_id", "id", "index"}
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    features = [
        col
        for col in numeric_cols
        if col != target_col and normalize_name(col) not in id_like and "site" not in normalize_name(col)
    ]
    if not features:
        print("Numeric columns:", numeric_cols)
        raise ValueError("No usable numeric feature columns found.")
    return features


def load_dataset(csv_path: Path) -> tuple[pd.DataFrame, list[str], str, np.ndarray, np.ndarray, np.ndarray]:
    df = pd.read_csv(csv_path)
    print("Loaded:", csv_path)
    print("Shape:", df.shape)
    print("Columns:", df.columns.tolist())
    print("First rows:")
    print(df.head())

    target_col = identify_target_column(df)
    feature_cols = identify_feature_columns(df, target_col)
    print(f"Feature columns ({len(feature_cols)}): {feature_cols}")
    print(f"Target column: {target_col}")

    X_raw = df[feature_cols].to_numpy(dtype=float)
    y = df[target_col].to_numpy(dtype=float)
    # Raw PCA / original-notebook-style mode:
    # selection, training, and evaluation all use pc1...pc10 directly.
    # Do not apply StandardScaler here.
    X = X_raw.copy()
    print("StandardScaler used: No; using raw pc1...pc10 coordinates.")
    return df, feature_cols, target_col, X_raw, X, y


def stable_sigmoid(z: np.ndarray) -> np.ndarray:
    """Numerically stable sigmoid for occupation probabilities."""
    out = np.empty_like(z, dtype=float)
    positive = z >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-z[positive]))
    exp_z = np.exp(z[~positive])
    out[~positive] = exp_z / (1.0 + exp_z)
    return out


def compute_occupation(deltaE: np.ndarray, X_c: float, T: float = DEFAULT_T) -> np.ndarray:
    """Fermi-Dirac-like GB occupation probability for each site.

    P_i = 1 / (1 + ((1 - X_c) / X_c) * exp(DeltaE_i / (kB*T)))
    """
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
    """Solve the mass-balance equation for the matrix concentration X_c."""
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
        # Fallback bisection keeps the script runnable if scipy is unavailable.
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
    """Return (X_c, X_GB) for a segregation energy spectrum."""
    X_c = solve_Xc(deltaE, X_tot=X_tot, f_GB=f_GB, T=T)
    occupation = compute_occupation(deltaE, X_c, T=T)
    return X_c, float(np.mean(occupation))


def select_random(X: np.ndarray, n_labels: int, random_state: int) -> np.ndarray:
    if n_labels > len(X):
        raise ValueError("n_labels cannot exceed number of samples")
    rng = np.random.default_rng(random_state)
    return rng.choice(len(X), size=n_labels, replace=False)


def select_kmeans(X: np.ndarray, n_labels: int, random_state: int) -> np.ndarray:
    if n_labels > len(X):
        raise ValueError("n_labels cannot exceed number of samples")

    kmeans = KMeans(n_clusters=n_labels, random_state=random_state, n_init=10)
    kmeans.fit(X)
    nearest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, X)

    selected: list[int] = []
    seen = set()
    for idx in nearest.tolist():
        if idx not in seen:
            selected.append(idx)
            seen.add(idx)

    # Rare duplicate-centroid-nearest cases are filled with unselected points
    # nearest to any centroid. This keeps the budget exact and deterministic.
    if len(selected) < n_labels:
        min_dist = kmeans.transform(X).min(axis=1)
        for idx in np.argsort(min_dist):
            idx = int(idx)
            if idx not in seen:
                selected.append(idx)
                seen.add(idx)
            if len(selected) == n_labels:
                break

    # Defensive fallback, in case the distance fill is ever insufficient.
    if len(selected) < n_labels:
        rng = np.random.default_rng(random_state)
        remaining = np.array([idx for idx in range(len(X)) if idx not in seen], dtype=int)
        fill = rng.choice(remaining, size=n_labels - len(selected), replace=False)
        selected.extend(fill.tolist())

    return np.array(selected[:n_labels], dtype=int)


def bootstrap_uncertainty(
    X: np.ndarray,
    y: np.ndarray,
    labeled: np.ndarray,
    unlabeled: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Fit bootstrap LinearRegression models on labeled points only."""
    _, sigma = bootstrap_predictions(X, y, labeled, unlabeled, n_bootstrap, rng)
    return sigma


def bootstrap_predictions(
    X: np.ndarray,
    y: np.ndarray,
    labeled: np.ndarray,
    unlabeled: np.ndarray,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Return bootstrap mean and standard deviation for unlabeled points.

    All bootstrap models are fit only on resampled currently labeled points.
    The true y values for unlabeled points are never used in acquisition.
    """
    n_features = X.shape[1]
    coefs = np.empty((n_bootstrap, n_features), dtype=float)
    intercepts = np.empty(n_bootstrap, dtype=float)

    for b in range(n_bootstrap):
        sample = rng.choice(labeled, size=len(labeled), replace=True)
        model = LinearRegression()
        model.fit(X[sample], y[sample])
        coefs[b] = model.coef_
        intercepts[b] = model.intercept_

    # Batch prediction is equivalent to calling each model, but much faster.
    predictions = X[unlabeled] @ coefs.T + intercepts
    return predictions.mean(axis=1), predictions.std(axis=1)


def select_active_bootstrap(
    X: np.ndarray,
    y: np.ndarray,
    initial_n: int = 50,
    batch_size: int = 20,
    final_n: int = 100,
    n_bootstrap: int = 30,
    random_state: int = 0,
) -> np.ndarray:
    if final_n > len(X):
        raise ValueError("final_n cannot exceed number of samples")
    if initial_n <= 0 or batch_size <= 0:
        raise ValueError("initial_n and batch_size must be positive")

    rng = np.random.default_rng(random_state)
    start_n = min(initial_n, final_n)
    labeled = rng.choice(len(X), size=start_n, replace=False)
    labeled_set = set(labeled.tolist())

    while len(labeled) < final_n:
        unlabeled = np.array([idx for idx in range(len(X)) if idx not in labeled_set], dtype=int)
        uncertainty = bootstrap_uncertainty(X, y, labeled, unlabeled, n_bootstrap, rng)
        acquire_n = min(batch_size, final_n - len(labeled))
        new_indices = unlabeled[np.argsort(uncertainty)[-acquire_n:]]
        labeled = np.concatenate([labeled, new_indices])
        labeled_set.update(new_indices.tolist())

    return labeled.astype(int)


def select_tail_aware_active(
    X: np.ndarray,
    y: np.ndarray,
    initial_n: int = 50,
    batch_size: int = 20,
    final_n: int = 100,
    n_bootstrap: int = 30,
    lambda_uncertainty: float = 1.0,
    random_state: int = 0,
) -> np.ndarray:
    """Tail-aware bootstrap active learning.

    Acquisition score for each unlabeled site:

        score = -mu + lambda_uncertainty * sigma

    where mu is the bootstrap mean predicted segregation energy and sigma is
    the bootstrap prediction standard deviation. The -mu term favors sites
    predicted to have more negative, segregation-relevant energies; the sigma
    term preserves uncertainty seeking. Only currently labeled y values are
    used to train the bootstrap models.
    """
    if final_n > len(X):
        raise ValueError("final_n cannot exceed number of samples")
    if initial_n <= 0 or batch_size <= 0:
        raise ValueError("initial_n and batch_size must be positive")

    rng = np.random.default_rng(random_state)
    if final_n <= initial_n:
        return rng.choice(len(X), size=final_n, replace=False).astype(int)

    labeled = rng.choice(len(X), size=initial_n, replace=False)
    labeled_set = set(labeled.tolist())

    while len(labeled) < final_n:
        unlabeled = np.array([idx for idx in range(len(X)) if idx not in labeled_set], dtype=int)
        mu, sigma = bootstrap_predictions(X, y, labeled, unlabeled, n_bootstrap, rng)
        # More negative predicted energy is favored by -mu; uncertainty is
        # favored by the lambda-weighted sigma term.
        score = -mu + lambda_uncertainty * sigma
        acquire_n = min(batch_size, final_n - len(labeled))
        new_indices = unlabeled[np.argsort(score)[-acquire_n:]]
        labeled = np.concatenate([labeled, new_indices])
        labeled_set.update(new_indices.tolist())

    return labeled.astype(int)


def select_xgb_aware_active(
    X: np.ndarray,
    y: np.ndarray,
    initial_n: int = 50,
    batch_size: int = 20,
    final_n: int = 100,
    n_bootstrap: int = 30,
    X_tot: float = DEFAULT_X_TOT,
    f_GB: float = DEFAULT_F_GB,
    T: float = DEFAULT_T,
    random_state: int = 0,
) -> np.ndarray:
    """Occupation-sensitivity-aware bootstrap active learning.

    The acquisition targets sites where energy uncertainty is expected to have
    the largest effect on predicted GB occupation. Bootstrap models produce a
    mean spectrum mu and uncertainty sigma. The mean spectrum is used to solve
    X_c from mass balance, and each unlabeled point receives:

        score = sigma_i * P_i * (1 - P_i)

    P_i * (1 - P_i) is largest near the occupation transition. True y values
    of unlabeled sites are not used during acquisition.
    """
    if final_n > len(X):
        raise ValueError("final_n cannot exceed number of samples")
    if initial_n <= 0 or batch_size <= 0:
        raise ValueError("initial_n and batch_size must be positive")

    rng = np.random.default_rng(random_state)
    if final_n <= initial_n:
        return rng.choice(len(X), size=final_n, replace=False).astype(int)

    labeled = rng.choice(len(X), size=initial_n, replace=False)
    labeled_set = set(labeled.tolist())
    all_indices = np.arange(len(X), dtype=int)

    while len(labeled) < final_n:
        unlabeled = np.array([idx for idx in range(len(X)) if idx not in labeled_set], dtype=int)
        mu_all, sigma_all = bootstrap_predictions(X, y, labeled, all_indices, n_bootstrap, rng)
        X_c = solve_Xc(mu_all, X_tot=X_tot, f_GB=f_GB, T=T)
        occupation = compute_occupation(mu_all[unlabeled], X_c, T=T)
        sensitivity = occupation * (1.0 - occupation)
        score = sigma_all[unlabeled] * sensitivity
        acquire_n = min(batch_size, final_n - len(labeled))
        new_indices = unlabeled[np.argsort(score)[-acquire_n:]]
        labeled = np.concatenate([labeled, new_indices])
        labeled_set.update(new_indices.tolist())

    return labeled.astype(int)


def evaluate_selection(
    X: np.ndarray,
    y: np.ndarray,
    selected: np.ndarray,
    XGB_true: float | None = None,
    X_tot: float = DEFAULT_X_TOT,
    f_GB: float = DEFAULT_F_GB,
    T: float = DEFAULT_T,
) -> tuple[Evaluation, np.ndarray]:
    model = LinearRegression()
    model.fit(X[selected], y[selected])
    y_pred = model.predict(X)
    residual = y_pred - y
    tail_cutoff = np.percentile(y, 10)
    tail_mask = y <= tail_cutoff
    if XGB_true is None:
        _, XGB_true = compute_xgb(y, X_tot=X_tot, f_GB=f_GB, T=T)
    _, XGB_pred = compute_xgb(y_pred, X_tot=X_tot, f_GB=f_GB, T=T)
    evaluation = Evaluation(
        mae=float(np.mean(np.abs(residual))),
        rmse=float(np.sqrt(np.mean(residual**2))),
        r2=float(r2_score(y, y_pred)),
        tail_mae=float(np.mean(np.abs(residual[tail_mask]))),
        xgb_true=float(XGB_true),
        xgb_pred=float(XGB_pred),
        xgb_abs_error=float(abs(XGB_pred - XGB_true)),
    )
    return evaluation, y_pred


def selected_indices_frame(df: pd.DataFrame, selected: np.ndarray) -> pd.DataFrame:
    out = pd.DataFrame({"site_index": selected})
    if "row_id" in df.columns:
        out["row_id"] = df.iloc[selected]["row_id"].to_numpy()
    if "site_id_if_verified" in df.columns:
        out["site_id_if_verified"] = df.iloc[selected]["site_id_if_verified"].to_numpy()
    return out


def prediction_frame(df: pd.DataFrame, y: np.ndarray, y_pred: np.ndarray) -> pd.DataFrame:
    out = pd.DataFrame(
        {
            "site_index": np.arange(len(y)),
            "true_deltaE": y,
            "predicted_deltaE": y_pred,
            "residual": y_pred - y,
        }
    )
    if "row_id" in df.columns:
        out.insert(1, "row_id", df["row_id"].to_numpy())
    if "site_id_if_verified" in df.columns:
        out.insert(2, "site_id_if_verified", df["site_id_if_verified"].to_numpy())
    return out


def add_result(
    rows: list[dict],
    budget: int,
    method: str,
    repeat: int,
    seed: int,
    selected: np.ndarray,
    evaluation: Evaluation,
) -> None:
    rows.append(
        {
            "budget": budget,
            "method": method,
            "repeat": repeat,
            "seed": seed,
            "n_selected": len(selected),
            "mae": evaluation.mae,
            "rmse": evaluation.rmse,
            "r2": evaluation.r2,
            "tail_mae": evaluation.tail_mae,
            "xgb_true": evaluation.xgb_true,
            "xgb_pred": evaluation.xgb_pred,
            "xgb_abs_error": evaluation.xgb_abs_error,
        }
    )


def run_learning_curves(
    X: np.ndarray,
    y: np.ndarray,
    budgets: list[int],
    random_repeats: int,
    active_repeats: int,
    n_bootstrap: int,
    lambda_uncertainty: float,
    X_tot: float,
    f_GB: float,
    T: float,
) -> pd.DataFrame:
    rows: list[dict] = []
    _, XGB_true = compute_xgb(y, X_tot=X_tot, f_GB=f_GB, T=T)
    print(f"True X_GB at T={T:g} K, X_tot={X_tot:g}, f_GB={f_GB:g}: {XGB_true:.8f}")

    for budget in budgets:
        print(f"\nBudget {budget}")

        for repeat in range(random_repeats):
            seed = 1000 + budget * 100 + repeat
            selected = select_random(X, budget, seed)
            evaluation, _ = evaluate_selection(X, y, selected, XGB_true, X_tot, f_GB, T)
            add_result(rows, budget, "random", repeat, seed, selected, evaluation)
        random_mae = [row["mae"] for row in rows if row["budget"] == budget and row["method"] == "random"]
        print(f"  random MAE: {np.mean(random_mae):.3f} +/- {np.std(random_mae, ddof=1):.3f}")

        selected = select_kmeans(X, budget, random_state=42)
        evaluation, _ = evaluate_selection(X, y, selected, XGB_true, X_tot, f_GB, T)
        add_result(rows, budget, "kmeans_raw_pca", 0, 42, selected, evaluation)
        print(f"  k-means MAE: {evaluation.mae:.3f}")

        for repeat in range(active_repeats):
            seed = 2000 + budget * 100 + repeat
            if budget < 50:
                selected = select_random(X, budget, seed)
            else:
                selected = select_active_bootstrap(
                    X,
                    y,
                    initial_n=50,
                    batch_size=20,
                    final_n=budget,
                    n_bootstrap=n_bootstrap,
                    random_state=seed,
                )
            evaluation, _ = evaluate_selection(X, y, selected, XGB_true, X_tot, f_GB, T)
            add_result(rows, budget, "active_bootstrap_raw_pca", repeat, seed, selected, evaluation)
        active_mae = [
            row["mae"] for row in rows if row["budget"] == budget and row["method"] == "active_bootstrap_raw_pca"
        ]
        print(f"  active bootstrap MAE: {np.mean(active_mae):.3f} +/- {np.std(active_mae, ddof=1):.3f}")

        for repeat in range(active_repeats):
            seed = 3000 + budget * 100 + repeat
            if budget < 50:
                selected = select_random(X, budget, seed)
            else:
                selected = select_tail_aware_active(
                    X,
                    y,
                    initial_n=50,
                    batch_size=20,
                    final_n=budget,
                    n_bootstrap=n_bootstrap,
                    lambda_uncertainty=lambda_uncertainty,
                    random_state=seed,
                )
            evaluation, _ = evaluate_selection(X, y, selected, XGB_true, X_tot, f_GB, T)
            add_result(rows, budget, "tail_aware_active_raw_pca", repeat, seed, selected, evaluation)
        tail_active_mae = [
            row["mae"] for row in rows if row["budget"] == budget and row["method"] == "tail_aware_active_raw_pca"
        ]
        print(
            "  tail-aware active MAE: "
            f"{np.mean(tail_active_mae):.3f} +/- {np.std(tail_active_mae, ddof=1):.3f}"
        )

        for repeat in range(active_repeats):
            seed = 4000 + budget * 100 + repeat
            if budget < 50:
                selected = select_random(X, budget, seed)
            else:
                selected = select_xgb_aware_active(
                    X,
                    y,
                    initial_n=50,
                    batch_size=20,
                    final_n=budget,
                    n_bootstrap=n_bootstrap,
                    X_tot=X_tot,
                    f_GB=f_GB,
                    T=T,
                    random_state=seed,
                )
            evaluation, _ = evaluate_selection(X, y, selected, XGB_true, X_tot, f_GB, T)
            add_result(rows, budget, "xgb_aware_active_raw_pca", repeat, seed, selected, evaluation)
        xgb_active_mae = [
            row["mae"] for row in rows if row["budget"] == budget and row["method"] == "xgb_aware_active_raw_pca"
        ]
        xgb_active_err = [
            row["xgb_abs_error"] for row in rows if row["budget"] == budget and row["method"] == "xgb_aware_active_raw_pca"
        ]
        print(
            "  XGB-aware active MAE: "
            f"{np.mean(xgb_active_mae):.3f} +/- {np.std(xgb_active_mae, ddof=1):.3f}; "
            f"XGB abs error: {np.mean(xgb_active_err):.5f} +/- {np.std(xgb_active_err, ddof=1):.5f}"
        )

    return pd.DataFrame(rows)


def summarize_metric(results: pd.DataFrame, metric: str) -> pd.DataFrame:
    summary = (
        results.groupby(["budget", "method"], as_index=False)
        .agg(
            mean=(metric, "mean"),
            std=(metric, lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else 0.0),
            n_runs=(metric, "size"),
        )
        .sort_values(["budget", "method"])
    )
    return summary


def plot_learning_curve(
    results: pd.DataFrame,
    metric: str,
    ylabel: str,
    title: str,
    output: Path,
) -> None:
    summary = summarize_metric(results, metric)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    colors = {
        "random": "tab:blue",
        "kmeans_raw_pca": "tab:orange",
        "active_bootstrap_raw_pca": "tab:green",
        "tail_aware_active_raw_pca": "tab:red",
        "xgb_aware_active_raw_pca": "tab:purple",
    }
    labels = {
        "random": "random",
        "kmeans_raw_pca": "kmeans_raw_pca",
        "active_bootstrap_raw_pca": "active_bootstrap_raw_pca",
        "tail_aware_active_raw_pca": "tail_aware_active_raw_pca",
        "xgb_aware_active_raw_pca": "xgb_aware_active_raw_pca",
    }

    for method in ["random", "kmeans_raw_pca", "active_bootstrap_raw_pca", "tail_aware_active_raw_pca", "xgb_aware_active_raw_pca"]:
        data = summary[summary["method"] == method]
        x = data["budget"].to_numpy()
        y = data["mean"].to_numpy()
        yerr = data["std"].to_numpy()
        ax.plot(x, y, marker="o", color=colors[method], label=labels[method])
        if method != "kmeans_raw_pca":
            ax.fill_between(x, y - yerr, y + yerr, color=colors[method], alpha=0.16, linewidth=0)

    ax.set_xlabel("Number of labeled sites")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_selected_points(
    X_raw: np.ndarray,
    feature_cols: list[str],
    selected_kmeans: np.ndarray,
    selected_active: np.ndarray,
    selected_tail_aware: np.ndarray | None,
    selected_xgb_aware: np.ndarray | None,
    output: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.8))
    ax.scatter(X_raw[:, 0], X_raw[:, 1], s=3, alpha=0.14, color="0.45", rasterized=True, label="All sites")
    ax.scatter(
        X_raw[selected_kmeans, 0],
        X_raw[selected_kmeans, 1],
        s=34,
        marker="x",
        linewidths=1.2,
        color="tab:orange",
        label="kmeans_raw_pca, n=100",
    )
    ax.scatter(
        X_raw[selected_active, 0],
        X_raw[selected_active, 1],
        s=22,
        marker="o",
        facecolors="none",
        edgecolors="tab:green",
        linewidths=1.0,
        label="active_bootstrap_raw_pca, n=100",
    )
    if selected_tail_aware is not None:
        ax.scatter(
            X_raw[selected_tail_aware, 0],
            X_raw[selected_tail_aware, 1],
            s=26,
            marker="^",
            facecolors="none",
            edgecolors="tab:red",
            linewidths=1.0,
            label="tail_aware_active_raw_pca, n=100",
        )
    if selected_xgb_aware is not None:
        ax.scatter(
            X_raw[selected_xgb_aware, 0],
            X_raw[selected_xgb_aware, 1],
            s=28,
            marker="s",
            facecolors="none",
            edgecolors="tab:purple",
            linewidths=1.0,
            label="xgb_aware_active_raw_pca, n=100",
        )
    ax.set_xlabel(feature_cols[0])
    ax.set_ylabel(feature_cols[1])
    ax.set_title("Selected Al-Mg GB training sites in raw PCA space")
    ax.grid(True, alpha=0.2)
    ax.legend(frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_predicted_vs_true(
    y: np.ndarray,
    pred_kmeans: np.ndarray,
    pred_active: np.ndarray,
    pred_tail_aware: np.ndarray | None,
    pred_xgb_aware: np.ndarray | None,
    mae_kmeans: float,
    mae_active: float,
    mae_tail_aware: float | None,
    mae_xgb_aware: float | None,
    output: Path,
) -> None:
    n_panels = 2 + int(pred_tail_aware is not None) + int(pred_xgb_aware is not None)
    fig, axes = plt.subplots(1, n_panels, figsize=(5.2 * n_panels, 4.7), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    panels = [
        (
            axes[0],
            pred_kmeans,
            f"kmeans_raw_pca, 100 labels\nRepresentative run; MAE = {mae_kmeans:.2f} kJ/mol",
            "tab:orange",
        ),
        (
            axes[1],
            pred_active,
            f"active_bootstrap_raw_pca, 100 labels\nRepresentative run; MAE = {mae_active:.2f} kJ/mol",
            "tab:green",
        ),
    ]
    pred_arrays = [pred_kmeans, pred_active]
    if pred_tail_aware is not None:
        pred_arrays.append(pred_tail_aware)
        panels.append(
            (
                axes[2],
                pred_tail_aware,
                f"tail_aware_active_raw_pca, 100 labels\nRepresentative run; MAE = {mae_tail_aware:.2f} kJ/mol",
                "tab:red",
            )
        )
    if pred_xgb_aware is not None:
        pred_arrays.append(pred_xgb_aware)
        panels.append(
            (
                axes[len(panels)],
                pred_xgb_aware,
                f"xgb_aware_active_raw_pca, 100 labels\nRepresentative run; MAE = {mae_xgb_aware:.2f} kJ/mol",
                "tab:purple",
            )
        )
    lo = float(min([y.min(), *[pred.min() for pred in pred_arrays]]))
    hi = float(max([y.max(), *[pred.max() for pred in pred_arrays]]))
    pad = 0.04 * (hi - lo)
    lo -= pad
    hi += pad
    for ax, pred, title, color in panels:
        ax.scatter(y, pred, s=4, alpha=0.18, color=color, rasterized=True)
        ax.plot([lo, hi], [lo, hi], color="black", linestyle="--", linewidth=1.0)
        ax.set_title(title)
        ax.set_xlabel("True DeltaE (kJ/mol)")
        ax.grid(True, alpha=0.22)
    axes[0].set_ylabel("Predicted DeltaE (kJ/mol)")
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)


def print_final_summary(results: pd.DataFrame) -> None:
    at100 = results[results["budget"] == 100]
    means = at100.groupby("method", as_index=False).agg(
        mae=("mae", "mean"),
        tail_mae=("tail_mae", "mean"),
        xgb_abs_error=("xgb_abs_error", "mean"),
    )
    best_mae = means.loc[means["mae"].idxmin()]
    best_tail = means.loc[means["tail_mae"].idxmin()]
    best_xgb = means.loc[means["xgb_abs_error"].idxmin()]

    kmeans_mae = float(means.loc[means["method"] == "kmeans_raw_pca", "mae"].iloc[0])
    active_mae = float(means.loc[means["method"] == "active_bootstrap_raw_pca", "mae"].iloc[0])
    tail_active_mae = float(means.loc[means["method"] == "tail_aware_active_raw_pca", "mae"].iloc[0])
    xgb_active_mae = float(means.loc[means["method"] == "xgb_aware_active_raw_pca", "mae"].iloc[0])
    kmeans_tail = float(means.loc[means["method"] == "kmeans_raw_pca", "tail_mae"].iloc[0])
    active_tail = float(means.loc[means["method"] == "active_bootstrap_raw_pca", "tail_mae"].iloc[0])
    tail_active_tail = float(means.loc[means["method"] == "tail_aware_active_raw_pca", "tail_mae"].iloc[0])
    xgb_active_tail = float(means.loc[means["method"] == "xgb_aware_active_raw_pca", "tail_mae"].iloc[0])
    kmeans_xgb = float(means.loc[means["method"] == "kmeans_raw_pca", "xgb_abs_error"].iloc[0])
    active_xgb = float(means.loc[means["method"] == "active_bootstrap_raw_pca", "xgb_abs_error"].iloc[0])
    tail_active_xgb = float(means.loc[means["method"] == "tail_aware_active_raw_pca", "xgb_abs_error"].iloc[0])
    xgb_active_xgb = float(means.loc[means["method"] == "xgb_aware_active_raw_pca", "xgb_abs_error"].iloc[0])

    print("\nConcise summary")
    print("----------------")
    print(f"Best method by MAE at 100 labels: {best_mae['method']} ({best_mae['mae']:.3f} kJ/mol)")
    print(
        f"Best method by tail MAE at 100 labels: {best_tail['method']} "
        f"({best_tail['tail_mae']:.3f} kJ/mol)"
    )
    print(
        f"Best method by XGB_abs_error at 100 labels: {best_xgb['method']} "
        f"({best_xgb['xgb_abs_error']:.6f})"
    )
    print(
        "kmeans_raw_pca 100-label results: "
        f"MAE={kmeans_mae:.3f} kJ/mol, tail_MAE={kmeans_tail:.3f} kJ/mol, "
        f"XGB_abs_error={kmeans_xgb:.6f}"
    )
    print(f"active_bootstrap_raw_pca beats kmeans_raw_pca by MAE: {active_mae < kmeans_mae}")
    print(f"active_bootstrap_raw_pca beats kmeans_raw_pca by tail MAE: {active_tail < kmeans_tail}")
    print(
        "active_bootstrap_raw_pca beats kmeans_raw_pca by XGB_abs_error: "
        f"{active_xgb < kmeans_xgb}"
    )
    print(f"tail_aware_active_raw_pca beats kmeans_raw_pca by MAE: {tail_active_mae < kmeans_mae}")
    print(
        "tail_aware_active_raw_pca beats kmeans_raw_pca by tail MAE: "
        f"{tail_active_tail < kmeans_tail}"
    )
    print(
        "tail_aware_active_raw_pca beats kmeans_raw_pca by XGB_abs_error: "
        f"{tail_active_xgb < kmeans_xgb}"
    )
    print(f"xgb_aware_active_raw_pca beats kmeans_raw_pca by MAE: {xgb_active_mae < kmeans_mae}")
    print(
        "xgb_aware_active_raw_pca beats kmeans_raw_pca by tail MAE: "
        f"{xgb_active_tail < kmeans_tail}"
    )
    print(f"xgb_aware_active_raw_pca beats kmeans_raw_pca by XGB_abs_error: {xgb_active_xgb < kmeans_xgb}")
    print(
        "xgb_aware_active_raw_pca beats active_bootstrap_raw_pca by XGB_abs_error: "
        f"{xgb_active_xgb < active_xgb}"
    )

    print("\nScientific interpretation")
    print("-------------------------")
    print("kmeans_raw_pca is an unsupervised raw-PCA feature-space coverage baseline.")
    print("active_bootstrap_raw_pca selects sites with high model uncertainty in raw PCA space.")
    print("tail_aware_active_raw_pca selects sites that are predicted to be low-energy and uncertain.")
    print("xgb_aware_active_raw_pca selects sites where energy uncertainty most affects GB occupation.")
    print("Global MAE measures average energy prediction quality over all grain-boundary sites.")
    print("Tail MAE matters because the lowest-energy sites dominate segregation tendency.")
    print("XGB_abs_error measures the error in predicted equilibrium GB solute concentration.")
    if xgb_active_xgb < kmeans_xgb or xgb_active_xgb < active_xgb:
        print("The new method only counts as better where the numerical results above show an improvement.")
    else:
        print("The new method does not beat the stated baselines at 100 labels on XGB_abs_error.")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        type=Path,
        default=PACKAGE_ROOT / "data" / "al_mg_pca_deltaE_verified.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PACKAGE_ROOT / "outputs" / "main_comparison",
    )
    parser.add_argument("--figures-dir", type=Path, default=None)
    parser.add_argument("--random-repeats", type=int, default=20)
    parser.add_argument("--active-repeats", type=int, default=10)
    parser.add_argument("--n-bootstrap", type=int, default=30)
    parser.add_argument("--lambda-uncertainty", type=float, default=1.0)
    parser.add_argument("--x-tot", type=float, default=DEFAULT_X_TOT)
    parser.add_argument("--f-gb", type=float, default=DEFAULT_F_GB)
    parser.add_argument("--temperature", type=float, default=DEFAULT_T)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.figures_dir is None:
        args.figures_dir = args.output_dir / "figures"
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    df, feature_cols, target_col, X_raw, X, y = load_dataset(args.csv)

    results = run_learning_curves(
        X,
        y,
        DEFAULT_BUDGETS,
        random_repeats=args.random_repeats,
        active_repeats=args.active_repeats,
        n_bootstrap=args.n_bootstrap,
        lambda_uncertainty=args.lambda_uncertainty,
        X_tot=args.x_tot,
        f_GB=args.f_gb,
        T=args.temperature,
    )
    results.to_csv(args.output_dir / "selection_results.csv", index=False)
    summary = (
        results.groupby(["budget", "method"], as_index=False)
        .agg(
            mae_mean=("mae", "mean"),
            mae_std=("mae", lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else 0.0),
            tail_mae_mean=("tail_mae", "mean"),
            tail_mae_std=("tail_mae", lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else 0.0),
            xgb_abs_error_mean=("xgb_abs_error", "mean"),
            xgb_abs_error_std=(
                "xgb_abs_error",
                lambda x: float(np.std(x, ddof=1)) if len(x) > 1 else 0.0,
            ),
            rmse_mean=("rmse", "mean"),
            r2_mean=("r2", "mean"),
            n_runs=("mae", "size"),
        )
    )
    summary.to_csv(args.output_dir / "selection_results_summary.csv", index=False)
    _, XGB_true = compute_xgb(y, X_tot=args.x_tot, f_GB=args.f_gb, T=args.temperature)

    selected_random_100 = select_random(X, 100, random_state=100100)
    selected_kmeans_100 = select_kmeans(X, 100, random_state=42)
    selected_active_100 = select_active_bootstrap(
        X,
        y,
        initial_n=50,
        batch_size=20,
        final_n=100,
        n_bootstrap=args.n_bootstrap,
        random_state=100,
    )
    selected_tail_aware_100 = select_tail_aware_active(
        X,
        y,
        initial_n=50,
        batch_size=20,
        final_n=100,
        n_bootstrap=args.n_bootstrap,
        lambda_uncertainty=args.lambda_uncertainty,
        random_state=100,
    )
    selected_xgb_aware_100 = select_xgb_aware_active(
        X,
        y,
        initial_n=50,
        batch_size=20,
        final_n=100,
        n_bootstrap=args.n_bootstrap,
        X_tot=args.x_tot,
        f_GB=args.f_gb,
        T=args.temperature,
        random_state=100,
    )

    selected_indices_frame(df, selected_random_100).to_csv(
        args.output_dir / "selected_indices_random_100.csv", index=False
    )
    selected_indices_frame(df, selected_kmeans_100).to_csv(
        args.output_dir / "selected_indices_kmeans_raw_pca_100.csv", index=False
    )
    selected_indices_frame(df, selected_active_100).to_csv(
        args.output_dir / "selected_indices_active_bootstrap_raw_pca_100.csv", index=False
    )
    selected_indices_frame(df, selected_tail_aware_100).to_csv(
        args.output_dir / "selected_indices_tail_aware_active_raw_pca_100.csv", index=False
    )
    selected_indices_frame(df, selected_xgb_aware_100).to_csv(
        args.output_dir / "selected_indices_xgb_aware_active_raw_pca_100.csv", index=False
    )

    eval_kmeans_100, pred_kmeans_100 = evaluate_selection(
        X, y, selected_kmeans_100, XGB_true, args.x_tot, args.f_gb, args.temperature
    )
    eval_active_100, pred_active_100 = evaluate_selection(
        X, y, selected_active_100, XGB_true, args.x_tot, args.f_gb, args.temperature
    )
    eval_tail_aware_100, pred_tail_aware_100 = evaluate_selection(
        X, y, selected_tail_aware_100, XGB_true, args.x_tot, args.f_gb, args.temperature
    )
    eval_xgb_aware_100, pred_xgb_aware_100 = evaluate_selection(
        X, y, selected_xgb_aware_100, XGB_true, args.x_tot, args.f_gb, args.temperature
    )
    prediction_frame(df, y, pred_kmeans_100).to_csv(
        args.output_dir / "predictions_kmeans_raw_pca_100.csv", index=False
    )
    prediction_frame(df, y, pred_active_100).to_csv(
        args.output_dir / "predictions_active_bootstrap_raw_pca_100.csv", index=False
    )
    prediction_frame(df, y, pred_tail_aware_100).to_csv(
        args.output_dir / "predictions_tail_aware_active_raw_pca_100.csv", index=False
    )
    prediction_frame(df, y, pred_xgb_aware_100).to_csv(
        args.output_dir / "predictions_xgb_aware_active_raw_pca_100.csv", index=False
    )

    title = "Training-site selection in raw PCA space for Al-Mg GB segregation energy prediction"
    plot_learning_curve(
        results,
        "mae",
        "MAE in DeltaE prediction (kJ/mol)",
        title,
        args.figures_dir / "mae_vs_labels.png",
    )
    plot_learning_curve(
        results,
        "tail_mae",
        "Tail MAE, lowest-energy 10% (kJ/mol)",
        title,
        args.figures_dir / "tail_mae_vs_labels.png",
    )
    plot_learning_curve(
        results,
        "xgb_abs_error",
        "Absolute error in X_GB",
        title,
        args.figures_dir / "xgb_error_vs_labels.png",
    )
    plot_selected_points(
        X_raw,
        feature_cols,
        selected_kmeans_100,
        selected_active_100,
        selected_tail_aware_100,
        selected_xgb_aware_100,
        args.figures_dir / "pca_selected_points.png",
    )
    plot_selected_points(
        X_raw,
        feature_cols,
        selected_kmeans_100,
        selected_active_100,
        selected_tail_aware_100,
        selected_xgb_aware_100,
        args.figures_dir / "pca_selected_points_all_methods.png",
    )
    plot_predicted_vs_true(
        y,
        pred_kmeans_100,
        pred_active_100,
        None,
        None,
        eval_kmeans_100.mae,
        eval_active_100.mae,
        None,
        None,
        args.figures_dir / "predicted_vs_true_100labels.png",
    )
    plot_predicted_vs_true(
        y,
        pred_kmeans_100,
        pred_active_100,
        pred_tail_aware_100,
        pred_xgb_aware_100,
        eval_kmeans_100.mae,
        eval_active_100.mae,
        eval_tail_aware_100.mae,
        eval_xgb_aware_100.mae,
        args.figures_dir / "predicted_vs_true_100labels_all_methods.png",
    )

    print("\nSaved outputs:")
    print(f"  tables and per-site CSVs: {args.output_dir}")
    print(f"  {args.figures_dir}/mae_vs_labels.png")
    print(f"  {args.figures_dir}/tail_mae_vs_labels.png")
    print(f"  {args.figures_dir}/xgb_error_vs_labels.png")
    print(f"  {args.figures_dir}/pca_selected_points.png")
    print(f"  {args.figures_dir}/pca_selected_points_all_methods.png")
    print(f"  {args.figures_dir}/predicted_vs_true_100labels.png")
    print(f"  {args.figures_dir}/predicted_vs_true_100labels_all_methods.png")

    print_final_summary(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
