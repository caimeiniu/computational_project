#!/usr/bin/env python3
"""
All-data 10D raw PCA + LinearRegression baseline for Al-Mg GB segregation.

This is an isolated baseline script. It writes only inside:
    project/point_selection/outputs/all_data_baseline/

It does not rerun any selection strategies and does not modify the main results.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PACKAGE_ROOT / "outputs" / "all_data_baseline"
DATA_PATH = PACKAGE_ROOT / "data" / "al_mg_pca_deltaE_verified.csv"
SUMMARY_PATH = PACKAGE_ROOT / "results" / "reference" / "selection_results_summary.csv"

FEATURE_COLS = [f"pc{i}" for i in range(1, 11)]
TARGET_COL = "deltaE_kJmol"
KB_KJMOL_K = 0.008314462618
DEFAULT_T = 600.0
DEFAULT_X_TOT = 0.05
DEFAULT_F_GB = 0.10


def stable_sigmoid(z: np.ndarray) -> np.ndarray:
    """Numerically stable sigmoid for occupation probabilities."""
    out = np.empty_like(z, dtype=float)
    positive = z >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-z[positive]))
    exp_z = np.exp(z[~positive])
    out[~positive] = exp_z / (1.0 + exp_z)
    return out


def compute_occupation(deltaE: np.ndarray, X_c: float, T: float = DEFAULT_T) -> np.ndarray:
    """Compute Fermi-Dirac-like GB occupation probabilities for each site."""
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
    """Solve mass balance for matrix concentration X_c."""
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
) -> tuple[float, float, np.ndarray]:
    """Return X_GB, X_c, and occupation probabilities for a spectrum."""
    X_c = solve_Xc(deltaE, X_tot=X_tot, f_GB=f_GB, T=T)
    p = compute_occupation(deltaE, X_c, T=T)
    return float(np.mean(p)), float(X_c), p


def load_data() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    df = pd.read_csv(DATA_PATH)
    missing = [col for col in FEATURE_COLS + [TARGET_COL] if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    if df[FEATURE_COLS + [TARGET_COL]].isna().any().any():
        raise ValueError("Input contains NaN values in feature or target columns.")
    X = df[FEATURE_COLS].to_numpy(dtype=float)
    y = df[TARGET_COL].to_numpy(dtype=float)
    return df, X, y


def evaluate_all_data_model(X: np.ndarray, y: np.ndarray) -> tuple[dict[str, float], np.ndarray]:
    """Train LinearRegression on all rows and evaluate on all rows."""
    model = LinearRegression()
    model.fit(X, y)
    y_pred = model.predict(X)
    residual = y_pred - y
    tail_cutoff = np.percentile(y, 10)
    tail_mask = y <= tail_cutoff
    xgb_true, x_c_true, _ = compute_xgb(y)
    xgb_pred, x_c_pred, _ = compute_xgb(y_pred)
    metrics = {
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "r2": float(r2_score(y, y_pred)),
        "tail_mae": float(np.mean(np.abs(residual[tail_mask]))),
        "xgb_true": xgb_true,
        "xgb_pred": xgb_pred,
        "xgb_abs_error": float(abs(xgb_pred - xgb_true)),
        "x_c_true": x_c_true,
        "x_c_pred": x_c_pred,
    }
    return metrics, y_pred


def load_kmeans_100_summary() -> dict[str, float]:
    summary = pd.read_csv(SUMMARY_PATH)
    mask = (summary["method"] == "kmeans_raw_pca") & (summary["budget"] == 100)
    if not mask.any():
        raise ValueError("Could not find method=kmeans_raw_pca, budget=100 in selection_results_summary.csv")
    row = summary.loc[mask].iloc[0]
    return {
        "mae": float(row["mae_mean"]),
        "tail_mae": float(row["tail_mae_mean"]),
        "xgb_abs_error": float(row["xgb_abs_error_mean"]),
        "rmse": float(row["rmse_mean"]),
        "r2": float(row["r2_mean"]),
    }


def save_predictions(df: pd.DataFrame, y: np.ndarray, y_pred: np.ndarray) -> None:
    out = pd.DataFrame(
        {
            "true_deltaE": y,
            "predicted_deltaE": y_pred,
            "residual": y_pred - y,
        }
    )
    if "row_id" in df.columns:
        out.insert(0, "row_id", df["row_id"].to_numpy())
    else:
        out.insert(0, "row_id", np.arange(len(df)))
    if "site_id_if_verified" in df.columns:
        out.insert(1, "site_id_if_verified", df["site_id_if_verified"].to_numpy())
    else:
        out.insert(1, "site_id_if_verified", np.nan)
    out.to_csv(OUTPUT_DIR / "predictions_all_data_10pc_linear_regression.csv", index=False)


def save_metric_csv(metrics: dict[str, float], n_rows: int) -> None:
    row = {
        "model_name": "all_data_10pc_linear_regression",
        "n_training_sites": n_rows,
        "n_evaluation_sites": n_rows,
        "feature_space": "raw pc1...pc10; no StandardScaler; no new PCA",
        "model": "sklearn.linear_model.LinearRegression",
        "mae": metrics["mae"],
        "rmse": metrics["rmse"],
        "r2": metrics["r2"],
        "tail_mae": metrics["tail_mae"],
        "xgb_true": metrics["xgb_true"],
        "xgb_pred": metrics["xgb_pred"],
        "xgb_abs_error": metrics["xgb_abs_error"],
    }
    pd.DataFrame([row]).to_csv(
        OUTPUT_DIR / "all_data_10pc_linear_regression_baseline.csv", index=False
    )


def save_plot(y: np.ndarray, y_pred: np.ndarray, metrics: dict[str, float]) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.4), dpi=180)
    ax.scatter(y, y_pred, s=3, alpha=0.22, color="#1f77b4", linewidths=0, rasterized=True)
    low = float(min(np.min(y), np.min(y_pred)))
    high = float(max(np.max(y), np.max(y_pred)))
    pad = 0.04 * (high - low)
    lims = [low - pad, high + pad]
    ax.plot(lims, lims, color="#111827", linestyle="--", linewidth=1.2, label="y = x")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("true DeltaE (kJ/mol)")
    ax.set_ylabel("predicted DeltaE (kJ/mol)")
    ax.set_title(
        "All-data 10D raw PCA LinearRegression baseline\n"
        f"MAE = {metrics['mae']:.3f} kJ/mol, R2 = {metrics['r2']:.3f}"
    )
    ax.legend(frameon=False)
    ax.grid(True, color="#d8dee9", linewidth=0.6, alpha=0.8)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "predicted_vs_true_all_data_10pc_lr.png")
    plt.close(fig)


def write_report(metrics: dict[str, float], kmeans: dict[str, float], n_rows: int) -> None:
    mae_gap = kmeans["mae"] - metrics["mae"]
    tail_gap = kmeans["tail_mae"] - metrics["tail_mae"]
    xgb_gap = kmeans["xgb_abs_error"] - metrics["xgb_abs_error"]

    close_text = (
        "The 100-label kmeans_raw_pca result is still noticeably above the all-data "
        "linear-regression baseline in global MAE, so label budget and training-site "
        "selection remain important sources of error."
        if mae_gap > 0.5
        else "The 100-label kmeans_raw_pca result is fairly close to the all-data "
        "linear-regression baseline in global MAE, suggesting the current 10D PCA + "
        "LinearRegression model ceiling may already be a major limitation."
    )
    tail_text = (
        "Tail MAE should be interpreted separately: the all-data fit is the best in-sample "
        "linear fit available to this feature/model setup, but low-energy tail errors can "
        "remain even when all labels are used."
    )
    xgb_text = (
        "The all-data fit does not necessarily minimize XGB_abs_error, because XGB is a nonlinear "
        "functional of the whole predicted energy spectrum. In this run, 100-label k-means has a "
        "slightly smaller XGB_abs_error than the all-data linear fit, even though the all-data fit "
        "has lower global and tail MAE."
        if xgb_gap < 0
        else "The all-data fit also improves XGB_abs_error relative to 100-label k-means in this run, "
        "but XGB remains a nonlinear spectrum-level metric and should be interpreted alongside MAE."
    )

    report = f"""# All-Data 10D Raw PCA LinearRegression Baseline

## A. What Was Done

- All `{n_rows:,}` GB sites were used for training.
- Raw `pc1` ... `pc10` were used as features.
- The model was `sklearn.linear_model.LinearRegression()`.
- No `StandardScaler` was used.
- No new PCA was fitted.
- The same `{n_rows:,}` sites were used for evaluation.

This is an in-sample all-data baseline. It is not an independent generalization test; it is a ceiling/lower-bound reference for the current dataset and feature/model setup.

## B. Why This Baseline Matters

This is the all-data baseline for the current 10D raw PCA + LinearRegression setup. It estimates the best achievable in-sample performance of this feature/model combination on the current dataset. It helps determine whether the 100-label `kmeans_raw_pca` result is mainly limited by label budget / selection, or whether the 10D PCA + LinearRegression model ceiling is already limiting performance.

## C. Results

| Metric | Value |
|---|---:|
| global MAE | {metrics['mae']:.6f} kJ/mol |
| RMSE | {metrics['rmse']:.6f} kJ/mol |
| R2 | {metrics['r2']:.6f} |
| tail MAE, lowest-energy 10% | {metrics['tail_mae']:.6f} kJ/mol |
| XGB_true | {metrics['xgb_true']:.8f} |
| XGB_pred | {metrics['xgb_pred']:.8f} |
| XGB_abs_error | {metrics['xgb_abs_error']:.8f} |
| Xc_true | {metrics['x_c_true']:.8f} |
| Xc_pred | {metrics['x_c_pred']:.8f} |

## D. Comparison With Current 100-Label `kmeans_raw_pca`

| Metric | 100-label kmeans_raw_pca | all-data 10D LR baseline | gap: kmeans100 - all-data |
|---|---:|---:|---:|
| global MAE | {kmeans['mae']:.6f} | {metrics['mae']:.6f} | {mae_gap:.6f} |
| tail MAE | {kmeans['tail_mae']:.6f} | {metrics['tail_mae']:.6f} | {tail_gap:.6f} |
| XGB_abs_error | {kmeans['xgb_abs_error']:.8f} | {metrics['xgb_abs_error']:.8f} | {xgb_gap:.8f} |

## E. Interpretation

{close_text}

{tail_text}

{xgb_text}

This baseline should be used cautiously: because it trains and evaluates on all sites, it does not measure out-of-sample prediction. Instead, it tells us what the current 10D raw PCA + LinearRegression model can achieve when label scarcity is removed.

## F. Recommended Presentation Wording

1. “As a ceiling reference for the current feature/model setup, we trained LinearRegression on all 82,645 raw-PCA sites and evaluated in sample.”
2. “The gap between 100-label k-means and the all-data fit separates label-budget effects from the intrinsic 10D PCA + linear-model ceiling.”
3. “This all-data baseline is not an independent test set result; it is a lower-bound reference for achievable reconstruction error on this dataset.”
"""
    (OUTPUT_DIR / "all_data_10pc_lr_baseline_report.md").write_text(report, encoding="utf-8")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df, X, y = load_data()
    metrics, y_pred = evaluate_all_data_model(X, y)
    kmeans = load_kmeans_100_summary()

    save_metric_csv(metrics, len(df))
    save_predictions(df, y, y_pred)
    save_plot(y, y_pred, metrics)
    write_report(metrics, kmeans, len(df))

    print("dataset shape:", df.shape)
    print("features used:", FEATURE_COLS)
    print("target used:", TARGET_COL)
    print("StandardScaler used: No")
    print(f"all-data MAE: {metrics['mae']:.6f} kJ/mol")
    print(f"all-data tail MAE: {metrics['tail_mae']:.6f} kJ/mol")
    print(f"all-data XGB_abs_error: {metrics['xgb_abs_error']:.8f}")
    print(f"100-label kmeans_raw_pca MAE: {kmeans['mae']:.6f} kJ/mol")
    print(f"mae_gap: {kmeans['mae'] - metrics['mae']:.6f}")
    print(f"tail_mae_gap: {kmeans['tail_mae'] - metrics['tail_mae']:.6f}")
    print(f"xgb_error_gap: {kmeans['xgb_abs_error'] - metrics['xgb_abs_error']:.8f}")
    print("report path:", OUTPUT_DIR / "all_data_10pc_lr_baseline_report.md")


if __name__ == "__main__":
    main()
