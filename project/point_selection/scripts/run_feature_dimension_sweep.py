#!/usr/bin/env python3
"""
Feature baseline sweep: all-data LinearRegression on 10/20/50/100 raw SOAP PCA features.

All outputs are written inside:
    project/point_selection/outputs/feature_dimension_sweep/

This script does not modify main selection results, existing figures, PPT files, or
the all_data_baseline folder.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PACKAGE_ROOT / "outputs" / "feature_dimension_sweep"
DATA_PATH = PACKAGE_ROOT / "data" / "al_mg_pca_deltaE_verified.csv"
RESULTS_PATH = OUTPUT_DIR / "pca_10_to_100_baseline_results.csv"
REPORT_PATH = OUTPUT_DIR / "pca_10_to_100_baseline_report.md"
PREDICTIONS_DIR = OUTPUT_DIR / "predictions"
FIGURES_DIR = OUTPUT_DIR / "figures"

PCA_DIMS = [10, 20, 50, 100]
FEATURE_COLS_10D = [f"pc{i}" for i in range(1, 11)]
TARGET_COL = "deltaE_kJmol"
KB_KJMOL_K = 0.008314462618
DEFAULT_T = 600.0
DEFAULT_X_TOT = 0.05
DEFAULT_F_GB = 0.10


def find_soap_matrix(explicit_path: Path | None = None) -> Path:
    if explicit_path is not None:
        path = explicit_path.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"SOAP matrix not found: {path}")
        return path

    candidates = sorted(PACKAGE_ROOT.rglob("GB_SOAP_Al_Mg.npy"))
    if not candidates:
        raise FileNotFoundError(
            "Could not find GB_SOAP_Al_Mg.npy under point_selection. "
            "Download it from Zenodo record 4107058 and pass --soap /path/to/GB_SOAP_Al_Mg.npy."
        )
    candidates.sort(
        key=lambda p: (
            "machine_learning_notebook" not in str(p),
            len(str(p)),
            str(p),
        )
    )
    return candidates[0]


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


def compute_xgb(deltaE: np.ndarray) -> tuple[float, float]:
    X_c = solve_Xc(deltaE)
    p = compute_occupation(deltaE, X_c)
    return float(np.mean(p)), X_c


def load_inputs(soap_arg: Path | None = None) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, Path]:
    df = pd.read_csv(DATA_PATH)
    missing = [col for col in ["row_id", "site_id_if_verified", TARGET_COL] + FEATURE_COLS_10D if col not in df.columns]
    if missing:
        raise ValueError(f"Verified CSV is missing required columns: {missing}")
    if df[[TARGET_COL] + FEATURE_COLS_10D].isna().any().any():
        raise ValueError("Verified CSV contains NaNs in target or 10D PCA columns.")

    soap_path = find_soap_matrix(soap_arg)
    soap = np.load(soap_path)
    if soap.ndim != 2:
        raise ValueError(f"SOAP matrix must be 2D, got shape {soap.shape}")
    if soap.shape[0] != len(df):
        raise ValueError(
            "SOAP row count does not match verified CSV row count. "
            f"SOAP shape={soap.shape}, verified rows={len(df)}. "
            "Do not continue without an explicit alignment map."
        )
    if soap.shape[1] < max(PCA_DIMS):
        raise ValueError(f"SOAP feature dimension {soap.shape[1]} is smaller than requested PCA dim {max(PCA_DIMS)}")

    print(f"SOAP matrix path: {soap_path}")
    print(f"SOAP matrix shape: {soap.shape}, dtype={soap.dtype}")
    print(f"Verified CSV shape: {df.shape}")
    print("Alignment check: SOAP rows == verified CSV rows; existing PCA audit confirms 10D exact match.")
    return df, soap, df[TARGET_COL].to_numpy(dtype=float), soap_path


def verify_10d_against_csv(x_pca_10: np.ndarray, df: pd.DataFrame) -> dict[str, float | str]:
    csv_pcs = df[FEATURE_COLS_10D].to_numpy(dtype=float)
    signs = []
    aligned = x_pca_10.copy()
    correlations = []
    for i in range(10):
        corr = float(np.corrcoef(aligned[:, i], csv_pcs[:, i])[0, 1])
        sign = 1.0 if corr >= 0 else -1.0
        aligned[:, i] *= sign
        signs.append(int(sign))
        correlations.append(abs(corr))
    max_abs_diff = float(np.max(np.abs(aligned - csv_pcs)))
    rmse = float(np.sqrt(np.mean((aligned - csv_pcs) ** 2)))
    return {
        "signs": ",".join(map(str, signs)),
        "mean_abs_corr": float(np.mean(correlations)),
        "max_abs_diff_after_sign_alignment": max_abs_diff,
        "rmse_after_sign_alignment": rmse,
    }


def fit_pca_and_evaluate(
    soap: np.ndarray,
    y: np.ndarray,
    n_components: int,
    xgb_true: float,
) -> tuple[dict[str, float | int | str], np.ndarray, np.ndarray, PCA, float]:
    start = time.perf_counter()
    pca = PCA(n_components=n_components, svd_solver="full")
    X_pca = pca.fit_transform(soap)
    model = LinearRegression()
    model.fit(X_pca, y)
    y_pred = model.predict(X_pca)
    residual = y_pred - y
    tail_cutoff = np.percentile(y, 10)
    tail_mask = y <= tail_cutoff
    xgb_pred, _ = compute_xgb(y_pred)
    runtime = time.perf_counter() - start
    row = {
        "n_components": n_components,
        "cumulative_explained_variance": float(np.sum(pca.explained_variance_ratio_)),
        "n_training_sites": int(len(y)),
        "n_evaluation_sites": int(len(y)),
        "model": "sklearn.linear_model.LinearRegression",
        "standard_scaler_used": "No",
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "r2": float(r2_score(y, y_pred)),
        "tail_mae": float(np.mean(np.abs(residual[tail_mask]))),
        "xgb_true": float(xgb_true),
        "xgb_pred": float(xgb_pred),
        "xgb_abs_error": float(abs(xgb_pred - xgb_true)),
        "runtime_seconds": float(runtime),
    }
    return row, X_pca, y_pred, pca, runtime


def save_predictions(df: pd.DataFrame, y: np.ndarray, y_pred: np.ndarray, n_components: int) -> None:
    out = pd.DataFrame(
        {
            "row_id": df["row_id"].to_numpy(),
            "site_id_if_verified": df["site_id_if_verified"].to_numpy(),
            "true_deltaE": y,
            "predicted_deltaE": y_pred,
            "residual": y_pred - y,
        }
    )
    out.to_csv(PREDICTIONS_DIR / f"predictions_all_data_{n_components}pc_lr.csv", index=False)


def make_single_plot(results: pd.DataFrame, y_col: str, ylabel: str, title: str, filename: str) -> Path:
    fig, ax = plt.subplots(figsize=(6.8, 4.9), dpi=300)
    ax.plot(results["n_components"], results[y_col], marker="o", linewidth=2.0, color="#2563eb")
    ax.set_xlabel("Number of PCA components")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(PCA_DIMS)
    ax.grid(True, color="#d8dee9", linewidth=0.7, alpha=0.85)
    fig.tight_layout()
    out = FIGURES_DIR / filename
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out


def make_figures(results: pd.DataFrame) -> list[Path]:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    paths = [
        make_single_plot(
            results,
            "mae",
            "Global MAE in DeltaE prediction (kJ/mol)",
            "PCA components vs global MAE",
            "pca_components_vs_global_mae.png",
        ),
        make_single_plot(
            results,
            "tail_mae",
            "Tail MAE, lowest-energy 10% (kJ/mol)",
            "PCA components vs low-energy tail MAE",
            "pca_components_vs_tail_mae.png",
        ),
        make_single_plot(
            results,
            "xgb_abs_error",
            "XGB_abs_error",
            "PCA components vs XGB error",
            "pca_components_vs_xgb_error.png",
        ),
        make_single_plot(
            results,
            "cumulative_explained_variance",
            "Cumulative explained variance ratio",
            "PCA explained variance from 10 to 100 components",
            "pca_explained_variance_10_to_100.png",
        ),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), dpi=300, constrained_layout=True)
    specs = [
        ("mae", "Global MAE (kJ/mol)", "Global MAE"),
        ("tail_mae", "Tail MAE (kJ/mol)", "Tail MAE"),
        ("xgb_abs_error", "XGB_abs_error", "XGB error"),
    ]
    for ax, (col, ylabel, title) in zip(axes, specs):
        ax.plot(results["n_components"], results[col], marker="o", linewidth=2.0, color="#2563eb")
        ax.set_xlabel("PCA components")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xticks(PCA_DIMS)
        ax.grid(True, color="#d8dee9", linewidth=0.7, alpha=0.85)
    fig.suptitle("All-data LinearRegression baseline vs PCA feature dimension", fontsize=15, y=1.04)
    summary_path = FIGURES_DIR / "pca_components_summary.png"
    fig.savefig(summary_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    paths.append(summary_path)
    return paths


def write_report(
    results: pd.DataFrame,
    soap_path: Path,
    soap_shape: tuple[int, int],
    csv_shape: tuple[int, int],
    alignment_check: dict[str, float | str],
) -> None:
    row10 = results[results["n_components"] == 10].iloc[0]
    row100 = results[results["n_components"] == 100].iloc[0]
    mae_improvement = float(row10["mae"] - row100["mae"])
    tail_improvement = float(row10["tail_mae"] - row100["tail_mae"])
    xgb_change = float(row100["xgb_abs_error"] - row10["xgb_abs_error"])
    ev_increase = float(row100["cumulative_explained_variance"] - row10["cumulative_explained_variance"])
    substantial = mae_improvement > 0.1 or tail_improvement > 0.1

    table_lines = [
        "| n_components | cumulative explained variance | global MAE | tail MAE | XGB_abs_error | R2 |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in results.iterrows():
        table_lines.append(
            f"| {int(r['n_components'])} | {r['cumulative_explained_variance']:.9f} | "
            f"{r['mae']:.6f} | {r['tail_mae']:.6f} | {r['xgb_abs_error']:.8f} | {r['r2']:.6f} |"
        )

    interpretation = (
        "Increasing from 10 to 100 PCs substantially lowers at least one MAE metric, which suggests that 10D PCA discards some target-relevant information for this all-data linear baseline."
        if substantial
        else "Increasing from 10 to 100 PCs does not substantially lower the MAE metrics, which suggests that 10D PCA is not the dominant feature bottleneck for this all-data linear baseline."
    )

    report = f"""# PCA Feature-Dimension Baseline: 10D To 100D

## A. What Was Done

- Used aligned `GB_SOAP_Al_Mg.npy`.
- SOAP path: `{soap_path}`.
- SOAP shape: `{soap_shape}`.
- Verified target CSV: `{DATA_PATH}`.
- Verified CSV shape: `{csv_shape}`.
- Fit PCA directly on SOAP with `n_components = 10, 20, 50, 100`.
- Used `PCA(..., svd_solver="full")`.
- Used no `StandardScaler`.
- Used no whitening.
- Trained all-data `LinearRegression()` for each PCA dimension.
- Evaluated predictions on all aligned GB sites.

Alignment check:

- SOAP row count equals verified CSV row count: `{soap_shape[0]} == {csv_shape[0]}`.
- Existing PCA audit classified the 10D PCA coordinates as an exact match.
- This script rechecked reproduced 10D PCA against verified `pc1...pc10` after sign alignment:
  - signs: `{alignment_check['signs']}`
  - mean absolute PC correlation: `{alignment_check['mean_abs_corr']:.12f}`
  - max abs difference after sign alignment: `{alignment_check['max_abs_diff_after_sign_alignment']:.12g}`
  - RMSE after sign alignment: `{alignment_check['rmse_after_sign_alignment']:.12g}`

## B. Why This Matters

PCA retains variance in SOAP feature space, but feature variance is not the same as target relevance. Extra PCs beyond 10 have low variance, but they may still contain information relevant to segregation energy, especially for rare local environments or the low-energy tail. This experiment tests whether the original 10D PCA representation is a feature bottleneck for an all-data linear model.

## C. Results Table

{chr(10).join(table_lines)}

## D. Compare 10D With 100D

- global MAE improvement, 10D - 100D: `{mae_improvement:.6f}` kJ/mol
- tail MAE improvement, 10D - 100D: `{tail_improvement:.6f}` kJ/mol
- XGB_abs_error change, 100D - 10D: `{xgb_change:.8f}`
- explained variance increase, 100D - 10D: `{ev_increase:.9f}`

## E. Bottleneck Interpretation

1. Does increasing PCA components from 10 to 100 lower global MAE? `{mae_improvement > 0}`.
2. Does it lower tail MAE? `{tail_improvement > 0}`.
3. Does it lower XGB_abs_error? `{xgb_change < 0}`.
4. Is 10D PCA likely a feature bottleneck? {interpretation}
5. Is the improvement large enough to matter compared with selection-strategy differences? Compare the MAE improvement above with the selection-strategy gaps in the extended-budget experiment; if it is small, selection/model-class effects may dominate, while a large drop would indicate feature dimensionality matters.

## F. Important Caution

- This is an all-data in-sample `LinearRegression` feature baseline.
- It is not an independent generalization benchmark.
- More PCA variance does not necessarily mean more predictive information.
- `LinearRegression()` minimizes squared error, not MAE, tail MAE, or XGB_abs_error.

## G. Recommended Presentation Wording

The 10D raw-PCA representation was tested against 20D, 50D, and 100D PCA features using the same all-data `LinearRegression` model. This checks whether low-variance SOAP PCs contain target-relevant information that the original accelerated model discards. Because this is an in-sample all-data feature baseline, it should be interpreted as a diagnostic of representation capacity rather than a generalization benchmark. If the 100D baseline does not improve much over 10D, the remaining error is more likely due to the linear model class or target metric than to PCA dimensionality alone.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--soap",
        type=Path,
        default=None,
        help="Path to GB_SOAP_Al_Mg.npy from Zenodo record 4107058.",
    )
    args = parser.parse_args()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PREDICTIONS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df, soap, y, soap_path = load_inputs(args.soap)
    xgb_true, _ = compute_xgb(y)
    rows = []
    alignment_check: dict[str, float | str] | None = None
    completed = []

    for n_components in PCA_DIMS:
        print(f"\nRunning PCA dimension {n_components}")
        row, X_pca, y_pred, pca, runtime = fit_pca_and_evaluate(soap, y, n_components, xgb_true)
        rows.append(row)
        completed.append(n_components)
        print(
            f"  done in {runtime:.1f}s | EV={row['cumulative_explained_variance']:.9f}, "
            f"MAE={row['mae']:.6f}, tail={row['tail_mae']:.6f}, XGBerr={row['xgb_abs_error']:.8f}"
        )
        if n_components == 10:
            alignment_check = verify_10d_against_csv(X_pca, df)
            print(
                "  10D alignment check: "
                f"mean_abs_corr={alignment_check['mean_abs_corr']:.12f}, "
                f"max_abs_diff={alignment_check['max_abs_diff_after_sign_alignment']:.3g}"
            )
            save_predictions(df, y, y_pred, n_components=10)
        if n_components == 100:
            save_predictions(df, y, y_pred, n_components=100)

    results = pd.DataFrame(rows)
    results.to_csv(RESULTS_PATH, index=False)
    fig_paths = make_figures(results)
    if alignment_check is None:
        raise RuntimeError("10D alignment check was not computed.")
    write_report(results, soap_path, tuple(soap.shape), tuple(df.shape), alignment_check)

    row10 = results[results["n_components"] == 10].iloc[0]
    row100 = results[results["n_components"] == 100].iloc[0]
    mae_improvement = float(row10["mae"] - row100["mae"])
    tail_improvement = float(row10["tail_mae"] - row100["tail_mae"])
    xgb_change = float(row100["xgb_abs_error"] - row10["xgb_abs_error"])
    substantially_lowers = mae_improvement > 0.1 or tail_improvement > 0.1

    print("\nConsole summary")
    print("---------------")
    print(f"SOAP matrix path: {soap_path}")
    print(f"SOAP matrix shape: {soap.shape}")
    print(f"Verified CSV shape: {df.shape}")
    print(f"PCA dimensions completed: {completed}")
    print(f"10D cumulative explained variance: {row10['cumulative_explained_variance']:.9f}")
    print(f"100D cumulative explained variance: {row100['cumulative_explained_variance']:.9f}")
    print(
        f"10D MAE / tail MAE / XGB error: "
        f"{row10['mae']:.6f} / {row10['tail_mae']:.6f} / {row10['xgb_abs_error']:.8f}"
    )
    print(
        f"100D MAE / tail MAE / XGB error: "
        f"{row100['mae']:.6f} / {row100['tail_mae']:.6f} / {row100['xgb_abs_error']:.8f}"
    )
    print(f"100D substantially lowers baseline: {substantially_lowers}")
    print(f"10D -> 100D global MAE improvement: {mae_improvement:.6f}")
    print(f"10D -> 100D tail MAE improvement: {tail_improvement:.6f}")
    print(f"10D -> 100D XGB_abs_error change: {xgb_change:.8f}")
    print(f"Report path: {REPORT_PATH}")
    print("Figures:")
    for path in fig_paths:
        print(f"  {path}")


if __name__ == "__main__":
    main()
