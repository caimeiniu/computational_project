#!/usr/bin/env python3
"""Export Delta E results to PNG graphs and an Excel-friendly workbook.

The script reads the CSV produced by analyze_cuni_deltae.py and writes:
  - histogram comparison PNG
  - ranked Delta E PNG
  - summary CSV
  - .xlsx workbook when pandas/openpyxl or xlsxwriter is available
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ks_2samp, kstest, skewnorm


WAGIH_NI_CU = {
    "mu_kj_mol": -2.0,
    "sigma_kj_mol": 8.0,
    "alpha": 0.0,
    "range_min_kj_mol": -30.0,
    "range_max_kj_mol": 30.0,
}


def read_results(path: Path) -> list[dict[str, str]]:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError("Reading Excel input requires pandas/openpyxl.") from exc
        data = pd.read_excel(path, sheet_name="DeltaE_Data")
        return data.to_dict(orient="records")

    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def numeric_column(rows: list[dict[str, str]], key: str) -> np.ndarray:
    return np.array([float(row[key]) for row in rows], dtype=float)


def summary_rows(deltae: np.ndarray, wagih_raw: np.ndarray | None = None) -> list[tuple[str, float | int | str]]:
    ours_alpha, ours_mu, ours_sigma = skewnorm.fit(deltae)
    outside = np.mean(
        (deltae < WAGIH_NI_CU["range_min_kj_mol"])
        | (deltae > WAGIH_NI_CU["range_max_kj_mol"])
    )
    rows: list[tuple[str, float | int | str]] = [
        ("n", len(deltae)),
        ("mean_kj_mol", float(np.mean(deltae))),
        ("std_kj_mol", float(np.std(deltae, ddof=1))),
        ("min_kj_mol", float(np.min(deltae))),
        ("max_kj_mol", float(np.max(deltae))),
        ("median_kj_mol", float(np.median(deltae))),
        ("q05_kj_mol", float(np.quantile(deltae, 0.05))),
        ("q95_kj_mol", float(np.quantile(deltae, 0.95))),
        ("ours_fit_mu_kj_mol", float(ours_mu)),
        ("ours_fit_sigma_kj_mol", float(ours_sigma)),
        ("ours_fit_alpha", float(ours_alpha)),
        ("wagih_mu_kj_mol", WAGIH_NI_CU["mu_kj_mol"]),
        ("wagih_sigma_kj_mol", WAGIH_NI_CU["sigma_kj_mol"]),
        ("wagih_range_min_kj_mol", WAGIH_NI_CU["range_min_kj_mol"]),
        ("wagih_range_max_kj_mol", WAGIH_NI_CU["range_max_kj_mol"]),
        ("fraction_outside_wagih_range", float(outside)),
    ]
    if wagih_raw is not None:
        wagih_alpha, wagih_mu, wagih_sigma = skewnorm.fit(wagih_raw)
        ks = ks_2samp(deltae, wagih_raw)
        rows.extend(
            [
                ("comparison_type", "raw_wagih_two_sample"),
                ("wagih_raw_n", len(wagih_raw)),
                ("wagih_raw_mean_kj_mol", float(np.mean(wagih_raw))),
                ("wagih_raw_std_kj_mol", float(np.std(wagih_raw, ddof=1))),
                ("wagih_raw_min_kj_mol", float(np.min(wagih_raw))),
                ("wagih_raw_max_kj_mol", float(np.max(wagih_raw))),
                ("wagih_raw_fit_mu_kj_mol", float(wagih_mu)),
                ("wagih_raw_fit_sigma_kj_mol", float(wagih_sigma)),
                ("wagih_raw_fit_alpha", float(wagih_alpha)),
                ("ks_d_vs_wagih_raw", float(ks.statistic)),
                ("ks_p_vs_wagih_raw", float(ks.pvalue)),
            ]
        )
    else:
        ks = kstest(
            deltae,
            lambda x: skewnorm.cdf(
                x,
                WAGIH_NI_CU["alpha"],
                loc=WAGIH_NI_CU["mu_kj_mol"],
                scale=WAGIH_NI_CU["sigma_kj_mol"],
            ),
        )
        rows.extend(
            [
                ("comparison_type", "approx_wagih_fit"),
                ("ks_d_vs_wagih_approx", float(ks.statistic)),
                ("ks_p_vs_wagih_approx", float(ks.pvalue)),
            ]
        )
    return rows


def write_summary_csv(path: Path, rows: list[tuple[str, float | int | str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def plot_histogram(path: Path, deltae: np.ndarray, title: str, wagih_raw: np.ndarray | None = None) -> None:
    if wagih_raw is not None:
        x_min = min(deltae.min(), wagih_raw.min()) - 5.0
        x_max = max(deltae.max(), wagih_raw.max()) + 5.0
    else:
        x_min = min(deltae.min(), WAGIH_NI_CU["range_min_kj_mol"]) - 5.0
        x_max = max(deltae.max(), WAGIH_NI_CU["range_max_kj_mol"]) + 5.0
    x = np.linspace(x_min, x_max, 500)
    ours_alpha, ours_mu, ours_sigma = skewnorm.fit(deltae)
    if wagih_raw is not None:
        wagih_alpha, wagih_mu, wagih_sigma = skewnorm.fit(wagih_raw)
        ks = ks_2samp(deltae, wagih_raw)
    else:
        wagih_alpha = WAGIH_NI_CU["alpha"]
        wagih_mu = WAGIH_NI_CU["mu_kj_mol"]
        wagih_sigma = WAGIH_NI_CU["sigma_kj_mol"]
        ks = kstest(
            deltae,
            lambda values: skewnorm.cdf(values, wagih_alpha, loc=wagih_mu, scale=wagih_sigma),
        )
    mean = np.mean(deltae)
    std = np.std(deltae, ddof=1)

    fig, ax = plt.subplots(figsize=(9.2, 5.8), constrained_layout=True)
    if wagih_raw is not None:
        ax.hist(
            wagih_raw,
            bins=50,
            density=True,
            alpha=0.48,
            color="#73C995",
            edgecolor="white",
            linewidth=0.35,
            label=f"Wagih (n={len(wagih_raw)})",
        )
    ax.hist(
        deltae,
        bins=35,
        density=True,
        alpha=0.58,
        color="#E98B91",
        edgecolor="white",
        linewidth=0.5,
        label=f"Ours production (n={len(deltae)})",
    )
    ax.plot(
        x,
        skewnorm.pdf(
            x,
            wagih_alpha,
            loc=wagih_mu,
            scale=wagih_sigma,
        ),
        color="#157A35",
        linewidth=2.2,
        label=f"Wagih fit  μ={wagih_mu:+.1f} σ={wagih_sigma:.1f} α={wagih_alpha:+.2f}",
    )
    ax.plot(
        x,
        skewnorm.pdf(x, ours_alpha, loc=ours_mu, scale=ours_sigma),
        color="#A51D2D",
        linewidth=2.2,
        linestyle="--",
        label=f"Ours fit  μ={ours_mu:+.1f}  σ={ours_sigma:.1f}  α={ours_alpha:+.2f}",
    )
    if wagih_raw is None:
        ax.axvspan(
            WAGIH_NI_CU["range_min_kj_mol"],
            WAGIH_NI_CU["range_max_kj_mol"],
            color="#111111",
            alpha=0.07,
            label="Wagih (n=N/A)",
        )
    ax.axvline(mean, color="#A51D2D", linewidth=1.3, alpha=0.65, linestyle=":")
    ax.set_title(f"{title} (KS D={ks.statistic:.3f}, p={ks.pvalue:.3g})", fontsize=15)
    ax.set_xlabel(r"$\Delta E_{\rm seg}$ (kJ/mol)")
    ax.set_ylabel("Probability density")
    handles, labels = ax.get_legend_handles_labels()
    wagih_label = next(label for label in labels if label.startswith("Wagih (n="))
    order = [
        labels.index(wagih_label),
        labels.index(f"Ours production (n={len(deltae)})"),
        next(i for i, label in enumerate(labels) if label.startswith("Wagih fit")),
        next(i for i, label in enumerate(labels) if label.startswith("Ours fit")),
    ]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], frameon=True, fontsize=9.2, loc="upper left")
    ax.grid(axis="y", alpha=0.2)

    stat_text = "\n".join(
        [
            f"n = {len(deltae)}",
            f"mean = {mean:+.2f} kJ/mol",
            f"std = {std:.2f} kJ/mol",
            f"min/max = {deltae.min():+.1f} / {deltae.max():+.1f}",
            f"outside Wagih range = {np.mean((deltae < WAGIH_NI_CU['range_min_kj_mol']) | (deltae > WAGIH_NI_CU['range_max_kj_mol'])):.1%}",
        ]
    )
    ax.text(
        0.985,
        0.965,
        stat_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#CCCCCC", "alpha": 0.9},
    )
    fig.savefig(path, dpi=240)
    plt.close(fig)


def plot_ranked(path: Path, deltae: np.ndarray, title: str) -> None:
    ranked = np.sort(deltae)
    fig, ax = plt.subplots(figsize=(7.4, 4.6), constrained_layout=True)
    ax.plot(np.arange(1, len(ranked) + 1), ranked, color="#4C78A8", linewidth=2.0)
    ax.axhline(WAGIH_NI_CU["range_min_kj_mol"], color="#111111", alpha=0.35, linestyle=":")
    ax.axhline(WAGIH_NI_CU["range_max_kj_mol"], color="#111111", alpha=0.35, linestyle=":")
    ax.axhline(0.0, color="#E45756", alpha=0.65, linewidth=1.2)
    ax.set_title(title)
    ax.set_xlabel("Ranked sampled site")
    ax.set_ylabel("Delta E_seg (kJ/mol)")
    ax.grid(alpha=0.2)
    fig.savefig(path, dpi=240)
    plt.close(fig)


def write_xlsx(path: Path, rows: list[dict[str, str]], summary: list[tuple[str, float | int | str]]) -> bool:
    try:
        import pandas as pd
    except ImportError:
        return False

    data = pd.DataFrame(rows)
    for column in data.columns:
        data[column] = pd.to_numeric(data[column], errors="ignore")
    summary_df = pd.DataFrame(summary, columns=["metric", "value"])
    wagih_df = pd.DataFrame([WAGIH_NI_CU])

    try:
        with pd.ExcelWriter(path, engine="openpyxl") as writer:
            summary_df.to_excel(writer, sheet_name="Summary", index=False)
            data.to_excel(writer, sheet_name="DeltaE_Data", index=False)
            wagih_df.to_excel(writer, sheet_name="Wagih_Reference", index=False)
        return True
    except Exception:
        try:
            with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
                summary_df.to_excel(writer, sheet_name="Summary", index=False)
                data.to_excel(writer, sheet_name="DeltaE_Data", index=False)
                wagih_df.to_excel(writer, sheet_name="Wagih_Reference", index=False)
            return True
        except Exception:
            return False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="cuni/nicu Delta E result CSV or Excel workbook.")
    parser.add_argument("--wagih-results", type=Path, help="Optional raw Wagih Delta E CSV or Excel workbook.")
    parser.add_argument("--outdir", type=Path, default=Path("nicu_export"))
    parser.add_argument("--prefix", default="nicu_3d")
    parser.add_argument("--alloy-label", default="Ni(Cu)", help="Label for plot titles, e.g. 'Au(Pt)' or 'Ni(Cu)'.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = read_results(args.results)
    if not rows:
        raise ValueError(f"No rows found in {args.results}")
    deltae = numeric_column(rows, "delta_e_kj_mol")
    wagih_raw = None
    if args.wagih_results:
        wagih_rows = read_results(args.wagih_results)
        wagih_raw = numeric_column(wagih_rows, "delta_e_kj_mol")
    summary = summary_rows(deltae, wagih_raw=wagih_raw)

    summary_csv = args.outdir / f"{args.prefix}_summary.csv"
    hist_png = args.outdir / f"{args.prefix}_deltae_histogram_vs_wagih.png"
    ranked_png = args.outdir / f"{args.prefix}_ranked_deltae.png"
    xlsx = args.outdir / f"{args.prefix}_deltae_analysis.xlsx"

    write_summary_csv(summary_csv, summary)
    plot_histogram(hist_png, deltae, f"{args.alloy_label} Delta E Spectrum vs Wagih", wagih_raw=wagih_raw)
    plot_ranked(ranked_png, deltae, f"{args.alloy_label} Ranked Segregation Energies")
    wrote_xlsx = write_xlsx(xlsx, rows, summary)

    print(f"Wrote {summary_csv}")
    print(f"Wrote {hist_png}")
    print(f"Wrote {ranked_png}")
    if wrote_xlsx:
        print(f"Wrote {xlsx}")
    else:
        print("Could not write .xlsx because pandas plus openpyxl/xlsxwriter is unavailable.")
        print("The summary CSV and original results CSV can still be opened in Excel.")


if __name__ == "__main__":
    main()
