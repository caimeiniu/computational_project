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
from scipy.stats import skewnorm


WAGIH_NI_CU = {
    "mu_kj_mol": -2.0,
    "sigma_kj_mol": 8.0,
    "alpha": 0.0,
    "range_min_kj_mol": -30.0,
    "range_max_kj_mol": 30.0,
}


def read_results(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def numeric_column(rows: list[dict[str, str]], key: str) -> np.ndarray:
    return np.array([float(row[key]) for row in rows], dtype=float)


def summary_rows(deltae: np.ndarray) -> list[tuple[str, float | int | str]]:
    outside = np.mean(
        (deltae < WAGIH_NI_CU["range_min_kj_mol"])
        | (deltae > WAGIH_NI_CU["range_max_kj_mol"])
    )
    return [
        ("n", len(deltae)),
        ("mean_kj_mol", float(np.mean(deltae))),
        ("std_kj_mol", float(np.std(deltae, ddof=1))),
        ("min_kj_mol", float(np.min(deltae))),
        ("max_kj_mol", float(np.max(deltae))),
        ("median_kj_mol", float(np.median(deltae))),
        ("q05_kj_mol", float(np.quantile(deltae, 0.05))),
        ("q95_kj_mol", float(np.quantile(deltae, 0.95))),
        ("wagih_mu_kj_mol", WAGIH_NI_CU["mu_kj_mol"]),
        ("wagih_sigma_kj_mol", WAGIH_NI_CU["sigma_kj_mol"]),
        ("wagih_range_min_kj_mol", WAGIH_NI_CU["range_min_kj_mol"]),
        ("wagih_range_max_kj_mol", WAGIH_NI_CU["range_max_kj_mol"]),
        ("fraction_outside_wagih_range", float(outside)),
    ]


def write_summary_csv(path: Path, rows: list[tuple[str, float | int | str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def plot_histogram(path: Path, deltae: np.ndarray, title: str) -> None:
    x_min = min(deltae.min(), WAGIH_NI_CU["range_min_kj_mol"]) - 5.0
    x_max = max(deltae.max(), WAGIH_NI_CU["range_max_kj_mol"]) + 5.0
    x = np.linspace(x_min, x_max, 500)

    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    ax.hist(deltae, bins=35, density=True, alpha=0.72, color="#4C78A8", label="Our Ni(Cu), 500 sites")
    ax.plot(
        x,
        skewnorm.pdf(
            x,
            WAGIH_NI_CU["alpha"],
            loc=WAGIH_NI_CU["mu_kj_mol"],
            scale=WAGIH_NI_CU["sigma_kj_mol"],
        ),
        color="#111111",
        linewidth=2.2,
        label="Wagih Ni(Cu) approx. fit",
    )
    ax.axvspan(
        WAGIH_NI_CU["range_min_kj_mol"],
        WAGIH_NI_CU["range_max_kj_mol"],
        color="#111111",
        alpha=0.07,
        label="Wagih approx. range",
    )
    ax.axvline(np.mean(deltae), color="#E45756", linewidth=1.7, linestyle="--", label="Our mean")
    ax.set_title(title)
    ax.set_xlabel("Delta E_seg (kJ/mol)")
    ax.set_ylabel("Probability density")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", alpha=0.2)
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
    parser.add_argument("--results", type=Path, required=True, help="cuni/nicu Delta E result CSV.")
    parser.add_argument("--outdir", type=Path, default=Path("nicu_export"))
    parser.add_argument("--prefix", default="nicu_3d")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = read_results(args.results)
    if not rows:
        raise ValueError(f"No rows found in {args.results}")
    deltae = numeric_column(rows, "delta_e_kj_mol")
    summary = summary_rows(deltae)

    summary_csv = args.outdir / f"{args.prefix}_summary.csv"
    hist_png = args.outdir / f"{args.prefix}_deltae_histogram_vs_wagih.png"
    ranked_png = args.outdir / f"{args.prefix}_ranked_deltae.png"
    xlsx = args.outdir / f"{args.prefix}_deltae_analysis.xlsx"

    write_summary_csv(summary_csv, summary)
    plot_histogram(hist_png, deltae, "Ni(Cu) Delta E Spectrum vs Wagih")
    plot_ranked(ranked_png, deltae, "Ni(Cu) Ranked Segregation Energies")
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
