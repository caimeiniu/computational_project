#!/usr/bin/env python3
"""Collect Cu-Ni ΔE results and compare with Wagih/Larsen/Schuh Ni(Cu)."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skewnorm


EV_TO_KJ_MOL = 96.4853321233
FINAL_PE_RE = re.compile(r"FINAL_PE_EV\s+(\S+)\s+(\d+)\s+([-+0-9.eE]+)")
WAGIH_NI_CU = {
    "mu_kj_mol": -2.0,
    "sigma_kj_mol": 8.0,
    "alpha": 0.0,
    "range_kj_mol": (-30.0, 30.0),
}


def read_sites(path: Path) -> dict[int, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {int(row["atom_id"]): row for row in csv.DictReader(handle)}


def parse_logs(jobs_dir: Path) -> list[dict[str, object]]:
    rows = []
    for log_path in sorted(jobs_dir.glob("log_*.lammps")):
        text = log_path.read_text(encoding="utf-8", errors="replace")
        matches = FINAL_PE_RE.findall(text)
        if not matches:
            continue
        role, atom_id, pe = matches[-1]
        rows.append(
            {
                "role": role,
                "atom_id": int(atom_id),
                "pe_eV": float(pe),
                "log": str(log_path),
            }
        )
    return rows


def compute_deltae(rows: list[dict[str, object]], sites: dict[int, dict[str, str]]) -> list[dict[str, object]]:
    bulk = [row for row in rows if row["role"] == "bulk_reference"]
    if len(bulk) < 1:
        raise ValueError("Expected at least one bulk_reference log, found 0")
    bulk_energies = np.array([float(row["pe_eV"]) for row in bulk], dtype=float)
    bulk_pe = float(np.mean(bulk_energies))
    bulk_std = float(np.std(bulk_energies, ddof=1)) if len(bulk_energies) > 1 else 0.0
    output = []
    for row in rows:
        if row["role"] != "gb_sample":
            continue
        atom_id = int(row["atom_id"])
        site = sites.get(atom_id, {})
        delta_e = (float(row["pe_eV"]) - bulk_pe) * EV_TO_KJ_MOL
        output.append(
            {
                "atom_id": atom_id,
                "delta_e_kj_mol": delta_e,
                "pe_eV": row["pe_eV"],
                "bulk_reference_mean_pe_eV": bulk_pe,
                "bulk_reference_std_pe_eV": bulk_std,
                "bulk_reference_count": len(bulk_energies),
                "x_A": site.get("x_A", ""),
                "y_A": site.get("y_A", ""),
                "z_A": site.get("z_A", ""),
                "first_neighbor_count": site.get("first_neighbor_count", ""),
            }
        )
    return output


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "atom_id",
        "delta_e_kj_mol",
        "pe_eV",
        "bulk_reference_mean_pe_eV",
        "bulk_reference_std_pe_eV",
        "bulk_reference_count",
        "x_A",
        "y_A",
        "z_A",
        "first_neighbor_count",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_comparison(path: Path, deltae: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    x_min = min(deltae.min(initial=-30.0), WAGIH_NI_CU["range_kj_mol"][0]) - 5.0
    x_max = max(deltae.max(initial=30.0), WAGIH_NI_CU["range_kj_mol"][1]) + 5.0
    x = np.linspace(x_min, x_max, 500)
    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    ax.hist(deltae, bins=35, density=True, alpha=0.65, label="Our 3D Cu-Ni sample")
    ax.plot(
        x,
        skewnorm.pdf(
            x,
            WAGIH_NI_CU["alpha"],
            loc=WAGIH_NI_CU["mu_kj_mol"],
            scale=WAGIH_NI_CU["sigma_kj_mol"],
        ),
        color="black",
        linewidth=2.0,
        label="Wagih Ni(Cu) approx. fit",
    )
    ax.axvspan(*WAGIH_NI_CU["range_kj_mol"], color="0.2", alpha=0.08, label="Wagih approx. range")
    ax.set_xlabel("Delta E_seg (kJ/mol)")
    ax.set_ylabel("Probability density")
    ax.legend(frameon=False)
    ax.set_title("3D Cu-Ni Segregation-Energy Spectrum")
    fig.savefig(path, dpi=220)


def summarize(deltae: np.ndarray, results: list[dict[str, object]]) -> str:
    loc, scale = float(np.mean(deltae)), float(np.std(deltae, ddof=1)) if len(deltae) > 1 else 0.0
    outside = np.mean((deltae < WAGIH_NI_CU["range_kj_mol"][0]) | (deltae > WAGIH_NI_CU["range_kj_mol"][1]))
    bulk_count = int(results[0].get("bulk_reference_count", 1)) if results else 0
    bulk_std_ev = float(results[0].get("bulk_reference_std_pe_eV", 0.0)) if results else 0.0
    return "\n".join(
        [
            f"n = {len(deltae)}",
            f"bulk references = {bulk_count}",
            f"bulk reference std = {bulk_std_ev * EV_TO_KJ_MOL:.3f} kJ/mol",
            f"mean = {loc:.3f} kJ/mol",
            f"std = {scale:.3f} kJ/mol",
            f"min/max = {deltae.min():.3f} / {deltae.max():.3f} kJ/mol",
            f"Wagih Ni(Cu) approx: mu={WAGIH_NI_CU['mu_kj_mol']} kJ/mol, sigma={WAGIH_NI_CU['sigma_kj_mol']} kJ/mol, range={WAGIH_NI_CU['range_kj_mol']} kJ/mol",
            f"fraction outside Wagih approx range = {outside:.2%}",
        ]
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sites", type=Path, default=Path("data/cuni_3d/cuni_3d_deltae_sites.csv"))
    parser.add_argument("--jobs-dir", type=Path, default=Path("data/cuni_3d/lammps/deltae_jobs"))
    parser.add_argument("--output", type=Path, default=Path("data/cuni_3d/cuni_3d_deltae_results.csv"))
    parser.add_argument("--figure", type=Path, default=Path("data/cuni_3d/cuni_3d_deltae_vs_wagih.png"))
    parser.add_argument("--summary", type=Path, default=Path("data/cuni_3d/cuni_3d_deltae_summary.txt"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sites = read_sites(args.sites)
    parsed = parse_logs(args.jobs_dir)
    results = compute_deltae(parsed, sites)
    if not results:
        raise ValueError(f"No GB sample results found in {args.jobs_dir}")
    write_csv(args.output, results)
    deltae = np.array([float(row["delta_e_kj_mol"]) for row in results])
    plot_comparison(args.figure, deltae)
    summary = summarize(deltae, results)
    args.summary.write_text(summary + "\n", encoding="utf-8")
    print(summary)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.figure}")


if __name__ == "__main__":
    main()
