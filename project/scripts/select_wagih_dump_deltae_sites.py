#!/usr/bin/env python3
"""Select substitution sites directly from a Wagih/Zenodo segregation dump."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


ENERGY_COLUMNS = (
    "seg_kJ_per_mol",
    "delta_e_kj_mol",
    "deltaE",
    "DeltaE",
    "dEseg",
    "Eseg",
    "segregation_energy",
)


def read_dump(path: Path, energy_column: str | None) -> tuple[list[dict[str, float]], np.ndarray, list[str], str]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    bounds: list[tuple[float, float]] = []
    columns: list[str] | None = None
    atoms_start: int | None = None
    for i, line in enumerate(lines):
        if line.startswith("ITEM: BOX BOUNDS"):
            bounds = [(float(lines[i + j].split()[0]), float(lines[i + j].split()[1])) for j in range(1, 4)]
        elif line.startswith("ITEM: ATOMS"):
            columns = line.split()[2:]
            atoms_start = i + 1
            break
    if columns is None or atoms_start is None:
        raise ValueError(f"No ITEM: ATOMS section found in {path}")
    if len(bounds) != 3:
        raise ValueError(f"No 3D BOX BOUNDS section found in {path}")

    selected_energy = energy_column
    if selected_energy is None:
        selected_energy = next((name for name in ENERGY_COLUMNS if name in columns), None)
    if selected_energy is None or selected_energy not in columns:
        raise ValueError(f"Could not find energy column in dump columns: {columns}")

    id_idx = columns.index("id")
    energy_idx = columns.index(selected_energy)
    coord_names = []
    for axis in "xyz":
        if axis in columns:
            coord_names.append(axis)
        elif f"{axis}s" in columns:
            coord_names.append(f"{axis}s")
        else:
            raise ValueError(f"Could not find {axis} or {axis}s coordinate column in {columns}")
    coord_idx = [columns.index(name) for name in coord_names]
    lo = np.array([item[0] for item in bounds], dtype=float)
    hi = np.array([item[1] for item in bounds], dtype=float)
    box = hi - lo

    rows: list[dict[str, float]] = []
    for line in lines[atoms_start:]:
        if line.startswith("ITEM:"):
            break
        parts = line.split()
        if len(parts) <= max(id_idx, energy_idx, *coord_idx):
            continue
        xyz = []
        for value, name, axis_lo, axis_box in zip((float(parts[idx]) for idx in coord_idx), coord_names, lo, box):
            xyz.append(axis_lo + value * axis_box if name.endswith("s") else value)
        rows.append(
            {
                "atom_id": int(float(parts[id_idx])),
                "x_A": xyz[0],
                "y_A": xyz[1],
                "z_A": xyz[2],
                "wagih_delta_e_kj_mol": float(parts[energy_idx]),
            }
        )
    return rows, box, columns, selected_energy


def select_rows(
    rows: list[dict[str, float]],
    box: np.ndarray,
    samples: int,
    bulk_references: int,
    seed: int,
    bulk_min_distance: float,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    energy = np.array([row["wagih_delta_e_kj_mol"] for row in rows], dtype=float)
    positions = np.array([[row["x_A"], row["y_A"], row["z_A"]] for row in rows], dtype=float)
    positions = positions % box

    gb_indices = np.flatnonzero(energy != 0.0)
    bulk_indices = np.flatnonzero(energy == 0.0)
    if len(gb_indices) == 0:
        raise ValueError("No nonzero Wagih segregation-energy sites found.")
    if len(bulk_indices) == 0:
        raise ValueError("No zero-energy bulk-reference candidates found.")

    chosen_gb = rng.choice(gb_indices, size=min(samples, len(gb_indices)), replace=False)
    gb_tree = cKDTree(positions[gb_indices], boxsize=box)
    distance_to_gb, _ = gb_tree.query(positions[bulk_indices], k=1)
    eligible_bulk_mask = distance_to_gb >= bulk_min_distance
    eligible_bulk = bulk_indices[eligible_bulk_mask]
    eligible_distances = distance_to_gb[eligible_bulk_mask]
    if len(eligible_bulk) < bulk_references:
        raise ValueError(
            f"Only {len(eligible_bulk)} zero-energy bulk candidates are at least "
            f"{bulk_min_distance:.3f} A from a nonzero Wagih site; requested {bulk_references}. "
            "Try a smaller --bulk-min-distance."
        )
    chosen_bulk = eligible_bulk[np.argsort(eligible_distances)[::-1]][:bulk_references]
    bulk_distance_by_idx = {int(idx): float(distance_to_gb[pos]) for pos, idx in enumerate(bulk_indices)}

    selected: list[dict[str, object]] = []
    for rank, idx in enumerate(chosen_bulk, start=1):
        row = rows[int(idx)]
        selected.append(
            {
                "role": "bulk_reference",
                "atom_id": int(row["atom_id"]),
                "x_A": row["x_A"],
                "y_A": row["y_A"],
                "z_A": row["z_A"],
                "is_gb": 0,
                "first_neighbor_count": "",
                "distance_to_nearest_gb_A": bulk_distance_by_idx[int(idx)],
                "bulk_distance_rank": rank,
                "wagih_delta_e_kj_mol": row["wagih_delta_e_kj_mol"],
            }
        )
    for idx in chosen_gb:
        row = rows[int(idx)]
        selected.append(
            {
                "role": "gb_sample",
                "atom_id": int(row["atom_id"]),
                "x_A": row["x_A"],
                "y_A": row["y_A"],
                "z_A": row["z_A"],
                "is_gb": 1,
                "first_neighbor_count": "",
                "distance_to_nearest_gb_A": 0.0,
                "bulk_distance_rank": "",
                "wagih_delta_e_kj_mol": row["wagih_delta_e_kj_mol"],
            }
        )
    return selected


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "role",
        "atom_id",
        "x_A",
        "y_A",
        "z_A",
        "is_gb",
        "first_neighbor_count",
        "distance_to_nearest_gb_A",
        "bulk_distance_rank",
        "wagih_delta_e_kj_mol",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--energy-column")
    parser.add_argument("--samples", type=int, default=1000)
    parser.add_argument("--bulk-references", type=int, default=50)
    parser.add_argument(
        "--bulk-min-distance",
        type=float,
        default=8.0,
        help="Minimum distance in Angstrom from a bulk reference to the nearest nonzero Wagih site.",
    )
    parser.add_argument("--seed", type=int, default=500)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows, box, columns, selected_energy = read_dump(args.dump, args.energy_column)
    selected = select_rows(rows, box, args.samples, args.bulk_references, args.seed, args.bulk_min_distance)
    write_csv(args.output, selected)
    gb_count = sum(1 for row in rows if row["wagih_delta_e_kj_mol"] != 0.0)
    bulk_count = len(rows) - gb_count
    print(f"Dump columns: {' '.join(columns)}")
    print(f"Used energy column: {selected_energy}")
    print(f"Nonzero Wagih GB-like sites: {gb_count}")
    print(f"Zero-energy bulk candidates: {bulk_count}")
    print(f"Bulk minimum distance: {args.bulk_min_distance} A")
    print(f"Wrote {len(selected)} selected sites to {args.output}")


if __name__ == "__main__":
    main()
