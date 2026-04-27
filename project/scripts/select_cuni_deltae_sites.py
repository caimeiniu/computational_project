#!/usr/bin/env python3
"""Select GB and bulk sites for 3D Cu-Ni segregation-energy calculations."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def read_lammps_atoms(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8").splitlines()
    lo = np.zeros(3, dtype=float)
    hi = np.zeros(3, dtype=float)
    for line in lines:
        if "xlo xhi" in line:
            parts = line.split()
            lo[0], hi[0] = float(parts[0]), float(parts[1])
        elif "ylo yhi" in line:
            parts = line.split()
            lo[1], hi[1] = float(parts[0]), float(parts[1])
        elif "zlo zhi" in line:
            parts = line.split()
            lo[2], hi[2] = float(parts[0]), float(parts[1])
    box = hi - lo
    if np.any(box <= 0.0):
        raise ValueError(f"Could not read orthorhombic box bounds from {path}")

    start = None
    for idx, line in enumerate(lines):
        if line.strip().startswith("Atoms"):
            start = idx + 2
            break
    if start is None:
        raise ValueError(f"No Atoms section found in {path}")

    ids: list[int] = []
    positions: list[list[float]] = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if len(parts) < 5:
            continue
        try:
            atom_id = int(parts[0])
            xyz = [float(parts[2]), float(parts[3]), float(parts[4])]
        except ValueError:
            break
        ids.append(atom_id)
        positions.append(xyz)
    wrapped = (np.array(positions, dtype=float) - lo) % box
    return np.array(ids, dtype=int), wrapped, box


def classify_by_neighbor_count(
    positions: np.ndarray, box: np.ndarray, nearest_neighbor: float, cutoff_factor: float
) -> tuple[np.ndarray, np.ndarray]:
    tree = cKDTree(positions, boxsize=box)
    neighbor_lists = tree.query_ball_point(positions, r=cutoff_factor * nearest_neighbor)
    counts = np.array([len(items) - 1 for items in neighbor_lists], dtype=int)
    is_bulk = counts == 12
    is_gb = ~is_bulk
    return is_gb, counts


def select_sites(
    ids: np.ndarray,
    positions: np.ndarray,
    is_gb: np.ndarray,
    neighbor_counts: np.ndarray,
    samples: int,
    bulk_references: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    gb_indices = np.flatnonzero(is_gb)
    bulk_indices = np.flatnonzero(~is_gb)
    if len(gb_indices) == 0:
        raise ValueError("No GB-like sites found. Try a larger neighbor cutoff.")
    if len(bulk_indices) == 0:
        raise ValueError("No bulk-like sites found. Try a smaller neighbor cutoff.")

    chosen_gb = rng.choice(gb_indices, size=min(samples, len(gb_indices)), replace=False)
    center = positions.mean(axis=0)
    bulk_dist2 = np.sum((positions[bulk_indices] - center) ** 2, axis=1)
    sorted_bulk = bulk_indices[np.argsort(bulk_dist2)]
    chosen_bulk = sorted_bulk[: min(bulk_references, len(sorted_bulk))]

    rows = []
    for idx in chosen_bulk:
        rows.append(["bulk_reference", ids[idx], *positions[idx], 0, neighbor_counts[idx]])
    for idx in chosen_gb:
        rows.append(["gb_sample", ids[idx], *positions[idx], 1, neighbor_counts[idx]])
    return np.array(rows, dtype=object)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/cuni_3d/cuni_3d_polycrystal.lammps"))
    parser.add_argument("--samples", type=int, default=500)
    parser.add_argument("--bulk-references", type=int, default=20)
    parser.add_argument("--lattice", type=float, default=3.615)
    parser.add_argument("--cutoff-factor", type=float, default=1.25, help="First-neighbor cutoff relative to a/sqrt(2).")
    parser.add_argument("--seed", type=int, default=500)
    parser.add_argument("--output", type=Path, default=Path("data/cuni_3d/cuni_3d_deltae_sites.csv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ids, positions, box = read_lammps_atoms(args.data)
    nearest_neighbor = args.lattice / np.sqrt(2.0)
    is_gb, neighbor_counts = classify_by_neighbor_count(positions, box, nearest_neighbor, args.cutoff_factor)
    rows = select_sites(ids, positions, is_gb, neighbor_counts, args.samples, args.bulk_references, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        args.output,
        rows,
        delimiter=",",
        header="role,atom_id,x_A,y_A,z_A,is_gb,first_neighbor_count",
        comments="",
        fmt="%s",
    )
    print(f"Atoms: {len(ids)}")
    print(f"GB-like sites: {int(is_gb.sum())} ({is_gb.mean():.2%})")
    print(f"Bulk-like sites: {int((~is_gb).sum())}")
    print(f"Wrote {args.samples} GB samples plus {min(args.bulk_references, int((~is_gb).sum()))} bulk references to {args.output}")


if __name__ == "__main__":
    main()
