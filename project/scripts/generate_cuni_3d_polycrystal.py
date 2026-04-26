#!/usr/bin/env python3
"""Generate a 3D Voronoi FCC Cu polycrystal as a LAMMPS data file.

This is a lightweight Atomsk-style fallback for the prototype Cu-Ni workflow:
random grain centers, random 3D orientations, Voronoi assignment with periodic
distances, and close-pair cleanup at grain boundaries.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation


FCC_BASIS = np.array(
    [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
    ],
    dtype=float,
)

ELEMENT_MASS = {
    "Cu": 63.546,
    "Ni": 58.693,
}

ELEMENT_LATTICE_A = {
    "Cu": 3.615,
    "Ni": 3.520,
}


@dataclass(frozen=True)
class Grain:
    center: np.ndarray
    rotation: Rotation


def periodic_delta(points: np.ndarray, center: np.ndarray, box: float) -> np.ndarray:
    delta = points - center
    return delta - box * np.round(delta / box)


def nearest_grain_ids(points: np.ndarray, centers: np.ndarray, box: float) -> np.ndarray:
    deltas = points[:, None, :] - centers[None, :, :]
    deltas -= box * np.round(deltas / box)
    dist2 = np.einsum("ijk,ijk->ij", deltas, deltas)
    return np.argmin(dist2, axis=1)


def grain_lattice_candidates(grain: Grain, box: float, lattice: float) -> np.ndarray:
    """Return one grain's rotated FCC lattice candidates covering the full box."""

    half_diagonal = 0.5 * np.sqrt(3.0) * box
    n = int(np.ceil(half_diagonal / lattice)) + 3
    grid = np.arange(-n, n + 1, dtype=float)
    ijk = np.stack(np.meshgrid(grid, grid, grid, indexing="ij"), axis=-1).reshape(-1, 3)
    lattice_points = (ijk[:, None, :] + FCC_BASIS[None, :, :]).reshape(-1, 3) * lattice
    rotated = grain.rotation.apply(lattice_points)
    return (rotated + grain.center) % box


def remove_close_pairs(
    positions: np.ndarray, grain_ids: np.ndarray, box: float, min_distance: float
) -> tuple[np.ndarray, np.ndarray]:
    """Remove one atom from each too-close pair, preferring balanced grain retention."""

    tree = cKDTree(positions, boxsize=box)
    pairs = tree.query_pairs(r=min_distance, output_type="ndarray")
    if len(pairs) == 0:
        return positions, grain_ids

    grain_counts = np.bincount(grain_ids)
    remove: set[int] = set()
    for i, j in pairs:
        if i in remove or j in remove:
            continue
        gi = grain_ids[i]
        gj = grain_ids[j]
        if grain_counts[gi] > grain_counts[gj]:
            victim = i
        elif grain_counts[gj] > grain_counts[gi]:
            victim = j
        else:
            victim = max(i, j)
        remove.add(int(victim))
        grain_counts[grain_ids[victim]] -= 1

    keep = np.ones(len(positions), dtype=bool)
    keep[list(remove)] = False
    return positions[keep], grain_ids[keep]


def generate_polycrystal(
    box_nm: float,
    grains: int,
    lattice: float,
    seed: int,
    min_distance_factor: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    box = box_nm * 10.0
    rng = np.random.default_rng(seed)
    centers = rng.random((grains, 3)) * box
    rotations = Rotation.random(grains, random_state=rng)
    grain_defs = [Grain(centers[i], rotations[i]) for i in range(grains)]

    all_positions: list[np.ndarray] = []
    all_grain_ids: list[np.ndarray] = []
    for grain_id, grain in enumerate(grain_defs):
        candidates = grain_lattice_candidates(grain, box, lattice)
        owners = nearest_grain_ids(candidates, centers, box)
        owned = candidates[owners == grain_id]
        all_positions.append(owned)
        all_grain_ids.append(np.full(len(owned), grain_id, dtype=int))

    positions = np.concatenate(all_positions, axis=0)
    grain_ids = np.concatenate(all_grain_ids, axis=0)
    positions, grain_ids = remove_close_pairs(
        positions, grain_ids, box, min_distance_factor * lattice / np.sqrt(2.0)
    )
    order = np.lexsort((positions[:, 2], positions[:, 1], positions[:, 0]))
    return positions[order], grain_ids[order], centers


def write_lammps_data(path: Path, positions: np.ndarray, box: float, solvent: str, solute: str) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"3D Voronoi FCC {solvent} polycrystal for {solvent}({solute})\n\n")
        handle.write(f"{len(positions)} atoms\n\n")
        handle.write("2 atom types\n\n")
        handle.write(f"0.0 {box:.10f} xlo xhi\n")
        handle.write(f"0.0 {box:.10f} ylo yhi\n")
        handle.write(f"0.0 {box:.10f} zlo zhi\n\n")
        handle.write("Masses\n\n")
        handle.write(f"1 {ELEMENT_MASS[solvent]}\n")
        handle.write(f"2 {ELEMENT_MASS[solute]}\n\n")
        handle.write("Atoms # atomic\n\n")
        for atom_id, (x, y, z) in enumerate(positions, start=1):
            handle.write(f"{atom_id} 1 {x:.10f} {y:.10f} {z:.10f}\n")


def write_metadata(path: Path, positions: np.ndarray, grain_ids: np.ndarray, centers: np.ndarray) -> None:
    header = "atom_id,x_A,y_A,z_A,grain_id"
    rows = np.column_stack((np.arange(1, len(positions) + 1), positions, grain_ids))
    np.savetxt(path, rows, delimiter=",", header=header, comments="", fmt=["%d", "%.10f", "%.10f", "%.10f", "%d"])

    centers_path = path.with_name(path.stem + "_grain_centers.csv")
    np.savetxt(
        centers_path,
        np.column_stack((np.arange(len(centers)), centers)),
        delimiter=",",
        header="grain_id,x_A,y_A,z_A",
        comments="",
        fmt=["%d", "%.10f", "%.10f", "%.10f"],
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--box-nm", type=float, default=10.0, help="Cubic box length in nm.")
    parser.add_argument("--grains", type=int, default=8, help="Number of Voronoi grains.")
    parser.add_argument("--solvent", choices=sorted(ELEMENT_MASS), default="Cu", help="Matrix/solvent element, atom type 1.")
    parser.add_argument("--solute", choices=sorted(ELEMENT_MASS), default="Ni", help="Substitutional solute element, atom type 2.")
    parser.add_argument("--lattice", type=float, default=None, help="FCC lattice parameter in Angstrom. Defaults to the solvent value.")
    parser.add_argument("--seed", type=int, default=20260425)
    parser.add_argument("--min-distance-factor", type=float, default=0.72, help="Close-pair cutoff relative to FCC nearest-neighbor distance.")
    parser.add_argument("--output", type=Path, default=Path("data/cuni_3d/cuni_3d_polycrystal.lammps"))
    parser.add_argument("--metadata", type=Path, default=Path("data/cuni_3d/cuni_3d_polycrystal_sites.csv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.solvent == args.solute:
        raise ValueError("--solvent and --solute must be different")
    if args.lattice is None:
        args.lattice = ELEMENT_LATTICE_A[args.solvent]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.parent.mkdir(parents=True, exist_ok=True)
    positions, grain_ids, centers = generate_polycrystal(
        box_nm=args.box_nm,
        grains=args.grains,
        lattice=args.lattice,
        seed=args.seed,
        min_distance_factor=args.min_distance_factor,
    )
    box = args.box_nm * 10.0
    write_lammps_data(args.output, positions, box, args.solvent, args.solute)
    write_metadata(args.metadata, positions, grain_ids, centers)
    print(f"Wrote {len(positions)} atoms to {args.output}")
    print(f"Wrote site metadata to {args.metadata}")


if __name__ == "__main__":
    main()
