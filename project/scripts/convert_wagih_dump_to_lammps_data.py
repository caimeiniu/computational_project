#!/usr/bin/env python3
"""Convert a Wagih/Zenodo LAMMPS dump structure to a LAMMPS data file."""

from __future__ import annotations

import argparse
from pathlib import Path


ELEMENT_MASS = {
    "Au": 196.96657,
    "Cu": 63.546,
    "Ni": 58.693,
    "Pt": 195.084,
}


def read_dump(path: Path) -> tuple[list[tuple[int, float, float, float]], list[tuple[float, float]], list[str]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    bounds: list[tuple[float, float]] = []
    columns: list[str] | None = None
    atoms_start: int | None = None

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("ITEM: BOX BOUNDS"):
            bounds = []
            for j in range(1, 4):
                parts = lines[i + j].split()
                bounds.append((float(parts[0]), float(parts[1])))
            i += 4
            continue
        if line.startswith("ITEM: ATOMS"):
            columns = line.split()[2:]
            atoms_start = i + 1
            break
        i += 1

    if columns is None or atoms_start is None:
        raise ValueError(f"No ITEM: ATOMS section found in {path}")
    if len(bounds) != 3:
        raise ValueError(f"No 3D BOX BOUNDS section found in {path}")

    id_idx = columns.index("id") if "id" in columns else None
    coord_names = []
    for axis in "xyz":
        if axis in columns:
            coord_names.append(axis)
        elif f"{axis}s" in columns:
            coord_names.append(f"{axis}s")
        else:
            raise ValueError(f"Could not find {axis} or {axis}s coordinate column in {columns}")
    coord_idx = [columns.index(name) for name in coord_names]

    atoms: list[tuple[int, float, float, float]] = []
    for offset, line in enumerate(lines[atoms_start:], start=1):
        if line.startswith("ITEM:"):
            break
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        atom_id = int(float(parts[id_idx])) if id_idx is not None else offset
        xyz = []
        for value, name, (lo, hi) in zip((float(parts[idx]) for idx in coord_idx), coord_names, bounds):
            if name.endswith("s"):
                xyz.append(lo + value * (hi - lo))
            else:
                xyz.append(value)
        atoms.append((atom_id, xyz[0], xyz[1], xyz[2]))

    atoms.sort(key=lambda row: row[0])
    return atoms, bounds, columns


def write_data(path: Path, atoms: list[tuple[int, float, float, float]], bounds: list[tuple[float, float]], solvent: str, solute: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"Wagih dump structure converted for {solvent}({solute}) substitution jobs\n\n")
        handle.write(f"{len(atoms)} atoms\n\n")
        handle.write("2 atom types\n\n")
        for axis, (lo, hi) in zip("xyz", bounds):
            handle.write(f"{lo:.10f} {hi:.10f} {axis}lo {axis}hi\n")
        handle.write("\nMasses\n\n")
        handle.write(f"1 {ELEMENT_MASS[solvent]}\n")
        handle.write(f"2 {ELEMENT_MASS[solute]}\n\n")
        handle.write("Atoms # atomic\n\n")
        for atom_id, x, y, z in atoms:
            handle.write(f"{atom_id} 1 {x:.10f} {y:.10f} {z:.10f}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--solvent", choices=sorted(ELEMENT_MASS), default="Au")
    parser.add_argument("--solute", choices=sorted(ELEMENT_MASS), default="Pt")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.solvent == args.solute:
        raise ValueError("--solvent and --solute must be different")
    atoms, bounds, columns = read_dump(args.dump)
    write_data(args.output, atoms, bounds, args.solvent, args.solute)
    print(f"Dump columns: {' '.join(columns)}")
    print(f"Wrote {len(atoms)} atoms to {args.output}")


if __name__ == "__main__":
    main()
