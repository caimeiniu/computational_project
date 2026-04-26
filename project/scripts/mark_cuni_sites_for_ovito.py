#!/usr/bin/env python3
"""Mark sampled Cu-Ni Delta E sites as separate atom types for OVITO."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_site_types(path: Path) -> dict[int, int]:
    site_types = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            atom_id = int(row["atom_id"])
            site_types[atom_id] = 3 if row["role"] == "bulk_reference" else 2
    return site_types


def rewrite_lammps_types(input_path: Path, output_path: Path, site_types: dict[int, int]) -> None:
    lines = input_path.read_text(encoding="utf-8").splitlines()
    output = []
    in_atoms = False
    atom_section_started = False

    for line in lines:
        stripped = line.strip()
        if stripped.endswith("atom types"):
            output.append("3 atom types")
            continue
        if stripped.startswith("Atoms"):
            in_atoms = True
            atom_section_started = False
            output.append(line)
            continue
        if in_atoms:
            if not stripped:
                output.append(line)
                if atom_section_started:
                    in_atoms = False
                continue
            atom_section_started = True
            parts = line.split()
            if len(parts) >= 5:
                atom_id = int(parts[0])
                parts[1] = str(site_types.get(atom_id, 1))
                output.append(" ".join(parts))
                continue
        output.append(line)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(output) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/cuni_3d/lammps/cuni_3d_relaxed.lammps"))
    parser.add_argument("--sites", type=Path, default=Path("data/cuni_3d/cuni_3d_deltae_sites.csv"))
    parser.add_argument("--output", type=Path, default=Path("data/cuni_3d/cuni_3d_sites_for_ovito.lammps"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    site_types = read_site_types(args.sites)
    rewrite_lammps_types(args.data, args.output, site_types)
    print(f"Wrote OVITO-marked structure to {args.output}")
    print("Type 1 = Cu matrix, type 2 = sampled GB sites, type 3 = bulk reference")


if __name__ == "__main__":
    main()
