#!/usr/bin/env python3
"""Rewrite Pt(Au) LAMMPS data types for OVITO GB/solute visualization.

Output type convention:
  1 = Pt bulk
  2 = Au bulk
  3 = Pt GB
  4 = Au GB
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def _is_atom_line(line: str) -> bool:
    parts = line.split()
    if len(parts) < 5:
        return False
    try:
        int(parts[0])
        int(parts[1])
    except ValueError:
        return False
    return True


def rewrite_data(input_path: Path, output_path: Path, gb_mask: np.ndarray) -> dict[str, float]:
    lines = input_path.read_text().splitlines(keepends=True)
    out_lines: list[str] = []
    in_atoms = False
    saw_atoms = False
    counts = {1: 0, 2: 0, 3: 0, 4: 0}

    for line in lines:
        stripped = line.strip()
        parts = stripped.split()

        if len(parts) >= 3 and parts[1] == "atom" and parts[2] == "types":
            out_lines.append("4 atom types\n")
            continue

        if stripped == "Atoms" or stripped.startswith("Atoms "):
            in_atoms = True
            saw_atoms = True
            out_lines.append(line)
            continue

        if in_atoms and _is_atom_line(line):
            cols = line.split()
            atom_id = int(cols[0])
            original_type = int(cols[1])
            if atom_id < 1 or atom_id > len(gb_mask):
                raise ValueError(f"Atom id {atom_id} outside GB mask length {len(gb_mask)}")
            is_gb = bool(gb_mask[atom_id - 1])
            if original_type == 1 and not is_gb:
                new_type = 1
            elif original_type == 2 and not is_gb:
                new_type = 2
            elif original_type == 1 and is_gb:
                new_type = 3
            elif original_type == 2 and is_gb:
                new_type = 4
            else:
                raise ValueError(f"Unexpected original atom type {original_type} at id {atom_id}")
            cols[1] = str(new_type)
            counts[new_type] += 1
            out_lines.append(" ".join(cols) + "\n")
            continue

        out_lines.append(line)

    if not saw_atoms:
        raise ValueError(f"Could not find Atoms section in {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(out_lines))

    n_total = sum(counts.values())
    n_gb = counts[3] + counts[4]
    n_au = counts[2] + counts[4]
    return {
        "n_total": n_total,
        "n_gb": n_gb,
        "n_au": n_au,
        "x_total": n_au / n_total if n_total else float("nan"),
        "x_gb": counts[4] / n_gb if n_gb else float("nan"),
        "type1_pt_bulk": counts[1],
        "type2_au_bulk": counts[2],
        "type3_pt_gb": counts[3],
        "type4_au_gb": counts[4],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input LAMMPS data file")
    parser.add_argument("--gb-mask", required=True, help="Boolean .npy mask, id-1 indexed")
    parser.add_argument("--output", required=True, help="Output OVITO-friendly LAMMPS data file")
    args = parser.parse_args()

    gb_mask = np.asarray(np.load(args.gb_mask), dtype=bool)
    stats = rewrite_data(Path(args.input), Path(args.output), gb_mask)
    print(f"Wrote {args.output}")
    print(
        "type map: 1=Pt bulk, 2=Au bulk, 3=Pt GB, 4=Au GB; "
        f"X_total={stats['x_total']:.6f}, X_GB={stats['x_gb']:.6f}"
    )


if __name__ == "__main__":
    main()
