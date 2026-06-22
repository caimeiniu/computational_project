#!/usr/bin/env python3
"""Create a Pt(Au) LAMMPS data file with Au preferentially seeded on GB atoms."""

from __future__ import annotations

import argparse
import json
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


def _count_atoms(lines: list[str]) -> int:
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[1] == "atoms":
            return int(parts[0])
    raise ValueError("Could not find '<N> atoms' header line")


def _seed_types(
    *,
    n_atoms: int,
    gb_mask: np.ndarray,
    xc: float,
    seed: int,
    solvent_type: int,
    solute_type: int,
) -> tuple[np.ndarray, dict]:
    if not (0.0 <= xc <= 1.0):
        raise ValueError(f"--xc must be in [0,1], got {xc}")
    if gb_mask.shape != (n_atoms,):
        raise ValueError(f"GB mask length {gb_mask.size} does not match atom count {n_atoms}")

    rng = np.random.default_rng(seed)
    n_solute = int(round(xc * n_atoms))
    gb_ids = np.where(gb_mask)[0]
    bulk_ids = np.where(~gb_mask)[0]

    types = np.full(n_atoms, solvent_type, dtype=np.int16)
    if n_solute <= len(gb_ids):
        chosen_gb = rng.choice(gb_ids, size=n_solute, replace=False)
        chosen_bulk = np.array([], dtype=np.int64)
    else:
        chosen_gb = gb_ids
        n_bulk_needed = n_solute - len(gb_ids)
        if n_bulk_needed > len(bulk_ids):
            raise ValueError("Requested more solute atoms than total atoms")
        chosen_bulk = rng.choice(bulk_ids, size=n_bulk_needed, replace=False)

    types[chosen_gb] = solute_type
    types[chosen_bulk] = solute_type
    meta = {
        "xc_target": xc,
        "seed": seed,
        "n_atoms": n_atoms,
        "n_solute": n_solute,
        "x_total_actual": n_solute / n_atoms,
        "n_gb": int(gb_mask.sum()),
        "n_gb_solute": int((types[gb_mask] == solute_type).sum()),
        "x_gb_initial": float((types[gb_mask] == solute_type).mean()),
        "n_bulk": int((~gb_mask).sum()),
        "n_bulk_solute": int((types[~gb_mask] == solute_type).sum()),
        "x_bulk_initial": float((types[~gb_mask] == solute_type).mean()),
    }
    return types, meta


def rewrite_data_types(in_path: Path, out_path: Path, types: np.ndarray) -> None:
    lines = in_path.read_text().splitlines(keepends=True)
    out_lines: list[str] = []
    in_atoms = False
    saw_atoms = False

    for line in lines:
        stripped = line.strip()
        if stripped == "Atoms" or stripped.startswith("Atoms "):
            in_atoms = True
            saw_atoms = True
            out_lines.append(line)
            continue

        if in_atoms and _is_atom_line(line):
            parts = line.split()
            atom_id = int(parts[0])
            if atom_id < 1 or atom_id > len(types):
                raise ValueError(f"Atom id {atom_id} outside 1..{len(types)}")
            parts[1] = str(int(types[atom_id - 1]))
            out_lines.append(" ".join(parts) + "\n")
            continue

        out_lines.append(line)

    if not saw_atoms:
        raise ValueError("Could not find Atoms section")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(out_lines))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input all-Pt annealed LAMMPS data file")
    parser.add_argument("--gb-mask", required=True, help="Boolean .npy GB mask, id-1 indexed")
    parser.add_argument("--xc", type=float, required=True, help="Target total Au fraction")
    parser.add_argument("--seed", type=int, default=20260531)
    parser.add_argument("--out", required=True, help="Output seeded LAMMPS data file")
    parser.add_argument("--meta-json", help="Optional JSON metadata output")
    parser.add_argument("--solvent-type", type=int, default=1)
    parser.add_argument("--solute-type", type=int, default=2)
    args = parser.parse_args()

    in_path = Path(args.input)
    lines = in_path.read_text().splitlines()
    n_atoms = _count_atoms(lines)
    gb_mask = np.asarray(np.load(args.gb_mask), dtype=bool)
    types, meta = _seed_types(
        n_atoms=n_atoms,
        gb_mask=gb_mask,
        xc=args.xc,
        seed=args.seed,
        solvent_type=args.solvent_type,
        solute_type=args.solute_type,
    )
    rewrite_data_types(in_path, Path(args.out), types)
    if args.meta_json:
        Path(args.meta_json).write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print(
        "GB_SEED_OK "
        f"out={args.out} X_total={meta['x_total_actual']:.6f} "
        f"X_GB_initial={meta['x_gb_initial']:.6f} "
        f"X_bulk_initial={meta['x_bulk_initial']:.6f}"
    )


if __name__ == "__main__":
    main()
