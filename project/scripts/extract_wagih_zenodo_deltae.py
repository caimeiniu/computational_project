#!/usr/bin/env python3
"""Extract raw Wagih Zenodo segregation energies from a LAMMPS dump file.

Use --list-candidates on the extracted Zenodo directory to find likely Ni(Cu)
files, then pass the chosen dump file to this script. The parser expects a
LAMMPS text dump with an "ITEM: ATOMS ..." header.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


ENERGY_COLUMN_CANDIDATES = (
    "delta_e_kj_mol",
    "deltaE",
    "DeltaE",
    "dE",
    "dEseg",
    "Eseg",
    "seg_energy",
    "segregation_energy",
    "c_deltaE",
    "c_eng",
)


def list_candidates(root: Path, solvent: str, solute: str) -> None:
    patterns = [
        f"*{solvent}*{solute}*",
        f"*{solvent.lower()}*{solute.lower()}*",
        f"*{solvent}_{solute}*",
        f"*{solute}_in_{solvent}*",
    ]
    seen: set[Path] = set()
    for pattern in patterns:
        for path in root.rglob(pattern):
            if path.is_file() and path not in seen:
                seen.add(path)
                print(path)


def find_atoms_header(path: Path) -> tuple[list[str], int]:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line.startswith("ITEM: ATOMS"):
                return line.split()[2:], line_number + 1
    raise ValueError(f"No 'ITEM: ATOMS' header found in {path}")


def choose_energy_column(columns: list[str], requested: str | None) -> str:
    if requested:
        if requested not in columns:
            raise ValueError(f"Requested column {requested!r} not in dump columns: {columns}")
        return requested
    for candidate in ENERGY_COLUMN_CANDIDATES:
        if candidate in columns:
            return candidate
    raise ValueError(
        "Could not auto-detect energy column. "
        f"Columns are: {columns}. Re-run with --energy-column COLUMN."
    )


def read_dump_values(path: Path, energy_column: str | None, drop_zero: bool) -> tuple[list[dict[str, float]], str, list[str]]:
    columns, start_line = find_atoms_header(path)
    selected_column = choose_energy_column(columns, energy_column)
    energy_idx = columns.index(selected_column)
    id_idx = columns.index("id") if "id" in columns else None
    rows: list[dict[str, float]] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for _ in range(start_line):
            next(handle)
        for line in handle:
            if line.startswith("ITEM:"):
                break
            parts = line.split()
            if len(parts) <= energy_idx:
                continue
            atom_id = int(float(parts[id_idx])) if id_idx is not None else len(rows) + 1
            delta_e = float(parts[energy_idx])
            if drop_zero and delta_e == 0.0:
                continue
            rows.append({"atom_id": atom_id, "delta_e_kj_mol": delta_e})
    return rows, selected_column, columns


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["atom_id", "delta_e_kj_mol"])
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--list-candidates", type=Path, help="Extracted Zenodo directory to search.")
    parser.add_argument("--solvent", default="Ni")
    parser.add_argument("--solute", default="Cu")
    parser.add_argument("--dump", type=Path, help="Chosen Wagih LAMMPS dump file.")
    parser.add_argument("--energy-column", help="Energy column name from ITEM: ATOMS header.")
    parser.add_argument("--drop-zero", action="store_true", help="Drop exact zero values, often non-GB placeholders in Wagih dumps.")
    parser.add_argument("--output", type=Path, default=Path("wagih_nicu_raw_deltae.csv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.list_candidates:
        list_candidates(args.list_candidates, args.solvent, args.solute)
        return
    if not args.dump:
        raise SystemExit("Provide --dump FILE, or use --list-candidates DIR first.")
    rows, selected_column, columns = read_dump_values(args.dump, args.energy_column, args.drop_zero)
    write_csv(args.output, rows)
    print(f"Dump columns: {' '.join(columns)}")
    print(f"Used energy column: {selected_column}")
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
