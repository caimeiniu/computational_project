#!/usr/bin/env python3
"""Build a pre-segregated Al(Mg) data file: place all N_Mg = X_c · N_total
Mg atoms at randomly-chosen GB sites (no bulk Mg).

Used as the second initial condition for the X_c=5e-2 equilibration check
(canonical-FD predicts X_GB^∞ ≈ 0.228; this IC starts at X_GB(0) ≈ ceiling
0.267 → HMC must pull X_GB *down* toward 0.228 from above. If random-IC
HMC also climbs from below to 0.228, the two curves bracket the
equilibrium and equilibration is verified.)

Inputs:
  --in-data   annealed pure-Al lmp data (atom_style atomic, type 1 only)
  --gb-mask   bool npy of length N_total; True = GB site (atom id-1 indexed)
  --xc        total Mg fraction (used to set N_Mg = round(X_c · N_total))
  --seed      RNG seed for which GB sites get Mg
  --out-data  output lmp data file with type 2 (Mg) at chosen GB IDs

The output preserves all coordinates, box, and atom ordering of the input;
only the type column is rewritten for the chosen Mg IDs.
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path
import numpy as np


def rewrite_data_types(in_path: Path, out_path: Path,
                       mg_ids: np.ndarray, comment: str) -> dict:
    """Stream the LAMMPS data file, rewriting the type column for atom IDs
    in `mg_ids`. Returns a dict of audit counts.
    """
    mg_set = set(int(x) for x in mg_ids.tolist())
    section_re = re.compile(r"^\s*(Atoms|Velocities|Bonds|Angles|Masses|"
                             r"Pair Coeffs|Dihedrals|Impropers)\b")
    # State machine:
    #   in_atoms = False            outside Atoms section
    #   in_atoms = True, seen=False just after "Atoms # ..." header,
    #                                blank line(s) before the data block
    #   in_atoms = True, seen=True  inside the atom-line block; first blank
    #                                terminates the section
    in_atoms = False
    seen_atom = False
    n_changed = 0
    n_atoms_seen = 0

    with in_path.open() as fin, out_path.open("w") as fout:
        # rewrite top comment line for traceability
        first = fin.readline()
        fout.write(first.rstrip() + f"  [pre_segregate: {comment}]\n")

        for line in fin:
            stripped = line.strip()
            if section_re.match(stripped):
                in_atoms = stripped.startswith("Atoms")
                seen_atom = False
                fout.write(line)
                continue

            if in_atoms and stripped:
                parts = line.split()
                # atom_style atomic: id type x y z [ix iy iz]
                if len(parts) >= 5 and parts[0].lstrip("-").isdigit():
                    aid = int(parts[0])
                    atype = int(parts[1])
                    n_atoms_seen += 1
                    seen_atom = True
                    if aid in mg_set and atype != 2:
                        parts[1] = "2"
                        n_changed += 1
                        line = " ".join(parts) + "\n"
                fout.write(line)
            else:
                fout.write(line)
                # blank line AFTER atoms terminates the section; blank
                # lines BEFORE the first atom (header padding) do not.
                if in_atoms and not stripped and seen_atom:
                    in_atoms = False
                    seen_atom = False

    return {"n_atoms_seen": n_atoms_seen, "n_changed_to_Mg": n_changed,
            "n_requested": len(mg_set)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--in-data", type=Path, required=True)
    p.add_argument("--gb-mask", type=Path, required=True)
    p.add_argument("--xc", type=float, required=True,
                   help="total Mg fraction; N_Mg = round(X_c · N_total)")
    p.add_argument("--seed", type=int, default=20260426)
    p.add_argument("--out-data", type=Path, required=True)
    args = p.parse_args()

    gb_mask = np.load(args.gb_mask).astype(bool)
    n_total = int(gb_mask.size)
    n_gb = int(gb_mask.sum())
    n_mg = int(round(args.xc * n_total))
    if n_mg <= 0:
        sys.exit(f"X_c={args.xc} too small: N_Mg={n_mg}")

    # GB / bulk atom IDs are 1-indexed (LAMMPS convention).
    gb_ids = np.where(gb_mask)[0] + 1
    bulk_ids = np.where(~gb_mask)[0] + 1
    rng = np.random.default_rng(args.seed)

    if n_mg <= n_gb:
        # Standard preseg: fill `n_mg` random GB sites; X_GB(0) = n_mg/n_gb,
        # X_bulk(0) = 0. Bracket equilibrium from above (Mg over-saturated at GB).
        mg_ids = rng.choice(gb_ids, size=n_mg, replace=False)
        n_mg_gb, n_mg_bulk = n_mg, 0
        x_gb_init = n_mg / n_gb
        x_bulk_init = 0.0
    else:
        # Mixed preseg: GB ceiling is below total-Mg requirement. Fill ALL
        # GB sites with Mg; remaining Mg goes to random bulk sites. Starts
        # at X_GB(0)=1, X_bulk(0)=overflow/N_bulk — still above canon-FD,
        # so the descent direction is still GB→bulk.
        n_mg_gb = n_gb
        n_mg_bulk = n_mg - n_gb
        if n_mg_bulk > bulk_ids.size:
            sys.exit(f"X_c={args.xc} requires {n_mg} Mg but box has only "
                     f"{n_total} atoms.")
        mg_ids = np.concatenate([
            gb_ids,
            rng.choice(bulk_ids, size=n_mg_bulk, replace=False),
        ])
        x_gb_init = 1.0
        x_bulk_init = n_mg_bulk / bulk_ids.size

    mg_ids.sort()  # for reproducibility / human-readability of the diff

    comment = (f"X_c={args.xc:g}  N_Mg={n_mg}  "
               f"GB={n_mg_gb} bulk={n_mg_bulk}  "
               f"X_GB(0)={x_gb_init:.4f} X_bulk(0)={x_bulk_init:.4f}  "
               f"seed={args.seed}")
    audit = rewrite_data_types(args.in_data, args.out_data, mg_ids, comment)

    # sanity
    assert audit["n_atoms_seen"] == n_total, (
        f"atoms seen {audit['n_atoms_seen']} != mask length {n_total}")
    assert audit["n_changed_to_Mg"] == n_mg, (
        f"changed {audit['n_changed_to_Mg']} != requested {n_mg}")

    print(f"wrote {args.out_data}")
    print(f"  N_total = {n_total}  N_GB = {n_gb}  GB_frac = {n_gb/n_total:.4f}")
    print(f"  X_c = {args.xc}  =>  N_Mg = {n_mg}")
    print(f"  N_Mg(GB) = {n_mg_gb}  N_Mg(bulk) = {n_mg_bulk}")
    print(f"  X_GB(0) = {x_gb_init:.4f}   X_bulk(0) = {x_bulk_init:.4f}")
    print(f"  seed = {args.seed}  (deterministic site choice)")


if __name__ == "__main__":
    main()
