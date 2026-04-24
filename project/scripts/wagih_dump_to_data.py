"""Convert Wagih's Zenodo annealed dump0 to a LAMMPS data file aligned with
our sample_delta_e.py expectations (2 atom types, atoms id-sorted 1..N).

Also builds a gb_mask.npy of length N_atoms where mask[id-1]=True iff
atom_id is in Wagih's seg_energies_Al_Mg.txt site list.

Usage:
    python wagih_dump_to_data.py \\
        --dump  .../heated_minimized_Al_polycrystal.dump0 \\
        --seg   .../seg_energies_Al_Mg.txt \\
        --out-data  out/wagih_Al_200A.lmp \\
        --out-mask  out/wagih_gb_mask.npy
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def read_dump(path: Path):
    with open(path) as f:
        lines = f.readlines()
    # Header parse
    i = 0
    n_atoms = 0
    box = np.zeros((3, 2))
    data_start = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            n_atoms = int(lines[i + 1])
            i += 2
        elif line.startswith("ITEM: BOX BOUNDS"):
            for k in range(3):
                lo, hi = map(float, lines[i + 1 + k].split()[:2])
                box[k] = (lo, hi)
            i += 4
        elif line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]   # names after "ITEM: ATOMS"
            data_start = i + 1
            break
        else:
            i += 1
    data = np.loadtxt(lines[data_start:data_start + n_atoms])
    # Map columns
    col = {name: idx for idx, name in enumerate(cols)}
    x = data[:, col["x"]]
    y = data[:, col["y"]]
    z = data[:, col["z"]]
    ids = data[:, col["id"]].astype(np.int64)
    types = data[:, col["type"]].astype(np.int64)
    # sort by id so that row i has atom_id = i+1 (assuming 1..N)
    order = np.argsort(ids)
    return {
        "n_atoms": n_atoms, "box": box,
        "x": x[order], "y": y[order], "z": z[order],
        "id": ids[order], "type": types[order],
    }


def write_data(d: dict, path: Path, n_atom_types: int = 2):
    with open(path, "w") as f:
        f.write("LAMMPS data file from Wagih Zenodo dump0\n\n")
        f.write(f"{d['n_atoms']} atoms\n")
        f.write(f"{n_atom_types} atom types\n\n")
        for k, ax in enumerate("xyz"):
            f.write(f"{d['box'][k, 0]:.9f} {d['box'][k, 1]:.9f} {ax}lo {ax}hi\n")
        f.write("\nAtoms\n\n")
        for i in range(d["n_atoms"]):
            f.write(f"{int(d['id'][i])} {int(d['type'][i])} "
                    f"{d['x'][i]:.9f} {d['y'][i]:.9f} {d['z'][i]:.9f}\n")


def build_mask(seg_txt: Path, n_atoms: int, max_id: int) -> np.ndarray:
    mask = np.zeros(max_id, dtype=bool)
    with open(seg_txt) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            sid = int(parts[0])
            mask[sid - 1] = True   # LAMMPS 1-based → np 0-based
    return mask


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dump", required=True)
    p.add_argument("--seg", required=True)
    p.add_argument("--out-data", required=True)
    p.add_argument("--out-mask", required=True)
    args = p.parse_args()

    d = read_dump(Path(args.dump))
    print(f"read dump: N={d['n_atoms']}, "
          f"id range [{d['id'].min()}, {d['id'].max()}], "
          f"types unique {np.unique(d['type'])}")
    # expect ids 1..N contiguous after sort
    if not np.array_equal(d["id"], np.arange(1, d["n_atoms"] + 1)):
        missing = np.setdiff1d(np.arange(1, d["n_atoms"] + 1), d["id"])
        print(f"WARNING: ids not contiguous; {len(missing)} gaps")
    write_data(d, Path(args.out_data), n_atom_types=2)
    print(f"wrote data file {args.out_data}")

    max_id = int(d["id"].max())
    mask = build_mask(Path(args.seg), n_atoms=d["n_atoms"], max_id=max_id)
    print(f"GB mask: {mask.sum()} True / {mask.size} total "
          f"(f_gb = {mask.mean():.4f})")
    Path(args.out_mask).parent.mkdir(parents=True, exist_ok=True)
    np.save(args.out_mask, mask)
    print(f"wrote mask {args.out_mask}")


if __name__ == "__main__":
    main()
