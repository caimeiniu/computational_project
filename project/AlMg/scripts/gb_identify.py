"""Grain-boundary identification via LAMMPS CNA (common neighbor analysis).

Produces a per-atom boolean mask — `True` for a grain-boundary atom, `False`
for a bulk-crystal atom — following Wagih, Larsen & Schuh (Nat. Commun. 2020):
    "all atoms NOT identified as bulk crystal (FCC/BCC/HCP) are GB atoms".

Internally runs `lmp` once with `compute cna/atom`, which is LAMMPS'
**fixed-cutoff** CNA — not the adaptive CNA (a-CNA) of OVITO/Stukowski 2012.
For the bulk/GB binary classification this is equivalent at the accuracy we
need (relaxed polycrystal, small thermal expansion). If the project later
needs fine-grained GB-character analysis, switch to OVITO's
`CommonNeighborAnalysisModifier` or LAMMPS `compute ptm/atom` (Polyhedral
Template Matching).

No OVITO dependency. The LAMMPS binary must be on `PATH` (on Euler, `module
load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4` before calling
this script).

**Semantic note (FCC parent)**: HCP-labelled atoms inside an FCC polycrystal
are **stacking faults** (planar defects), not grain boundaries. Under the
strict "bulk = parent structure only" definition used here, stacking faults
get flagged as GB (mask=True). This is consistent with Wagih's treatment
(Methods §"GB site identification") but worth knowing if you see a non-trivial
HCP count in `info["cna_counts"]`.

Usage as a module (downstream Phase 3 / HMC analysis):
    from gb_identify import compute_gb_mask
    mask, info = compute_gb_mask(
        "poly_Al_100A_8g_annealed.lmp",
        parent_structure="fcc",
        lattice_a=4.05,
    )
    gb_ids = np.where(mask)[0] + 1   # LAMMPS atom ids are 1-based
    print(info["f_gb"], info["cna_counts"])

Usage from the shell (quick inspection, writes sidecar files):
    python gb_identify.py poly_Al_100A_8g_annealed.lmp \\
        --lattice-a 4.05 \\
        --out-mask gb_mask.npy \\
        --out-report gb_info.json \\
        --out-dump gb_cna.dump          # CNA-labelled dump for OVITO

The mask is ordered by LAMMPS atom id (dump is sort-id'd), so
`mask[i]` corresponds to the atom with id `i+1`. This assumes ids are
contiguous 1..N (true for the output of `generate_polycrystal.py` and of
any LAMMPS run that didn't delete atoms).

CNA cutoff (LAMMPS recommendation, set automatically from `lattice_a`):
    FCC : 0.854 * a
    BCC : 1.207 * a
    HCP : 1.207 * a
Override with `cna_cutoff=<angstroms>` if your structure is under noticeable
thermal expansion and the default shells no longer split cleanly.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np


# LAMMPS `compute cna/atom` output codes.
_CNA_LABELS = {0: "unknown", 1: "fcc", 2: "hcp", 3: "bcc", 4: "ico", 5: "other"}
_BULK_CNA_INT = {"fcc": 1, "bcc": 3, "hcp": 2}

# Midpoint-between-neighbor-shells cutoffs (LAMMPS docs).
_CUTOFF_FORMULAS = {
    "fcc": lambda a: 0.854 * a,
    "bcc": lambda a: 1.207 * a,
    "hcp": lambda a: 1.207 * a,
}


def compute_gb_mask(
    data_file: str | Path,
    *,
    parent_structure: str = "fcc",
    lattice_a: float = 4.05,
    cna_cutoff: float | None = None,
    lmp_binary: str = "lmp",
    out_dump_path: str | Path | None = None,
) -> tuple[np.ndarray, dict]:
    """Classify every atom in a LAMMPS data file as GB or bulk via CNA.

    Parameters
    ----------
    data_file : str or Path
        LAMMPS data file, typically the post-anneal relaxed polycrystal
        (e.g. `poly_Al_100A_8g_annealed.lmp`).
    parent_structure : {"fcc", "bcc", "hcp"}
        Which CNA class counts as bulk. All other classes — including
        "other" and "unknown" — are flagged as GB.
    lattice_a : float
        Lattice parameter in Å, used to set the default CNA cutoff
        (ignored if `cna_cutoff` is provided).
    cna_cutoff : float, optional
        Override the CNA neighbor cutoff (Å).
    lmp_binary : str
        Name/path of the LAMMPS executable. Default `"lmp"` assumes the
        Euler module is already loaded.
    out_dump_path : str or Path, optional
        If provided, copy the CNA LAMMPS dump (columns `id type x y z c_cna`)
        to this location for OVITO inspection. Includes coordinates so OVITO
        can render atoms coloured by CNA class without needing the source
        data file.

    Returns
    -------
    mask : (N,) np.ndarray[bool]
        `mask[i] == True` iff the atom with LAMMPS id `i+1` is a GB atom.
        Length equals the number of atoms in `data_file`.
    info : dict
        Human-readable summary with `n_atoms`, `n_gb`, `f_gb`,
        `cna_counts` (per-label histogram), `parent_structure`,
        `cna_cutoff_angstrom`, `data_file`.

    Raises
    ------
    FileNotFoundError
        `data_file` does not exist, or `lmp_binary` is not on `PATH`.
    RuntimeError
        LAMMPS exits non-zero (stderr tail is in the message).
    """
    data_file = Path(data_file).resolve()
    if not data_file.exists():
        raise FileNotFoundError(f"data file not found: {data_file}")
    if parent_structure not in _CUTOFF_FORMULAS:
        raise ValueError(
            f"unsupported parent_structure {parent_structure!r}; "
            f"choose from {sorted(_CUTOFF_FORMULAS)}"
        )
    if shutil.which(lmp_binary) is None:
        raise FileNotFoundError(
            f"{lmp_binary!r} is not on PATH. On Euler, first run:\n"
            f"  module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4"
        )

    if cna_cutoff is None:
        cna_cutoff = _CUTOFF_FORMULAS[parent_structure](lattice_a)

    with tempfile.TemporaryDirectory(prefix="gb_cna_") as raw_wd:
        wd = Path(raw_wd)
        deck = wd / "cna.lammps"
        dump = wd / "cna.dump"
        _write_cna_deck(deck, data_file, dump.name, cna_cutoff)

        proc = subprocess.run(
            [lmp_binary, "-in", str(deck), "-log", "cna.log"],
            cwd=wd,
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"LAMMPS CNA run failed (exit {proc.returncode}).\n"
                f"--- stderr tail ---\n{proc.stderr[-2000:]}\n"
                f"--- stdout tail ---\n{proc.stdout[-2000:]}"
            )

        cna_values = _parse_cna_dump(dump)

        if out_dump_path is not None:
            out_dump_path = Path(out_dump_path)
            out_dump_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy(dump, out_dump_path)

    bulk_int = _BULK_CNA_INT[parent_structure]
    mask = cna_values != bulk_int
    counts = {
        label: int((cna_values == code).sum())
        for code, label in _CNA_LABELS.items()
    }
    info = {
        "n_atoms": int(len(cna_values)),
        "n_gb": int(mask.sum()),
        "f_gb": float(mask.mean()),
        "cna_counts": counts,
        "parent_structure": parent_structure,
        "bulk_cna_int": bulk_int,
        "cna_cutoff_angstrom": float(cna_cutoff),
        "data_file": str(data_file),
    }
    return mask, info


def _write_cna_deck(deck_path: Path, data_file: Path, dump_name: str, cutoff: float) -> None:
    """Write the minimal LAMMPS deck used for CNA classification.

    `pair_style zero` gives LAMMPS a well-defined neighbor list without
    computing forces — CNA is purely geometric so no real interaction is
    needed. `run 0` executes only one force/neighbor build pass.
    """
    # LAMMPS requires `mass` for every declared atom type even with run 0.
    # CNA is geometric, so placeholder masses (1.0) are fine. We read the
    # type count from the data file so the deck adapts to 1 / 2 / N types.
    n_types = _count_atom_types(data_file)
    mass_lines = "\n".join(f"mass         {t} 1.0" for t in range(1, n_types + 1))
    deck_path.write_text(
        "units        metal\n"
        "atom_style   atomic\n"
        "boundary     p p p\n"
        f"read_data    {data_file}\n"
        f"{mass_lines}\n"
        "pair_style   zero 10.0\n"
        "pair_coeff   * *\n"
        f"compute      cna all cna/atom {cutoff}\n"
        # id + type + coords + cna so downstream OVITO can colour-by-CNA
        # directly from this dump without loading the data file again.
        f"dump         d all custom 1 {dump_name} id type x y z c_cna\n"
        "dump_modify  d sort id\n"
        "run          0\n"
    )


def _count_atom_types(data_file: Path) -> int:
    """Read the `N atom types` line out of a LAMMPS data-file header."""
    with open(data_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3 and parts[1] == "atom" and parts[2] == "types":
                return int(parts[0])
            if line.strip() == "Atoms":
                break
    raise ValueError(f"no 'N atom types' header line in {data_file}")


def _parse_cna_dump(dump_path: Path) -> np.ndarray:
    """Parse the CNA dump; return id-sorted CNA ints (last column)."""
    with open(dump_path) as f:
        lines = f.readlines()
    n_atoms = None
    start = None
    for i, line in enumerate(lines):
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            n_atoms = int(lines[i + 1])
        if line.startswith("ITEM: ATOMS"):
            start = i + 1
            break
    if start is None:
        raise ValueError(f"no 'ITEM: ATOMS' section in {dump_path}")
    arr = np.loadtxt(lines[start:], dtype=float)
    if n_atoms is not None and len(arr) != n_atoms:
        raise ValueError(
            f"dump truncated: header says {n_atoms} atoms, file has {len(arr)} rows"
        )
    # Last column is c_cna regardless of how many prepend columns were added.
    return arr[:, -1].astype(np.int8)


def _cli() -> None:
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("data_file", help="LAMMPS data file (post-anneal polycrystal)")
    p.add_argument(
        "--parent",
        choices=sorted(_CUTOFF_FORMULAS),
        default="fcc",
        help="parent crystal structure for bulk classification",
    )
    p.add_argument(
        "--lattice-a",
        type=float,
        default=4.05,
        help="lattice parameter [Å] for default CNA cutoff (default = Al). "
             "Ignored when --cutoff is given.",
    )
    p.add_argument(
        "--cutoff",
        type=float,
        default=None,
        help="CNA neighbor cutoff [Å]; overrides --lattice-a",
    )
    p.add_argument("--lmp", default="lmp", help="LAMMPS binary name/path")
    p.add_argument(
        "--out-mask", default=None, help="write bool GB mask as .npy here"
    )
    p.add_argument(
        "--out-report",
        default=None,
        help="write JSON with f_gb + CNA counts here",
    )
    p.add_argument(
        "--out-dump",
        default=None,
        help="copy CNA-labelled LAMMPS dump here (for OVITO)",
    )
    args = p.parse_args()

    mask, info = compute_gb_mask(
        args.data_file,
        parent_structure=args.parent,
        lattice_a=args.lattice_a,
        cna_cutoff=args.cutoff,
        lmp_binary=args.lmp,
        out_dump_path=args.out_dump,
    )

    if args.out_mask:
        Path(args.out_mask).parent.mkdir(parents=True, exist_ok=True)
        np.save(args.out_mask, mask)
    if args.out_report:
        Path(args.out_report).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_report).write_text(json.dumps(info, indent=2))

    print(f"n_atoms         = {info['n_atoms']}")
    print(f"n_gb            = {info['n_gb']}")
    print(f"f_gb            = {info['f_gb']:.4f}")
    print(f"cna_cutoff      = {info['cna_cutoff_angstrom']:.3f} Å")
    print(f"parent          = {info['parent_structure']}")
    print(f"cna_counts      = {info['cna_counts']}")


if __name__ == "__main__":
    _cli()
