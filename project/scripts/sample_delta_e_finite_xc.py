"""Phase 3 (extended) — per-site ΔE_eff spectrum on finite-X_c substrate.

Companion to ``sample_delta_e.py``: extends per-site Mg-insertion ΔE
measurement from a pure-Al substrate to an arbitrary Al-Mg substrate
that already carries Mg solute (e.g. an HMC-equilibrated finite-X_c
configuration). The resulting ΔE_eff(i; X_c) spectrum is the central
quantity for testing Wagih's independent-site assumption

    ΔE_i(X_c) = ΔE_i(0)    ∀ X_c

at finite background concentration. Compared to the parent script:

    1. Substrate may be any LAMMPS data file containing types {1, 2}
       (Al, Mg). The script reads atom types in addition to positions.
    2. GB site selection filters to **empty** GB sites (currently
       type==1) — already-Mg sites cannot accept another Mg.
    3. Bulk reference selection filters to type==1 atoms that are ALSO
       ≥ bulk_separation Å from any Mg atom (not just any GB atom),
       so the μ_Mg^bulk reference is not contaminated by neighbouring
       solute. The Mg-distance and GB-distance use the same
       bulk_separation value.
    4. Output metadata records X_c_actual = N_Mg / N_total and
       X_gb_actual = N_Mg^GB / N_GB of the substrate, plus the
       background_data_file path.

Everything else (per-site CG-relax loop, checkpoint via _run_meta.json
+ _results.csv, ΔE_eff = E_GB^Mg − mean(E_bulk^Mg) definition) is
unchanged from `sample_delta_e.py`. The bulk_e_mean reference is
recomputed per substrate (each X_c gets its own μ_Mg^bulk) — this
keeps the per-substrate spectrum self-consistent and isolates the
shift in segregation energies from any X_c-dependence of the bulk
chemical potential.

Usage:
    python sample_delta_e_finite_xc.py \\
        --background-data data/snapshots/hmc_T500_Xc0.10_preseg_final.lmp \\
        --gb-mask gb_mask.npy \\
        --potential /.../Al-Mg.eam.fs \\
        --n-gb 500 --n-bulk 5 --seed 42 \\
        --mpi-ranks 32 \\
        --xc-label 0.10 \\
        --work-dir /cluster/scratch/.../delta_e_xc0.10 \\
        --out-npz delta_e_results_xc0.10.npz \\
        --out-json delta_e_meta_xc0.10.json

Requires `lmp` on PATH. On Euler:
    module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


_DEFAULT_ELEMENTS = ("Al", "Mg")
_DEFAULT_MASSES = (26.9815, 24.3050)

_SITE_ENERGY_RE = re.compile(r"SITE_ENERGY\s+site_id=(\d+)\s+pe=(-?\d+\.\d+(?:[eE][+-]?\d+)?)")
_STOP_CRITERION_RE = re.compile(r"Stopping criterion\s*=\s*(.+)")

_CSV_FIELDS = ("idx", "site_id", "is_gb", "pe", "stop_reason", "wall_s")
_META_KEYS_TO_MATCH = (
    "background_data_file", "sample_seed", "n_gb_sites", "n_bulk_refs",
    "bulk_separation", "gb_sample_ids", "bulk_sample_ids",
    "cg_etol", "cg_ftol", "cg_maxiter", "cg_maxeval",
)


def compute_delta_e_spectrum(
    background_data_file: str | Path,
    gb_mask: np.ndarray | str | Path,
    *,
    potential_file: str | Path,
    n_gb_sites: int = 50,
    n_bulk_refs: int = 5,
    solute_type: int = 2,
    bulk_separation: float = 8.0,
    sample_seed: int = 0,
    cg_etol: float = 1.0e-25,
    cg_ftol: float = 1.0e-25,
    cg_maxiter: int = 50000,
    cg_maxeval: int = 5000000,
    lmp_binary: str = "lmp",
    mpi_launcher: list[str] | None = None,
    work_dir: str | Path | None = None,
    keep_per_site_files: bool = False,
    elements: tuple[str, str] = _DEFAULT_ELEMENTS,
    masses: tuple[float, float] = _DEFAULT_MASSES,
) -> dict:
    """Sample empty GB + clean bulk sites on a finite-X_c substrate, run
    per-site Al→Mg CG relaxation, return ΔE_eff.

    Differences from sample_delta_e.compute_delta_e_spectrum:
        - GB site pool is restricted to mask atoms with type==1 (empty GB).
        - Bulk reference pool is restricted to type==1 atoms that are
          ≥ bulk_separation Å from any GB atom AND any Mg atom.
        - Background X_c, X_gb are computed from substrate and stored in
          the result.
        - bulk_e_mean is recomputed per substrate → ΔE_eff is self-consistent
          to the test X_c.

    Parameters
    ----------
    background_data_file : str or Path
        LAMMPS data file (atom_style atomic) with two atom types:
        type 1 = Al, type 2 = Mg. May be pure Al (X_c=0) or any
        finite-X_c configuration (e.g. HMC-equilibrated snapshot).
    gb_mask : np.ndarray of bool or path to .npy
        Per-atom GB mask; length must equal atom count of the
        background substrate. The mask should reflect GB atoms in the
        underlying lattice geometry, not whether each atom is currently
        Al or Mg.
    potential_file : str or Path
        Absolute path to the Al-Mg EAM/fs file (e.g. NIST Mendelev 2009).
    n_gb_sites : int
        How many empty (type==1) GB atoms to sample.
    n_bulk_refs : int
        How many clean bulk reference atoms to sample and average.
    solute_type : int
        LAMMPS atom type for Mg (default 2 to match our data file convention).
    bulk_separation : float
        Minimum distance (Å) from a bulk reference candidate to the
        nearest GB atom AND the nearest Mg atom. Default 8 Å ≈ 2 NN
        shells in FCC Al.
    sample_seed : int
        RNG seed for site selection.
    cg_etol, cg_ftol, cg_maxiter, cg_maxeval
        LAMMPS `minimize` tolerances; defaults match `sample_delta_e.py`
        (1e-25 / 1e-25 / 50000 / 5000000).
    lmp_binary, mpi_launcher, work_dir, keep_per_site_files, elements, masses
        See `sample_delta_e.compute_delta_e_spectrum` for semantics.

    Returns
    -------
    result : dict
        Keys:
            gb_site_ids, gb_e_mg, gb_delta_e, gb_stop_reasons,
            bulk_ref_ids, bulk_e_mg, bulk_e_mean, bulk_e_std,
            bulk_stop_reasons,
            n_atoms, n_type1, n_type2, n_gb_total, n_gb_type2,
            x_c_actual            : float, n_type2 / n_atoms of substrate
            x_gb_actual           : float, n_gb_type2 / n_gb_total of substrate
            background_data_file  : str
            potential_file, sample_seed, bulk_separation,
            cg_etol, cg_ftol, cg_maxiter, cg_maxeval,
            wall_seconds_total, wall_seconds_per_site, n_sites_resumed.

    Raises
    ------
    FileNotFoundError, ValueError, RuntimeError
        On missing files, mask/atom-count mismatch, too few empty GB
        candidates, too few clean bulk Al candidates after the
        Mg-distance filter, or LAMMPS non-zero exit.
    """
    background_data_file = Path(background_data_file).resolve()
    potential_file = Path(potential_file).resolve()
    if not background_data_file.exists():
        raise FileNotFoundError(f"background data file not found: {background_data_file}")
    if not potential_file.exists():
        raise FileNotFoundError(f"potential file not found: {potential_file}")
    if shutil.which(lmp_binary) is None:
        raise FileNotFoundError(
            f"{lmp_binary!r} is not on PATH. On Euler:\n"
            f"  module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4"
        )

    if isinstance(gb_mask, (str, Path)):
        gb_mask = np.load(gb_mask)
    gb_mask = np.asarray(gb_mask, dtype=bool)

    # Read positions AND types so we can filter by Al/Mg.
    positions, types, box = _read_lammps_data_positions_types(background_data_file)
    n_atoms = positions.shape[0]
    if gb_mask.shape != (n_atoms,):
        raise ValueError(
            f"gb_mask length {gb_mask.shape} does not match n_atoms={n_atoms}"
        )

    # Background composition statistics ------------------------------------
    n_type1 = int((types == 1).sum())
    n_type2 = int((types == 2).sum())
    n_gb_total = int(gb_mask.sum())
    n_gb_type2 = int(((types == 2) & gb_mask).sum())
    x_c_actual = n_type2 / n_atoms
    x_gb_actual = n_gb_type2 / n_gb_total if n_gb_total > 0 else 0.0

    # Site selection --------------------------------------------------------
    rng = np.random.default_rng(sample_seed)

    # GB candidates: GB atoms that are currently Al (type==1).
    gb_ids_0based = np.where(gb_mask & (types == 1))[0]
    if len(gb_ids_0based) < n_gb_sites:
        raise ValueError(
            f"only {len(gb_ids_0based)} EMPTY (Al) GB atoms available "
            f"(of {n_gb_total} total GB; X_GB={x_gb_actual:.3f}); "
            f"requested {n_gb_sites}. Lower --n-gb or pick a less-saturated background."
        )
    gb_sampled = rng.choice(gb_ids_0based, size=n_gb_sites, replace=False)

    bulk_sampled = _sample_bulk_refs(
        positions=positions,
        types=types,
        box=box,
        gb_mask=gb_mask,
        n_bulk_refs=n_bulk_refs,
        bulk_separation=bulk_separation,
        rng=rng,
    )

    all_site_ids_0based = np.concatenate([gb_sampled, bulk_sampled])
    all_site_ids = all_site_ids_0based + 1        # LAMMPS ids are 1-based

    # Sanity: every selected site must currently be Al; setting type 2 on
    # an already-Mg site would yield a zero ΔE silently.
    selected_types = types[all_site_ids_0based]
    if not np.all(selected_types == 1):
        raise RuntimeError(
            f"internal bug: {(selected_types != 1).sum()} selected sites are "
            f"not type 1. Site selection filter must be wrong."
        )

    # Work directory --------------------------------------------------------
    if work_dir is None:
        work_dir = background_data_file.parent / "delta_e_run"
    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    mpi_launcher = list(mpi_launcher) if mpi_launcher else []

    # Checkpoint setup ------------------------------------------------------
    meta_expected = {
        "background_data_file": str(background_data_file),
        "sample_seed": sample_seed,
        "n_gb_sites": n_gb_sites,
        "n_bulk_refs": n_bulk_refs,
        "bulk_separation": bulk_separation,
        "gb_sample_ids": gb_sampled.tolist(),
        "bulk_sample_ids": bulk_sampled.tolist(),
    }
    existing, csv_path = _load_or_init_checkpoint(work_dir, meta_expected)

    # Per-site LAMMPS runs --------------------------------------------------
    energies = np.full(len(all_site_ids), np.nan, dtype=float)
    stop_reasons = np.array(["" for _ in all_site_ids], dtype=object)
    per_site_time = np.zeros(len(all_site_ids), dtype=float)
    n_resumed = 0
    t0_total = time.time()
    for k, site_id in enumerate(all_site_ids):
        sid = int(site_id)
        if sid in existing:
            energies[k] = existing[sid]["pe"]
            stop_reasons[k] = existing[sid]["stop_reason"]
            per_site_time[k] = existing[sid]["wall_s"]
            n_resumed += 1
            print(
                f"  [{k+1:4d}/{len(all_site_ids)}] site_id={sid:7d}  "
                f"pe={energies[k]:14.6f} eV  (cached, {stop_reasons[k]})",
                flush=True,
            )
            continue

        t0 = time.time()
        deck = work_dir / f"site_{sid}.lammps"
        logf = work_dir / f"site_{sid}.log"
        _write_site_deck(
            deck_path=deck,
            log_path=logf,
            data_file=background_data_file,
            potential_file=potential_file,
            site_id=sid,
            solute_type=solute_type,
            etol=cg_etol,
            ftol=cg_ftol,
            maxiter=cg_maxiter,
            maxeval=cg_maxeval,
            elements=elements,
            masses=masses,
        )
        cmd = [*mpi_launcher, lmp_binary, "-in", str(deck)]
        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
        if proc.returncode != 0:
            raise RuntimeError(
                f"LAMMPS failed at site id={sid} (k={k}, exit {proc.returncode}).\n"
                f"--- stderr tail ---\n{proc.stderr[-2000:]}\n"
                f"--- stdout tail ---\n{proc.stdout[-2000:]}"
            )
        pe, stop = _parse_site_result(proc.stdout, sid)
        energies[k] = pe
        stop_reasons[k] = stop
        per_site_time[k] = time.time() - t0
        _append_result_row(
            csv_path,
            {
                "idx": k,
                "site_id": sid,
                "is_gb": int(k < n_gb_sites),
                "pe": pe,
                "stop_reason": stop,
                "wall_s": per_site_time[k],
            },
        )
        if not keep_per_site_files:
            deck.unlink(missing_ok=True)
            logf.unlink(missing_ok=True)
        print(
            f"  [{k+1:4d}/{len(all_site_ids)}] site_id={sid:7d}  "
            f"pe={energies[k]:14.6f} eV  ({per_site_time[k]:.1f} s, {stop})",
            flush=True,
        )

    wall_total = time.time() - t0_total

    # Split GB and bulk results --------------------------------------------
    gb_e_mg = energies[:n_gb_sites]
    bulk_e_mg = energies[n_gb_sites:]
    bulk_e_mean = float(np.mean(bulk_e_mg))
    bulk_e_std = float(np.std(bulk_e_mg, ddof=1)) if n_bulk_refs > 1 else 0.0
    gb_delta_e = gb_e_mg - bulk_e_mean

    return {
        "gb_site_ids": all_site_ids[:n_gb_sites].astype(np.int64),
        "gb_e_mg": gb_e_mg,
        "gb_delta_e": gb_delta_e,
        "gb_stop_reasons": np.asarray(stop_reasons[:n_gb_sites], dtype=str),
        "bulk_ref_ids": all_site_ids[n_gb_sites:].astype(np.int64),
        "bulk_e_mg": bulk_e_mg,
        "bulk_e_mean": bulk_e_mean,
        "bulk_e_std": bulk_e_std,
        "bulk_stop_reasons": np.asarray(stop_reasons[n_gb_sites:], dtype=str),
        "n_atoms": n_atoms,
        "n_type1": n_type1,
        "n_type2": n_type2,
        "n_gb_total": n_gb_total,
        "n_gb_type2": n_gb_type2,
        "x_c_actual": x_c_actual,
        "x_gb_actual": x_gb_actual,
        "potential_file": str(potential_file),
        "background_data_file": str(background_data_file),
        "sample_seed": sample_seed,
        "bulk_separation": bulk_separation,
        "cg_etol": cg_etol,
        "cg_ftol": cg_ftol,
        "cg_maxiter": cg_maxiter,
        "cg_maxeval": cg_maxeval,
        "wall_seconds_total": wall_total,
        "wall_seconds_per_site": per_site_time.tolist(),
        "n_sites_resumed": n_resumed,
    }


def _sample_bulk_refs(
    *,
    positions: np.ndarray,
    types: np.ndarray,
    box: np.ndarray,
    gb_mask: np.ndarray,
    n_bulk_refs: int,
    bulk_separation: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return 0-based indices of n_bulk_refs Al atoms that are
    (a) not GB, (b) currently type 1 (Al), (c) ≥ bulk_separation Å from
    any GB atom, (d) ≥ bulk_separation Å from any Mg atom.

    Constraint (d) keeps the μ_Mg^bulk reference clean at finite X_c —
    a bulk Al sitting next to existing Mg solute would carry chemical
    contamination from solute-solute interactions.
    """
    candidate_mask = (~gb_mask) & (types == 1)
    candidate_idx = np.where(candidate_mask)[0]
    if len(candidate_idx) == 0:
        raise ValueError(
            "no bulk (~gb) Al (type==1) atoms in substrate; cannot place bulk reference"
        )

    candidate_positions = positions[candidate_idx]

    # Distance to any GB atom
    gb_positions = positions[gb_mask]
    if len(gb_positions) > 0:
        gb_tree = cKDTree(gb_positions, boxsize=box)
        gb_dists, _ = gb_tree.query(candidate_positions, k=1)
        keep = gb_dists >= bulk_separation
    else:
        keep = np.ones(len(candidate_idx), dtype=bool)

    # Distance to any Mg atom (new for finite-X_c substrate)
    mg_positions = positions[types == 2]
    if len(mg_positions) > 0:
        mg_tree = cKDTree(mg_positions, boxsize=box)
        mg_dists, _ = mg_tree.query(candidate_positions, k=1)
        keep = keep & (mg_dists >= bulk_separation)

    candidates = candidate_idx[keep]
    if len(candidates) < n_bulk_refs:
        raise ValueError(
            f"only {len(candidates)} bulk Al candidates ≥ {bulk_separation} Å "
            f"from any GB atom AND any Mg atom; need {n_bulk_refs}. "
            f"Lower --bulk-separation or pick a less-saturated substrate."
        )
    return rng.choice(candidates, size=n_bulk_refs, replace=False)


def _read_lammps_data_positions_types(
    data_file: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Minimal parser for atom_style atomic data files; returns positions,
    types, and box.

    Returns
    -------
    positions : (N, 3) float
        Atom positions, wrapped into [0, L).
    types : (N,) int
        Per-atom integer type (1 or 2 for Al-Mg). id-sorted.
    box_lengths : (3,) float

    Atoms section is `id type x y z [ix iy iz]` (atom_style atomic).
    """
    with open(data_file) as f:
        text = f.read()
    n_atoms = None
    lo = {"x": None, "y": None, "z": None}
    hi = {"x": None, "y": None, "z": None}
    atoms_section_start = None
    lines = text.splitlines()
    for i, line in enumerate(lines):
        s = line.split()
        if len(s) >= 2 and s[1] == "atoms":
            n_atoms = int(s[0])
        elif len(s) >= 4 and s[2] in ("xlo", "ylo", "zlo") and s[3] in ("xhi", "yhi", "zhi"):
            axis = s[2][0]
            lo[axis] = float(s[0])
            hi[axis] = float(s[1])
        elif line.strip() == "Atoms" or line.strip().startswith("Atoms "):
            atoms_section_start = i + 2     # skip section header + blank line
            break
    if n_atoms is None or atoms_section_start is None:
        raise ValueError(f"couldn't parse LAMMPS data header in {data_file}")
    box_lengths = np.array([hi["x"] - lo["x"], hi["y"] - lo["y"], hi["z"] - lo["z"]])

    rows = []
    for raw in lines[atoms_section_start : atoms_section_start + n_atoms]:
        parts = raw.split()
        if len(parts) < 5:
            raise ValueError(f"short Atoms line: {raw!r}")
        rows.append((
            int(parts[0]),       # id
            int(parts[1]),       # type
            float(parts[2]),     # x
            float(parts[3]),     # y
            float(parts[4]),     # z
        ))
    rows.sort(key=lambda r: r[0])    # id-sorted
    if len(rows) != n_atoms:
        raise ValueError(f"Atoms section has {len(rows)} rows, expected {n_atoms}")
    positions = np.array([[x, y, z] for _, _, x, y, z in rows])
    types = np.array([t for _, t, _, _, _ in rows], dtype=np.int64)
    positions = positions - np.array([lo["x"], lo["y"], lo["z"]])
    positions = np.mod(positions, box_lengths)
    return positions, types, box_lengths


def _write_site_deck(
    *,
    deck_path: Path,
    log_path: Path,
    data_file: Path,
    potential_file: Path,
    site_id: int,
    solute_type: int,
    etol: float,
    ftol: float,
    maxiter: int,
    maxeval: int,
    elements: tuple[str, str] = _DEFAULT_ELEMENTS,
    masses: tuple[float, float] = _DEFAULT_MASSES,
) -> None:
    """Write the LAMMPS deck that relaxes one Al→Mg substitution at site_id.

    Substrate is read as-is — types of all OTHER atoms remain whatever
    they were in the data file, so for a finite-X_c substrate the existing
    Mg distribution is preserved across the relaxation.
    """
    el1, el2 = elements
    m1, m2 = masses
    deck_path.write_text(
        f"log {log_path}\n"
        "units        metal\n"
        "atom_style   atomic\n"
        "boundary     p p p\n"
        f"read_data    {data_file}\n"
        f"mass         1 {m1}\n"
        f"mass         2 {m2}\n"
        "pair_style   eam/fs\n"
        f"pair_coeff   * * {potential_file} {el1} {el2}\n"
        "neighbor     1.0 bin\n"
        "neigh_modify every 10 delay 0 check yes\n"
        f"set atom {site_id} type {solute_type}\n"
        "thermo       50\n"
        "thermo_style custom step pe fnorm\n"
        "min_style    cg\n"
        f"minimize     {etol} {ftol} {maxiter} {maxeval}\n"
        "variable     E equal pe\n"
        f'print        "SITE_ENERGY site_id={site_id} pe=${{E}}"\n'
    )


def _parse_site_result(stdout: str, expected_site_id: int) -> tuple[float, str]:
    """Return (pe, stop_reason). Raises on missing or mismatched SITE_ENERGY."""
    match = _SITE_ENERGY_RE.search(stdout)
    if match is None:
        raise RuntimeError(
            f"no SITE_ENERGY line for site {expected_site_id} in LAMMPS stdout.\n"
            f"--- stdout tail ---\n{stdout[-2000:]}"
        )
    found_id = int(match.group(1))
    if found_id != expected_site_id:
        raise RuntimeError(
            f"SITE_ENERGY id mismatch: expected {expected_site_id}, got {found_id}"
        )
    pe = float(match.group(2))
    stop_match = _STOP_CRITERION_RE.search(stdout)
    stop_reason = stop_match.group(1).strip() if stop_match else "unknown"
    return pe, stop_reason


def _load_or_init_checkpoint(
    work_dir: Path, meta_expected: dict
) -> tuple[dict[int, dict], Path]:
    """Same protocol as the parent script — fresh: write meta, return empty.
    Resume: verify meta matches; load completed per-site rows."""
    meta_file = work_dir / "_run_meta.json"
    csv_file = work_dir / "_results.csv"

    if meta_file.exists():
        saved = json.loads(meta_file.read_text())
        for key in _META_KEYS_TO_MATCH:
            if saved.get(key) != meta_expected.get(key):
                raise ValueError(
                    f"checkpoint in {meta_file} mismatches current params:\n"
                    f"  key={key!r}\n"
                    f"  saved    = {saved.get(key)!r}\n"
                    f"  current  = {meta_expected.get(key)!r}\n"
                    f"Pick a different --work-dir or delete {work_dir!s}."
                )
    else:
        meta_file.write_text(json.dumps(meta_expected, indent=2))

    existing: dict[int, dict] = {}
    if csv_file.exists():
        with open(csv_file) as f:
            for row in csv.DictReader(f):
                existing[int(row["site_id"])] = {
                    "pe": float(row["pe"]),
                    "stop_reason": row["stop_reason"],
                    "wall_s": float(row["wall_s"]),
                }
    return existing, csv_file


def _append_result_row(csv_path: Path, row: dict) -> None:
    is_new = not csv_path.exists()
    with open(csv_path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=_CSV_FIELDS)
        if is_new:
            w.writeheader()
        w.writerow(row)


def _save_results(result: dict, npz_path: Path, json_path: Path | None = None) -> None:
    """Persist numerics to .npz and metadata to .json."""
    npz_keys = (
        "gb_site_ids", "gb_e_mg", "gb_delta_e", "gb_stop_reasons",
        "bulk_ref_ids", "bulk_e_mg", "bulk_stop_reasons",
    )
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        npz_path,
        **{k: result[k] for k in npz_keys},
        bulk_e_mean=np.array(result["bulk_e_mean"]),
        bulk_e_std=np.array(result["bulk_e_std"]),
        x_c_actual=np.array(result["x_c_actual"]),
        x_gb_actual=np.array(result["x_gb_actual"]),
    )
    if json_path is not None:
        meta = {k: v for k, v in result.items() if k not in npz_keys}
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(meta, indent=2, default=float))


def _cli() -> None:
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--background-data", required=True,
                   help="LAMMPS data file (Al-Mg, types {1,2}); may be pure Al or finite-X_c")
    p.add_argument("--gb-mask", required=True, help=".npy bool mask from gb_identify.py")
    p.add_argument("--potential", required=True, help="Al-Mg EAM/fs file")
    p.add_argument("--n-gb", type=int, default=50, help="number of empty GB sites to sample")
    p.add_argument("--n-bulk", type=int, default=5, help="number of clean bulk ref sites")
    p.add_argument("--seed", type=int, default=0, help="sample_seed for site selection")
    p.add_argument("--solute-type", type=int, default=2, help="LAMMPS type for solute")
    p.add_argument("--elements", default="Al Mg",
                   help='Element symbols for type 1 / type 2 in the EAM file '
                        '(e.g. "Cu Ni"). Order must match the EAM/fs file header.')
    p.add_argument("--masses", default="26.9815 24.305",
                   help='Atomic masses (g/mol) for type 1 / type 2.')
    p.add_argument("--bulk-separation", type=float, default=8.0,
                   help="min Å from any GB atom AND any Mg atom for a bulk reference")
    p.add_argument("--xc-label", default=None,
                   help="optional label for stdout/metadata (e.g. '0.10'); does not affect "
                        "computation, only stored in result for traceability")
    p.add_argument("--etol", type=float, default=1.0e-25)
    p.add_argument("--ftol", type=float, default=1.0e-25)
    p.add_argument("--maxiter", type=int, default=50000)
    p.add_argument("--maxeval", type=int, default=5000000)
    p.add_argument("--lmp", default="lmp", help="LAMMPS binary")
    p.add_argument("--mpi-ranks", type=int, default=1,
                   help="MPI ranks per LAMMPS invocation (0 = no launcher)")
    p.add_argument("--mpi-cmd", default=None,
                   help="launcher binary; auto-detect srun vs mpirun if omitted")
    p.add_argument("--work-dir", default=None,
                   help="per-site deck/log dir (default: <substrate dir>/delta_e_run)")
    p.add_argument("--keep-per-site-files", action="store_true",
                   help="keep individual site decks + logs after run")
    p.add_argument("--out-npz", required=True, help="output .npz path")
    p.add_argument("--out-json", default=None, help="output metadata .json")
    args = p.parse_args()

    if args.mpi_ranks <= 1:
        launcher = None
    else:
        cmd = args.mpi_cmd
        if cmd is None:
            cmd = "srun" if "SLURM_JOB_ID" in os.environ else "mpirun"
        if cmd == "srun":
            launcher = ["srun", "-n", str(args.mpi_ranks)]
        elif cmd == "mpirun":
            launcher = ["mpirun", "-np", str(args.mpi_ranks)]
        else:
            launcher = [cmd, "-np", str(args.mpi_ranks)]

    el_tuple = tuple(args.elements.split())
    if len(el_tuple) != 2:
        raise SystemExit(f"--elements must be 2 symbols, got {el_tuple!r}")
    mass_tuple = tuple(float(x) for x in args.masses.split())
    if len(mass_tuple) != 2:
        raise SystemExit(f"--masses must be 2 floats, got {mass_tuple!r}")

    result = compute_delta_e_spectrum(
        background_data_file=args.background_data,
        gb_mask=args.gb_mask,
        potential_file=args.potential,
        n_gb_sites=args.n_gb,
        n_bulk_refs=args.n_bulk,
        solute_type=args.solute_type,
        bulk_separation=args.bulk_separation,
        sample_seed=args.seed,
        cg_etol=args.etol,
        cg_ftol=args.ftol,
        cg_maxiter=args.maxiter,
        cg_maxeval=args.maxeval,
        lmp_binary=args.lmp,
        mpi_launcher=launcher,
        work_dir=args.work_dir,
        keep_per_site_files=args.keep_per_site_files,
        elements=el_tuple,
        masses=mass_tuple,
    )

    if args.xc_label is not None:
        result["xc_label"] = args.xc_label

    _save_results(result, Path(args.out_npz), Path(args.out_json) if args.out_json else None)

    dE = result["gb_delta_e"]
    all_stops = np.concatenate([result["gb_stop_reasons"], result["bulk_stop_reasons"]])
    bad = [r for r in all_stops if r not in ("energy tolerance", "force tolerance")]
    label_str = f" (xc_label={args.xc_label})" if args.xc_label else ""
    print(
        "\n"
        f"background      = {result['background_data_file']}{label_str}\n"
        f"X_c_actual      = {result['x_c_actual']:.4f}   "
        f"(n_type2={result['n_type2']}, n_atoms={result['n_atoms']})\n"
        f"X_gb_actual     = {result['x_gb_actual']:.4f}   "
        f"(n_gb_type2={result['n_gb_type2']}, n_gb={result['n_gb_total']})\n"
        f"n_gb_sites       = {len(dE)}\n"
        f"n_sites_resumed  = {result['n_sites_resumed']}/{len(all_stops)}\n"
        f"CG stop reasons  = {dict((r, int((all_stops==r).sum())) for r in set(all_stops.tolist()))}\n"
        + (f"  ⚠ {len(bad)} site(s) did NOT reach energy/force tolerance\n" if bad else "")
        + f"bulk_e_mean      = {result['bulk_e_mean']:.4f} eV   "
        f"(std {result['bulk_e_std']:.4f}, n={len(result['bulk_e_mg'])})\n"
        f"ΔE_eff [eV]      min={dE.min():+.4f}  max={dE.max():+.4f}  "
        f"mean={dE.mean():+.4f}  median={np.median(dE):+.4f}\n"
        f"ΔE_eff [kJ/mol]  min={dE.min()*96.485:+.2f}  max={dE.max()*96.485:+.2f}  "
        f"mean={dE.mean()*96.485:+.2f}\n"
        f"wall_time_total  = {result['wall_seconds_total']:.1f} s\n",
        flush=True,
    )


if __name__ == "__main__":
    _cli()
