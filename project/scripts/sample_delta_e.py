"""Phase 3 — per-site ΔE_seg spectrum via single-substitution LAMMPS relaxations.

Given a relaxed pure-Al polycrystal (post-anneal LAMMPS data file), a per-atom
GB boolean mask (from `gb_identify.py`), and the Al-Mg EAM potential, this
driver:

    1. Samples `n_gb_sites` GB atoms uniformly at random from `mask == True`.
    2. Samples `n_bulk_refs` bulk reference atoms that are at least
       `bulk_separation` Å from the nearest GB atom (interior of a grain).
    3. For each sampled site, writes a small LAMMPS deck that:
           read_data <annealed polycrystal>
           pair_style eam/fs ; pair_coeff * * <potfile> Al Mg
           set atom <site_id> type 2     # Al -> Mg
           minimize <etol> <ftol> ...
           print "SITE_ENERGY site_id=<id> pe=<pe>"
       runs LAMMPS, and parses the relaxed total energy.
    4. Computes
           ΔE_seg(i) = E_GB^Mg(i) − mean(E_bulk^Mg)
       for each GB site, and reports the mean / spread of the bulk reference.

The result is a reproducible {site_ids, ΔE, metadata} bundle ready to feed
into a skew-normal fit and compare against Wagih's Al(Mg) Fig 2 parameters.

One LAMMPS process per site (overhead ~2 s start-up × 55 sites ≈ 2 min total,
negligible next to the CG cost). This keeps state cleanly isolated — no
cross-contamination between relaxations — at the price of some wasted init.

Usage as a module:
    import numpy as np
    from sample_delta_e import compute_delta_e_spectrum
    result = compute_delta_e_spectrum(
        "poly_Al_100A_8g_annealed.lmp",
        gb_mask=np.load("gb_mask.npy"),
        potential_file="/.../Al-Mg.eam.fs",
        n_gb_sites=50, n_bulk_refs=5, sample_seed=42,
    )
    # result["gb_delta_e"], result["gb_site_ids"], result["bulk_e_mean"], ...

Usage from the shell:
    python sample_delta_e.py \\
        --annealed poly_Al_100A_8g_annealed.lmp \\
        --gb-mask gb_mask.npy \\
        --potential /.../Al-Mg.eam.fs \\
        --n-gb 50 --n-bulk 5 --seed 42 \\
        --mpi-ranks 16 \\
        --work-dir /cluster/scratch/.../delta_e_run \\
        --out-npz delta_e_results.npz \\
        --out-json delta_e_meta.json

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


# Default masses (g/mol) for Al(1) / Mg(2) in an atom_style atomic data file.
_MASSES = {1: 26.9815, 2: 24.3050}

_SITE_ENERGY_RE = re.compile(r"SITE_ENERGY\s+site_id=(\d+)\s+pe=(-?\d+\.\d+(?:[eE][+-]?\d+)?)")
_STOP_CRITERION_RE = re.compile(r"Stopping criterion\s*=\s*(.+)")

_CSV_FIELDS = ("idx", "site_id", "is_gb", "pe", "stop_reason", "wall_s")
_META_KEYS_TO_MATCH = (
    "annealed_file", "sample_seed", "n_gb_sites", "n_bulk_refs",
    "bulk_separation", "gb_sample_ids", "bulk_sample_ids",
)


def compute_delta_e_spectrum(
    annealed_data_file: str | Path,
    gb_mask: np.ndarray | str | Path,
    *,
    potential_file: str | Path,
    n_gb_sites: int = 50,
    n_bulk_refs: int = 5,
    solute_type: int = 2,
    bulk_separation: float = 8.0,
    sample_seed: int = 0,
    cg_etol: float = 1.0e-8,
    cg_ftol: float = 1.0e-10,
    cg_maxiter: int = 5000,
    cg_maxeval: int = 50000,
    lmp_binary: str = "lmp",
    mpi_launcher: list[str] | None = None,
    work_dir: str | Path | None = None,
    keep_per_site_files: bool = False,
) -> dict:
    """Sample GB + bulk sites, run per-site Al→Mg CG relaxation, return ΔE.

    Parameters
    ----------
    annealed_data_file : str or Path
        Post-anneal pure-Al polycrystal (LAMMPS data file, atom_style atomic,
        two atom types declared with all atoms type 1).
    gb_mask : np.ndarray of bool or path to .npy
        Per-atom GB mask (True = GB atom); length must equal atom count.
    potential_file : str or Path
        Absolute path to the Al-Mg EAM/fs file (e.g. NIST Mendelev 2009).
    n_gb_sites : int
        How many GB atoms to sample.
    n_bulk_refs : int
        How many bulk reference atoms to sample and average.
    solute_type : int
        LAMMPS atom type for Mg (default 2 to match our data file convention).
    bulk_separation : float
        Minimum distance (Å) from a bulk ref atom to the nearest GB atom.
        Default 8 Å ≈ 2 NN shells in FCC Al; chosen to avoid the GB's
        local relaxation field contaminating the bulk reference.
    sample_seed : int
        RNG seed for site selection. Independent from structure/solute/swap/
        velocity seeds elsewhere in the pipeline.
    cg_etol, cg_ftol, cg_maxiter, cg_maxeval
        LAMMPS `minimize` tolerances. Tighter than the anneal deck (ΔE can
        be < 1 kJ/mol = 0.01 eV so we want sub-meV numerics).
    lmp_binary : str
        LAMMPS executable on PATH.
    mpi_launcher : list[str], optional
        Command prefix for MPI launch, e.g. ``["srun", "-n", "16"]`` inside
        a SLURM job or ``["mpirun", "-np", "16"]`` on a login node. None
        (default) runs LAMMPS without MPI launcher (serial).
    work_dir : str or Path, optional
        Persistent work directory for per-site decks + logs. Defaults to
        a subdirectory of the annealed file's parent, auto-created.
    keep_per_site_files : bool
        If True, preserve each site's deck + log file; otherwise delete
        after energy is extracted.

    Returns
    -------
    result : dict
        Keys:
            gb_site_ids        (n_gb_sites,) int, 1-based LAMMPS atom ids
            gb_e_mg            (n_gb_sites,) float, relaxed PE [eV]
            gb_delta_e         (n_gb_sites,) float, ΔE_seg [eV]
            bulk_ref_ids       (n_bulk_refs,) int
            bulk_e_mg          (n_bulk_refs,) float
            bulk_e_mean        float
            bulk_e_std         float
            n_atoms            int
            potential_file     str
            annealed_file      str
            sample_seed        int
            bulk_separation    float
            wall_seconds_total float
            wall_seconds_per_site list[float]

    Raises
    ------
    FileNotFoundError, ValueError, RuntimeError
        On missing files, mask/atom-count mismatch, too few bulk candidates,
        or LAMMPS non-zero exit.
    """
    annealed_data_file = Path(annealed_data_file).resolve()
    potential_file = Path(potential_file).resolve()
    if not annealed_data_file.exists():
        raise FileNotFoundError(f"annealed file not found: {annealed_data_file}")
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

    # Read positions so we can enforce the bulk-ref distance constraint.
    positions, box = _read_lammps_data_positions(annealed_data_file)
    n_atoms = positions.shape[0]
    if gb_mask.shape != (n_atoms,):
        raise ValueError(
            f"gb_mask length {gb_mask.shape} does not match n_atoms={n_atoms}"
        )

    # Site selection --------------------------------------------------------
    rng = np.random.default_rng(sample_seed)

    gb_ids_0based = np.where(gb_mask)[0]
    if len(gb_ids_0based) < n_gb_sites:
        raise ValueError(
            f"only {len(gb_ids_0based)} GB atoms available; requested {n_gb_sites}"
        )
    gb_sampled = rng.choice(gb_ids_0based, size=n_gb_sites, replace=False)

    bulk_sampled = _sample_bulk_refs(
        positions=positions,
        box=box,
        gb_mask=gb_mask,
        n_bulk_refs=n_bulk_refs,
        bulk_separation=bulk_separation,
        rng=rng,
    )

    all_site_ids_0based = np.concatenate([gb_sampled, bulk_sampled])
    all_site_ids = all_site_ids_0based + 1        # LAMMPS ids are 1-based

    # Work directory --------------------------------------------------------
    if work_dir is None:
        work_dir = annealed_data_file.parent / "delta_e_run"
    work_dir = Path(work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    mpi_launcher = list(mpi_launcher) if mpi_launcher else []

    # Checkpoint setup ------------------------------------------------------
    # Per-site results stream to `_results.csv` after each LAMMPS run. If this
    # function is rerun with matching params against the same work_dir, sites
    # already completed are skipped and read back from the csv. Prevents
    # losing 40/55 sites of compute when a SLURM job hits time limit.
    meta_expected = {
        "annealed_file": str(annealed_data_file),
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
            data_file=annealed_data_file,
            potential_file=potential_file,
            site_id=sid,
            solute_type=solute_type,
            etol=cg_etol,
            ftol=cg_ftol,
            maxiter=cg_maxiter,
            maxeval=cg_maxeval,
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
        "potential_file": str(potential_file),
        "annealed_file": str(annealed_data_file),
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
    box: np.ndarray,
    gb_mask: np.ndarray,
    n_bulk_refs: int,
    bulk_separation: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Return 0-based indices of n_bulk_refs atoms ≥ bulk_separation from any GB atom."""
    gb_positions = positions[gb_mask]
    if len(gb_positions) == 0:
        # No GB atoms — every interior atom qualifies. Just sample randomly.
        candidates = np.where(~gb_mask)[0]
    else:
        gb_tree = cKDTree(gb_positions, boxsize=box)
        bulk_idx = np.where(~gb_mask)[0]
        # KDTree.query with k=1 returns nearest-neighbor distance for each query.
        dists, _ = gb_tree.query(positions[bulk_idx], k=1)
        candidates = bulk_idx[dists >= bulk_separation]

    if len(candidates) < n_bulk_refs:
        raise ValueError(
            f"only {len(candidates)} bulk candidates ≥ {bulk_separation} Å from any "
            f"GB atom; need {n_bulk_refs}. Lower --bulk-separation or increase box size."
        )
    return rng.choice(candidates, size=n_bulk_refs, replace=False)


def _read_lammps_data_positions(data_file: Path) -> tuple[np.ndarray, np.ndarray]:
    """Minimal parser for atom_style atomic data files.

    Returns (positions (N,3) float, box_lengths (3,) float). Atoms are
    ordered by appearance in the Atoms section, which — for files written
    by our generator or LAMMPS' write_data — equals id-sorted 1..N order.
    """
    with open(data_file) as f:
        text = f.read()
    # Parse header: N atoms, box bounds.
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
            atoms_section_start = i + 2     # skip the section header + blank line
            break
    if n_atoms is None or atoms_section_start is None:
        raise ValueError(f"couldn't parse LAMMPS data header in {data_file}")
    box_lengths = np.array([hi["x"] - lo["x"], hi["y"] - lo["y"], hi["z"] - lo["z"]])

    # Atoms lines: id type x y z [ix iy iz]. We keep x,y,z from columns 2..4.
    rows = []
    for raw in lines[atoms_section_start : atoms_section_start + n_atoms]:
        parts = raw.split()
        if len(parts) < 5:
            raise ValueError(f"short Atoms line: {raw!r}")
        rows.append((int(parts[0]), float(parts[2]), float(parts[3]), float(parts[4])))
    rows.sort(key=lambda r: r[0])            # id-sorted
    positions = np.array([[x, y, z] for _, x, y, z in rows])
    if positions.shape[0] != n_atoms:
        raise ValueError(f"Atoms section has {positions.shape[0]} rows, expected {n_atoms}")
    # Wrap into [0, L) in case the file used out-of-box coordinates.
    positions = positions - np.array([lo["x"], lo["y"], lo["z"]])
    positions = np.mod(positions, box_lengths)
    return positions, box_lengths


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
) -> None:
    """Write the LAMMPS deck that relaxes one Al→Mg substitution."""
    deck_path.write_text(
        f"log {log_path}\n"
        "units        metal\n"
        "atom_style   atomic\n"
        "boundary     p p p\n"
        f"read_data    {data_file}\n"
        f"mass         1 {_MASSES[1]}\n"
        f"mass         2 {_MASSES[2]}\n"
        "pair_style   eam/fs\n"
        f"pair_coeff   * * {potential_file} Al Mg\n"
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
    # LAMMPS prints 'Stopping criterion = <reason>' once per `minimize`. Values
    # we care about: "energy tolerance" / "force tolerance" (good); "max
    # iterations" / "max force evaluations" (CG didn't fully converge —
    # still returns the last pe, but downstream should flag these);
    # "linesearch alpha is zero" (bad: stuck, pe may be unreliable).
    stop_match = _STOP_CRITERION_RE.search(stdout)
    stop_reason = stop_match.group(1).strip() if stop_match else "unknown"
    return pe, stop_reason


def _load_or_init_checkpoint(
    work_dir: Path, meta_expected: dict
) -> tuple[dict[int, dict], Path]:
    """Wire up the `_run_meta.json` / `_results.csv` checkpoint pair.

    Fresh run: write meta, return empty dict + csv path.
    Resume: verify meta matches current params (annealed file, seed, n_sites,
    sampled ids), load completed per-site rows into a dict keyed by site_id.
    Mismatch: raise so the user doesn't unknowingly merge incompatible data.
    """
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
    """Append one per-site result line to the checkpoint csv (write header once)."""
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
    p.add_argument("--annealed", required=True, help="annealed LAMMPS data file")
    p.add_argument("--gb-mask", required=True, help=".npy bool mask from gb_identify.py")
    p.add_argument("--potential", required=True, help="Al-Mg EAM/fs file")
    p.add_argument("--n-gb", type=int, default=50, help="number of GB sites")
    p.add_argument("--n-bulk", type=int, default=5, help="number of bulk ref sites")
    p.add_argument("--seed", type=int, default=0, help="sample_seed for site selection")
    p.add_argument("--solute-type", type=int, default=2, help="LAMMPS type for Mg")
    p.add_argument("--bulk-separation", type=float, default=8.0,
                   help="min Å from any GB atom for a bulk reference")
    p.add_argument("--etol", type=float, default=1.0e-8)
    p.add_argument("--ftol", type=float, default=1.0e-10)
    p.add_argument("--maxiter", type=int, default=5000)
    p.add_argument("--maxeval", type=int, default=50000)
    p.add_argument("--lmp", default="lmp", help="LAMMPS binary")
    p.add_argument("--mpi-ranks", type=int, default=1,
                   help="MPI ranks per LAMMPS invocation (0 = no launcher)")
    p.add_argument("--mpi-cmd", default=None,
                   help="launcher binary; auto-detect srun vs mpirun if omitted")
    p.add_argument("--work-dir", default=None,
                   help="per-site deck/log dir (default: <annealed dir>/delta_e_run)")
    p.add_argument("--keep-per-site-files", action="store_true",
                   help="keep individual site decks + logs after run")
    p.add_argument("--out-npz", required=True, help="output .npz path")
    p.add_argument("--out-json", default=None, help="output metadata .json")
    args = p.parse_args()

    # Auto-select launcher: srun if inside SLURM, else mpirun.
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

    result = compute_delta_e_spectrum(
        annealed_data_file=args.annealed,
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
    )

    _save_results(result, Path(args.out_npz), Path(args.out_json) if args.out_json else None)

    dE = result["gb_delta_e"]
    all_stops = np.concatenate([result["gb_stop_reasons"], result["bulk_stop_reasons"]])
    bad = [r for r in all_stops if r not in ("energy tolerance", "force tolerance")]
    print(
        "\n"
        f"n_gb_sites       = {len(dE)}\n"
        f"n_sites_resumed  = {result['n_sites_resumed']}/{len(all_stops)}\n"
        f"CG stop reasons  = {dict((r, int((all_stops==r).sum())) for r in set(all_stops.tolist()))}\n"
        + (f"  ⚠ {len(bad)} site(s) did NOT reach energy/force tolerance\n" if bad else "")
        + f"bulk_e_mean      = {result['bulk_e_mean']:.4f} eV   "
        f"(std {result['bulk_e_std']:.4f}, n={len(result['bulk_e_mg'])})\n"
        f"ΔE_seg [eV]      min={dE.min():+.4f}  max={dE.max():+.4f}  "
        f"mean={dE.mean():+.4f}  median={np.median(dE):+.4f}\n"
        f"ΔE_seg [kJ/mol]  min={dE.min()*96.485:+.2f}  max={dE.max()*96.485:+.2f}  "
        f"mean={dE.mean()*96.485:+.2f}\n"
        f"wall_time_total  = {result['wall_seconds_total']:.1f} s\n",
        flush=True,
    )


if __name__ == "__main__":
    _cli()
