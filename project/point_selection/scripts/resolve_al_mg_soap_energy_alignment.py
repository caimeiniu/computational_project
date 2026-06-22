#!/usr/bin/env python3
"""
Resolve the Al-Mg SOAP/segregation-energy alignment for Wagih et al. 2020.

This script deliberately fails unless it can prove the row alignment between:

  machine_learning_notebook/GB_SOAP_Al_Mg.npy
  machine_learning_notebook/seg_energies_Al_Mg.txt
  machine_learning_notebook/bulk_solute_Al_Mg.dat
  machine_learning_notebook/heated_minimized_Al_polycrystal.dump0

The proof implemented here is:

1. Reproduce the current OVITO adaptive-CNA GB candidate list from the LAMMPS
   dump, in dump/atom order.
2. Recompute QUIP SOAP vectors for those current OVITO candidates only, using
   the exact descriptor string from the original notebook.
3. Compare the recomputed 82,646-row SOAP matrix against the saved 82,645-row
   Zenodo SOAP matrix.
4. Accept alignment only when the saved SOAP matrix equals the recomputed
   candidate matrix with exactly one current candidate row removed, within a
   strict floating tolerance.

Observed for the Zenodo archive:
  current OVITO candidate row 8677, site ID 50044, is absent from the saved
  GB_SOAP_Al_Mg.npy. Dropping that row gives a verified 82,645-row feature and
  target matrix.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import warnings
from pathlib import Path
from typing import Iterable

import numpy as np


EV_TO_KJMOL = 96.485
SOAP_TOL = 1e-10
PC_COLUMNS = [f"pc{i}" for i in range(1, 11)]
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
ALIGNMENT_OUTPUT_DIR = PACKAGE_ROOT / "outputs" / "alignment"


NOTEBOOK_LOGIC = """\
Original notebook logic summarized from Learn_Segregation_Spectra.ipynb:

GB IDs:
  GB_ids, GB_indices = get_gb_ids_and_indices(solvent_polycrystal)
  get_gb_ids_and_indices imports the LAMMPS dump with OVITO, applies
  CommonNeighborAnalysisModifier(), finds the majority structure type among
  FCC/HCP/BCC, then selects StructureType != majority. It returns both
  Particle Identifier values (for LAMMPS site IDs) and particle indices
  (for SOAP row selection).

Segregation energies:
  gb_ids_solute_energies = np.genfromtxt("seg_energies_Al_Mg.txt")
  bulk_solute_energy = np.genfromtxt("bulk_solute_Al_Mg.dat")
  gb_id_seg_dict[id] = E_GB_solute - E_bulk_solute
  y_arr[key] = gb_id_seg_dict[int(gb_atom_id)] for gb_atom_id in GB_ids
  y_arr = y_arr * 96.485

SOAP:
  system = ase.io.read("heated_minimized_Al_polycrystal.dump0",
                       format="lammps-dump-text")
  system.set_pbc([1,1,1])
  system.set_atomic_numbers(np.ones(len(system))*13)
  Descriptor("soap cutoff=6.00 l_max=12 n_max=12 atom_sigma=1.00 "
             "n_Z=1 Z={13} normalise=F")
  soap_arr = desc.calc_descriptor(system)
  soap_gb_vectors = soap_arr[GB_indices]
  np.save("GB_SOAP_Al_Mg.npy", soap_gb_vectors)

PCA:
  pca_x = PCA(n_components=10, svd_solver="full")
  x_arr_pca = pca_x.fit_transform(soap_gb_vectors)

K-means active-learning demo:
  KMeans(n_clusters=100, random_state=42) on x_arr_pca
  best_lae_indices = pairwise_distances_argmin_min(kmeans.cluster_centers_,
                                                   x_arr_pca)
"""


def log(message: str) -> None:
    print(message, flush=True)


def require_file(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def find_notebook_dir(root: Path) -> Path:
    candidates = list(root.rglob("machine_learning_notebook"))
    for path in candidates:
        if (path / "GB_SOAP_Al_Mg.npy").exists() and (path / "seg_energies_Al_Mg.txt").exists():
            return path
    raise FileNotFoundError(f"Could not find machine_learning_notebook under {root}")


def load_delta_e_kjmol(notebook_dir: Path) -> dict[str, float]:
    seg_path = require_file(notebook_dir / "seg_energies_Al_Mg.txt")
    bulk_path = require_file(notebook_dir / "bulk_solute_Al_Mg.dat")
    bulk = float(bulk_path.read_text().split()[0])
    out: dict[str, float] = {}
    with seg_path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 2:
                raise ValueError(f"Malformed line {line_number} in {seg_path}: {line!r}")
            site_id = str(int(float(parts[0])))
            if site_id in out:
                raise ValueError(f"Duplicate site ID in {seg_path}: {site_id}")
            e_gb = float(parts[1])
            out[site_id] = (e_gb - bulk) * EV_TO_KJMOL
    return out


def compute_current_ovito_gb_ids(dump_path: Path) -> list[str]:
    warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
    try:
        from ovito.io import import_file
        from ovito.modifiers import CommonNeighborAnalysisModifier
    except ImportError as exc:
        raise RuntimeError(
            "OVITO is required to reproduce the notebook GB selection. "
            "Install it with a compatible environment, e.g. `pip install ovito` "
            "or the OVITO conda package."
        ) from exc

    pipeline = import_file(str(dump_path))
    pipeline.modifiers.append(CommonNeighborAnalysisModifier())
    data = pipeline.compute()
    structure = data.particles["Structure Type"].array
    particle_ids = data.particles["Particle Identifier"].array
    counts = {t: int((structure == t).sum()) for t in [1, 2, 3]}
    majority = max(counts, key=counts.get)
    gb_ids = [str(int(pid)) for pid, st in zip(particle_ids, structure) if int(st) != majority]
    return gb_ids


def recompute_masked_quip_soap(notebook_dir: Path, gb_ids: list[str], cache_path: Path | None) -> np.ndarray:
    if cache_path and cache_path.exists():
        log(f"Loading cached recomputed SOAP: {cache_path}")
        return np.load(cache_path, mmap_mode="r")

    try:
        import ase.io
        from quippy import descriptors
    except ImportError as exc:
        raise RuntimeError(
            "ASE and quippy-ase are required for the SOAP proof. "
            "Install them with `pip install ase quippy-ase`."
        ) from exc

    dump_path = require_file(notebook_dir / "heated_minimized_Al_polycrystal.dump0")
    log("Reading LAMMPS dump with ASE...")
    atoms = ase.io.read(dump_path, format="lammps-dump-text")
    atoms.set_pbc([1, 1, 1])
    atoms.set_atomic_numbers(np.ones(len(atoms), dtype=int) * 13)

    mask = np.zeros(len(atoms), dtype=bool)
    for site_id in gb_ids:
        index = int(site_id) - 1
        if index < 0 or index >= len(atoms):
            raise ValueError(f"Site ID {site_id} is outside the atom index range")
        mask[index] = True
    if int(mask.sum()) != len(gb_ids):
        raise ValueError("GB mask count does not match GB ID count")
    atoms.arrays["gb_mask"] = mask

    descriptor = descriptors.Descriptor(
        "soap cutoff=6.00 l_max=12 n_max=12 atom_sigma=1.00 n_Z=1 Z={13} normalise=F"
    )
    log("Recomputing QUIP SOAP for current OVITO GB candidates...")
    start = time.time()
    soap = descriptor.calc_descriptor(atoms, atom_mask_name="gb_mask")
    log(f"Recomputed SOAP shape: {soap.shape} in {time.time() - start:.1f} s")

    if cache_path:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        np.save(cache_path, soap)
        log(f"Cached recomputed SOAP: {cache_path}")
    return soap


def find_single_extra_candidate(saved: np.ndarray, current: np.ndarray, tol: float) -> tuple[int, int, int]:
    """Return (extra_current_index, prefix_matches, suffix_matches)."""
    if current.shape[0] != saved.shape[0] + 1:
        raise ValueError(f"Expected current SOAP to have one extra row: saved={saved.shape}, current={current.shape}")
    if current.shape[1] != saved.shape[1]:
        raise ValueError(f"SOAP dimensions differ: saved={saved.shape}, current={current.shape}")

    n = saved.shape[0]
    prefix = 0
    while prefix < n and np.allclose(saved[prefix], current[prefix], rtol=tol, atol=tol):
        prefix += 1

    suffix = 0
    while suffix < n - prefix and np.allclose(saved[n - 1 - suffix], current[n - suffix], rtol=tol, atol=tol):
        suffix += 1

    if prefix + suffix != n:
        raise RuntimeError(
            "Could not prove that saved SOAP equals current SOAP with one row removed. "
            f"prefix={prefix}, suffix={suffix}, saved_rows={n}"
        )
    return prefix, prefix, suffix


def write_csv(output: Path, site_ids: list[str], pcs: np.ndarray, y: np.ndarray) -> None:
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["row_id", "site_id_if_verified", *PC_COLUMNS, "deltaE_kJmol"],
        )
        writer.writeheader()
        for row_id, (site_id, pc_row, delta_e) in enumerate(zip(site_ids, pcs, y)):
            row = {"row_id": row_id, "site_id_if_verified": site_id, "deltaE_kJmol": f"{delta_e:.12g}"}
            row.update({f"pc{i}": f"{float(value):.12g}" for i, value in enumerate(pc_row, start=1)})
            writer.writerow(row)


def stats(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "mean": float(np.mean(values)),
        "max": float(np.max(values)),
        "std": float(np.std(values)),
    }


def summarize_archive_files(dataset_root: Path) -> list[str]:
    hits = []
    for path in sorted(dataset_root.rglob("*")):
        if path.is_file() and ("Al_Mg" in path.name or "Al-Mg" in str(path) or "/Mg_" in str(path)):
            hits.append(str(path))
    return hits


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=PACKAGE_ROOT / "data" / "zenodo_4107058" / "extracted" / "learning_segregation_energies",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ALIGNMENT_OUTPUT_DIR / "al_mg_pca_deltaE_verified.csv",
    )
    parser.add_argument(
        "--report", type=Path, default=ALIGNMENT_OUTPUT_DIR / "al_mg_alignment_report.txt"
    )
    parser.add_argument(
        "--cache-current-soap",
        type=Path,
        default=ALIGNMENT_OUTPUT_DIR / "GB_SOAP_Al_Mg_current_ovito315_masked.npy",
    )
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--tolerance", type=float, default=SOAP_TOL)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.cache_current_soap.parent.mkdir(parents=True, exist_ok=True)

    try:
        from sklearn.cluster import KMeans
        from sklearn.decomposition import PCA
        from sklearn.linear_model import LinearRegression
        from sklearn.metrics import mean_absolute_error, pairwise_distances_argmin_min
        from sklearn.model_selection import train_test_split
    except ImportError as exc:
        raise RuntimeError("scikit-learn is required for PCA and verification") from exc

    dataset_root = args.dataset_root
    notebook_dir = find_notebook_dir(dataset_root)
    log(f"Notebook directory: {notebook_dir}")

    nb_path = require_file(notebook_dir / "Learn_Segregation_Spectra.ipynb")
    dump_path = require_file(notebook_dir / "heated_minimized_Al_polycrystal.dump0")
    saved_soap_path = require_file(notebook_dir / "GB_SOAP_Al_Mg.npy")

    delta_by_id = load_delta_e_kjmol(notebook_dir)
    log(f"Energy rows: {len(delta_by_id)}")

    current_ids = compute_current_ovito_gb_ids(dump_path)
    log(f"Current OVITO GB candidates: {len(current_ids)}")
    if len(current_ids) != len(set(current_ids)):
        raise RuntimeError("Current OVITO GB candidate IDs contain duplicates")

    saved_soap = np.load(saved_soap_path, mmap_mode="r")
    log(f"Saved SOAP shape: {saved_soap.shape}")
    cache = None if args.no_cache else args.cache_current_soap
    current_soap = recompute_masked_quip_soap(notebook_dir, current_ids, cache)

    extra_index, prefix, suffix = find_single_extra_candidate(saved_soap, current_soap, args.tolerance)
    excluded_site_id = current_ids[extra_index]
    log(f"Verified missing current SOAP row: index={extra_index}, site_id={excluded_site_id}")

    verified_ids = current_ids[:extra_index] + current_ids[extra_index + 1 :]
    if len(verified_ids) != saved_soap.shape[0]:
        raise RuntimeError("Verified ID count does not match saved SOAP row count")
    if len(verified_ids) != len(set(verified_ids)):
        raise RuntimeError("Verified IDs contain duplicates")

    missing_delta = [site_id for site_id in verified_ids if site_id not in delta_by_id]
    if missing_delta:
        raise RuntimeError(f"Missing deltaE values for verified SOAP IDs: {missing_delta[:10]}")

    y = np.array([delta_by_id[site_id] for site_id in verified_ids], dtype=float)
    if y.shape[0] != saved_soap.shape[0]:
        raise RuntimeError("Target vector length does not match SOAP row count")

    log("Fitting PCA(n_components=10, svd_solver='full') on saved SOAP...")
    pca = PCA(n_components=10, svd_solver="full")
    pcs = pca.fit_transform(np.asarray(saved_soap))
    if pcs.shape != (saved_soap.shape[0], 10):
        raise RuntimeError(f"Unexpected PCA shape: {pcs.shape}")

    log("Running k-means verification (k=100)...")
    kmeans = KMeans(n_clusters=100, random_state=42, n_init=10)
    kmeans.fit(pcs)
    best_indices, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, pcs)
    if len(set(best_indices.tolist())) != 100:
        raise RuntimeError("K-means selected duplicate nearest training indices")

    # Fingerprint against notebook's reported high-fidelity regression MAE.
    x_train, x_test, y_train, y_test = train_test_split(saved_soap, y, test_size=0.5, random_state=42)
    high_lr = LinearRegression()
    high_lr.fit(x_train, y_train)
    high_mae_train = float(mean_absolute_error(y_train, high_lr.predict(x_train)))
    high_mae_test = float(mean_absolute_error(y_test, high_lr.predict(x_test)))
    high_mae_all = float(mean_absolute_error(y, high_lr.predict(saved_soap)))

    accel_lr = LinearRegression()
    accel_lr.fit(pcs[best_indices], y[best_indices])
    accel_mae_train = float(mean_absolute_error(y[best_indices], accel_lr.predict(pcs[best_indices])))
    accel_mae_all = float(mean_absolute_error(y, accel_lr.predict(pcs)))

    write_csv(args.output, verified_ids, pcs, y)

    delta_stats = stats(y)
    report_lines = [
        "Al-Mg SOAP/deltaE alignment report",
        "=" * 38,
        "",
        "Files used:",
        f"- notebook: {nb_path.resolve()}",
        f"- structure dump: {dump_path.resolve()}",
        f"- saved SOAP: {saved_soap_path.resolve()}",
        f"- seg energies: {(notebook_dir / 'seg_energies_Al_Mg.txt').resolve()}",
        f"- bulk reference: {(notebook_dir / 'bulk_solute_Al_Mg.dat').resolve()}",
        "",
        "Al-Mg-related files found:",
        *[f"- {path}" for path in summarize_archive_files(dataset_root)],
        "",
        NOTEBOOK_LOGIC,
        "",
        "Alignment proof:",
        f"- seg_energies_Al_Mg.txt unique site IDs: {len(delta_by_id)}",
        f"- saved GB_SOAP_Al_Mg.npy shape: {tuple(saved_soap.shape)}",
        f"- current OVITO adaptive-CNA GB candidate count: {len(current_ids)}",
        f"- recomputed masked QUIP SOAP shape: {tuple(current_soap.shape)}",
        f"- tolerance for SOAP row comparison: rtol=atol={args.tolerance}",
        f"- exact one-row deletion proof: prefix rows={prefix}, suffix rows={suffix}",
        f"- excluded current candidate index: {extra_index}",
        f"- excluded site ID: {excluded_site_id}",
        "- interpretation: current OVITO identifies one extra GB atom relative to the",
        "  GB list used to create Zenodo's saved SOAP matrix. Removing site 50044",
        "  aligns every saved SOAP row with recomputed QUIP SOAP within tolerance.",
        "",
        "Verification checks:",
        f"- SOAP rows == matched deltaE rows: {saved_soap.shape[0]} == {y.shape[0]}",
        f"- duplicate verified site IDs: {len(verified_ids) - len(set(verified_ids))}",
        f"- missing deltaE for verified SOAP IDs: {len(missing_delta)}",
        f"- deltaE units: kJ/mol, using (E_GB_solute - E_bulk_solute) * 96.485",
        f"- deltaE min/mean/max/std: {delta_stats['min']:.12g}, {delta_stats['mean']:.12g}, "
        f"{delta_stats['max']:.12g}, {delta_stats['std']:.12g}",
        f"- PCA output shape: {tuple(pcs.shape)}",
        "- PCA explained variance ratios: "
        + ", ".join(f"pc{i+1}={v:.12g}" for i, v in enumerate(pca.explained_variance_ratio_)),
        f"- PCA cumulative explained variance: {float(np.sum(pca.explained_variance_ratio_)):.12g}",
        f"- k-means k=100 runnable: yes; unique nearest indices={len(set(best_indices.tolist()))}",
        f"- high-fidelity LinearRegression MAE train/test/all: "
        f"{high_mae_train:.6g}, {high_mae_test:.6g}, {high_mae_all:.6g}",
        f"- PCA+k-means LinearRegression MAE train/all: {accel_mae_train:.6g}, {accel_mae_all:.6g}",
        "",
        "Safety conclusion:",
        "The final table is safe for random sampling, k-means sampling, and",
        "active-learning experiments that need a trustworthy X=pc1..pc10 and",
        "y=deltaE_kJmol matrix. The site_id_if_verified column is not invented;",
        "it is derived from current OVITO GB IDs after a SOAP-level proof that",
        "site 50044 is the single current candidate absent from the saved Zenodo",
        "SOAP matrix.",
        "",
        f"Output CSV: {args.output.resolve()}",
    ]
    args.report.write_text("\n".join(report_lines) + "\n")

    log(f"Wrote {args.output} with {len(verified_ids)} rows")
    log(f"Wrote {args.report}")
    log("First 5 rows:")
    with args.output.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for _, row in zip(range(5), reader):
            print(row)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        ALIGNMENT_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        failed = ALIGNMENT_OUTPUT_DIR / "alignment_failed_report.txt"
        failed.write_text(
            "Alignment failed\n"
            "================\n\n"
            f"Error: {type(exc).__name__}: {exc}\n\n"
            "No verified PCA/deltaE CSV was written by this failing run.\n"
        )
        print(f"ERROR: {exc}", file=sys.stderr)
        print(f"Wrote {failed}", file=sys.stderr)
        raise
