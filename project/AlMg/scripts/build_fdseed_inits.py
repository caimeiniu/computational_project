#!/usr/bin/env python3
"""Build Fermi-Dirac-seeded initial conditions for the Hybrid-Monte-Carlo
(HMC) sweep at multiple total Mg fractions X_c.

Each output places  round(X_GB^FD(T, X_c) * N_GB)  Mg atoms at random
grain-boundary (GB) sites and the remaining Mg at random bulk sites, so
the HMC trajectory STARTS exactly at the canonical Fermi-Dirac prediction.

Why this initial condition (rather than the standard "fill GB first, spill
to bulk" preseg used by pre_segregate.py with no --xgb-init):

  - Standard preseg sets X_GB(0) ~ X_c / GB_fraction, which for X_c = 0.20
    or 0.30 saturates the GB completely (X_GB(0) = 1.0). The HMC trace
    must descend a long way (~0.3 in X_GB units) before reaching the FD
    prediction, and a 24-h job cannot finish that descent. The resulting
    "X_GB(end of run) is an upper bound on the equilibrium X_GB" still
    sits well above the FD curve — it does NOT prove the dilute-limit
    Fermi-Dirac prediction over-estimates the true equilibrium.

  - FD-seeded init sets X_GB(0) = X_GB^FD exactly. Any net descent in the
    HMC trace is direct evidence that the equilibrium X_GB is BELOW the
    FD prediction (the central project claim). Even a partially
    equilibrated trace gives a meaningful one-sided bound on
    X_GB(equilibrium) - X_GB^FD < 0.

Outputs:
  - poly_AlMg_200A_fdseed_T{T}K_Xc{xc}.lmp     (one per X_c)
  - fdseed_T{T}K_manifest.json                  (list of files + values)

Usage:
  python3 scripts/build_fdseed_inits.py \
      --in-data /cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g_annealed.lmp \
      --gb-mask data/snapshots/gb_mask_200A.npy \
      --delta-e-npz data/snapshots/delta_e_results_n500_200A_tight.npz \
      --out-dir data/snapshots \
      --T 500.0 \
      --xc-list 0.10,0.15,0.20,0.30
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

# Reuse the project's canonical (closed-box) Fermi-Dirac solver.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fermi_dirac_predict import x_gb_canonical


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--in-data", type=Path, required=True,
                   help="annealed pure-Al data file (LAMMPS atom_style atomic, "
                        "type 1 only)")
    p.add_argument("--gb-mask", type=Path, required=True,
                   help="boolean GB mask (.npy), True = GB site")
    p.add_argument("--delta-e-npz", type=Path, required=True,
                   help="segregation energy NPZ from sample_delta_e.py "
                        "(used to compute X_GB^FD)")
    p.add_argument("--T", type=float, default=500.0,
                   help="temperature in K (default 500)")
    p.add_argument("--xc-list", default="0.10,0.15,0.20,0.30",
                   help="comma-separated total Mg fractions")
    p.add_argument("--out-dir", type=Path, required=True,
                   help="directory to write the .lmp files and manifest")
    p.add_argument("--seed", type=int, default=20260505,
                   help="RNG seed for site selection (deterministic)")
    p.add_argument("--pre-segregate-py", type=Path,
                   default=Path(__file__).resolve().parent / "pre_segregate.py",
                   help="path to pre_segregate.py")
    args = p.parse_args()

    # Validate inputs early so we fail fast.
    for path in (args.in_data, args.gb_mask, args.delta_e_npz,
                 args.pre_segregate_py):
        if not path.exists():
            sys.exit(f"input not found: {path}")

    mask = np.load(args.gb_mask).astype(bool)
    n_total = int(mask.size)
    n_gb = int(mask.sum())
    delta_e = np.asarray(np.load(args.delta_e_npz)["gb_delta_e"])

    print(f"box:  N_total = {n_total}   N_GB = {n_gb}   "
          f"GB_fraction = {n_gb / n_total:.4f}")
    print(f"spectrum:  n = {delta_e.size}   "
          f"mean = {delta_e.mean()*96.485:+.2f} kJ/mol   "
          f"frac(ΔE<0) = {(delta_e < 0).mean():.3f}")
    print(f"T = {args.T:g} K   seed = {args.seed}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "T_K": args.T,
        "n_total": n_total,
        "n_gb": n_gb,
        "gb_fraction": n_gb / n_total,
        "delta_e_npz": str(args.delta_e_npz),
        "in_data": str(args.in_data),
        "seed": args.seed,
        "entries": [],
    }

    for xc_str in args.xc_list.split(","):
        xc = float(xc_str)
        x_gb_fd, x_bulk_fd = x_gb_canonical(delta_e, args.T, xc, n_gb, n_total)
        out_lmp = args.out_dir / (
            f"poly_AlMg_200A_fdseed_T{args.T:g}K_Xc{xc:g}.lmp"
        )

        cmd = [
            sys.executable, str(args.pre_segregate_py),
            "--in-data", str(args.in_data),
            "--gb-mask", str(args.gb_mask),
            "--xc", str(xc),
            "--xgb-init", f"{x_gb_fd:.6f}",
            "--seed", str(args.seed),
            "--out-data", str(out_lmp),
        ]
        print(f"\n[X_c = {xc:g}]  target X_GB(0) = X_GB^FD = {x_gb_fd:.4f}   "
              f"target X_bulk(0) = X_bulk^FD = {x_bulk_fd:.4f}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.stdout:
            print(result.stdout.rstrip())
        if result.returncode != 0:
            print("--- pre_segregate.py STDERR ---", file=sys.stderr)
            print(result.stderr, file=sys.stderr)
            sys.exit(f"pre_segregate.py failed for X_c = {xc}")

        manifest["entries"].append({
            "X_c": xc,
            "X_GB_FD": x_gb_fd,
            "X_bulk_FD": x_bulk_fd,
            "out_file": str(out_lmp),
        })

    manifest_path = args.out_dir / f"fdseed_T{args.T:g}K_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"\nWrote {len(manifest['entries'])} initial conditions.")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
