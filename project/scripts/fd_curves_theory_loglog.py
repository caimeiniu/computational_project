#!/usr/bin/env python3
"""Theory-only reference figure for the master-figure layout: four
predicted curves vs X_c at T=500 K on log-log axes, NO HMC points.

The four curves:
  - GC-FD: textbook Fermi-Dirac (Wagih eq. 2 with bulk reservoir at X_c).
    Diverges at large X_c because the formula assumes an infinite bulk
    reservoir to back-fill GB. Educationally important: this is what
    Wagih's framework strictly predicts in the thermodynamic limit.
  - canon-FD: mass-conserving form. Solves self-consistently for the
    post-segregation bulk composition X_b given a closed simulation
    box. Bound below by no-segregation diagonal, above by ceiling.
  - ceiling: closed-box maximum X_GB = min(1, X_c · N_total / N_GB).
    Hit when ALL Mg sits in GB and bulk is pure Al; an absolute upper
    bound regardless of T.
  - diagonal: X_GB = X_c (no segregation; thermal scrambling at high T).

Interpretation: GC-FD ≈ canon-FD at low X_c (dilute limit, reservoir
unchanged). At X_c above ~0.05 they diverge because the closed-box
runs out of bulk Mg to feed GB. canon-FD asymptotes to ceiling once X_c
saturates the GB site reservoir; GC-FD continues climbing toward 1.
The reader needs to see this divergence to understand why the
HMC/measurement comparison in the headline figure uses canon-FD as the
target, not GC-FD.

Output:
    output/fd_curves_theory_T500_loglog.{png,json}
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fermi_dirac_predict import (load_ours, load_wagih, x_gb_curve,
                                  x_gb_canonical_curve)

OUT = Path("/cluster/home/cainiu/Computational_modeling/project/output")
SPECTRUM_NPZ = Path("/cluster/scratch/cainiu/production_AlMg_200A/"
                    "delta_e_results_n500_200A_tight.npz")
GB_MASK_NPY = Path("/cluster/scratch/cainiu/production_AlMg_200A/"
                   "gb_mask_200A.npy")
WAGIH_SEG_TXT = Path("/cluster/scratch/cainiu/wagih_zenodo/"
                     "learning_segregation_energies/machine_learning_notebook/"
                     "seg_energies_Al_Mg.txt")
WAGIH_BULK_DAT = Path("/cluster/scratch/cainiu/wagih_zenodo/"
                      "learning_segregation_energies/machine_learning_notebook/"
                      "bulk_solute_Al_Mg.dat")
N_TOTAL_WAGIH = 483_425
N_GB_WAGIH    =  82_646


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--temp", type=float, default=500.0)
    ap.add_argument("--xc-min", type=float, default=1e-3)
    ap.add_argument("--xc-max", type=float, default=0.5)
    ap.add_argument("--out-prefix", default="fd_curves_theory_T500_loglog")
    ap.add_argument("--no-wagih", action="store_true",
                    help="omit the Wagih spectrum overlay (use only ours).")
    args = ap.parse_args()

    gb_mask = np.load(GB_MASK_NPY).astype(bool)
    N_TOTAL = int(gb_mask.size)
    N_GB = int(gb_mask.sum())
    GB_FRAC = N_GB / N_TOTAL

    dE = load_ours(SPECTRUM_NPZ)
    dE_wagih = (None if args.no_wagih
                else load_wagih(WAGIH_SEG_TXT, WAGIH_BULK_DAT))

    xc = np.geomspace(args.xc_min, args.xc_max, 200)
    T = args.temp
    gc_curve         = x_gb_curve(dE, T, xc)
    canon_ours, _    = x_gb_canonical_curve(dE, T, xc, N_GB, N_TOTAL)
    if dE_wagih is not None:
        canon_wagih, _ = x_gb_canonical_curve(dE_wagih, T, xc,
                                              N_GB_WAGIH, N_TOTAL_WAGIH)
    else:
        canon_wagih = None
    ceiling = np.minimum(1.0, xc * N_TOTAL / N_GB)

    # --- figure ---
    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    ax.plot(xc, gc_curve, "-",  color="C0", lw=2.0,
            label=r"GC-FD  $X_{GB}^{\mathrm{FD}}$  (Wagih eq. 2)")
    ax.plot(xc, canon_ours, "-", color="C2", lw=2.0,
            label=r"canon-FD (ours, mass-conserving)")
    if canon_wagih is not None:
        ax.plot(xc, canon_wagih, "--", color="C2", lw=1.4, alpha=0.85,
                label=r"canon-FD (Wagih spectrum)")
    ax.plot(xc, ceiling, ":", color="k", lw=1.4,
            label=r"closed-box ceiling  $X_c\cdot N_{\mathrm{tot}}/N_{\mathrm{GB}}$")
    ax.plot(xc, xc, "--", color="0.5", lw=1.0,
            label=r"diagonal  $X_{GB}=X_c$  (no segregation)")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(args.xc_min, args.xc_max)
    ax.set_ylim(args.xc_min, 1.05)
    ax.set_xlabel(r"bulk Mg fraction  $X_c$")
    ax.set_ylabel(r"$X_{GB}$")
    ax.set_title(f"Al(Mg) GB occupancy: theory predictions at "
                 fr"T={T:g} K  ($N_{{\mathrm{{GB}}}}/N_{{\mathrm{{tot}}}}"
                 fr"={GB_FRAC:.3f}$)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.95)

    fig.tight_layout()
    out_png = OUT / f"{args.out_prefix}.png"
    fig.savefig(out_png, dpi=140)
    print(f"wrote {out_png}")

    # --- json (for traceability / re-plot) ---
    summary = {
        "T": T, "xc": xc.tolist(),
        "geometry_ours":  {"N_total": N_TOTAL, "N_GB": N_GB,
                           "GB_frac": GB_FRAC,
                           "spectrum_n": int(dE.size),
                           "spectrum_mean_kjmol": float(dE.mean()*96.485)},
        "geometry_wagih": ({"N_total": N_TOTAL_WAGIH, "N_GB": N_GB_WAGIH,
                            "spectrum_n": int(dE_wagih.size),
                            "spectrum_mean_kjmol":
                                float(dE_wagih.mean()*96.485)}
                           if dE_wagih is not None else None),
        "GC_FD":           gc_curve.tolist(),
        "canon_FD_ours":   canon_ours.tolist(),
        "canon_FD_wagih":  (canon_wagih.tolist()
                            if canon_wagih is not None else None),
        "ceiling":         ceiling.tolist(),
        "diagonal":        xc.tolist(),
    }
    out_json = OUT / f"{args.out_prefix}.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"wrote {out_json}")


if __name__ == "__main__":
    main()
