#!/usr/bin/env python3
"""Headline figure for the (T=500 K) row of the master report figure.

Plots HMC X_GB measurements against three reference curves:
  - canon-FD (mass-conserving Wagih eq. 2; apples-to-apples target)
  - GC-FD    (textbook Wagih eq. 2; thermodynamic-limit target)
  - ceiling  = min(1, X_c · N_total / N_GB)

The HMC point set is configurable; default = the 5-point set planned in
the master figure memory (verify-preseg X_c=5e-2 + sweep X_c ∈
{0.10, 0.15, 0.20, 0.30}). Points loaded from `output/*_xgb.json` files
written by `scripts/hmc_xgb_timeseries.py`. Equilibrated vs
non-equilibrium (kinetic-floor) points may be passed in two separate
lists so the figure draws them with distinct markers.

Output (linear-x, linear-y, X_c ∈ [0, 0.4]):
  output/hmc_vs_fd_T500_5pt.png   panel (d) candidate
  output/hmc_vs_fd_T500_5pt.json  per-point HMC + canon-FD + gap
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fermi_dirac_predict import (load_ours, x_gb, x_gb_curve,
                                  x_gb_canonical_curve)

OUT = Path("/cluster/home/cainiu/Computational_modeling/project/output")
SPECTRUM_NPZ = Path("/cluster/scratch/cainiu/production_AlMg_200A/"
                    "delta_e_results_n500_200A_tight.npz")
GB_MASK_NPY = Path("/cluster/scratch/cainiu/production_AlMg_200A/"
                   "gb_mask_200A.npy")


def load_hmc_pt(json_path: Path) -> dict:
    d = json.load(open(json_path))
    return {
        "Xc": d["X_c"],
        "X_HMC": d["x_gb"]["mean"],
        "lo": d["x_gb"]["ci95_lo"],
        "hi": d["x_gb"]["ci95_hi"],
        "first": d["series"]["x_gb"][0],
        "last": d["series"]["x_gb"][-1],
        "n_frames": d["x_gb"]["n_frames_total"],
        "accepts": d["swap"]["total_accepts"],
        "attempts": d["swap"]["total_attempts"],
        "PE_drift_eV": d["thermo_prod"]["PE_drift_eV"],
        "source": str(json_path),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--equilibrated", nargs="*", default=[
        "hmc_T500_Xc5e-2_verify-preseg_xgb.json"],
        help="JSONs of equilibrated HMC measurements (drawn solid).")
    ap.add_argument("--kinetic-floor", nargs="*", default=[
        "hmc_T500_Xc0.10_xgb.json",
        "hmc_T500_Xc0.15_xgb.json",
        "hmc_T500_Xc0.20_xgb.json",
        "hmc_T500_Xc0.30_xgb_partial.json"],
        help="JSONs of non-equilibrium HMC measurements (drawn open).")
    ap.add_argument("--out-prefix", default="hmc_vs_fd_T500_5pt")
    ap.add_argument("--temp", type=float, default=500.0)
    ap.add_argument("--xc-max", type=float, default=0.40,
                    help="x-axis upper bound (linear scale)")
    args = ap.parse_args()

    eq_paths = [OUT / f for f in args.equilibrated]
    kf_paths = [OUT / f for f in args.kinetic_floor]

    # ---- box geometry + spectrum ----
    gb_mask = np.load(GB_MASK_NPY).astype(bool)
    N_TOTAL = int(gb_mask.size)
    N_GB = int(gb_mask.sum())
    GB_FRAC = N_GB / N_TOTAL
    print(f"box: N_total={N_TOTAL}, N_GB={N_GB}, GB_frac={GB_FRAC:.4f}")

    dE = load_ours(SPECTRUM_NPZ)
    print(f"spectrum: n={dE.size}, mean={dE.mean()*96.485:+.2f} kJ/mol")

    # ---- reference curves ----
    T = args.temp
    xc_grid = np.linspace(1e-4, args.xc_max, 400)
    x_gb_gc = x_gb_curve(dE, T, xc_grid)
    x_gb_canon, x_bulk_canon = x_gb_canonical_curve(dE, T, xc_grid,
                                                     N_GB, N_TOTAL)
    ceiling = np.minimum(1.0, xc_grid * N_TOTAL / N_GB)

    # ---- HMC points (equilibrated + kinetic floor) ----
    eq_pts = [load_hmc_pt(p) for p in eq_paths]
    kf_pts = [load_hmc_pt(p) for p in kf_paths]

    # canonical FD predictions at each point's X_c
    for pt in eq_pts + kf_pts:
        pt["X_FD_canon"] = x_gb_canonical_curve(
            dE, T, np.array([pt["Xc"]]), N_GB, N_TOTAL)[0][0]
        pt["X_FD_GC"] = x_gb(dE, T, pt["Xc"])
        pt["ceiling"] = min(1.0, pt["Xc"] * N_TOTAL / N_GB)
        pt["gap_HMC_minus_canon"] = pt["X_HMC"] - pt["X_FD_canon"]

    # ---- JSON dump for traceability ----
    summary = {
        "T": T,
        "geometry": {"N_total": N_TOTAL, "N_GB": N_GB, "GB_frac": GB_FRAC,
                     "spectrum_n": int(dE.size),
                     "spectrum_mean_kjmol": float(dE.mean() * 96.485)},
        "curves": {
            "X_c": xc_grid.tolist(),
            "X_GB_canon": x_gb_canon.tolist(),
            "X_GB_GC":    x_gb_gc.tolist(),
            "X_GB_ceiling": ceiling.tolist(),
        },
        "equilibrated": eq_pts,
        "kinetic_floor": kf_pts,
    }
    out_json = OUT / f"{args.out_prefix}.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"wrote {out_json}")

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(7.5, 5.4))

    ax.plot(xc_grid, x_gb_gc, "-", color="C0", lw=2,
            label=r"GC-FD  $X_{\mathrm{GB}}^{\mathrm{FD}}$ (Wagih eq. 2)")
    ax.plot(xc_grid, x_gb_canon, "-", color="C2", lw=2,
            label=r"canon-FD  $X_{\mathrm{GB}}^{\mathrm{FD}}$ (mass-conserving)")
    ax.plot(xc_grid, ceiling, ":", color="k", lw=1.0,
            label=r"closed-box ceiling $X_c \cdot N_\mathrm{tot}/N_\mathrm{GB}$")
    ax.plot(xc_grid, xc_grid, "--", color="gray", lw=0.9,
            label=r"$X_{\mathrm{GB}} = X_c$ (no segregation)")

    # equilibrated HMC: filled red circles + errorbars
    if eq_pts:
        x_eq = np.array([p["Xc"] for p in eq_pts])
        y_eq = np.array([p["X_HMC"] for p in eq_pts])
        lo_eq = np.array([p["lo"] for p in eq_pts])
        hi_eq = np.array([p["hi"] for p in eq_pts])
        ax.errorbar(x_eq, y_eq,
                    yerr=[y_eq - lo_eq, hi_eq - y_eq],
                    fmt="o", color="C3", ms=10, capsize=4, zorder=6,
                    label=r"$X_{\mathrm{GB}}^{\mathrm{HMC}}$ equilibrated")

    # kinetic-floor: open red squares (non-equilibrium random-IC)
    if kf_pts:
        x_kf = np.array([p["Xc"] for p in kf_pts])
        y_kf = np.array([p["X_HMC"] for p in kf_pts])
        lo_kf = np.array([p["lo"] for p in kf_pts])
        hi_kf = np.array([p["hi"] for p in kf_pts])
        ax.errorbar(x_kf, y_kf,
                    yerr=[y_kf - lo_kf, hi_kf - y_kf],
                    fmt="s", mfc="white", color="C3", ms=9, capsize=4,
                    mew=1.6, zorder=5,
                    label=r"$X_{\mathrm{GB}}^{\mathrm{HMC}}$ random-IC (kinetic floor)")

    ax.set_xlabel(r"total Mg fraction  $X_c$")
    ax.set_ylabel(r"$X_{\mathrm{GB}}$")
    ax.set_title(f"Al(Mg) GB occupancy at T={T:g} K  "
                 fr"($N_\mathrm{{GB}}/N_\mathrm{{tot}} = {GB_FRAC:.3f}$)")
    ax.set_xlim(0, args.xc_max)
    ax.set_ylim(0, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.5, framealpha=0.95)

    fig.tight_layout()
    out_png = OUT / f"{args.out_prefix}.png"
    fig.savefig(out_png, dpi=140)
    print(f"wrote {out_png}")

    # ---- console table ----
    print()
    hdr = f"{'X_c':>7} {'X_HMC':>9} {'GC-FD':>9} {'canon-FD':>9} {'ceiling':>9} {'gap':>9}  source"
    print(hdr)
    print("-" * len(hdr))
    print("(equilibrated)")
    for pt in eq_pts:
        print(f"{pt['Xc']:>7.3f} {pt['X_HMC']:>9.4f} {pt['X_FD_GC']:>9.4f} "
              f"{pt['X_FD_canon']:>9.4f} {pt['ceiling']:>9.4f} "
              f"{pt['gap_HMC_minus_canon']:>+9.4f}  "
              f"{Path(pt['source']).name}")
    print("(kinetic floor / non-equilibrium)")
    for pt in kf_pts:
        print(f"{pt['Xc']:>7.3f} {pt['X_HMC']:>9.4f} {pt['X_FD_GC']:>9.4f} "
              f"{pt['X_FD_canon']:>9.4f} {pt['ceiling']:>9.4f} "
              f"{pt['gap_HMC_minus_canon']:>+9.4f}  "
              f"{Path(pt['source']).name}")


if __name__ == "__main__":
    main()
