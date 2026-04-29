#!/usr/bin/env python3
"""Headline figure for the (T=500 K) row of the master report figure.

Plots HMC X_GB measurements against four reference curves and a shaded
canon-FD reference band:
  - canon-FD (ours)  — mass-conserving Wagih eq. 2 on our n=500 spectrum
  - canon-FD (Wagih) — same eq. on Wagih's n=82,646 LAMMPS-truth spectrum
    The ours/Wagih band visualises *reference uncertainty* — at low X_c
    the band is wide (~0.013 in X_GB) because the deep-tail reference
    differs; at high X_c it collapses (~0.002).
  - GC-FD    — textbook Wagih eq. 2 (thermodynamic-limit target)
  - ceiling  = min(1, X_c · N_total / N_GB)

HMC measurements are drawn in three classes (panel-d markers per
2026-04-28 evening plan):
  ● equilibrated bracket          — solid red circle + CI bar
  ▽ preseg upper bound (descending) — open red down-triangle
  □ random-IC kinetic floor (legacy) — open red square

Output (linear-x, linear-y, X_c ∈ [0, 0.4]):
  output/panel_d_draft_T500.png   panel (d) draft
  output/panel_d_draft_T500.json  per-point HMC + canon-FD + gap
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
from fermi_dirac_predict import (load_ours, load_wagih, x_gb, x_gb_curve,
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
# Wagih's polycrystal counts (his 20 nm box, 16 grains; 2026-04-28 audit).
N_TOTAL_WAGIH = 483_425
N_GB_WAGIH    =  82_646


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
        help="JSONs of equilibrated HMC measurements (drawn solid ●).")
    ap.add_argument("--upper-bound", nargs="*", default=[
        "hmc_T500_Xc0.10_preseg.json",
        "hmc_T500_Xc0.15_preseg.json",
        "hmc_T500_Xc0.20_preseg.json",
        "hmc_T500_Xc0.30_preseg.json"],
        help="JSONs of preseg-IC trajectories that have not equilibrated; "
             "use the trajectory END VALUE as a one-sided upper bound (▽).")
    ap.add_argument("--kinetic-floor", nargs="*", default=[],
        help="(legacy) JSONs of random-IC kinetic-floor measurements (□).")
    ap.add_argument("--out-prefix", default="hmc_vs_canonfd_T500")
    ap.add_argument("--temp", type=float, default=500.0)
    ap.add_argument("--xc-max", type=float, default=0.40,
                    help="x-axis upper bound (linear scale)")
    ap.add_argument("--no-wagih", action="store_true",
                    help="skip the canon-Wagih reference band overlay.")
    args = ap.parse_args()

    eq_paths = [OUT / f for f in args.equilibrated]
    ub_paths = [OUT / f for f in args.upper_bound]
    kf_paths = [OUT / f for f in args.kinetic_floor]

    # ---- box geometry + spectrum ----
    gb_mask = np.load(GB_MASK_NPY).astype(bool)
    N_TOTAL = int(gb_mask.size)
    N_GB = int(gb_mask.sum())
    GB_FRAC = N_GB / N_TOTAL
    print(f"ours box: N_total={N_TOTAL}, N_GB={N_GB}, GB_frac={GB_FRAC:.4f}")

    dE = load_ours(SPECTRUM_NPZ)
    print(f"ours spectrum: n={dE.size}, mean={dE.mean()*96.485:+.2f} kJ/mol")

    dE_wagih = None
    if not args.no_wagih:
        try:
            dE_wagih = load_wagih(WAGIH_SEG_TXT, WAGIH_BULK_DAT)
            print(f"Wagih spectrum: n={dE_wagih.size}, "
                  f"mean={dE_wagih.mean()*96.485:+.2f} kJ/mol  "
                  f"(box: N_total={N_TOTAL_WAGIH}, N_GB={N_GB_WAGIH})")
        except FileNotFoundError as e:
            print(f"  (Wagih spectrum unavailable, skipping band: {e})")

    # ---- reference curves ----
    T = args.temp
    xc_grid = np.linspace(1e-4, args.xc_max, 400)
    x_gb_gc = x_gb_curve(dE, T, xc_grid)
    x_gb_canon_ours, _ = x_gb_canonical_curve(dE, T, xc_grid,
                                               N_GB, N_TOTAL)
    if dE_wagih is not None:
        x_gb_canon_wagih, _ = x_gb_canonical_curve(
            dE_wagih, T, xc_grid, N_GB_WAGIH, N_TOTAL_WAGIH)
    else:
        x_gb_canon_wagih = None
    ceiling = np.minimum(1.0, xc_grid * N_TOTAL / N_GB)

    # ---- HMC points ----
    eq_pts = [load_hmc_pt(p) for p in eq_paths]
    ub_pts = [load_hmc_pt(p) for p in ub_paths]
    kf_pts = [load_hmc_pt(p) for p in kf_paths]

    # For upper-bound markers we report the trajectory LAST value (one-sided
    # upper bound on the equilibrium), not the mean of the descending series.
    for pt in ub_pts:
        pt["X_HMC_mean_descending"] = pt["X_HMC"]
        pt["X_HMC"] = pt["last"]
        pt["lo"] = pt["last"]   # one-sided: no lower CI from this single trajectory
        pt["hi"] = pt["last"]

    # canonical FD predictions at each point's X_c (using ours spectrum)
    for pt in eq_pts + ub_pts + kf_pts:
        pt["X_FD_canon_ours"] = x_gb_canonical_curve(
            dE, T, np.array([pt["Xc"]]), N_GB, N_TOTAL)[0][0]
        if dE_wagih is not None:
            pt["X_FD_canon_wagih"] = x_gb_canonical_curve(
                dE_wagih, T, np.array([pt["Xc"]]), N_GB_WAGIH, N_TOTAL_WAGIH)[0][0]
        else:
            pt["X_FD_canon_wagih"] = None
        pt["X_FD_GC"] = x_gb(dE, T, pt["Xc"])
        pt["ceiling"] = min(1.0, pt["Xc"] * N_TOTAL / N_GB)
        pt["gap_HMC_minus_canon_ours"] = pt["X_HMC"] - pt["X_FD_canon_ours"]
        pt["gap_HMC_minus_canon_wagih"] = (
            pt["X_HMC"] - pt["X_FD_canon_wagih"]
            if pt["X_FD_canon_wagih"] is not None else None)

    # ---- JSON dump for traceability ----
    summary = {
        "T": T,
        "geometry_ours":  {"N_total": N_TOTAL, "N_GB": N_GB,
                           "GB_frac": GB_FRAC,
                           "spectrum_n": int(dE.size),
                           "spectrum_mean_kjmol": float(dE.mean() * 96.485)},
        "geometry_wagih": ({"N_total": N_TOTAL_WAGIH, "N_GB": N_GB_WAGIH,
                            "spectrum_n": int(dE_wagih.size),
                            "spectrum_mean_kjmol":
                                float(dE_wagih.mean() * 96.485)}
                           if dE_wagih is not None else None),
        "curves": {
            "X_c": xc_grid.tolist(),
            "X_GB_canon_ours":  x_gb_canon_ours.tolist(),
            "X_GB_canon_wagih": (x_gb_canon_wagih.tolist()
                                 if x_gb_canon_wagih is not None else None),
            "X_GB_GC":          x_gb_gc.tolist(),
            "X_GB_ceiling":     ceiling.tolist(),
        },
        "equilibrated": eq_pts,
        "upper_bound":  ub_pts,
        "kinetic_floor": kf_pts,
    }
    out_json = OUT / f"{args.out_prefix}.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"wrote {out_json}")

    # ---- plot ----
    fig, ax = plt.subplots(figsize=(7.8, 5.6))

    # Reference band: shaded between canon-ours and canon-Wagih.
    if x_gb_canon_wagih is not None:
        lo_band = np.minimum(x_gb_canon_ours, x_gb_canon_wagih)
        hi_band = np.maximum(x_gb_canon_ours, x_gb_canon_wagih)
        ax.fill_between(xc_grid, lo_band, hi_band, color="C2", alpha=0.18,
                        zorder=1, label=r"canon-FD reference band "
                                         r"(ours $\leftrightarrow$ Wagih)")
        ax.plot(xc_grid, x_gb_canon_wagih, "-", color="C2", lw=1.2,
                alpha=0.9, label=r"canon-FD (Wagih, $n=82{,}646$)")
    ax.plot(xc_grid, x_gb_canon_ours, "-", color="C2", lw=2,
            label=r"canon-FD (ours, $n=500$)")
    ax.plot(xc_grid, x_gb_gc, "-", color="C0", lw=1.5,
            label=r"GC-FD  $X_{\mathrm{GB}}^{\mathrm{FD}}$")
    ax.plot(xc_grid, ceiling, ":", color="k", lw=1.0,
            label=r"closed-box ceiling")
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
                    label=r"$X_{\mathrm{GB}}^{\mathrm{HMC}}$ "
                          r"equilibrated bracket")

    # upper-bound HMC: open red down-triangles, descending arrows
    if ub_pts:
        x_ub = np.array([p["Xc"] for p in ub_pts])
        y_ub = np.array([p["X_HMC"] for p in ub_pts])  # = trajectory last
        ax.plot(x_ub, y_ub, "v", color="C3", mfc="white", mew=1.8, ms=11,
                zorder=5,
                label=r"$X_{\mathrm{GB}}^{\mathrm{HMC}}$ "
                      r"preseg upper bound (still descending)")
        # descent arrows from upper-bound to its canon-ours target — visual hint
        for p in ub_pts:
            ax.annotate("", xy=(p["Xc"], p["X_FD_canon_ours"]),
                        xytext=(p["Xc"], p["X_HMC"]),
                        arrowprops=dict(arrowstyle="->", color="C3",
                                        alpha=0.4, lw=1.2),
                        zorder=4)

    # kinetic-floor: open red squares (non-equilibrium random-IC, legacy)
    if kf_pts:
        x_kf = np.array([p["Xc"] for p in kf_pts])
        y_kf = np.array([p["X_HMC"] for p in kf_pts])
        lo_kf = np.array([p["lo"] for p in kf_pts])
        hi_kf = np.array([p["hi"] for p in kf_pts])
        ax.errorbar(x_kf, y_kf,
                    yerr=[y_kf - lo_kf, hi_kf - y_kf],
                    fmt="s", mfc="white", color="C3", ms=9, capsize=4,
                    mew=1.6, zorder=5,
                    label=r"$X_{\mathrm{GB}}^{\mathrm{HMC}}$ "
                          r"random-IC kinetic floor")

    ax.set_xlabel(r"total Mg fraction  $X_c$")
    ax.set_ylabel(r"$X_{\mathrm{GB}}$")
    ax.set_title(f"Al(Mg) GB occupancy at T={T:g} K  "
                 fr"($N_\mathrm{{GB}}/N_\mathrm{{tot}} = {GB_FRAC:.3f}$)")
    ax.set_xlim(0, args.xc_max)
    ax.set_ylim(0, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8.0, framealpha=0.95)

    fig.tight_layout()
    out_png = OUT / f"{args.out_prefix}.png"
    fig.savefig(out_png, dpi=140)
    print(f"wrote {out_png}")

    # ---- console table ----
    print()
    hdr = (f"{'X_c':>7} {'X_HMC':>9} {'canon-O':>9} {'canon-W':>9} "
           f"{'gap-O':>9} {'gap-W':>9}  source")
    print(hdr)
    print("-" * len(hdr))
    def _row(pts, tag):
        if not pts:
            return
        print(f"({tag})")
        for pt in pts:
            cw = pt['X_FD_canon_wagih']
            gw = pt['gap_HMC_minus_canon_wagih']
            print(f"{pt['Xc']:>7.3f} {pt['X_HMC']:>9.4f} "
                  f"{pt['X_FD_canon_ours']:>9.4f} "
                  f"{(cw if cw is not None else float('nan')):>9.4f} "
                  f"{pt['gap_HMC_minus_canon_ours']:>+9.4f} "
                  f"{(gw if gw is not None else float('nan')):>+9.4f}  "
                  f"{Path(pt['source']).name}")
    _row(eq_pts, "equilibrated")
    _row(ub_pts, "preseg upper bound")
    _row(kf_pts, "kinetic floor / random-IC")


if __name__ == "__main__":
    main()
