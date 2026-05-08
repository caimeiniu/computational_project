#!/usr/bin/env python3
"""T-axis subpanel: dilute-limit-breakdown robustness across temperature.

Companion to ``canonical_fd_compare_5pt.py`` which sweeps X_c at fixed T=500 K.
This script sweeps T at fixed X_c (default X_c=0.10, the central point of the
broken band) and produces a TWO-PANEL figure:

  Left  (sub-A, "equilibrium" view): canon-FD(T) curves overlaid with HMC
        measurements at each T. Two claims testable here:
          (1) does CI95 exclude FD at every T?           → T-robust break
          (2) is the gap monotonic / steeper at low T?    → enthalpic sig
  Right (sub-B, "kinetic" view): X_GB(t) trajectories starting from
        X_GB^FD seed for each T, plotted on a common time axis. Shape
        comparison sidesteps the kinetic caveat about post-burnin means
        being upper bounds at low T (CHANGELOG 2026-05-07 afternoon).

This is a copy + rename of canonical_fd_compare_5pt.py per the
``no-in-place-script-edits`` rule (feedback memory) — the canonical X_c-axis
script is preserved unchanged.

Usage::

    python3 scripts/canonical_fd_compare_t_axis.py \
        --jsons hmc_T300_Xc0.10_fdseed.json hmc_T500_Xc0.10_fdseed.json \
        --x-c 0.10 \
        --out-prefix panel_d_T_axis_X_c0.10_2pt_draft

Drop --jsons / --x-c when task C lands (T=700, X_c=0.10) to extend to 3 pts.
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
from fermi_dirac_predict import (load_ours, load_wagih, x_gb,
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


def load_hmc_pt(json_path: Path) -> dict:
    d = json.load(open(json_path))
    return {
        "T":         float(d["T"]),
        "Xc":        float(d["X_c"]),
        "X_HMC":     d["x_gb"]["mean"],
        "lo":        d["x_gb"]["ci95_lo"],
        "hi":        d["x_gb"]["ci95_hi"],
        "first":     d["series"]["x_gb"][0],
        "last":      d["series"]["x_gb"][-1],
        "timesteps": d["series"]["timestep"],
        "x_gb_t":    d["series"]["x_gb"],
        "n_frames":  d["x_gb"]["n_frames_total"],
        "imbalance_postburnin": (d.get("swap_decomposition", {})
                                  .get("post_burnin", {})
                                  .get("imbalance_signed")),
        "PE_drift_eV": d["thermo_prod"]["PE_drift_eV"],
        "source":    str(json_path),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jsons", nargs="+", required=True,
                    help="HMC fdseed JSONs (one per T, all at the same X_c).")
    ap.add_argument("--x-c", type=float, required=True,
                    help="Fixed X_c at which the T sweep was run.")
    ap.add_argument("--t-min", type=float, default=200.0,
                    help="lower T bound for FD reference curves (K)")
    ap.add_argument("--t-max", type=float, default=1000.0,
                    help="upper T bound for FD reference curves (K)")
    ap.add_argument("--out-prefix", default="panel_d_T_axis_draft")
    ap.add_argument("--no-wagih", action="store_true")
    ap.add_argument("--draft-banner", default="",
                    help="optional text overlaid in figure (e.g. "
                         "'draft: T=700 pending task C')")
    args = ap.parse_args()

    pts = [load_hmc_pt(OUT / f) for f in args.jsons]
    # sanity: every point at the requested X_c
    for p in pts:
        if abs(p["Xc"] - args.x_c) > 1e-6:
            raise SystemExit(f"X_c mismatch: {p['source']} has X_c={p['Xc']}, "
                             f"expected {args.x_c}")
    pts.sort(key=lambda p: p["T"])

    # ---- box geometry + spectrum ----
    gb_mask = np.load(GB_MASK_NPY).astype(bool)
    N_TOTAL = int(gb_mask.size)
    N_GB = int(gb_mask.sum())
    GB_FRAC = N_GB / N_TOTAL

    dE = load_ours(SPECTRUM_NPZ)
    dE_wagih = None
    if not args.no_wagih:
        try:
            dE_wagih = load_wagih(WAGIH_SEG_TXT, WAGIH_BULK_DAT)
        except FileNotFoundError as e:
            print(f"  (Wagih spectrum unavailable, skipping band: {e})")

    # ---- reference curves: X_GB vs T at fixed X_c ----
    T_grid = np.linspace(args.t_min, args.t_max, 80)
    Xc_arr = np.array([args.x_c])

    fd_canon_ours = np.empty_like(T_grid)
    fd_canon_wagih = (np.empty_like(T_grid) if dE_wagih is not None else None)
    fd_gc = np.empty_like(T_grid)
    for i, T in enumerate(T_grid):
        fd_canon_ours[i] = x_gb_canonical_curve(dE, T, Xc_arr, N_GB,
                                                 N_TOTAL)[0][0]
        if fd_canon_wagih is not None:
            fd_canon_wagih[i] = x_gb_canonical_curve(
                dE_wagih, T, Xc_arr, N_GB_WAGIH, N_TOTAL_WAGIH)[0][0]
        fd_gc[i] = x_gb(dE, T, args.x_c)
    ceiling_const = min(1.0, args.x_c * N_TOTAL / N_GB)

    # canonical FD prediction at each HMC point's T (for gap)
    for p in pts:
        p["X_FD_canon_ours"] = float(x_gb_canonical_curve(
            dE, p["T"], Xc_arr, N_GB, N_TOTAL)[0][0])
        if dE_wagih is not None:
            p["X_FD_canon_wagih"] = float(x_gb_canonical_curve(
                dE_wagih, p["T"], Xc_arr, N_GB_WAGIH, N_TOTAL_WAGIH)[0][0])
        else:
            p["X_FD_canon_wagih"] = None
        p["X_FD_GC"]      = float(x_gb(dE, p["T"], args.x_c))
        p["ceiling"]      = ceiling_const
        p["gap_HMC_minus_canon_ours"]  = p["X_HMC"] - p["X_FD_canon_ours"]
        p["gap_HMC_minus_canon_wagih"] = (
            p["X_HMC"] - p["X_FD_canon_wagih"]
            if p["X_FD_canon_wagih"] is not None else None)

    # ---- JSON dump ----
    summary = {
        "X_c_fixed":       args.x_c,
        "geometry_ours":  {"N_total": N_TOTAL, "N_GB": N_GB,
                           "GB_frac": GB_FRAC,
                           "spectrum_n": int(dE.size)},
        "curves_T_axis":  {
            "T_K":              T_grid.tolist(),
            "X_GB_canon_ours":  fd_canon_ours.tolist(),
            "X_GB_canon_wagih": (fd_canon_wagih.tolist()
                                 if fd_canon_wagih is not None else None),
            "X_GB_GC":          fd_gc.tolist(),
            "X_GB_ceiling":     [ceiling_const] * len(T_grid),
        },
        "hmc_points":     pts,
    }
    out_json = OUT / f"{args.out_prefix}.json"
    out_json.write_text(json.dumps(summary, indent=2, default=str))
    print(f"wrote {out_json}")

    # ---- plot ----
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.4))

    # ---- LEFT panel: equilibrium-view T-axis ----
    if fd_canon_wagih is not None:
        lo_band = np.minimum(fd_canon_ours, fd_canon_wagih)
        hi_band = np.maximum(fd_canon_ours, fd_canon_wagih)
        axL.fill_between(T_grid, lo_band, hi_band, color="C2", alpha=0.18,
                         label="canon-FD reference band")
        axL.plot(T_grid, fd_canon_wagih, "-", color="C2", lw=1.2, alpha=0.85,
                 label=r"canon-FD (Wagih, $n=82{,}646$)")
    axL.plot(T_grid, fd_canon_ours, "-", color="C2", lw=2,
             label=r"canon-FD (ours, $n=500$)")
    axL.plot(T_grid, fd_gc, "-", color="C0", lw=1.5,
             label=r"GC-FD")
    axL.axhline(ceiling_const, ls=":", color="k", lw=1.0,
                label="closed-box ceiling")

    if pts:
        TT  = np.array([p["T"] for p in pts])
        YY  = np.array([p["X_HMC"] for p in pts])
        LO  = np.array([p["lo"] for p in pts])
        HI  = np.array([p["hi"] for p in pts])
        axL.errorbar(TT, YY, yerr=[YY - LO, HI - YY],
                     fmt="o", color="C3", ms=10, capsize=4, zorder=6,
                     label=r"$X_{\mathrm{GB}}^{\mathrm{HMC}}$ post-burnin")

    axL.set_xlabel("temperature  T  (K)")
    axL.set_ylabel(r"$X_{\mathrm{GB}}$")
    axL.set_title(f"Equilibrium view — $X_c = {args.x_c:g}$")
    axL.set_xlim(args.t_min, args.t_max)
    axL.set_ylim(0, max(ceiling_const + 0.05, float(YY.max() if pts else 0) + 0.05,
                        float(fd_canon_ours.max()) + 0.05))
    axL.grid(True, alpha=0.3)
    axL.legend(loc="upper right", fontsize=8.5, framealpha=0.95)

    # ---- RIGHT panel: kinetic-view trajectories ----
    cmap = plt.get_cmap("coolwarm")
    Ts_arr = np.array([p["T"] for p in pts]) if pts else np.array([])
    if Ts_arr.size:
        T_lo, T_hi = Ts_arr.min(), Ts_arr.max()
    else:
        T_lo = T_hi = 500.0
    def color_for_T(T):
        if T_hi == T_lo:
            return cmap(0.5)
        return cmap((T - T_lo) / (T_hi - T_lo))

    for p in pts:
        ts = np.asarray(p["timesteps"], dtype=float)
        t_ps = (ts - ts[0]) * 1e-3   # 1000 LAMMPS steps × dt=1 fs = 1 ps
        c = color_for_T(p["T"])
        axR.plot(t_ps, p["x_gb_t"], "-", color=c, lw=1.4,
                 label=fr"T = {p['T']:.0f} K  (HMC, fdseed)")
        axR.axhline(p["X_FD_canon_ours"], ls="--", color=c, lw=1.0, alpha=0.7,
                    label=fr"canon-FD T={p['T']:.0f} K = {p['X_FD_canon_ours']:.3f}")

    axR.set_xlabel("time after fdseed start  (ps)")
    axR.set_ylabel(r"$X_{\mathrm{GB}}(t)$")
    axR.set_title(f"Kinetic view — $X_c = {args.x_c:g}$")
    axR.grid(True, alpha=0.3)
    axR.legend(loc="upper right", fontsize=8.5, framealpha=0.95)

    if args.draft_banner:
        fig.suptitle(args.draft_banner, fontsize=11, color="C3", y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
    else:
        fig.tight_layout()
    out_png = OUT / f"{args.out_prefix}.png"
    fig.savefig(out_png, dpi=140)
    print(f"wrote {out_png}")

    # ---- console table ----
    print()
    hdr = (f"{'T (K)':>8} {'X_HMC':>9} {'CI95-lo':>9} {'CI95-hi':>9} "
           f"{'canon-O':>9} {'gap-O':>9} {'imbalance':>10}  source")
    print(hdr)
    print("-" * len(hdr))
    for p in pts:
        imb = p["imbalance_postburnin"]
        imb_s = f"{imb:>+10.3f}" if isinstance(imb, (int, float)) else f"{'n/a':>10s}"
        print(f"{p['T']:>8.0f} {p['X_HMC']:>9.4f} {p['lo']:>9.4f} "
              f"{p['hi']:>9.4f} {p['X_FD_canon_ours']:>9.4f} "
              f"{p['gap_HMC_minus_canon_ours']:>+9.4f} {imb_s}  "
              f"{Path(p['source']).name}")


if __name__ == "__main__":
    main()
