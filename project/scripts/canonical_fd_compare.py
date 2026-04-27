#!/usr/bin/env python3
"""Compare HMC X_GB to BOTH grand-canonical and canonical FD predictions.

Background: Wagih's FD eq. (P_i = 1/(1+(1-X_c)/X_c · exp(ΔE_i/kT))) treats
X_c as a *bulk* mole fraction held fixed by an infinite reservoir
(grand canonical). Our HMC is closed-box (canonical, total Mg fixed),
so the apples-to-apples comparison must self-consistently solve for the
post-segregation bulk fraction X_bulk under mass conservation:

    X_bulk · N_bulk + <P_i(T, X_bulk)> · N_GB = X_c_total · N_total

Today's three HMC points all sit at X_c_total ≤ 0.05 in our 200³ Å box
(N_GB / N_total = 0.187), where the grand-canonical FD prediction
exceeds X_c_total · N_total / N_GB and is therefore physically
unreachable in a closed box. Canonical FD is the legitimate target.

Outputs:
  output/canonical_fd_T500.json    canonical X_GB, X_bulk per X_c_total
  output/hmc_vs_fd_T500_canonical.png   3-curve plot: GC-FD, canon-FD, HMC
  output/hmc_vs_canonical_fd.csv   per-point comparison + gap
"""
from __future__ import annotations
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
GB_MASK_NPY = Path("/cluster/scratch/cainiu/production_AlMg_200A/gb_mask_200A.npy")

# ---------------------- box geometry from the actual mask ------------------
gb_mask = np.load(GB_MASK_NPY).astype(bool)
N_TOTAL = int(gb_mask.size)
N_GB = int(gb_mask.sum())
GB_FRAC = N_GB / N_TOTAL
print(f"box: N_total={N_TOTAL}, N_GB={N_GB}, GB_frac={GB_FRAC:.4f}")

# ---------------------- load ΔE spectrum -----------------------------------
dE = load_ours(SPECTRUM_NPZ)  # eV, N=500 sample
print(f"ΔE: n={dE.size}, mean={dE.mean()*96.485:+.2f} kJ/mol, "
      f"frac<0={float((dE<0).mean()):.3f}")

# ---------------------- evaluation grid ------------------------------------
T = 500.0
xc_grid = np.logspace(-5, np.log10(0.5), 200)

x_gb_gc = x_gb_curve(dE, T, xc_grid)
x_gb_canon, x_bulk_canon = x_gb_canonical_curve(dE, T, xc_grid, N_GB, N_TOTAL)

# physical ceiling: all Mg piled into GB
ceiling = np.minimum(1.0, xc_grid * N_TOTAL / N_GB)

# ---------------------- HMC measurements -----------------------------------
hmc_files = [
    OUT / "hmc_T500_Xc5e-4_xgb.json",
    OUT / "hmc_T500_Xc5e-3_xgb.json",
    OUT / "hmc_T500_Xc5e-2_xgb.json",
]
hmc_pts = []
for p in hmc_files:
    d = json.load(open(p))
    hmc_pts.append({
        "Xc": d["X_c"],
        "X_HMC": d["x_gb"]["mean"],
        "lo": d["x_gb"]["ci95_lo"],
        "hi": d["x_gb"]["ci95_hi"],
        "first": d["series"]["x_gb"][0],
        "last": d["series"]["x_gb"][-1],
        "accepts": d["swap"]["total_accepts"],
    })

# canonical FD evaluated AT each HMC X_c
hmc_pred = []
for pt in hmc_pts:
    xc = pt["Xc"]
    gc = x_gb(dE, T, xc)
    canon_curve_one, _ = x_gb_canonical_curve(dE, T, np.array([xc]),
                                              N_GB, N_TOTAL)
    canon = float(canon_curve_one[0])
    pt["X_FD_GC"] = gc
    pt["X_FD_canon"] = canon
    pt["ceiling"] = min(1.0, xc * N_TOTAL / N_GB)
    pt["gap_HMC_minus_canon"] = pt["X_HMC"] - canon
    hmc_pred.append(pt)

# ---------------------- save canonical curve JSON --------------------------
canon_out = {
    "T": T,
    "geometry": {"N_total": N_TOTAL, "N_GB": N_GB, "GB_frac": GB_FRAC,
                 "spectrum_n": int(dE.size),
                 "spectrum_mean_kjmol": float(dE.mean() * 96.485)},
    "X_c_total": xc_grid.tolist(),
    "X_GB_grand_canonical": x_gb_gc.tolist(),
    "X_GB_canonical": x_gb_canon.tolist(),
    "X_bulk_canonical": x_bulk_canon.tolist(),
    "X_GB_ceiling": ceiling.tolist(),
}
(OUT / "canonical_fd_T500.json").write_text(json.dumps(canon_out, indent=2))
print(f"wrote {OUT/'canonical_fd_T500.json'}")

# ---------------------- per-point CSV --------------------------------------
csv_path = OUT / "hmc_vs_canonical_fd.csv"
with csv_path.open("w") as f:
    f.write("Xc_total,X_HMC,X_HMC_lo,X_HMC_hi,"
            "X_FD_grand_canonical,X_FD_canonical,ceiling,gap_HMC_canon,"
            "X_HMC_first_frame,X_HMC_last_frame,accepts\n")
    for pt in hmc_pred:
        f.write(f"{pt['Xc']:.4e},{pt['X_HMC']:.5f},{pt['lo']:.5f},{pt['hi']:.5f},"
                f"{pt['X_FD_GC']:.5f},{pt['X_FD_canon']:.5f},{pt['ceiling']:.5f},"
                f"{pt['gap_HMC_minus_canon']:+.5f},"
                f"{pt['first']:.5f},{pt['last']:.5f},{pt['accepts']}\n")
print(f"wrote {csv_path}")

# ---------------------- plot -----------------------------------------------
xs_hmc = np.array([p["Xc"] for p in hmc_pred])
y_hmc = np.array([p["X_HMC"] for p in hmc_pred])
y_lo = np.array([p["lo"] for p in hmc_pred])
y_hi = np.array([p["hi"] for p in hmc_pred])

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

ax = axes[0]
ax.plot(xc_grid, x_gb_gc, "-", color="C0", lw=2,
        label=r"$X_{\mathrm{GB}}^{\mathrm{FD}}$ grand-canonical (Wagih eq.)")
ax.plot(xc_grid, x_gb_canon, "-", color="C2", lw=2,
        label=r"$X_{\mathrm{GB}}^{\mathrm{FD}}$ canonical (mass-conserving)")
ax.plot(xc_grid, ceiling, ":", color="k", lw=1,
        label=r"closed-box ceiling $X_c \cdot N_\mathrm{tot}/N_\mathrm{GB}$")
ax.plot(xc_grid, xc_grid, "--", color="gray", lw=0.7,
        label=r"$X_{\mathrm{GB}} = X_c$ (no segregation)")
ax.errorbar(xs_hmc, y_hmc,
            yerr=[y_hmc - y_lo, y_hi - y_hmc],
            fmt="o", color="C3", ms=9, capsize=4, zorder=5,
            label=r"$X_{\mathrm{GB}}^{\mathrm{HMC}}$ (50 ps PROD)")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"total Mg fraction $X_c$ (closed box)")
ax.set_ylabel(r"$X_{\mathrm{GB}}$")
ax.set_title(f"Al(Mg) GB occupancy at T={T:g} K — N_GB/N_tot = {GB_FRAC:.3f}")
ax.legend(loc="lower right", fontsize=8)
ax.grid(True, which="both", alpha=0.3)
ax.set_xlim(1e-5, 0.5)
ax.set_ylim(1e-5, 1.5)

ax = axes[1]
ax.semilogx(xc_grid, x_gb_gc / xc_grid, "-", color="C0", lw=2,
            label="grand-canonical FD")
ax.semilogx(xc_grid, x_gb_canon / xc_grid, "-", color="C2", lw=2,
            label="canonical FD")
ax.semilogx(xc_grid, ceiling / xc_grid, ":", color="k", lw=1,
            label="ceiling = N_tot/N_GB")
ax.errorbar(xs_hmc, y_hmc / xs_hmc, yerr=[(y_hmc-y_lo)/xs_hmc, (y_hi-y_hmc)/xs_hmc],
            fmt="o", color="C3", ms=9, capsize=4, zorder=5, label="HMC")
ax.axhline(1, color="gray", lw=0.6, ls=":")
ax.set_yscale("log")
ax.set_xlabel(r"total Mg fraction $X_c$")
ax.set_ylabel(r"enrichment $X_{\mathrm{GB}} / X_c$")
ax.set_title("Enrichment factor — canonical FD bounded by 1/GB_frac")
ax.legend(loc="best", fontsize=8)
ax.grid(True, which="both", alpha=0.3)

fig.tight_layout()
fig.savefig(OUT / "hmc_vs_fd_T500_canonical.png", dpi=130)
print(f"wrote {OUT/'hmc_vs_fd_T500_canonical.png'}")

# ---------------------- console table --------------------------------------
print()
print(f"{'X_c':>10} {'X_HMC':>10} {'X_FD_GC':>10} {'X_FD_can':>10} "
      f"{'ceiling':>10} {'gap_canon':>11}")
for pt in hmc_pred:
    print(f"{pt['Xc']:>10.1e} {pt['X_HMC']:>10.4f} {pt['X_FD_GC']:>10.4f} "
          f"{pt['X_FD_canon']:>10.4f} {pt['ceiling']:>10.4f} "
          f"{pt['gap_HMC_minus_canon']:>+11.4f}")
