#!/usr/bin/env python3
"""Audit: canon-FD reference curve from OUR n=500 spectrum vs Wagih's
n=82,635 Zenodo spectrum (same alloy = Al-Mg, same potential = Mendelev
2009).

Both spectra have already been argued statistically indistinguishable by
KS test at p > 0.5 (project bar, memory `reference_ks_test.md`). The
audit asks the next question: even if KS-indistinguishable, do they
predict different X_GB(X_c) under the canonical (mass-conserving) FD
formula? FD integrates exp(-ΔE/kT) which is exquisitely sensitive to
the deep-negative tail (5% quantile and lower) — a regime where KS
test has weak power.

Outputs:
  output/audit_canon_fd_curves.json  — canon-FD curves on both spectra
  output/audit_canon_fd_curves.png   — figure for the report
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fermi_dirac_predict import x_gb_canonical_curve, x_gb_curve, load_wagih

OUT  = Path("/cluster/home/cainiu/Computational_modeling/project/output")
OURS_NPZ = Path("/cluster/scratch/cainiu/production_AlMg_200A/"
                "delta_e_results_n500_200A_tight.npz")
# LAMMPS ground-truth ΔE_seg for every GB site (n=82,646). The companion
# `..._GB_segregation.dump` file in `accelerated_model/` is Wagih's
# ML-predicted spectrum (~4 kJ/mol MAE smoothing) — do NOT use it as the
# reference: it shifts the deep tail by ~0.5 kJ/mol and inflates the
# canon-FD difference by a few × 0.001 in X_GB.
WAGIH_SEG_TXT  = Path(
    "/cluster/scratch/cainiu/wagih_zenodo/learning_segregation_energies/"
    "machine_learning_notebook/seg_energies_Al_Mg.txt")
WAGIH_BULK_DAT = Path(
    "/cluster/scratch/cainiu/wagih_zenodo/learning_segregation_energies/"
    "machine_learning_notebook/bulk_solute_Al_Mg.dat")

KJ_PER_EV = 96.485
T = 500.0


def main():
    de_ours  = np.load(OURS_NPZ)["gb_delta_e"]                # eV, n=500
    de_wagih = load_wagih(WAGIH_SEG_TXT, WAGIH_BULK_DAT)      # eV, n=82,646 LAMMPS truth

    # Polycrystal counts (canon-FD needs these for mass conservation).
    N_total_ours, N_GB_ours  = 475_715, 89_042
    N_total_wag,  N_GB_wag   = 483_425, 82_646

    # X_c grid spans dilute → concentrated.
    xc_grid = np.geomspace(1e-3, 0.5, 60)

    canon_ours,  _ = x_gb_canonical_curve(de_ours,  T, xc_grid,
                                          N_GB_ours, N_total_ours)
    canon_wagih, _ = x_gb_canonical_curve(de_wagih, T, xc_grid,
                                          N_GB_wag, N_total_wag)

    # Tail diagnostics.
    def stats(x_eV):
        x = x_eV * KJ_PER_EV
        q = np.quantile(x, [0.01, 0.05, 0.10, 0.50, 0.90, 0.95, 0.99])
        return dict(n=int(len(x)), mean=float(x.mean()), std=float(x.std()),
                    q01=float(q[0]), q05=float(q[1]), q10=float(q[2]),
                    q50=float(q[3]), q90=float(q[4]), q95=float(q[5]),
                    q99=float(q[6]))

    spec_ours  = stats(de_ours)
    spec_wagih = stats(de_wagih)

    # HMC measurements with X_c attached
    hmc_pts = [
        ("X_c=0.05 (bracket equilibrated)", 0.05,  0.2375, 0.0050),
        ("X_c=0.10 (preseg upper bound)",   0.10,  0.3749, None),
        ("X_c=0.15 (preseg upper bound)",   0.15,  0.5888, None),
        ("X_c=0.20 (preseg upper bound)",   0.20,  0.7942, None),
    ]

    # --- write JSON ---
    out_json = OUT / "audit_canon_fd_curves.json"
    out_json.write_text(json.dumps({
        "T": T,
        "spectrum_ours":  spec_ours,
        "spectrum_wagih": spec_wagih,
        "xc_grid":   xc_grid.tolist(),
        "canon_ours":  canon_ours.tolist(),
        "canon_wagih": canon_wagih.tolist(),
        "canon_diff_W_minus_O": (canon_wagih - canon_ours).tolist(),
        "hmc_pts": hmc_pts,
    }, indent=2))
    print(f"wrote {out_json}")

    # --- figure ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4),
                             gridspec_kw=dict(width_ratios=[1.4, 1]))

    # Panel L: canon-FD curves overlay + HMC points
    ax = axes[0]
    ax.fill_between(xc_grid, canon_ours, canon_wagih,
                    color="0.85", label="canon-FD reference band")
    ax.plot(xc_grid, canon_ours,  "C0-",  lw=1.4,
            label=f"canon-FD (ours, n={spec_ours['n']})")
    ax.plot(xc_grid, canon_wagih, "C3-",  lw=1.4,
            label=f"canon-FD (Wagih, n={spec_wagih['n']})")
    ax.plot(xc_grid, xc_grid, "k:", lw=0.7, label=r"$X_{GB}=X_c$")

    # HMC points
    for label, xc, x_hmc, ci in hmc_pts:
        is_eq = ci is not None
        if is_eq:
            ax.errorbar([xc], [x_hmc], yerr=[ci],
                        fmt="o", color="k", ms=7, capsize=3,
                        label=label.split(" (")[0]
                              if "0.05" in label else None)
        else:
            ax.plot([xc], [x_hmc], "v", color="0.4", ms=7, mfc="none",
                    label="HMC upper bound" if "0.10" in label else None)

    ax.set_xscale("log"); ax.set_xlim(1e-3, 0.5)
    ax.set_xlabel(r"$X_c$ (bulk Mg fraction)")
    ax.set_ylabel(r"$X_{GB}$")
    ax.set_title(f"canon-FD: ours vs Wagih (T={T:.0f} K)")
    ax.legend(loc="lower right", fontsize=8)
    ax.grid(alpha=0.3, which="both")

    # Panel R: spectrum left-tail comparison (the FD-relevant part)
    ax = axes[1]
    ax.hist(de_ours  * KJ_PER_EV, bins=60, range=(-80, 50),
            density=True, alpha=0.55, color="C0",
            label=f"ours n={spec_ours['n']}")
    ax.hist(de_wagih * KJ_PER_EV, bins=120, range=(-80, 50),
            density=True, alpha=0.55, color="C3",
            label=f"Wagih n={spec_wagih['n']}")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel(r"$\Delta E_{seg}$  (kJ/mol)")
    ax.set_ylabel("PDF")
    ax.set_title("ΔE_seg histogram")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_png = OUT / "audit_canon_fd_curves.png"
    fig.savefig(out_png, dpi=160)
    print(f"wrote {out_png}")

    # --- summary table for stdout ---
    print()
    print(f"OURS  spectrum: mean={spec_ours['mean']:+.2f}  std={spec_ours['std']:.2f}  "
          f"q05={spec_ours['q05']:.1f}  q10={spec_ours['q10']:.1f}  kJ/mol  (n={spec_ours['n']})")
    print(f"WAGIH spectrum: mean={spec_wagih['mean']:+.2f}  std={spec_wagih['std']:.2f}  "
          f"q05={spec_wagih['q05']:.1f}  q10={spec_wagih['q10']:.1f}  kJ/mol  (n={spec_wagih['n']})")
    print()
    xc_print = [0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.40]
    print(f"{'X_c':>8} {'canon-ours':>12} {'canon-Wagih':>12} {'Δ(W−O)':>10}")
    for xc in xc_print:
        i = np.argmin(np.abs(xc_grid - xc))
        d = canon_wagih[i] - canon_ours[i]
        print(f"{xc:>8.3f} {canon_ours[i]:>12.4f} {canon_wagih[i]:>12.4f} {d:>+10.4f}")

if __name__ == "__main__":
    main()
