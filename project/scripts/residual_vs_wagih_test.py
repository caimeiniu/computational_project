"""Test: residual-vs-Wagih analysis on a single snapshot (X_c=0.075).

Replaces the [-30,-5] kJ/mol window slope from `solute_correlation_analysis.py`
panel 3 with a residual approach that does NOT silently confound ΔE with n_local:

    r_i = I_i^HMC - P_Wagih(ΔE_i, X_c, T)

E[r_i | n_local_i = n] = 0 in Wagih null, regardless of joint (ΔE, n_local)
geometry. Slope ∂r̄/∂n_local < 0 → site-level Mg-Mg repulsion.

Single X_c=0.075 to verify the curve is monotonic before extending to all six.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
from solute_correlation_analysis import read_lmp, T_K, KT_KJMOL, R_LOCAL, wagih_p_fd  # noqa

LABEL = "0.075_preseg"
SNAPSHOT = REPO / "data/snapshots/hmc_T500_Xc0.075_preseg_final.lmp"
X_C = 0.075
GB_MASK = REPO / "data/snapshots/gb_mask_200A.npy"
REFERENCE_NPZ = REPO / "data/snapshots/delta_e_results_n500_200A_tight.npz"
OUT_DIR = REPO / "output"


def main() -> None:
    print(f"--- Residual analysis: {LABEL}  (X_c={X_C}, T={T_K:.0f} K) ---")

    ref = np.load(REFERENCE_NPZ)
    ref_ids_1based = np.asarray(ref["gb_site_ids"], dtype=np.int64)
    ref_de_kjmol = np.asarray(ref["gb_delta_e"], dtype=float) * 96.485
    print(f"  reference sites: n={len(ref_ids_1based)}, "
          f"ΔE range [{ref_de_kjmol.min():.1f}, {ref_de_kjmol.max():.1f}] kJ/mol")

    positions, types, box = read_lmp(SNAPSHOT)
    is_mg = (types == 2)
    gb_mask = np.load(GB_MASK)
    n_gb_mg = int((is_mg & gb_mask).sum())
    n_gb = int(gb_mask.sum())
    x_gb = n_gb_mg / n_gb
    print(f"  X_GB={x_gb:.4f}  N_GB_Mg={n_gb_mg:,}")

    p_emp = is_mg[ref_ids_1based - 1].astype(int)

    mg_tree = cKDTree(positions[is_mg], boxsize=box)
    ref_positions = positions[ref_ids_1based - 1]
    nbr_counts = np.array(
        [len(nbrs) for nbrs in mg_tree.query_ball_point(ref_positions, R_LOCAL)]
    )
    n_local = nbr_counts - p_emp

    p_wagih = wagih_p_fd(ref_de_kjmol, X_C)
    residual = p_emp - p_wagih
    print(f"  p_Wagih: mean={p_wagih.mean():.3f}  range "
          f"[{p_wagih.min():.3f}, {p_wagih.max():.3f}]")
    print(f"  residual: mean={residual.mean():+.4f}  std={residual.std():.4f}")

    n_unique = np.arange(int(n_local.min()), int(n_local.max()) + 1)
    means, ses, counts, p_w_mean = [], [], [], []
    for n in n_unique:
        sel = n_local == n
        n_sel = int(sel.sum())
        if n_sel == 0:
            means.append(np.nan); ses.append(np.nan)
            counts.append(0); p_w_mean.append(np.nan); continue
        r_bar = float(residual[sel].mean())
        p_bar = float(p_emp[sel].mean())
        se = np.sqrt(p_bar * (1 - p_bar) / n_sel) if 0 < p_bar < 1 else 1.0 / max(n_sel, 1)
        means.append(r_bar); ses.append(1.96 * se)
        counts.append(n_sel); p_w_mean.append(float(p_wagih[sel].mean()))

    means = np.array(means); ses = np.array(ses)
    counts = np.array(counts); p_w_mean = np.array(p_w_mean)

    print()
    print(f"  {'n_local':>8s} {'n_sites':>8s} {'<p_W>':>8s} {'r̄':>10s} {'95% CI':>10s}")
    for n, m, s, c, pw in zip(n_unique, means, ses, counts, p_w_mean):
        if c == 0:
            continue
        print(f"  {n:>8d} {c:>8d} {pw:>8.3f} {m:>+10.4f} {s:>10.4f}")

    n_min_per_bin = 10
    valid = (counts >= n_min_per_bin) & ~np.isnan(means)
    coeff = None
    if valid.sum() >= 2:
        coeff = np.polyfit(n_unique[valid], means[valid], 1,
                           w=1.0 / np.maximum(ses[valid], 1e-3))
        print(f"\n  weighted linear fit (n_sel >= {n_min_per_bin}):  "
              f"slope = {coeff[0]:+.4f} per Mg-neighbour")
        print(f"                        intercept = {coeff[1]:+.4f}")

    r_global = float(residual.mean())
    n_used = int(counts[valid].sum())

    fig, ax = plt.subplots(figsize=(7.2, 4.7))
    ax.axhline(0, color="0.55", ls="--", lw=1.0)
    ax.axhline(r_global, color="#1f77b4", ls=":", lw=1.2)
    ax.errorbar(n_unique[valid], means[valid], yerr=ses[valid],
                fmt="o", color="#d62728", capsize=3, ms=8, lw=1.4,
                label=r"HMC residual")
    if coeff is not None:
        xfit = np.linspace(n_unique[valid].min(), n_unique[valid].max(), 50)
        ax.plot(xfit, np.polyval(coeff, xfit), "k--", lw=1.6, alpha=0.8,
                label=rf"slope $= {coeff[0]:+.3f}$ / Mg-nbr")

    x_label_text = n_unique[valid].max() + 0.15
    ax.text(x_label_text, 0.025, "Wagih baseline",
            fontsize=9.5, color="0.4", ha="right", va="bottom")
    ax.text(x_label_text, r_global - 0.015,
            rf"$\langle r \rangle = {r_global:+.2f}$",
            fontsize=9.5, color="#1f77b4", ha="right", va="top")

    ax.set_xlabel(rf"$n_{{Mg}}^{{local}}$  (Mg neighbours within "
                  rf"$r \leq {R_LOCAL:.0f}$ Å)", fontsize=12)
    ax.set_ylabel(
        r"residual  $\bar r = \langle I_{HMC} - P_{Wagih}(\Delta E_i) \rangle$",
        fontsize=12,
    )
    ax.set_title(
        rf"Site-level Mg–Mg repulsion via Wagih residual  "
        rf"($X_c={X_C}$, $T={T_K:.0f}$ K)",
        fontsize=12, pad=14,
    )
    ax.set_xlim(n_unique[valid].min() - 0.5, n_unique[valid].max() + 0.5)
    ax.set_ylim(-0.6, 0.5)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.94)
    fig.tight_layout()
    out_png = OUT_DIR / "residual_vs_wagih_test_Xc0.075.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {out_png}")

    out_json = OUT_DIR / "residual_vs_wagih_test_Xc0.075.json"
    out_json.write_text(json.dumps({
        "label": LABEL, "x_c": X_C, "x_gb": x_gb,
        "T_K": T_K, "kT_kjmol": KT_KJMOL, "r_local_angstrom": R_LOCAL,
        "n_unique": n_unique.tolist(), "n_per_bin": counts.tolist(),
        "p_wagih_mean_per_bin": p_w_mean.tolist(),
        "residual_mean": means.tolist(), "residual_err95": ses.tolist(),
        "slope_per_neighbour": float(coeff[0]) if coeff is not None else None,
        "intercept": float(coeff[1]) if coeff is not None else None,
        "n_sites_total": int(len(ref_ids_1based)),
    }, indent=2, default=float))
    print(f"  Saved: {out_json}")


if __name__ == "__main__":
    main()
