"""Mechanism Fig 3 (single panel, double-column paper figure):
Mg occupancy P_i vs n_Mg^local at fixed ΔE — site-level Mg-Mg
interaction signal.

At X_c=0.075, restrict to the favourable ΔE bin [-30, -15] kJ/mol
(n=112 sites). Within that window the per-site Wagih FD prediction is
held in [0.75, 0.99] (the gray band). Sites are grouped by integer
n_Mg^local (Mg neighbours within 5 Å) — one marker per integer, no
pre-binning; markers sit at integer x-positions (n_local is an atom
count, not a continuous variable). Empirical Mg-occupation fractions
are shown with Wald 95 % binomial CI.

Filtering: integers with fewer than N_SITES_MIN sites are dropped
from the plot (CIs span essentially [0, 1] and add no information);
they are still written to the JSON for full auditability. With
N_SITES_MIN=3 the figure shows n=1 through n=9 (9 markers).

Reading: low n_local (1-5) sits in P~0.33-0.57 with substantial CI
overlap (mostly below the Wagih band, but n=1 with 3 sites overlaps);
sharp drop at n_local=5→6 to P~0.13-0.15. ΔE is held within a 15
kJ/mol-wide window so the n_local dependence cannot be a confounded
ΔE effect → site-level Mg-Mg interaction signal.

Falsification check: same analysis on the neutral ΔE bin [-5, +5]
kJ/mol (n=128, Wagih predicts P ∈ [0.024, 0.213]) is computed and
emitted in the JSON for the caption — slope much weaker, consistent
with no signal where Wagih predicts little occupation.

No LAMMPS, no recomputation of ΔE — uses cached n=500 spectrum.

Re-running with the equilibrated snapshot (post-job 65430294): change
SNAPSHOT path; everything else stays.
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

XC = 0.075
SNAPSHOT = REPO / "AlMg/data/snapshots/hmc_T500_Xc0.075_preseg_final.lmp"
GB_MASK = REPO / "AlMg/data/snapshots/gb_mask_200A.npy"
REFERENCE_NPZ = REPO / "AlMg/data/snapshots/delta_e_results_n500_200A_tight.npz"
OUT_DIR = REPO / "output"
OUT_PREFIX = "03_mgmg_repulsion_fixed_dE"

# Favourable ΔE bin: narrow enough that Wagih FD prediction is roughly
# constant across it (P ∈ [0.75, 0.99] = 0.24 wide), wide enough to
# give ~22 sites per n_local sub-bin.
DE_FAV = (-30.0, -15.0)
# Neutral bin (control / falsification check); Mg-Mg repulsion has no
# headroom to register because Wagih's P is small throughout.
DE_NEU = (-5.0, +5.0)

# Plot every integer n_local — no pre-binning. Markers sit at integer
# x-positions (the only physically meaningful values: n_local is an atom
# count, not a continuous variable). Sub-bin merging tried earlier
# always introduced a half-integer marker plus the choice of where to
# cut, both of which raised "what physics drives this binning" questions.
# Per-integer keeps the underlying data visible.
N_LOCAL_MAX = 11

# Drop n_local values with fewer than this many sites in the favourable
# ΔE bin: at very small n_sites the Wald CI spans most of [0, 1] and
# adds visual noise without information. n_sites=1 gives a vacuous
# CI [0, 1]; n_sites=2 likewise. Setting this to 3 keeps the n=1
# (3 sites, CI [0.35, 1.00]) "Wagih at empty neighbourhood" anchor.
N_SITES_MIN = 3


def _wald_ci(p_hat: float, n: int) -> tuple[float, float]:
    if n == 0:
        return float("nan"), float("nan")
    if 0 < p_hat < 1:
        se = float(np.sqrt(p_hat * (1 - p_hat) / n))
    else:
        se = 1.0 / max(n, 1)
    return max(0.0, p_hat - 1.96 * se), min(1.0, p_hat + 1.96 * se)


def _bin_stats(p_arr: np.ndarray, nl_arr: np.ndarray, sel_de: np.ndarray):
    """Per-integer stats restricted to the ΔE selection.

    One row per integer n_local in [0, N_LOCAL_MAX]; rows with fewer
    than N_SITES_MIN sites are kept in the JSON (so the full dist is
    auditable) but flagged via `n_sites` so the plot can filter them.
    """
    p_in = p_arr[sel_de]
    nl_in = nl_arr[sel_de]
    rows = []
    for n_val in range(N_LOCAL_MAX + 1):
        sub = nl_in == n_val
        n_sub = int(sub.sum())
        p_hat = float(p_in[sub].mean()) if n_sub > 0 else float("nan")
        ci_lo, ci_hi = _wald_ci(p_hat, n_sub)
        rows.append({
            "n_local": n_val,
            "n_sites": n_sub,
            "p_hat": p_hat, "ci_lo": ci_lo, "ci_hi": ci_hi,
        })
    slope = float("nan")
    intercept = float("nan")
    if int(sel_de.sum()) > 5:
        coeff = np.polyfit(nl_in.astype(float), p_in.astype(float), 1)
        slope, intercept = float(coeff[0]), float(coeff[1])
    return rows, slope, intercept


def main() -> None:
    print(f"--- Mechanism Fig 3:  X_c={XC}, T={T_K:.0f} K, kT={KT_KJMOL:.2f} kJ/mol ---")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    ref = np.load(REFERENCE_NPZ)
    ref_ids_1based = np.asarray(ref["gb_site_ids"], dtype=np.int64)
    ref_de = np.asarray(ref["gb_delta_e"], dtype=float) * 96.485

    positions, types, box = read_lmp(SNAPSHOT)
    is_mg = (types == 2)
    p_emp = is_mg[ref_ids_1based - 1].astype(int)
    mg_tree = cKDTree(positions[is_mg], boxsize=box)
    ref_pos = positions[ref_ids_1based - 1]
    nbr_counts = np.array(
        [len(n) for n in mg_tree.query_ball_point(ref_pos, R_LOCAL)]
    )
    n_local = nbr_counts - p_emp  # exclude self if Mg

    sel_fav = (ref_de >= DE_FAV[0]) & (ref_de < DE_FAV[1])
    sel_neu = (ref_de >= DE_NEU[0]) & (ref_de < DE_NEU[1])

    # Wagih FD bands. wagih_p_fd is monotonically decreasing in ΔE.
    pw_fav_min = float(wagih_p_fd(np.array([DE_FAV[1]]), XC)[0])  # ΔE upper edge
    pw_fav_max = float(wagih_p_fd(np.array([DE_FAV[0]]), XC)[0])  # ΔE lower edge
    pw_fav_mean = float(wagih_p_fd(ref_de[sel_fav], XC).mean())
    pw_neu_min = float(wagih_p_fd(np.array([DE_NEU[1]]), XC)[0])
    pw_neu_max = float(wagih_p_fd(np.array([DE_NEU[0]]), XC)[0])

    rows_fav, slope_fav, intercept_fav = _bin_stats(p_emp, n_local, sel_fav)
    rows_neu, slope_neu, intercept_neu = _bin_stats(p_emp, n_local, sel_neu)

    n_fav = int(sel_fav.sum())
    n_neu = int(sel_neu.sum())

    print(f"  favorable ΔE [{DE_FAV[0]:+.0f}, {DE_FAV[1]:+.0f}] kJ/mol: "
          f"n={n_fav},  Wagih FD ∈ [{pw_fav_min:.3f}, {pw_fav_max:.3f}],  "
          f"slope={slope_fav:+.4f} per Mg-nbr")
    print(f"  neutral   ΔE [{DE_NEU[0]:+.0f}, {DE_NEU[1]:+.0f}] kJ/mol: "
          f"n={n_neu},  Wagih FD ∈ [{pw_neu_min:.3f}, {pw_neu_max:.3f}],  "
          f"slope={slope_neu:+.4f} per Mg-nbr")
    print()
    for r in rows_fav:
        if r["n_sites"] == 0:
            print(f"    fav n_local={r['n_local']:2d}: n=0  (no sites)")
            continue
        flag = "" if r["n_sites"] >= N_SITES_MIN else "  [excluded from plot, n<{}]".format(N_SITES_MIN)
        print(f"    fav n_local={r['n_local']:2d}: "
              f"n={r['n_sites']:3d}  P̂={r['p_hat']:.3f}  "
              f"CI=[{r['ci_lo']:.3f}, {r['ci_hi']:.3f}]{flag}")

    # ---- single-panel plot, sized for double-column paper figure ----
    fig, ax = plt.subplots(figsize=(5.2, 3.8))

    # Wagih band (favorable ΔE) — shaded; the band is the prediction
    # range, NOT a statistical confidence interval. The shaded region's
    # own boundaries are visible enough that explicit edge axhlines are
    # redundant — dropped per Fig-3 review pass 2 #4.
    ax.axhspan(pw_fav_min, pw_fav_max, color="0.78", alpha=0.45,
               zorder=0, label="Wagih FD prediction")

    # Empirical points at integer n_local — filter to n_sites >= N_SITES_MIN
    # so we don't pollute the panel with vacuous CIs from 1- or 2-site bins.
    # No connecting line: per-integer the data has real noise (small n_sites
    # at most integers) and a line would over-imply a smooth function.
    nl_int = np.array([r["n_local"] for r in rows_fav])
    p_hats = np.array([r["p_hat"] for r in rows_fav])
    ci_los = np.array([r["ci_lo"] for r in rows_fav])
    ci_his = np.array([r["ci_hi"] for r in rows_fav])
    n_sites_arr = np.array([r["n_sites"] for r in rows_fav])
    keep = (n_sites_arr >= N_SITES_MIN) & ~np.isnan(p_hats)
    ax.errorbar(
        nl_int[keep], p_hats[keep],
        yerr=[p_hats[keep] - ci_los[keep], ci_his[keep] - p_hats[keep]],
        fmt="o", color="#d62728", ms=5, capsize=2.5,
        elinewidth=0.8, mew=0.8,
        label="HMC empirical",
    )

    ax.set_xlabel(
        rf"$n_\mathrm{{Mg}}^\mathrm{{local}}$  (Mg neighbours within "
        rf"$r \leq {R_LOCAL:.0f}$ Å)",
        fontsize=10,
    )
    ax.set_ylabel(r"$P_i$  (probability site is Mg)", fontsize=10)
    ax.set_title(
        rf"Mg occupancy at $\Delta E_i \in [{DE_FAV[0]:.0f}, {DE_FAV[1]:.0f}]$ kJ/mol  "
        rf"($T={T_K:.0f}$ K, $X_c={XC}$)",
        fontsize=10.5, pad=16,
    )
    ax.set_xlim(-0.7, N_LOCAL_MAX + 0.7)
    ax.set_xticks(range(0, N_LOCAL_MAX + 1, 2))
    # Top headroom raised so the upper-right legend sits above the
    # Wagih band (which extends to P=0.99) instead of overlapping it.
    ax.set_ylim(-0.08, 1.22)
    ax.tick_params(labelsize=9)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
    ax.grid(alpha=0.25, lw=0.4)
    ax.legend(loc="upper right", fontsize=6.5, framealpha=0.92)

    fig.tight_layout(pad=1.4)
    out_png = OUT_DIR / f"{OUT_PREFIX}.png"
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {out_png}")

    out_json = OUT_DIR / f"{OUT_PREFIX}.json"
    out_json.write_text(json.dumps({
        "x_c": XC, "T_K": T_K, "kT_kjmol": KT_KJMOL,
        "r_local_angstrom": R_LOCAL,
        "snapshot": str(SNAPSHOT.relative_to(REPO)),
        "favorable_bin": {
            "delta_e_kjmol": list(DE_FAV),
            "n_sites": n_fav,
            "wagih_p_min": pw_fav_min,
            "wagih_p_max": pw_fav_max,
            "wagih_p_mean": pw_fav_mean,
            "linear_slope_per_neighbour": slope_fav,
            "linear_intercept": intercept_fav,
            "n_sites_min_for_plot": N_SITES_MIN,
            "per_integer_stats": rows_fav,
        },
        "neutral_bin_falsification": {
            "delta_e_kjmol": list(DE_NEU),
            "n_sites": n_neu,
            "wagih_p_min": pw_neu_min,
            "wagih_p_max": pw_neu_max,
            "linear_slope_per_neighbour": slope_neu,
            "linear_intercept": intercept_neu,
            "per_integer_stats": rows_neu,
        },
    }, indent=2, default=float))
    print(f"Saved: {out_json}")


if __name__ == "__main__":
    main()
