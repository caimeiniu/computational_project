"""Defense-prep figure replots: 3 audience-friendly mechanism panels.

Reads `output/solute_correlation_analysis.json` (already produced by
scripts/solute_correlation_analysis.py) and produces 3 separate PNGs
sized for slides:

    1. defense_MgMg_clustering.png     — g(r) clustering ratio, 3 curves
    2. defense_occupation_breakdown.png — P_i vs ΔE_i, 3 panels with
       favourable-energy band shaded
    3. defense_repulsion_summary.png   — slope vs X_c summary + zoom on
       steepest case

No LAMMPS, no recomputation. ~5 s on a login node.

Run:
    python scripts/replot_mechanism_for_defense.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO = Path(__file__).resolve().parent.parent
JSON_PATH = REPO / "output/solute_correlation_analysis.json"
OUT_DIR = REPO / "output"

REPRESENTATIVE = ["0.075_preseg", "0.15_preseg", "0.30_preseg"]
COLORS = {
    "0.075_preseg":    "#d62728",  # red — threshold breakdown
    "0.15_preseg":     "#1f77b4",  # blue — mid range
    "0.30_preseg":     "#2ca02c",  # green — saturation regime
    "0.10_multistart": "#7f7f7f",  # gray — kinetic-floor (UB) family
    "0.10_preseg":     "#9467bd",  # purple — additional preseg point
    "0.20_preseg":     "#ff7f0e",  # orange — additional preseg point
}


def _label(label: str, x_gb: float) -> str:
    x_c = float(label.split("_")[0])
    return rf"$X_c={x_c:.3f}$, $X_{{GB}}={x_gb:.3f}$"


def fig1_clustering(d: dict) -> None:
    r = np.array(d["g_MgMg_pair_correlation"]["r_axis"])
    fig, ax = plt.subplots(figsize=(8.0, 5.2))

    for label in REPRESENTATIVE:
        c = d["g_MgMg_pair_correlation"]["curves"][label]
        ratio = np.array(c["ratio"], dtype=float)
        ax.plot(r, ratio, color=COLORS[label], lw=2.2,
                label=_label(label, c["x_gb"]))

    ax.axhline(1.0, color="k", lw=0.8, ls="--", alpha=0.5)
    ax.text(24.5, 1.02, "uniform-random reference  ($g_{HMC}=g_{random}$)",
            fontsize=9, color="k", alpha=0.7, ha="right", va="bottom")

    ax.axvspan(3.0, 3.5, color="0.85", alpha=0.45, zorder=0)
    ax.text(3.25, 0.08, "1st NN shell\n(FCC Al)", fontsize=8.5,
            ha="center", va="bottom", alpha=0.8)

    ax.set_xlim(0, 25)
    ax.set_ylim(0, 1.55)
    ax.set_xlabel("pair separation  $r$  [Å]", fontsize=12)
    ax.set_ylabel(r"$g_{Mg-Mg}^{HMC}(r)\ /\ g_{Mg-Mg}^{random}(r)$",
                  fontsize=12)
    ax.set_title(
        "Mg–Mg spatial clustering on the GB (HMC vs uniform-random reference)\n"
        r"$> 1$: clustered    $< 1$: avoidant    "
        r"(non-random structure out to $r \sim 10$ Å)",
        fontsize=12.5,
    )
    ax.legend(loc="upper right", framealpha=0.92,
              title="HMC at $T=500$ K", title_fontsize=10, fontsize=10)
    ax.grid(alpha=0.3)

    fig.text(
        0.5, 0.005,
        "Aggregate spatial signal. The peak above $1$ is partly driven by "
        "the geometric proximity of deep-$\\Delta E$ binding sites; the "
        "$\\Delta E$-controlled (residual) interaction is in Fig. 3.",
        ha="center", fontsize=9, style="italic", color="0.4", wrap=True,
    )

    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(OUT_DIR / "defense_MgMg_clustering.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


def wagih_fd(de_kjmol: np.ndarray, x_c: float, kT_kjmol: float) -> np.ndarray:
    return 1.0 / (1.0 + ((1.0 - x_c) / x_c) * np.exp(de_kjmol / kT_kjmol))


def fig2_occupation_breakdown(d: dict) -> None:
    kT = d["kT_kjmol"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), sharey=True)

    for ax, label in zip(axes, REPRESENTATIVE):
        s = d["site_occupation_vs_energy"]["snapshots"][label]
        de = np.array(s["bin_centers_kjmol"])
        p = np.array(s["p_empirical"])
        lo = np.array(s["p_low95"])
        hi = np.array(s["p_high95"])
        x_c = s["x_c"]
        x_gb = s["x_gb"]

        ax.axvspan(de.min() - 5, 0, color="#a8d5a8", alpha=0.22, zorder=0)
        ax.axvline(0, color="#2a662a", lw=0.8, ls="-", alpha=0.5)

        de_grid = np.linspace(de.min(), de.max(), 200)
        p_wagih = wagih_fd(de_grid, x_c, kT)
        ax.plot(de_grid, p_wagih, "k-", lw=1.8,
                label="Wagih FD prediction")

        valid = ~np.isnan(p)
        ax.errorbar(
            de[valid], p[valid],
            yerr=[p[valid] - lo[valid], hi[valid] - p[valid]],
            fmt="o", color=COLORS[label], capsize=3,
            markersize=7, lw=1.2, label="HMC empirical",
        )

        favourable_mask = (de < 0) & valid
        if favourable_mask.any():
            gap = wagih_fd(de[favourable_mask], x_c, kT) - p[favourable_mask]
            i = int(np.argmax(gap))
            de_a = de[favourable_mask][i]
            p_a = p[favourable_mask][i]
            p_w_a = wagih_fd(de_a, x_c, kT)
            if p_w_a - p_a > 0.08:
                ax.annotate(
                    f"$\\Delta P_i \\approx {p_w_a - p_a:.2f}$",
                    xy=(de_a, p_a),
                    xytext=(de_a + 5, p_a - 0.20),
                    arrowprops=dict(arrowstyle="->", color="0.25", lw=1.2),
                    fontsize=10.5, color="0.15",
                    bbox=dict(boxstyle="round,pad=0.25",
                              fc="white", ec="0.6", lw=0.6),
                )

        ax.set_xlabel(r"$\Delta E_i$  [kJ/mol]   ($X_c=0$ reference)",
                      fontsize=11)
        ax.set_title(rf"$X_c={x_c:.3f}$  $\rightarrow$  "
                     rf"$X_{{GB}}={x_gb:.3f}$", fontsize=12.5)
        ax.set_ylim(-0.05, 1.08)
        ax.grid(alpha=0.25)
        ax.legend(loc="lower right", fontsize=9.5, framealpha=0.92)

    axes[0].set_ylabel(r"$P_i$ — probability site $i$ is Mg",
                       fontsize=11.5)
    fig.suptitle(
        "Where Wagih's site-independent FD breaks: empirical $P_i$ "
        "below sigmoid at favourable $\\Delta E_i$",
        y=1.02, fontsize=13,
    )
    fig.text(
        0.5, -0.02,
        "Black curve: Wagih FD  $P_i = 1\\,/\\,(1 + ((1-X_c)/X_c)\\exp(\\Delta E_i/kT))$  "
        "(independent-site prediction).   Coloured points: HMC equilibrium occupancy "
        "binned by $\\Delta E_i$, 95% CIs.   Green band:  $\\Delta E < 0$  "
        "(favourable binding).",
        ha="center", fontsize=9, style="italic", color="0.4",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_DIR / "defense_occupation_breakdown.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


def fig3_repulsion_summary(d: dict) -> None:
    snapshots = d["site_occupation_vs_density"]["snapshots"]
    de_window = d["site_occupation_vs_density"]["delta_e_window_kjmol"]
    r_local = d["site_occupation_vs_density"]["r_local_angstrom"]
    n_window = d["site_occupation_vs_density"]["n_sites_in_window"]

    preseg = ["0.075_preseg", "0.10_preseg", "0.15_preseg",
              "0.20_preseg", "0.30_preseg"]
    multistart = "0.10_multistart"

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 4.7))

    # === Left: slope vs X_c ===
    xs = np.array([snapshots[l]["x_c"] for l in preseg])
    ys = np.array([snapshots[l]["slope_per_neighbour"] for l in preseg])
    ax1.plot(xs, ys, "o-", color="#d62728", lw=2.2, ms=9,
             label="preseg HMC (equilibrating)")

    ms = snapshots[multistart]
    ax1.plot([ms["x_c"]], [ms["slope_per_neighbour"]],
             "s", color="0.45", ms=9, mfc="none", mew=1.8,
             label="multistart UB (kinetic-floor)")

    ax1.axhline(0, color="k", ls="--", lw=0.9, alpha=0.55)
    ax1.text(0.305, 0.004, "no interaction (slope = 0)",
             fontsize=9.5, color="0.35", ha="right", va="bottom")

    # Annotate steepest and saturation crossover
    s075 = snapshots["0.075_preseg"]["slope_per_neighbour"]
    s020 = snapshots["0.20_preseg"]["slope_per_neighbour"]
    ax1.annotate(
        f"slope = {s075:.3f}\n(steepest repulsion)",
        xy=(0.075, s075),
        xytext=(0.115, -0.108),
        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.2),
        fontsize=10.5, color="0.15", ha="left",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", lw=0.6),
    )
    ax1.annotate(
        "saturation regime\n(sites mostly full)",
        xy=(0.20, s020),
        xytext=(0.245, -0.045),
        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.2,
                        connectionstyle="arc3,rad=-0.25"),
        fontsize=10.5, color="0.15", ha="center",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", lw=0.6),
    )

    ax1.set_xlabel(r"$X_c$  (total Mg fraction)", fontsize=11.5)
    ax1.set_ylabel(
        r"slope  $\partial P_i / \partial n_{Mg}^{local}$  "
        r"[per Mg-neighbour]",
        fontsize=11.5,
    )
    ax1.set_title(
        f"Site-level Mg–Mg repulsion vs $X_c$\n"
        f"($\\Delta E \\in [{de_window[0]:.0f}, {de_window[1]:.0f}]$ kJ/mol,  "
        f"$r \\leq {r_local:.0f}$ Å,  $n={n_window}$ sites)",
        fontsize=12,
    )
    ax1.legend(loc="upper right", fontsize=9.5, framealpha=0.92)
    ax1.grid(alpha=0.3)
    ax1.set_xlim(0.04, 0.33)
    ax1.set_ylim(-0.13, 0.06)

    # === Right: zoom on X_c=0.075 (steepest signal) ===
    s = snapshots["0.075_preseg"]
    n_mg = np.array(s["n_mg_unique"])
    p = np.array(s["p_empirical"])
    err = np.array(s["p_err95"])
    counts = np.array(s["n_per_bin"])

    valid = (~np.isnan(p)) & (counts >= 3)
    ax2.errorbar(
        n_mg[valid], p[valid], yerr=err[valid],
        fmt="o", color="#d62728", capsize=3, ms=7, lw=1.3,
        label=f"empirical (n={n_window} sites in $\\Delta E$ window)",
    )

    if valid.sum() >= 2:
        coeff = np.polyfit(n_mg[valid], p[valid], 1,
                           w=1.0 / np.maximum(err[valid], 1e-3))
        xfit = np.linspace(n_mg.min(), n_mg.max(), 50)
        ax2.plot(xfit, np.polyval(coeff, xfit), "k--", lw=1.7,
                 label=rf"linear fit:  slope $= {coeff[0]:+.3f}$ / Mg-nbr")

    ax2.set_xlabel(rf"$n_{{Mg}}^{{local}}$  (Mg neighbours within "
                   rf"$r \leq {r_local:.0f}$ Å)", fontsize=11.5)
    ax2.set_ylabel(r"$P_i$  within  $\Delta E \in [-30, -5]$  kJ/mol",
                   fontsize=11.5)
    ax2.set_title(r"Detail:  $X_c=0.075$  —  direct repulsion evidence",
                  fontsize=12)
    ax2.set_ylim(-0.05, 1.08)
    ax2.grid(alpha=0.3)
    ax2.legend(loc="upper right", fontsize=10, framealpha=0.92)

    fig.suptitle(
        "Direct site-level repulsion: more Mg neighbours → lower occupation "
        "(at fixed $\\Delta E_i$)",
        y=1.02, fontsize=13,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_DIR / "defense_repulsion_summary.png", dpi=200,
                bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    d = json.loads(JSON_PATH.read_text())

    fig1_clustering(d)
    fig2_occupation_breakdown(d)
    fig3_repulsion_summary(d)

    print(f"Wrote 3 PNGs to {OUT_DIR}/:")
    print("    defense_MgMg_clustering.png")
    print("    defense_occupation_breakdown.png")
    print("    defense_repulsion_summary.png")


if __name__ == "__main__":
    main()
