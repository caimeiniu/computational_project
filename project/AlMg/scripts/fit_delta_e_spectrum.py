"""Fit the per-site ΔE_seg spectrum from sample_delta_e.py to a skew-normal.

Wagih uses the skew-normal PDF (Nat. Commun. 11:6376, eq. 1), labelled in
the paper with (μ, σ, α):
    F(ΔE; μ, σ, α) = (1/(σ√(2π))) exp(-(ΔE-μ)²/(2σ²))
                     × erfc(-α(ΔE-μ)/(σ√2))
Mathematically μ and σ are the *location* and *scale* of the skew-normal,
not the distribution's mean and std — but we keep Wagih's (μ, σ, α)
symbols verbatim so values copy directly between our outputs and the
paper/SI tables. scipy mapping: scipy.stats.skewnorm(loc=μ, scale=σ, a=α).

Reported alongside the fit: sample moments (min/max/mean/std/skew) and the
Wagih Al(Mg) `Mg^15` reference values (SI Fig. 3, potential = Mendelev 2009
= same Al-Mg.eam.fs used here) for apples-to-apples comparison.

Usage:
    python fit_delta_e_spectrum.py \\
        --npz /cluster/scratch/cainiu/prototype_AlMg_100A/delta_e_results_n500.npz \\
        --out-png /cluster/home/cainiu/Computational_modeling/project/output/delta_e_spectrum_n500.png \\
        --out-json /cluster/home/cainiu/Computational_modeling/project/output/delta_e_fit_n500.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats

EV_TO_KJMOL = 96.485

# Wagih Al(Mg) skew-normal parameters for the Mendelev 2009 potential
# (SI Ref [15], `Mg^15` panel in Supplementary Fig. 3 — same Al-Mg.eam.fs we use).
# Read from the SI figure panel; kJ/mol. Symbols match Wagih eq. 1 (μ, σ, α).
WAGIH_ALMG = {"mu": 9.0, "sigma": 23.0, "alpha": -2.3, "r2": 1.00,
              "range": (-60.0, 40.0),
              "source": "Wagih 2020 SI Fig. 3, Mg^15 (potential = Mendelev 2009)"}


def fit_skewnorm(delta_e_kjmol: np.ndarray) -> dict:
    a, loc, scale = stats.skewnorm.fit(delta_e_kjmol)
    return {"alpha": float(a), "mu": float(loc), "sigma": float(scale)}


def sample_moments(x: np.ndarray) -> dict:
    return {
        "n": int(x.size),
        "min": float(x.min()),
        "max": float(x.max()),
        "mean": float(x.mean()),
        "median": float(np.median(x)),
        "std": float(x.std(ddof=1)),
        "skew": float(stats.skew(x)),
    }


def plot_fit(delta_e_kjmol: np.ndarray, fit: dict, png_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.0, 4.0), dpi=120)
    ax.hist(delta_e_kjmol, bins=40, density=True, color="#4c72b0",
            alpha=0.65, edgecolor="white", label=f"ΔE_seg (n={delta_e_kjmol.size})")

    x = np.linspace(delta_e_kjmol.min() - 5, delta_e_kjmol.max() + 5, 400)
    pdf_fit = stats.skewnorm.pdf(x, a=fit["alpha"], loc=fit["mu"], scale=fit["sigma"])
    ax.plot(x, pdf_fit, "-", color="#c44e52", lw=2.0,
            label=f"skew-normal fit\nμ={fit['mu']:+.2f}, σ={fit['sigma']:.2f}, α={fit['alpha']:+.2f}")

    pdf_w = stats.skewnorm.pdf(x, a=WAGIH_ALMG["alpha"],
                               loc=WAGIH_ALMG["mu"], scale=WAGIH_ALMG["sigma"])
    ax.plot(x, pdf_w, "--", color="#55a868", lw=1.5,
            label=(f"Wagih Mg$^{{15}}$ (SI Fig. 3)\n"
                   f"μ={WAGIH_ALMG['mu']:+.1f}, σ={WAGIH_ALMG['sigma']:.1f}, "
                   f"α={WAGIH_ALMG['alpha']:+.1f}"))

    ax.axvline(0.0, color="k", lw=0.6, alpha=0.5)
    ax.set_xlabel("ΔE_seg  (kJ/mol)")
    ax.set_ylabel("probability density")
    ax.set_title("Al(Mg) per-site GB segregation spectrum")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--npz", required=True,
                   help=".npz from sample_delta_e.py (contains gb_delta_e in eV)")
    p.add_argument("--out-png", default=None, help="output histogram + fit plot")
    p.add_argument("--out-json", default=None, help="output fit-params JSON")
    args = p.parse_args()

    data = np.load(args.npz, allow_pickle=False)
    dE_ev = np.asarray(data["gb_delta_e"])
    dE = dE_ev * EV_TO_KJMOL

    moments = sample_moments(dE)
    fit = fit_skewnorm(dE)

    # Dump a human summary to stdout.
    print("Sample moments (kJ/mol):")
    for k, v in moments.items():
        print(f"  {k:<7s} = {v}")
    print("\nSkew-normal fit (kJ/mol, Wagih eq. 1 symbols):")
    print(f"  μ (location)  = {fit['mu']:+.3f}")
    print(f"  σ (scale)     = {fit['sigma']:.3f}")
    print(f"  α (shape)     = {fit['alpha']:+.3f}")

    print("\nWagih Al(Mg) reference (SI Fig. 3, Mg^15 = Mendelev 2009):")
    print(f"  μ = {WAGIH_ALMG['mu']:+.1f}, σ = {WAGIH_ALMG['sigma']:.1f}, "
          f"α = {WAGIH_ALMG['alpha']:+.2f}, R² = {WAGIH_ALMG['r2']:.2f}")
    print(f"  ΔE range        = {WAGIH_ALMG['range']}")

    if args.out_png:
        plot_fit(dE, fit, Path(args.out_png))
        print(f"\n→ wrote plot   {args.out_png}")
    if args.out_json:
        out = {"sample_moments_kjmol": moments,
               "skewnorm_fit_kjmol": fit,
               "wagih_alm_reference": WAGIH_ALMG,
               "source_npz": str(Path(args.npz).resolve())}
        Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_json).write_text(json.dumps(out, indent=2))
        print(f"→ wrote params {args.out_json}")


if __name__ == "__main__":
    main()
