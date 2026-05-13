"""Fit the per-site ΔE_seg spectrum from sample_delta_e_PtAu.py to a skew-normal.

Pt(Au) sibling of project/scripts/fit_delta_e_spectrum.py — kept as a
separate file so the canonical Al(Mg) version is untouched. Only three
substantive differences: title string, dashed Wagih reference values
(Pt(Au) instead of Al(Mg)), JSON key name.

Wagih uses the skew-normal PDF (Nat. Commun. 11:6376, eq. 1), labelled
in the paper with (μ, σ, α):
    F(ΔE; μ, σ, α) = (1/(σ√(2π))) exp(-(ΔE-μ)²/(2σ²))
                     × erfc(-α(ΔE-μ)/(σ√2))
Mathematically μ and σ are the *location* and *scale* of the skew-normal,
not the distribution's mean and std — but we keep Wagih's (μ, σ, α)
symbols verbatim so values copy directly between our outputs and the
paper/SI tables. scipy mapping: scipy.stats.skewnorm(loc=μ, scale=σ, a=α).

Reported alongside the fit: sample moments (min/max/mean/std/skew) and
the Wagih Pt(Au) reference values computed from the Zenodo dump
`Pt_Au_20nm_GB_segregation.dump` (potential = O'Brien 2017 = same
PtAu.eam.alloy used here, so apples-to-apples).

Usage:
    python fit_delta_e_spectrum_PtAu.py \\
        --npz /cluster/scratch/cainiu/prototype_PtAu_100A/delta_e_results_n500_PtAu_100A_tight.npz \\
        --out-png /cluster/home/cainiu/Computational_modeling/project/PtAu/output/delta_e_spectrum_PtAu_100A.png \\
        --out-json /cluster/home/cainiu/Computational_modeling/project/PtAu/output/delta_e_fit_PtAu_100A.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy import stats

EV_TO_KJMOL = 96.485

# Wagih Pt(Au) skew-normal parameters for the O'Brien 2017 EAM.
# Computed by scipy.stats.skewnorm.fit on Wagih's Pt_Au_20nm_GB_segregation.dump
# (n=97,440 GB sites, seg_kJ_per_mol column, accelerated_model database).
# Symbols match Wagih eq. 1 (μ, σ, α). Range reflects observed min/max in dump.
WAGIH_PTAU = {"mu": 3.65, "sigma": 11.92, "alpha": -1.42,
              "range": (-50.0, 30.0),
              "n_ref": 97440,
              "source": ("Wagih 2020 Zenodo accelerated_model "
                         "Pt_Au_20nm_GB_segregation.dump (O'Brien 2017 EAM)")}


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

    pdf_w = stats.skewnorm.pdf(x, a=WAGIH_PTAU["alpha"],
                               loc=WAGIH_PTAU["mu"], scale=WAGIH_PTAU["sigma"])
    ax.plot(x, pdf_w, "--", color="#55a868", lw=1.5,
            label=(f"Wagih Pt(Au) (Zenodo dump, n={WAGIH_PTAU['n_ref']})\n"
                   f"μ={WAGIH_PTAU['mu']:+.2f}, σ={WAGIH_PTAU['sigma']:.2f}, "
                   f"α={WAGIH_PTAU['alpha']:+.2f}"))

    ax.axvline(0.0, color="k", lw=0.6, alpha=0.5)
    ax.set_xlabel("ΔE_seg  (kJ/mol)")
    ax.set_ylabel("probability density")
    ax.set_title("Pt(Au) per-site GB segregation spectrum")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--npz", required=True,
                   help=".npz from sample_delta_e_PtAu.py (contains gb_delta_e in eV)")
    p.add_argument("--out-png", default=None, help="output histogram + fit plot")
    p.add_argument("--out-json", default=None, help="output fit-params JSON")
    args = p.parse_args()

    data = np.load(args.npz, allow_pickle=False)
    dE_ev = np.asarray(data["gb_delta_e"])
    dE = dE_ev * EV_TO_KJMOL

    moments = sample_moments(dE)
    fit = fit_skewnorm(dE)

    print("Sample moments (kJ/mol):")
    for k, v in moments.items():
        print(f"  {k:<7s} = {v}")
    print("\nSkew-normal fit (kJ/mol, Wagih eq. 1 symbols):")
    print(f"  μ (location)  = {fit['mu']:+.3f}")
    print(f"  σ (scale)     = {fit['sigma']:.3f}")
    print(f"  α (shape)     = {fit['alpha']:+.3f}")

    print("\nWagih Pt(Au) reference (Zenodo dump, O'Brien 2017 EAM):")
    print(f"  μ = {WAGIH_PTAU['mu']:+.2f}, σ = {WAGIH_PTAU['sigma']:.2f}, "
          f"α = {WAGIH_PTAU['alpha']:+.2f}  (n_ref={WAGIH_PTAU['n_ref']})")
    print(f"  ΔE range       = {WAGIH_PTAU['range']}")

    if args.out_png:
        plot_fit(dE, fit, Path(args.out_png))
        print(f"\n→ wrote plot   {args.out_png}")
    if args.out_json:
        out = {"sample_moments_kjmol": moments,
               "skewnorm_fit_kjmol": fit,
               "wagih_ptau_reference": WAGIH_PTAU,
               "source_npz": str(Path(args.npz).resolve())}
        Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out_json).write_text(json.dumps(out, indent=2))
        print(f"→ wrote params {args.out_json}")


if __name__ == "__main__":
    main()
