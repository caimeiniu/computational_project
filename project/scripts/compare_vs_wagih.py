"""Compare our N_GB=500 production Al(Mg) ΔE spectrum against Wagih's
~82k-site Zenodo dataset (same Mendelev 2009 Al-Mg.eam.fs potential).

Inputs:
  - ours:  .npz from sample_delta_e.py  (key `gb_delta_e` in eV,
           implicitly references our own E_bulk^Mg mean from the run)
  - wagih: learning_segregation_energies/machine_learning_notebook/
             seg_energies_Al_Mg.txt     (per-line `site_id E_GB^Mg` in eV)
           bulk_solute_Al_Mg.dat        (single value E_bulk^Mg in eV)

Outputs:
  stdout summary  + PNG overlay + JSON
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy import stats

EV_TO_KJMOL = 96.485


def load_ours(npz_path: Path) -> np.ndarray:
    d = np.load(npz_path)
    return np.asarray(d["gb_delta_e"]) * EV_TO_KJMOL  # already ΔE in eV → kJ/mol


def load_wagih(seg_txt: Path, bulk_dat: Path) -> np.ndarray:
    E_bulk_mg = float(open(bulk_dat).read().strip().split()[0])
    dE_eV = []
    with open(seg_txt) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            _sid, e_gb = parts
            dE_eV.append(float(e_gb) - E_bulk_mg)
    return np.array(dE_eV) * EV_TO_KJMOL


def moments(x: np.ndarray) -> dict:
    return {"n": int(x.size), "min": float(x.min()), "max": float(x.max()),
            "mean": float(x.mean()), "median": float(np.median(x)),
            "std": float(x.std(ddof=1)), "skew": float(stats.skew(x))}


def fit(x: np.ndarray) -> dict:
    a, loc, scale = stats.skewnorm.fit(x)
    return {"mu": float(loc), "sigma": float(scale), "alpha": float(a)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ours-npz", required=True)
    p.add_argument("--wagih-seg", required=True)
    p.add_argument("--wagih-bulk", required=True)
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    dE_ours = load_ours(Path(args.ours_npz))
    dE_wagih = load_wagih(Path(args.wagih_seg), Path(args.wagih_bulk))

    m_ours, m_wagih = moments(dE_ours), moments(dE_wagih)
    f_ours, f_wagih = fit(dE_ours), fit(dE_wagih)

    # KS two-sample test
    ks = stats.ks_2samp(dE_ours, dE_wagih)
    print("Sample sizes: ours n =", m_ours["n"], " Wagih n =", m_wagih["n"])
    print("\nOur sample moments (kJ/mol):", json.dumps(m_ours, indent=2))
    print("\nWagih sample moments (kJ/mol):", json.dumps(m_wagih, indent=2))
    print("\nOur skew-normal fit: μ={mu:+.2f} σ={sigma:.2f} α={alpha:+.2f}".format(**f_ours))
    print("Wagih skew-normal fit: μ={mu:+.2f} σ={sigma:.2f} α={alpha:+.2f}".format(**f_wagih))
    print(f"\nKS statistic: D = {ks.statistic:.4f},  p-value = {ks.pvalue:.4g}")

    # Plot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=130)
    bins = np.linspace(min(dE_ours.min(), dE_wagih.min()) - 2,
                        max(dE_ours.max(), dE_wagih.max()) + 2, 60)
    ax.hist(dE_wagih, bins=bins, density=True, color="#55a868", alpha=0.55,
            edgecolor="white", label=f"Wagih Zenodo (n={m_wagih['n']})")
    ax.hist(dE_ours, bins=bins, density=True, color="#c44e52", alpha=0.55,
            edgecolor="white", label=f"Ours production (n={m_ours['n']})")
    # overlay fits
    xs = np.linspace(bins[0], bins[-1], 400)
    ax.plot(xs, stats.skewnorm.pdf(xs, a=f_wagih["alpha"], loc=f_wagih["mu"],
            scale=f_wagih["sigma"]), "-", color="#2c662c", lw=1.8,
            label=(f"Wagih fit  μ={f_wagih['mu']:+.1f} σ={f_wagih['sigma']:.1f} "
                   f"α={f_wagih['alpha']:+.2f}"))
    ax.plot(xs, stats.skewnorm.pdf(xs, a=f_ours["alpha"], loc=f_ours["mu"],
            scale=f_ours["sigma"]), "--", color="#7a2222", lw=1.8,
            label=(f"Ours fit  μ={f_ours['mu']:+.1f} σ={f_ours['sigma']:.1f} "
                   f"α={f_ours['alpha']:+.2f}"))
    ax.axvline(0.0, color="k", lw=0.6, alpha=0.4)
    ax.set_xlabel("ΔE_seg  (kJ/mol)")
    ax.set_ylabel("probability density")
    ax.set_title(f"Al(Mg) ΔE_seg — ours vs Wagih Zenodo "
                 f"(KS D={ks.statistic:.3f}, p={ks.pvalue:.2g})")
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_png)
    plt.close(fig)

    # JSON
    out = {"ours": {"moments_kjmol": m_ours, "skewnorm_fit_kjmol": f_ours},
           "wagih": {"moments_kjmol": m_wagih, "skewnorm_fit_kjmol": f_wagih},
           "ks_two_sample": {"D": float(ks.statistic), "p": float(ks.pvalue)}}
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"\n→ wrote plot   {args.out_png}")
    print(f"→ wrote params {args.out_json}")


if __name__ == "__main__":
    main()
