"""Compare our N_GB=500 production Pt(Au) ΔE spectrum against Wagih's
accelerated_model database dump (same O'Brien 2017 PtAu.eam.alloy potential).

Adapted from project/scripts/compare_vs_wagih.py (Al-Mg). The Al-Mg path
loaded a pre-extracted (seg_energies_Al_Mg.txt + bulk_solute_Al_Mg.dat)
pair from machine_learning_notebook/. The Pt(Au) reference lives in the
accelerated_model dump at:

  segregation_spectra_database_accelerated_model/Pt/
    Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/
      Pt_Au_20nm_GB_segregation.dump

Format: LAMMPS dump with columns `id type x y z seg_kJ_per_mol`. The
`seg_kJ_per_mol` value is the per-site ΔE_seg already in kJ/mol. Bulk
atoms are stored with seg_kJ_per_mol = 0.0; GB atoms have non-zero
predictions from Wagih's accelerated SOAP-based model. The KS sample
is the non-zero rows.

Inputs:
  --ours-npz    .npz from sample_delta_e_PtAu.py  (key `gb_delta_e` in eV)
  --wagih-dump  the Pt_Au_20nm_GB_segregation.dump above

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
    return np.asarray(d["gb_delta_e"]) * EV_TO_KJMOL


def load_wagih_dump(dump_path: Path) -> np.ndarray:
    """Read Wagih accelerated_model dump and return per-GB-site ΔE in kJ/mol.

    Returns the non-zero `seg_kJ_per_mol` column (bulk atoms are stored
    as exactly 0.0 — they are not GB sites).
    """
    with open(dump_path) as f:
        lines = f.readlines()
    n_atoms = 0
    cols = None
    data_start = 0
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            n_atoms = int(lines[i + 1])
            i += 2
        elif line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            data_start = i + 1
            break
        else:
            i += 1
    if cols is None or "seg_kJ_per_mol" not in cols:
        raise ValueError(f"unexpected dump format; columns = {cols}")
    seg_idx = cols.index("seg_kJ_per_mol")
    data = np.loadtxt(lines[data_start:data_start + n_atoms])
    seg = data[:, seg_idx]
    gb = seg[seg != 0.0]
    return gb


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
    p.add_argument("--wagih-dump", required=True)
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    dE_ours = load_ours(Path(args.ours_npz))
    dE_wagih = load_wagih_dump(Path(args.wagih_dump))

    m_ours, m_wagih = moments(dE_ours), moments(dE_wagih)
    f_ours, f_wagih = fit(dE_ours), fit(dE_wagih)

    ks = stats.ks_2samp(dE_ours, dE_wagih)
    print("Sample sizes: ours n =", m_ours["n"], " Wagih n =", m_wagih["n"])
    print("\nOur sample moments (kJ/mol):", json.dumps(m_ours, indent=2))
    print("\nWagih sample moments (kJ/mol):", json.dumps(m_wagih, indent=2))
    print("\nOur skew-normal fit: μ={mu:+.2f} σ={sigma:.2f} α={alpha:+.2f}".format(**f_ours))
    print("Wagih skew-normal fit: μ={mu:+.2f} σ={sigma:.2f} α={alpha:+.2f}".format(**f_wagih))
    print(f"\nKS statistic: D = {ks.statistic:.4f},  p-value = {ks.pvalue:.4g}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.0, 4.5), dpi=130)
    bins = np.linspace(min(dE_ours.min(), dE_wagih.min()) - 2,
                       max(dE_ours.max(), dE_wagih.max()) + 2, 60)
    ax.hist(dE_wagih, bins=bins, density=True, color="#55a868", alpha=0.55,
            edgecolor="white", label=f"Wagih Zenodo (n={m_wagih['n']})")
    ax.hist(dE_ours, bins=bins, density=True, color="#c44e52", alpha=0.55,
            edgecolor="white", label=f"Ours polycrystal (n={m_ours['n']})")
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
    ax.set_xlabel(r"$\Delta E_\mathrm{seg}$  (kJ/mol)")
    ax.set_ylabel("probability density")
    ax.set_title(r"Pt(Au) $\Delta E_\mathrm{seg}$ — ours vs Wagih Zenodo "
                 + f"(KS D={ks.statistic:.3f}, p={ks.pvalue:.2g})", pad=14)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.95)
    fig.tight_layout()
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_png)
    plt.close(fig)

    out = {"ours": {"moments_kjmol": m_ours, "skewnorm_fit_kjmol": f_ours},
           "wagih": {"moments_kjmol": m_wagih, "skewnorm_fit_kjmol": f_wagih},
           "ks_two_sample": {"D": float(ks.statistic), "p": float(ks.pvalue)}}
    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"\n→ wrote plot   {args.out_png}")
    print(f"→ wrote params {args.out_json}")


if __name__ == "__main__":
    main()
