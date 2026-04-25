"""Pipeline-residual diagnostic: pair pe_ours(i) against Wagih's E_GB(i)
on the same atom_id, on the same annealed structure (Wagih's dump0).

Question answered: does our Phase-3 substitution + CG protocol reproduce
Wagih's per-site PE atom-for-atom on the SAME structure? A near-zero
residual means our pipeline is correct and any spectrum-level shift
between us and Wagih comes from structure-realization variance, not
pipeline differences.

Inputs:
  --ours        either:
                  *.csv  (partial, streamed from sample_delta_e.py),
                  *.npz  (final, keys `gb_e_mg` + `gb_site_ids`)
  --wagih-seg   learning_segregation_energies/.../seg_energies_Al_Mg.txt
  --out-png     scatter pe_ours vs E_GB_wagih + residual histogram
  --out-json    summary stats (mean, std, min, max in eV / meV / kJ/mol)
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np

EV_TO_KJMOL = 96.485
EV_TO_MEV   = 1000.0


def load_ours_csv(path: Path):
    sids, pes = [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if int(row["is_gb"]) != 1:
                continue
            sids.append(int(row["site_id"]))
            pes.append(float(row["pe"]))
    return np.array(sids, dtype=np.int64), np.array(pes, dtype=float)


def load_ours_npz(path: Path):
    d = np.load(path)
    sids = np.asarray(d["gb_site_ids"], dtype=np.int64)
    pes = np.asarray(d["gb_e_mg"], dtype=float)
    return sids, pes


def load_wagih_segenergies(path: Path):
    sids, egbs = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            sids.append(int(parts[0]))
            egbs.append(float(parts[1]))
    return np.array(sids, dtype=np.int64), np.array(egbs, dtype=float)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ours", required=True, help=".csv (streamed) or .npz (final)")
    p.add_argument("--wagih-seg", required=True)
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    ours_path = Path(args.ours)
    if ours_path.suffix == ".csv":
        sid_o, pe_o = load_ours_csv(ours_path)
        ours_kind = "csv (partial stream)"
    elif ours_path.suffix == ".npz":
        sid_o, pe_o = load_ours_npz(ours_path)
        ours_kind = "npz (final)"
    else:
        raise SystemExit(f"Unsupported --ours suffix: {ours_path.suffix}")

    sid_w, e_gb_w = load_wagih_segenergies(Path(args.wagih_seg))
    sid_to_egb = dict(zip(sid_w.tolist(), e_gb_w.tolist()))

    pe_w = np.array([sid_to_egb[int(s)] for s in sid_o], dtype=float)
    resid = pe_o - pe_w   # eV

    summary = {
        "ours_source": str(ours_path),
        "ours_kind": ours_kind,
        "wagih_seg": str(args.wagih_seg),
        "n_paired": int(sid_o.size),
        "residual_pe_eV": {
            "mean":   float(resid.mean()),
            "std":    float(resid.std(ddof=1)),
            "median": float(np.median(resid)),
            "min":    float(resid.min()),
            "max":    float(resid.max()),
        },
        "residual_pe_meV": {
            "mean":   float(resid.mean() * EV_TO_MEV),
            "std":    float(resid.std(ddof=1) * EV_TO_MEV),
            "min":    float(resid.min() * EV_TO_MEV),
            "max":    float(resid.max() * EV_TO_MEV),
        },
        "residual_pe_kJmol": {
            "mean":   float(resid.mean() * EV_TO_KJMOL),
            "std":    float(resid.std(ddof=1) * EV_TO_KJMOL),
            "min":    float(resid.min() * EV_TO_KJMOL),
            "max":    float(resid.max() * EV_TO_KJMOL),
        },
        "pearson_r": float(np.corrcoef(pe_o, pe_w)[0, 1]),
    }

    print(json.dumps(summary, indent=2))

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(summary, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), dpi=130)

    ax = axes[0]
    ax.scatter(pe_w, pe_o, s=8, alpha=0.5, color="C0")
    lo, hi = min(pe_w.min(), pe_o.min()), max(pe_w.max(), pe_o.max())
    ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="y = x")
    ax.set_xlabel("Wagih E_GB^Mg (eV)")
    ax.set_ylabel("Ours pe (eV)")
    ax.set_title(f"Paired per-site PE  (n = {sid_o.size})\n"
                 f"r = {summary['pearson_r']:.6f}")
    ax.legend(loc="upper left", fontsize=9)

    ax = axes[1]
    resid_meV = resid * EV_TO_MEV
    ax.hist(resid_meV, bins=40, color="C1", alpha=0.85)
    ax.axvline(resid_meV.mean(), color="k", lw=1.0,
               label=f"mean = {resid_meV.mean():+.2f} meV")
    ax.set_xlabel("Residual pe_ours − E_GB_wagih (meV)")
    ax.set_ylabel("count")
    ax.set_title(f"Residual distribution  "
                 f"(σ = {resid_meV.std(ddof=1):.2f} meV)")
    ax.legend(loc="upper right", fontsize=9)

    fig.tight_layout()
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_png)
    print(f"\nWrote {args.out_png}")
    print(f"Wrote {args.out_json}")


if __name__ == "__main__":
    main()
