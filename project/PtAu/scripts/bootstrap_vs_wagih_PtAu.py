"""Bootstrap CI for our N=500 Pt(Au) ΔE_seg statistics vs Wagih's
accelerated_model dump pool (same O'Brien 2017 PtAu.eam.alloy potential).

Adapted from project/scripts/bootstrap_vs_wagih.py (Al-Mg). Loads the
Wagih reference from the dump file's `seg_kJ_per_mol` column rather
than the pre-extracted (seg.txt + bulk.dat) pair used for Al-Mg.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy import stats

EV_TO_KJMOL = 96.485


def load_ours(npz_path: Path) -> np.ndarray:
    return np.asarray(np.load(npz_path)["gb_delta_e"]) * EV_TO_KJMOL


def load_wagih_dump(dump_path: Path) -> np.ndarray:
    with open(dump_path) as f:
        lines = f.readlines()
    n_atoms, cols, data_start = 0, None, 0
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
    return seg[seg != 0.0]


def six_stats(x: np.ndarray) -> tuple[float, float, float, float, float, float]:
    a, loc, scale = stats.skewnorm.fit(x)
    return (float(x.mean()), float(x.std(ddof=1)), float(stats.skew(x)),
            float(loc), float(scale), float(a))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ours-npz", required=True)
    p.add_argument("--wagih-dump", required=True)
    p.add_argument("--out-json", required=True)
    p.add_argument("--B", type=int, default=10000)
    p.add_argument("--n", type=int, default=500)
    p.add_argument("--seed", type=int, default=20260513)
    args = p.parse_args()

    dE_ours = load_ours(Path(args.ours_npz))
    pool = load_wagih_dump(Path(args.wagih_dump))
    print(f"Ours n={dE_ours.size}; Wagih pool n={pool.size}")

    obs = six_stats(dE_ours)
    rng = np.random.default_rng(args.seed)
    samples = np.empty((args.B, 6))
    for b in range(args.B):
        idx = rng.integers(0, pool.size, size=args.n)
        samples[b] = six_stats(pool[idx])
        if (b + 1) % 1000 == 0:
            print(f"  bootstrap {b+1}/{args.B}")

    names = ["sample_mean", "sample_std", "sample_skew",
             "skewnorm_mu", "skewnorm_sigma", "skewnorm_alpha"]
    out = {"n_per_sample": args.n, "B": args.B, "seed": args.seed,
           "observed_ours": dict(zip(names, obs)), "stats": {}}

    print(f"\n{'stat':<18}{'boot_mean':>12}{'boot_std':>10}"
          f"{'CI_2.5':>10}{'CI_97.5':>10}{'ours':>10}{'z':>8}{'pct':>8}{'inside?':>10}")
    for i, name in enumerate(names):
        col = samples[:, i]
        m, s = float(col.mean()), float(col.std(ddof=1))
        lo, hi = (float(np.percentile(col, 2.5)),
                  float(np.percentile(col, 97.5)))
        o = obs[i]
        z = (o - m) / s if s > 0 else float("nan")
        pct = float((col < o).mean() * 100.0)
        inside = lo <= o <= hi
        out["stats"][name] = {"boot_mean": m, "boot_std": s,
                              "CI95": [lo, hi], "ours": o,
                              "z": z, "percentile": pct, "inside95": inside}
        print(f"{name:<18}{m:>12.3f}{s:>10.3f}{lo:>10.3f}{hi:>10.3f}"
              f"{o:>10.3f}{z:>+8.2f}{pct:>8.1f}{'yes' if inside else 'NO':>10}")

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(out, indent=2))
    print(f"\n→ wrote {args.out_json}")


if __name__ == "__main__":
    main()
