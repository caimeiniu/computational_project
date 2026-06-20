"""Bootstrap CI for our N=500 ΔE_seg statistics vs Wagih's 82,646 pool.

Draws B sub-samples of size n=500 (with replacement) from Wagih's pool,
computes 6 statistics each (mean/std/skew + skewnorm μ/σ/α), and reports
where our observed statistics fall in the bootstrap distribution.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from scipy import stats

EV_TO_KJMOL = 96.485


def load_ours(npz_path: Path) -> np.ndarray:
    return np.asarray(np.load(npz_path)["gb_delta_e"]) * EV_TO_KJMOL


def load_wagih(seg_txt: Path, bulk_dat: Path) -> np.ndarray:
    e_bulk = float(open(bulk_dat).read().strip().split()[0])
    out = []
    with open(seg_txt) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                out.append(float(parts[1]) - e_bulk)
    return np.array(out) * EV_TO_KJMOL


def six_stats(x: np.ndarray) -> tuple[float, float, float, float, float, float]:
    a, loc, scale = stats.skewnorm.fit(x)
    return (float(x.mean()), float(x.std(ddof=1)), float(stats.skew(x)),
            float(loc), float(scale), float(a))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ours-npz", required=True)
    p.add_argument("--wagih-seg", required=True)
    p.add_argument("--wagih-bulk", required=True)
    p.add_argument("--out-json", required=True)
    p.add_argument("--B", type=int, default=10000)
    p.add_argument("--n", type=int, default=500)
    p.add_argument("--seed", type=int, default=20260425)
    args = p.parse_args()

    dE_ours = load_ours(Path(args.ours_npz))
    pool = load_wagih(Path(args.wagih_seg), Path(args.wagih_bulk))
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
