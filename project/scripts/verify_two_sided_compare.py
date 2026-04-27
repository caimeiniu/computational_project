#!/usr/bin/env python3
"""Two-sided IC equilibration verification overlay plot.

Reads two `*_xgb.json` files (random-IC and preseg-IC at the same T, X_c)
and produces a single PNG showing whether the two replicas have met at
the canonical-FD equilibrium target.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def split_half_drift(t, x, burnin_frac=0.2):
    n = len(t); k0 = int(burnin_frac * n)
    if n - k0 < 4:
        return None
    half = (n - k0) // 2
    return float(x[k0+half:].mean() - x[k0:k0+half].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rand-json", required=True, type=Path)
    ap.add_argument("--preseg-json", required=True, type=Path)
    ap.add_argument("--canon-fd", type=float, required=True,
                    help="canonical-FD equilibrium target X_GB")
    ap.add_argument("--ceiling", type=float, required=True,
                    help="closed-box ceiling X_c / f_gb")
    ap.add_argument("--xc", type=float, required=True)
    ap.add_argument("--temp", type=float, required=True)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    R = json.loads(args.rand_json.read_text())
    P = json.loads(args.preseg_json.read_text())

    tR  = np.array(R["series"]["timestep"]) * 1e-3
    xR  = np.array(R["series"]["x_gb"])
    tP  = np.array(P["series"]["timestep"]) * 1e-3
    xP  = np.array(P["series"]["x_gb"])

    drift_R = split_half_drift(tR, xR)
    drift_P = split_half_drift(tP, xP)

    fig, ax = plt.subplots(figsize=(8.5, 5.0))

    ax.plot(tR, xR, color="#1f77b4", lw=1.4,
            label=fr"random IC: $X_\mathrm{{GB}}={R['x_gb']['mean']:.4f}^{{+{R['x_gb']['ci95_hi']-R['x_gb']['mean']:.4f}}}_{{-{R['x_gb']['mean']-R['x_gb']['ci95_lo']:.4f}}}$  ($\Delta_{{1/2}}={drift_R:+.4f}$)")
    ax.fill_between(tR, R['x_gb']['ci95_lo'], R['x_gb']['ci95_hi'],
                    color="#1f77b4", alpha=0.10)

    ax.plot(tP, xP, color="#d62728", lw=1.4,
            label=fr"preseg IC: $X_\mathrm{{GB}}={P['x_gb']['mean']:.4f}^{{+{P['x_gb']['ci95_hi']-P['x_gb']['mean']:.4f}}}_{{-{P['x_gb']['mean']-P['x_gb']['ci95_lo']:.4f}}}$  ($\Delta_{{1/2}}={drift_P:+.4f}$)")
    ax.fill_between(tP, P['x_gb']['ci95_lo'], P['x_gb']['ci95_hi'],
                    color="#d62728", alpha=0.10)

    ax.axhline(args.canon_fd, color="black", ls="--", lw=1.2,
               label=fr"canon-FD target $={args.canon_fd:.3f}$")
    ax.axhline(args.ceiling, color="grey", ls=":", lw=1.0,
               label=fr"ceiling $X_c/f_\mathrm{{GB}}={args.ceiling:.3f}$")
    ax.axhline(args.xc, color="grey", ls=":", lw=1.0, alpha=0.5)
    ax.text(tR[-1], args.xc, f"  $X_c={args.xc:g}$", va="center",
            fontsize=8, color="grey")

    ax.set_xlabel("time (ps)")
    ax.set_ylabel(r"$X_\mathrm{GB}(t)$")
    ax.set_title(fr"Two-sided IC equilibration verification — T={args.temp:.0f} K, $X_c={args.xc:g}$"
                 f"\nrand wall-truncated at {tR[-1]:.0f} ps ({R['swap']['total_accepts']} swaps);"
                 f" preseg wall-truncated at {tP[-1]:.0f} ps ({P['swap']['total_accepts']} swaps)")
    ax.legend(loc="center right", fontsize=8.5)
    ax.set_ylim(0, max(args.ceiling*1.05, xP.max()*1.05))
    ax.grid(alpha=0.25)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    print(f"  rand    : {R['x_gb']['mean']:.4f} [{R['x_gb']['ci95_lo']:.4f}, {R['x_gb']['ci95_hi']:.4f}]"
          f"  half-half drift {drift_R:+.4f}")
    print(f"  preseg  : {P['x_gb']['mean']:.4f} [{P['x_gb']['ci95_lo']:.4f}, {P['x_gb']['ci95_hi']:.4f}]"
          f"  half-half drift {drift_P:+.4f}")
    print(f"  canon-FD: {args.canon_fd:.4f}")
    print(f"  bracket gap (preseg - rand): {P['x_gb']['mean'] - R['x_gb']['mean']:+.4f}")
    print(f"  preseg vs canon-FD gap     : {P['x_gb']['mean'] - args.canon_fd:+.4f}")
    print(f"  rand   vs canon-FD gap     : {R['x_gb']['mean'] - args.canon_fd:+.4f}")


if __name__ == "__main__":
    main()
