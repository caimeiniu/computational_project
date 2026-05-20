#!/usr/bin/env python3
"""Au-Au pair clustering on Pt(Au) GB compared with random GB occupancy."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from plot_hmc_site_diagnostics_PtAu import read_data_positions


def pair_distances(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    if positions.shape[0] < 2:
        return np.array([], dtype=float)
    dists: list[np.ndarray] = []
    for i in range(positions.shape[0] - 1):
        dr = positions[i + 1 :] - positions[i]
        dr -= box * np.round(dr / box)
        dists.append(np.sqrt(np.sum(dr * dr, axis=1)))
    return np.concatenate(dists)


def density_ratio_for_data(
    *,
    data_path: Path,
    gb_mask: np.ndarray,
    bins: np.ndarray,
    n_random: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    positions, types, box = read_data_positions(data_path)
    gb_ids = np.where(gb_mask)[0]
    gb_au_ids = gb_ids[types[gb_ids] == 2]
    x_gb = len(gb_au_ids) / len(gb_ids)

    hmc_hist, _ = np.histogram(pair_distances(positions[gb_au_ids], box), bins=bins)
    rng = np.random.default_rng(seed)
    rand_accum = np.zeros_like(hmc_hist, dtype=float)
    for _ in range(n_random):
        random_ids = rng.choice(gb_ids, size=len(gb_au_ids), replace=False)
        rand_hist, _ = np.histogram(pair_distances(positions[random_ids], box), bins=bins)
        rand_accum += rand_hist
    rand_mean = rand_accum / n_random
    ratio = np.divide(hmc_hist, rand_mean, out=np.full_like(rand_mean, np.nan), where=rand_mean > 0)
    centers = 0.5 * (bins[:-1] + bins[1:])
    return centers, ratio, x_gb, len(gb_au_ids) / positions.shape[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gb-mask", required=True)
    parser.add_argument("--out-png", required=True)
    parser.add_argument("--data", action="append", required=True, help="LAMMPS data file. Repeatable.")
    parser.add_argument("--label", action="append", required=True, help="Legend label. Repeatable.")
    parser.add_argument("--r-max", type=float, default=25.0)
    parser.add_argument("--dr", type=float, default=0.25)
    parser.add_argument("--n-random", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20260520)
    args = parser.parse_args()

    if len(args.data) != len(args.label):
        raise SystemExit("--data and --label must be supplied the same number of times")

    gb_mask = np.asarray(np.load(args.gb_mask), dtype=bool)
    bins = np.arange(0.0, args.r_max + args.dr, args.dr)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    colors = ["#e31a1c", "#1f78b4", "#33a02c", "#6a3d9a"]
    for i, (data, label) in enumerate(zip(args.data, args.label)):
        centers, ratio, x_gb, _ = density_ratio_for_data(
            data_path=Path(data),
            gb_mask=gb_mask,
            bins=bins,
            n_random=args.n_random,
            seed=args.seed + i,
        )
        ax.plot(centers, ratio, color=colors[i % len(colors)], lw=1.8, label=f"{label}, $X_{{GB}}={x_gb:.3f}$")

    ax.axhline(1.0, color="#777777", ls="--", lw=1.0)
    ax.axvspan(2.6, 3.1, color="#cccccc", alpha=0.45, label="1st NN shell")
    ax.set_xlabel(r"pair separation $r$ [Å]")
    ax.set_ylabel("Au-Au pair density / random reference")
    ax.set_title("Pt(Au) Au-Au spatial clustering on the GB")
    ax.set_xlim(0, args.r_max)
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()

    out = Path(args.out_png)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
