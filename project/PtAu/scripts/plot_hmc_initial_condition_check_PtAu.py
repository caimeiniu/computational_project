#!/usr/bin/env python3
"""Plot random-start vs GB-seeded HMC X_GB(t) for one Pt(Au) condition."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def _load_series(path: Path) -> tuple[np.ndarray, np.ndarray]:
    steps: list[float] = []
    xgb: list[float] = []
    with path.open() as fh:
        for row in csv.DictReader(fh):
            steps.append(float(row["step"]))
            xgb.append(float(row["x_gb"]))
    return np.asarray(steps), np.asarray(xgb)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--random-csv", required=True)
    parser.add_argument("--gbseed-csv", required=True)
    parser.add_argument("--out-png", required=True)
    parser.add_argument("--title", default="Pt(Au) initial-condition convergence")
    args = parser.parse_args()

    step_r, xgb_r = _load_series(Path(args.random_csv))
    step_g, xgb_g = _load_series(Path(args.gbseed_csv))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=160)
    ax.plot(step_r, xgb_r, color="#9b3d2e", lw=1.7, label="random start")
    ax.plot(step_g, xgb_g, color="#2f5d50", lw=1.7, label="GB-seeded start")

    tail_r = float(np.mean(xgb_r[-min(50, len(xgb_r)) :]))
    tail_g = float(np.mean(xgb_g[-min(50, len(xgb_g)) :]))
    ax.axhline(tail_r, color="#9b3d2e", ls=":", lw=1.0)
    ax.axhline(tail_g, color="#2f5d50", ls=":", lw=1.0)
    ax.text(
        0.98,
        0.06,
        f"tail means: random {tail_r:.4f}, GB-seeded {tail_g:.4f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
    )

    ax.set_xlabel("LAMMPS step")
    ax.set_ylabel(r"GB Au fraction $X_\mathrm{GB}$")
    ax.set_title(args.title)
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
