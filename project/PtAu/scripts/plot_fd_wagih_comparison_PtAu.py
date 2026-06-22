#!/usr/bin/env python3
"""Plot our Pt(Au) FD curves against Wagih's Pt(Au) FD baseline."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from fermi_dirac_predict_PtAu import (
    load_ours,
    load_wagih_dump,
    x_gb_canonical_curve,
    x_gb_curve,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ours-npz", required=True)
    parser.add_argument("--wagih-dump", required=True)
    parser.add_argument("--temp", type=float, default=700.0)
    parser.add_argument("--Xc-min", type=float, default=1e-4)
    parser.add_argument("--Xc-max", type=float, default=0.2)
    parser.add_argument("--Xc-points", type=int, default=180)
    parser.add_argument("--n-total", type=int, default=62096)
    parser.add_argument("--n-gb", type=int, default=23272)
    parser.add_argument("--out-png", required=True)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    dE_ours = load_ours(Path(args.ours_npz))
    dE_wagih = load_wagih_dump(Path(args.wagih_dump))
    x = np.logspace(np.log10(args.Xc_min), np.log10(args.Xc_max), args.Xc_points)

    ours_reservoir = x_gb_curve(dE_ours, args.temp, x)
    wagih_reservoir = x_gb_curve(dE_wagih, args.temp, x)
    ours_closed, ours_closed_bulk = x_gb_canonical_curve(
        dE_ours,
        args.temp,
        x,
        n_total=args.n_total,
        n_gb=args.n_gb,
    )

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "X",
                "ours_reservoir_FD",
                "wagih_reservoir_FD",
                "ours_closed_box_FD",
                "ours_closed_box_bulk_X",
            ]
        )
        writer.writerows(zip(x, ours_reservoir, wagih_reservoir, ours_closed, ours_closed_bulk))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7.0, 4.8), dpi=160)
    ax.plot(x, wagih_reservoir, color="#4c6f9f", lw=2.0, ls="--", label="Wagih reservoir FD")
    ax.plot(x, ours_reservoir, color="#2f5d50", lw=2.0, label="ours reservoir FD")
    ax.plot(x, ours_closed, color="#9b3d2e", lw=2.0, label="ours closed-box FD")
    ax.plot([args.Xc_min, args.Xc_max], [args.Xc_min, args.Xc_max], ":", color="#777777", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlabel(r"composition variable $X$")
    ax.set_ylabel(r"FD GB Au fraction $X_\mathrm{GB}^\mathrm{FD}$")
    ax.set_title(f"Pt(Au) FD comparison at {args.temp:g} K")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()

    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)

    print(f"Wrote {out_png}")
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
