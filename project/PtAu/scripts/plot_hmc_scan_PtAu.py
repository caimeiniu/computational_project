#!/usr/bin/env python3
"""Plot Pt(Au) HMC scan against closed-box Fermi-Dirac baseline."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from fermi_dirac_predict_PtAu import load_wagih_dump, x_gb_curve


def _load_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open() as fh:
        for row in csv.DictReader(fh):
            parsed: dict[str, float | str] = {}
            for key, value in row.items():
                parsed[key] = value if key == "summary" else float(value)
            rows.append(parsed)
    rows.sort(key=lambda row: float(row["X_total_target"]))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-csv", required=True, help="CSV from summarize_hmc_scan_PtAu.py")
    parser.add_argument("--out-png", required=True)
    parser.add_argument("--out-csv", help="Optional compact plotting table.")
    parser.add_argument("--title", default="Pt(Au) HMC vs closed-box FD at 700 K")
    parser.add_argument("--wagih-dump", help="Optional Wagih Pt(Au) dump for reservoir FD overlay.")
    parser.add_argument("--temp", type=float, default=700.0)
    parser.add_argument(
        "--slow-mixing-x",
        type=float,
        default=0.10,
        help="Mark this concentration as a slow-mixing/stress point.",
    )
    args = parser.parse_args()

    rows = _load_rows(Path(args.scan_csv))
    x = np.array([row["X_total_target"] for row in rows], dtype=float)
    hmc = np.array([row["X_GB_HMC"] for row in rows], dtype=float)
    hmc_std = np.array([row["X_GB_HMC_std"] for row in rows], dtype=float)
    fd = np.array([row["X_GB_FD_closed"] for row in rows], dtype=float)
    delta = hmc - fd
    ratio = hmc / fd
    slow = np.isclose(x, args.slow_mixing_x, rtol=0.0, atol=5e-5)
    main = ~slow
    wagih_x = None
    wagih_fd = None
    wagih_fd_at_hmc = None
    if args.wagih_dump:
        dE_wagih = load_wagih_dump(Path(args.wagih_dump))
        wagih_x = np.logspace(np.log10(max(float(x.min()), 1e-5)), np.log10(float(x.max())), 160)
        wagih_fd = x_gb_curve(dE_wagih, args.temp, wagih_x)
        wagih_fd_at_hmc = x_gb_curve(dE_wagih, args.temp, x)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax, ax_delta) = plt.subplots(
        2,
        1,
        figsize=(6.8, 6.4),
        dpi=160,
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )

    if wagih_x is not None and wagih_fd is not None:
        ax.plot(
            wagih_x,
            wagih_fd,
            color="#4c6f9f",
            lw=1.8,
            ls="--",
            label="Wagih reservoir FD",
        )
    ax.plot(x, fd, color="#2f5d50", lw=2.0, label="closed-box FD")
    ax.errorbar(
        x[main],
        hmc[main],
        yerr=hmc_std[main],
        color="#9b3d2e",
        marker="o",
        ms=4.5,
        lw=1.8,
        capsize=2.5,
        label="HMC",
    )
    if np.any(slow):
        ax.errorbar(
            x[slow],
            hmc[slow],
            yerr=hmc_std[slow],
            color="#777777",
            marker="o",
            mfc="white",
            mec="#777777",
            ms=5.5,
            lw=0.0,
            capsize=2.5,
            label="HMC slow-mixing check",
        )
        ax.annotate(
            "slow mixing",
            xy=(x[slow][0], hmc[slow][0]),
            xytext=(-58, -18),
            textcoords="offset points",
            arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "#666666"},
            fontsize=8,
            color="#555555",
        )
    ax.plot(x, x, ":", color="#777777", lw=1.0, label=r"$X_\mathrm{GB}=X_\mathrm{total}$")
    ax.set_ylabel(r"GB Au fraction $X_\mathrm{GB}$")
    ax.set_title(args.title)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, loc="upper left")

    ax_delta.axhline(0.0, color="#555555", lw=1.0)
    ax_delta.plot(x[main], delta[main], color="#5b4b8a", marker="s", ms=4.0, lw=1.5)
    if np.any(slow):
        ax_delta.plot(
            x[slow],
            delta[slow],
            color="#777777",
            marker="s",
            mfc="white",
            mec="#777777",
            ms=5.0,
            lw=0.0,
        )
    ax_delta.axvspan(0.015, 0.02, color="#d8b45a", alpha=0.22, lw=0)
    ax_delta.text(
        0.0175,
        max(delta) * 0.82,
        "onset",
        ha="center",
        va="center",
        fontsize=8,
        color="#5d4b16",
    )
    ax_delta.set_xlabel(r"total Au fraction $X_\mathrm{total}$")
    ax_delta.set_ylabel(r"HMC $-$ FD")
    ax_delta.grid(True, alpha=0.25)

    fig.tight_layout()
    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)

    if args.out_csv:
        out_csv = Path(args.out_csv)
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        with out_csv.open("w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=[
                    "X_total",
                    "X_GB_HMC",
                    "X_GB_FD_closed",
                    "X_GB_Wagih_reservoir_FD",
                    "HMC_minus_FD_closed",
                    "HMC_over_FD_closed",
                ],
            )
            writer.writeheader()
            if wagih_fd_at_hmc is None:
                wagih_fd_iter = [""] * len(x)
            else:
                wagih_fd_iter = [float(v) for v in wagih_fd_at_hmc]
            for xi, hi, fi, wi, di, ri in zip(x, hmc, fd, wagih_fd_iter, delta, ratio):
                wagih_at_x = ""
                if args.wagih_dump:
                    wagih_at_x = wi
                writer.writerow(
                    {
                        "X_total": xi,
                        "X_GB_HMC": hi,
                        "X_GB_FD_closed": fi,
                        "X_GB_Wagih_reservoir_FD": wagih_at_x,
                        "HMC_minus_FD_closed": di,
                        "HMC_over_FD_closed": ri,
                    }
                )

    print(f"Wrote {out_png}")
    if args.out_csv:
        print(f"Wrote {args.out_csv}")
    print("Breakdown summary: agreement through X=0.01, onset at X=0.015-0.02, clear by X=0.03.")


if __name__ == "__main__":
    main()
