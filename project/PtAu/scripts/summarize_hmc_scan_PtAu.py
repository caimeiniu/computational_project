#!/usr/bin/env python3
"""Summarize Pt(Au) HMC scan JSON files into one CSV table."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

from fermi_dirac_predict_PtAu import load_ours, x_gb as x_gb_reservoir, x_gb_canonical


def _infer_xc(path: Path, data: dict) -> float:
    if "xc" in data:
        return float(data["xc"])
    match = re.search(r"Xc([0-9]+(?:\.[0-9]+)?)", path.name)
    if not match:
        raise ValueError(f"Could not infer Xc from {path}")
    return float(match.group(1))


def _infer_temp(path: Path, data: dict) -> float:
    if "temp_K" in data:
        return float(data["temp_K"])
    match = re.search(r"T([0-9]+(?:\.[0-9]+)?)", path.name)
    if not match:
        raise ValueError(f"Could not infer T from {path}")
    return float(match.group(1))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "summaries",
        nargs="+",
        help="HMC *_xgb_summary.json files.",
    )
    parser.add_argument("--deltae-npz", help="Optional Pt(Au) ΔE NPZ for FD comparison.")
    parser.add_argument("--n-total", type=int, default=62096)
    parser.add_argument("--n-gb", type=int, default=23272)
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    dE = load_ours(Path(args.deltae_npz)) if args.deltae_npz else None

    rows = []
    for raw in args.summaries:
        path = Path(raw)
        data = json.loads(path.read_text())
        fd_pred = float(data.get("fd_pred", math.nan))
        x_gb_hmc = float(data["tail_mean_x_gb"])
        temp = _infer_temp(path, data)
        x_total = _infer_xc(path, data)
        fd_reservoir = math.nan
        fd_closed = math.nan
        fd_bulk = math.nan
        if dE is not None:
            fd_reservoir = x_gb_reservoir(dE, temp, x_total)
            fd_closed, fd_bulk = x_gb_canonical(
                dE,
                temp,
                x_total,
                n_total=args.n_total,
                n_gb=args.n_gb,
            )
            fd_pred = fd_closed
        row = {
            "summary": str(path),
            "T_K": temp,
            "X_total_target": x_total,
            "X_total_actual": float(data.get("final_x_total", data.get("xc", math.nan))),
            "X_GB_HMC": x_gb_hmc,
            "X_GB_HMC_std": float(data.get("tail_std_x_gb", math.nan)),
            "X_bulk_HMC": float(data.get("tail_mean_x_bulk", math.nan)),
            "X_bulk_HMC_std": float(data.get("tail_std_x_bulk", math.nan)),
            "X_GB_FD_reservoir": fd_reservoir,
            "X_GB_FD_closed": fd_pred,
            "X_bulk_FD_closed": fd_bulk,
            "HMC_minus_FD_closed": x_gb_hmc - fd_pred if math.isfinite(fd_pred) else math.nan,
            "HMC_over_FD_closed": x_gb_hmc / fd_pred if fd_pred > 0 else math.nan,
            "final_step": int(data.get("final_step", -1)),
            "n_frames": int(data.get("n_frames", -1)),
            "hmc_acceptance": float(data.get("hmc_acceptance_last_log", math.nan)),
        }
        rows.append(row)

    rows.sort(key=lambda r: (r["T_K"], r["X_total_target"]))
    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {out}")
    for row in rows:
        print(
            f"T={row['T_K']:g} K X={row['X_total_target']:.4g}: "
            f"X_GB={row['X_GB_HMC']:.4f}, "
            f"FD_closed={row['X_GB_FD_closed']:.4f}, "
            f"delta={row['HMC_minus_FD_closed']:+.4f}, "
            f"acc={row['hmc_acceptance']:.3f}"
        )


if __name__ == "__main__":
    main()
