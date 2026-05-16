#!/usr/bin/env python3
"""Summarize Pt(Au) HMC scan JSON files into one CSV table."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


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
    parser.add_argument("--out-csv", required=True)
    args = parser.parse_args()

    rows = []
    for raw in args.summaries:
        path = Path(raw)
        data = json.loads(path.read_text())
        fd_pred = float(data.get("fd_pred", math.nan))
        x_gb = float(data["tail_mean_x_gb"])
        row = {
            "summary": str(path),
            "T_K": _infer_temp(path, data),
            "X_total_target": _infer_xc(path, data),
            "X_total_actual": float(data.get("final_x_total", data.get("xc", math.nan))),
            "X_GB_HMC": x_gb,
            "X_GB_HMC_std": float(data.get("tail_std_x_gb", math.nan)),
            "X_bulk_HMC": float(data.get("tail_mean_x_bulk", math.nan)),
            "X_bulk_HMC_std": float(data.get("tail_std_x_bulk", math.nan)),
            "X_GB_FD": fd_pred,
            "HMC_minus_FD": x_gb - fd_pred if math.isfinite(fd_pred) else math.nan,
            "HMC_over_FD": x_gb / fd_pred if fd_pred > 0 else math.nan,
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
            f"X_bulk={row['X_bulk_HMC']:.4f}, "
            f"acc={row['hmc_acceptance']:.3f}"
        )


if __name__ == "__main__":
    main()
