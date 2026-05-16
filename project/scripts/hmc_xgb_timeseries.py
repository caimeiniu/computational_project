#!/usr/bin/env python3
"""Compute GB and bulk solute fractions from one or more HMC dump streams."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Post-process LAMMPS atom/swap dumps into X_GB(t)."
    )
    parser.add_argument(
        "--stub",
        action="append",
        required=True,
        help="Output stub without .dump/.log. May be supplied more than once.",
    )
    parser.add_argument("--gb-mask", required=True, help="Boolean .npy mask, id-1 indexed.")
    parser.add_argument("--xc", type=float, required=True, help="Nominal global solute fraction.")
    parser.add_argument("--temp", type=float, required=True, help="Run temperature in K.")
    parser.add_argument("--fd-pred", type=float, default=math.nan, help="FD X_GB prediction.")
    parser.add_argument("--solute-type", type=int, default=2, help="LAMMPS type id for solute.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix for csv/json/png.")
    parser.add_argument(
        "--tail-frames",
        type=int,
        default=50,
        help="Number of final frames used for plateau averages.",
    )
    return parser.parse_args()


def _read_numeric_log_rows(log_path: Path) -> list[dict[str, float]]:
    if not log_path.exists():
        return []

    rows: list[dict[str, float]] = []
    header: list[str] | None = None
    for raw in log_path.read_text(errors="replace").splitlines():
        parts = raw.split()
        if not parts:
            continue
        if parts[0] == "Step" and "Temp" in parts:
            header = parts
            continue
        if header is None or len(parts) != len(header):
            continue
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        rows.append(dict(zip(header, vals)))
    return rows


def dump_frames(dump_path: Path) -> Iterable[tuple[int, np.ndarray, np.ndarray]]:
    """Yield (step, ids, types) from a simple LAMMPS custom dump."""
    with dump_path.open("r", errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line:
                return
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            step = int(fh.readline().strip())

            marker = fh.readline()
            if not marker.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"{dump_path}: expected atom count after step {step}")
            n_atoms = int(fh.readline().strip())

            marker = fh.readline()
            if not marker.startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"{dump_path}: expected box bounds after step {step}")
            for _ in range(3):
                fh.readline()

            marker = fh.readline().split()
            if marker[:2] != ["ITEM:", "ATOMS"]:
                raise ValueError(f"{dump_path}: expected atom table after step {step}")
            columns = marker[2:]
            try:
                id_col = columns.index("id")
                type_col = columns.index("type")
            except ValueError as exc:
                raise ValueError(f"{dump_path}: dump must contain id and type columns") from exc

            ids = np.empty(n_atoms, dtype=np.int64)
            types = np.empty(n_atoms, dtype=np.int16)
            for i in range(n_atoms):
                parts = fh.readline().split()
                ids[i] = int(parts[id_col])
                types[i] = int(parts[type_col])
            yield step, ids, types


def summarize_dump(stub: Path, gb_mask: np.ndarray, solute_type: int) -> list[dict[str, float]]:
    dump_path = Path(str(stub) + ".dump")
    if not dump_path.exists():
        raise FileNotFoundError(f"Missing dump: {dump_path}")

    rows: list[dict[str, float]] = []
    gb_mask = np.asarray(gb_mask, dtype=bool)
    for frame_i, (step, ids, types) in enumerate(dump_frames(dump_path)):
        zero_based = ids - 1
        if zero_based.min() < 0 or zero_based.max() >= len(gb_mask):
            raise ValueError(
                f"{dump_path}: atom ids exceed gb mask length {len(gb_mask)} at step {step}"
            )

        is_gb = gb_mask[zero_based]
        is_solute = types == solute_type
        n_total = int(types.size)
        n_gb = int(is_gb.sum())
        n_solute = int(is_solute.sum())
        n_gb_solute = int((is_gb & is_solute).sum())
        n_bulk = n_total - n_gb
        n_bulk_solute = n_solute - n_gb_solute

        rows.append(
            {
                "stub": str(stub),
                "frame": frame_i,
                "step": int(step),
                "n_total": n_total,
                "n_solute": n_solute,
                "x_total": n_solute / n_total,
                "n_gb": n_gb,
                "n_gb_solute": n_gb_solute,
                "x_gb": n_gb_solute / n_gb if n_gb else math.nan,
                "n_bulk": n_bulk,
                "n_bulk_solute": n_bulk_solute,
                "x_bulk": n_bulk_solute / n_bulk if n_bulk else math.nan,
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "stub",
        "frame",
        "step",
        "n_total",
        "n_solute",
        "x_total",
        "n_gb",
        "n_gb_solute",
        "x_gb",
        "n_bulk",
        "n_bulk_solute",
        "x_bulk",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def maybe_plot(path: Path, rows: list[dict[str, float]], fd_pred: float) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    steps = np.array([row["step"] for row in rows], dtype=float)
    x_gb = np.array([row["x_gb"] for row in rows], dtype=float)
    x_bulk = np.array([row["x_bulk"] for row in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(steps, x_gb, label=r"$X_{GB}$", lw=1.8)
    ax.plot(steps, x_bulk, label=r"$X_{bulk}$", lw=1.2, alpha=0.75)
    if math.isfinite(fd_pred):
        ax.axhline(fd_pred, color="black", ls="--", lw=1.1, label=rf"FD {fd_pred:.4f}")
    ax.set_xlabel("LAMMPS step")
    ax.set_ylabel("Au fraction")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    gb_mask = np.load(args.gb_mask)
    stubs = [Path(s) for s in args.stub]

    rows: list[dict[str, float]] = []
    log_rows: list[dict[str, float]] = []
    for stub in stubs:
        rows.extend(summarize_dump(stub, gb_mask, args.solute_type))
        log_rows.extend(_read_numeric_log_rows(Path(str(stub) + ".log")))

    rows.sort(key=lambda r: (r["step"], r["stub"], r["frame"]))
    if not rows:
        raise SystemExit("No dump frames found.")

    out_prefix = Path(args.out_prefix)
    csv_path = out_prefix.with_name(out_prefix.name + "_xgb_timeseries.csv")
    json_path = out_prefix.with_name(out_prefix.name + "_xgb_summary.json")
    png_path = out_prefix.with_name(out_prefix.name + "_xgb_timeseries.png")

    write_csv(csv_path, rows)

    tail = rows[-min(args.tail_frames, len(rows)) :]
    x_gb_tail = np.array([r["x_gb"] for r in tail], dtype=float)
    x_bulk_tail = np.array([r["x_bulk"] for r in tail], dtype=float)
    final = rows[-1]

    summary = {
        "stubs": [str(s) for s in stubs],
        "gb_mask": str(Path(args.gb_mask)),
        "temp_K": args.temp,
        "xc": args.xc,
        "fd_pred": args.fd_pred,
        "n_frames": len(rows),
        "tail_frames": len(tail),
        "final_step": final["step"],
        "final_x_total": final["x_total"],
        "final_x_gb": final["x_gb"],
        "final_x_bulk": final["x_bulk"],
        "tail_mean_x_gb": float(np.nanmean(x_gb_tail)),
        "tail_std_x_gb": float(np.nanstd(x_gb_tail, ddof=1)) if len(x_gb_tail) > 1 else 0.0,
        "tail_mean_x_bulk": float(np.nanmean(x_bulk_tail)),
        "tail_std_x_bulk": float(np.nanstd(x_bulk_tail, ddof=1)) if len(x_bulk_tail) > 1 else 0.0,
        "tail_minus_fd": float(np.nanmean(x_gb_tail) - args.fd_pred)
        if math.isfinite(args.fd_pred)
        else math.nan,
    }

    if log_rows and "f_hmc[1]" in log_rows[-1] and "f_hmc[2]" in log_rows[-1]:
        attempts = log_rows[-1]["f_hmc[1]"]
        accepts = log_rows[-1]["f_hmc[2]"]
        summary["hmc_attempts_last_log"] = attempts
        summary["hmc_accepts_last_log"] = accepts
        summary["hmc_acceptance_last_log"] = accepts / attempts if attempts else math.nan

    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    maybe_plot(png_path, rows, args.fd_pred)

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    if png_path.exists():
        print(f"Wrote {png_path}")
    print(
        "Tail X_GB = "
        f"{summary['tail_mean_x_gb']:.4f} +/- {summary['tail_std_x_gb']:.4f}; "
        f"final X_GB = {summary['final_x_gb']:.4f}"
    )


if __name__ == "__main__":
    main()
