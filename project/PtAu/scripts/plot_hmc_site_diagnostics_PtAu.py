#!/usr/bin/env python3
"""Site-resolved Pt(Au) HMC diagnostics similar to Cainiu's Al(Mg) plots."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from fermi_dirac_predict_PtAu import KB_EV_PER_K, x_gb_canonical

EV_TO_KJMOL = 96.485


def load_deltae(npz_path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.load(npz_path)
    return np.asarray(data["gb_site_ids"], dtype=np.int64), np.asarray(data["gb_delta_e"], dtype=float)


def read_data_positions(data_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lines = data_path.read_text().splitlines()
    n_atoms = None
    bounds: list[tuple[float, float]] = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[1] == "atoms":
            n_atoms = int(parts[0])
        elif len(parts) >= 4 and parts[2:4] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
            bounds.append((float(parts[0]), float(parts[1])))
    if n_atoms is None or len(bounds) != 3:
        raise ValueError(f"Could not parse atom count/bounds from {data_path}")

    ids = np.empty(n_atoms, dtype=np.int64)
    types = np.empty(n_atoms, dtype=np.int16)
    pos = np.empty((n_atoms, 3), dtype=float)
    in_atoms = False
    count = 0
    for line in lines:
        stripped = line.strip()
        if stripped == "Atoms" or stripped.startswith("Atoms "):
            in_atoms = True
            continue
        if not in_atoms:
            continue
        parts = stripped.split()
        if len(parts) < 5:
            continue
        try:
            atom_id = int(parts[0])
            atom_type = int(parts[1])
            xyz = [float(parts[2]), float(parts[3]), float(parts[4])]
        except ValueError:
            continue
        ids[count] = atom_id
        types[count] = atom_type
        pos[count] = xyz
        count += 1
    if count != n_atoms:
        raise ValueError(f"Atoms section has {count} rows, expected {n_atoms}")
    order = np.argsort(ids)
    box = np.array([hi - lo for lo, hi in bounds], dtype=float)
    return pos[order], types[order], box


def dump_type_frames(dump_path: Path, n_atoms: int) -> np.ndarray:
    frames: list[np.ndarray] = []
    with dump_path.open(errors="replace") as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            fh.readline()
            if not fh.readline().startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError(f"Bad dump format in {dump_path}")
            n = int(fh.readline())
            if n != n_atoms:
                raise ValueError(f"Dump atom count {n} != data atom count {n_atoms}")
            if not fh.readline().startswith("ITEM: BOX BOUNDS"):
                raise ValueError(f"Bad box section in {dump_path}")
            for _ in range(3):
                fh.readline()
            header = fh.readline().split()
            cols = header[2:]
            id_col = cols.index("id")
            type_col = cols.index("type")
            types = np.empty(n_atoms, dtype=np.int16)
            for _ in range(n_atoms):
                parts = fh.readline().split()
                atom_id = int(parts[id_col])
                types[atom_id - 1] = int(parts[type_col])
            frames.append(types)
    if not frames:
        raise ValueError(f"No frames found in {dump_path}")
    return np.stack(frames, axis=0)


def fd_prob(dE_eV: np.ndarray, temp: float, x_bulk: float) -> np.ndarray:
    pref = (1.0 - x_bulk) / x_bulk
    expo = np.clip(dE_eV / (KB_EV_PER_K * temp), -700.0, 700.0)
    return 1.0 / (1.0 + pref * np.exp(expo))


def neighbor_lists_for_sites(
    positions: np.ndarray,
    box: np.ndarray,
    site_indices: np.ndarray,
    cutoff: float,
) -> dict[int, np.ndarray]:
    neigh: dict[int, np.ndarray] = {}
    for i in site_indices:
        p = positions[i]
        dr = positions - p
        dr -= box * np.round(dr / box)
        dist2 = np.sum(dr * dr, axis=1)
        mask = (dist2 <= cutoff * cutoff) & (dist2 > 1.0e-12)
        neigh[int(i)] = np.where(mask)[0]
    return neigh


def binned_means(x: np.ndarray, y: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centers: list[float] = []
    means: list[float] = []
    errs: list[float] = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (x >= lo) & (x < hi)
        if not np.any(mask):
            continue
        vals = y[mask]
        centers.append(0.5 * (lo + hi))
        means.append(float(np.mean(vals)))
        errs.append(float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0)
    return np.array(centers), np.array(means), np.array(errs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", required=True, help="LAMMPS HMC dump")
    parser.add_argument("--data", required=True, help="LAMMPS data file with coordinates")
    parser.add_argument("--deltae-npz", required=True)
    parser.add_argument("--gb-mask", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--temp", type=float, default=700.0)
    parser.add_argument("--xc", type=float, required=True)
    parser.add_argument("--energy-bin", nargs=2, type=float, default=[-30.0, -15.0], help="kJ/mol")
    parser.add_argument(
        "--energy-xlim",
        nargs=2,
        type=float,
        default=None,
        help="Optional x-axis limits for the P_i vs Delta E_i plot, in kJ/mol.",
    )
    parser.add_argument("--neighbor-cutoff", type=float, default=5.0)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    site_ids, dE = load_deltae(Path(args.deltae_npz))
    dE_kj = dE * EV_TO_KJMOL
    positions, _, box = read_data_positions(Path(args.data))
    frames = dump_type_frames(Path(args.dump), positions.shape[0])
    is_solute = frames == 2
    p_occ = is_solute[:, site_ids - 1].mean(axis=0)
    p_err = np.sqrt(np.maximum(p_occ * (1.0 - p_occ), 0.0) / frames.shape[0])

    summary = json.loads(Path(args.summary_json).read_text())
    x_bulk = float(summary["tail_mean_x_bulk"])
    fd_site = fd_prob(dE, args.temp, x_bulk)
    gb_mask = np.load(args.gb_mask)
    fd_xgb_closed, fd_xbulk_closed = x_gb_canonical(
        dE,
        args.temp,
        args.xc,
        n_total=positions.shape[0],
        n_gb=int(gb_mask.sum()),
    )

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    with (out_prefix.with_name(out_prefix.name + "_site_occupancy.csv")).open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["site_id", "deltaE_kJ_mol", "P_HMC", "P_HMC_err", "P_FD_at_HMC_bulk"])
        writer.writerows(zip(site_ids, dE_kj, p_occ, p_err, fd_site))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    bins = np.linspace(math.floor(dE_kj.min() / 5) * 5, math.ceil(dE_kj.max() / 5) * 5, 14)
    xc_bin, p_bin, e_bin = binned_means(dE_kj, p_occ, bins)
    dE_grid = np.linspace(dE_kj.min(), dE_kj.max(), 400)
    fd_grid = fd_prob(dE_grid / EV_TO_KJMOL, args.temp, x_bulk)

    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    ax.plot(dE_grid, fd_grid, color="black", lw=2.0, label="closed-box FD at HMC bulk")
    ax.errorbar(xc_bin, p_bin, yerr=e_bin, fmt="o", color="#e31a1c", capsize=3, label="HMC measurement")
    ax.axvspan(min(args.energy_bin), max(args.energy_bin), color="#d8ead8", alpha=0.45)
    ax.axvline(0.0, color="#76a476", lw=1.0)
    if args.energy_xlim is not None:
        ax.set_xlim(min(args.energy_xlim), max(args.energy_xlim))
    ax.set_xlabel(r"$\Delta E_i$  [kJ/mol]")
    ax.set_ylabel(r"Au occupation probability $P_i$")
    ax.set_title(
        rf"Per-site Au occupation vs $\Delta E_i$ at $X={args.xc:g}$ "
        rf"($X_{{GB}}={summary['tail_mean_x_gb']:.3f}$)"
    )
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_prefix.with_name(out_prefix.name + "_Pi_vs_deltaE.png"))
    plt.close(fig)

    neigh = neighbor_lists_for_sites(positions, box, site_ids - 1, args.neighbor_cutoff)
    selected = (dE_kj >= args.energy_bin[0]) & (dE_kj <= args.energy_bin[1])
    local_counts: list[float] = []
    selected_p: list[float] = []
    for idx, atom_id in enumerate(site_ids):
        if not selected[idx]:
            continue
        neighbors = neigh[int(atom_id - 1)]
        local_counts.append(float(is_solute[:, neighbors].sum(axis=1).mean()))
        selected_p.append(float(p_occ[idx]))
    local_counts_arr = np.array(local_counts)
    selected_p_arr = np.array(selected_p)
    if local_counts_arr.size == 0:
        raise ValueError(
            "No sampled GB sites fall inside the requested energy bin "
            f"[{args.energy_bin[0]}, {args.energy_bin[1]}] kJ/mol"
        )
    count_bins = np.arange(-0.5, max(1.5, math.ceil(local_counts_arr.max()) + 1.5), 1.0)
    cx, cy, ce = binned_means(np.rint(local_counts_arr), selected_p_arr, count_bins)

    fd_lo = float(fd_prob(np.array([args.energy_bin[0] / EV_TO_KJMOL]), args.temp, x_bulk)[0])
    fd_hi = float(fd_prob(np.array([args.energy_bin[1] / EV_TO_KJMOL]), args.temp, x_bulk)[0])
    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=160)
    ax.axhspan(min(fd_lo, fd_hi), max(fd_lo, fd_hi), color="#cccccc", alpha=0.55, label="FD prediction band")
    ax.errorbar(cx, cy, yerr=ce, fmt="o", color="#e31a1c", capsize=3, label="HMC empirical")
    ax.set_xlabel(rf"local Au neighbours within $r \leq {args.neighbor_cutoff:g}$ Å")
    ax.set_ylabel(r"Au occupation probability $P_i$")
    ax.set_title(
        rf"Au occupation at $\Delta E_i \in [{args.energy_bin[0]:g},{args.energy_bin[1]:g}]$ kJ/mol "
        rf"($T={args.temp:g}$ K, $X={args.xc:g}$)"
    )
    ax.set_ylim(-0.05, 1.2)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_prefix.with_name(out_prefix.name + "_Pi_vs_local_Au.png"))
    plt.close(fig)

    print(f"Wrote {out_prefix.with_name(out_prefix.name + '_Pi_vs_deltaE.png')}")
    print(f"Wrote {out_prefix.with_name(out_prefix.name + '_Pi_vs_local_Au.png')}")
    print(f"Closed-box FD from sampled sites: X_GB={fd_xgb_closed:.4f}, X_bulk={fd_xbulk_closed:.4f}")


if __name__ == "__main__":
    main()
