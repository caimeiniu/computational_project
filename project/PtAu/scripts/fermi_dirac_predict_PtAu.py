"""Fermi-Dirac (FD) GB-segregation predictor (Wagih 2020 eq. 2) — Pt(Au) sibling.

Pt(Au) copy of project/scripts/fermi_dirac_predict.py — kept as a separate
file so the canonical Al(Mg) version is untouched. Only three substantive
differences: title string, "100 Å" instead of "200 Å" in the legend, and
docstring framed for Pt(Au).

Given a set of per-site segregation energies {ΔE_i}_{i=1..N_GB} from
sample_delta_e_PtAu.py, evaluate the dilute-limit prediction

    P_i(T, X_c) = 1 / (1 + ((1 - X_c)/X_c) · exp(ΔE_i / kT))
    X_GB^FD(T, X_c) = (1/N_GB) · Σ_i P_i(T, X_c)

on a (T, X_c) grid. Optionally overlay the same calculation done on
Wagih's Pt(Au) ΔE pool (97,440 sites from Pt_Au_20nm_GB_segregation.dump)
as a reference curve.

The standard form above is *grand-canonical*: X_c is the bulk solute
mole fraction, held fixed by an implicit infinite reservoir. For a
closed simulation box (canonical / total-solute-conserved) this can
predict X_GB unreachable in finite-size: X_GB^FD × N_GB > X_c × N_total.
The canonical version (`x_gb_canonical`) self-consistently solves for
the post-segregation bulk fraction X_bulk under mass conservation, so
the comparison to a closed-box HMC is apples-to-apples.

Sign convention: ΔE_i = E_GB^solute - E_bulk^solute (eV). Negative
means the GB site is energetically preferred → P → 1 (segregates).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

EV_TO_KJMOL = 96.485
KJMOL_TO_EV = 1.0 / EV_TO_KJMOL
KB_EV_PER_K = 8.617333262e-5  # Boltzmann constant in eV/K


def fermi_dirac(dE_eV: np.ndarray, T: float, X_c: float) -> np.ndarray:
    if X_c <= 0.0:
        return np.zeros_like(dE_eV)
    if X_c >= 1.0:
        return np.ones_like(dE_eV)
    kT = KB_EV_PER_K * T
    expo = np.clip(dE_eV / kT, -700.0, 700.0)
    pref = (1.0 - X_c) / X_c
    return 1.0 / (1.0 + pref * np.exp(expo))


def x_gb(dE_eV: np.ndarray, T: float, X_c: float) -> float:
    return float(fermi_dirac(dE_eV, T, X_c).mean())


def x_gb_curve(dE_eV: np.ndarray, T: float, X_c_grid: np.ndarray) -> np.ndarray:
    return np.array([x_gb(dE_eV, T, x) for x in X_c_grid])


def load_ours(npz_path: Path) -> np.ndarray:
    return np.asarray(np.load(npz_path)["gb_delta_e"])  # eV


def load_wagih_dump(dump_path: Path) -> np.ndarray:
    """Load Wagih Pt(Au) per-site ΔE_seg from accelerated_model dump.

    Same adapter as compare_vs_wagih_PtAu.py / bootstrap_vs_wagih_PtAu.py:
    reads the `seg_kJ_per_mol` column, drops bulk atoms (seg == 0),
    returns ΔE in eV (FD predictor's native unit).
    """
    with open(dump_path) as f:
        lines = f.readlines()
    n_atoms, cols, data_start = 0, None, 0
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            n_atoms = int(lines[i + 1])
            i += 2
        elif line.startswith("ITEM: ATOMS"):
            cols = line.split()[2:]
            data_start = i + 1
            break
        else:
            i += 1
    if cols is None or "seg_kJ_per_mol" not in cols:
        raise ValueError(f"unexpected dump format; columns = {cols}")
    seg_idx = cols.index("seg_kJ_per_mol")
    data = np.loadtxt(lines[data_start:data_start + n_atoms])
    seg_kjmol = data[:, seg_idx]
    seg_kjmol = seg_kjmol[seg_kjmol != 0.0]
    return seg_kjmol * KJMOL_TO_EV


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ours-npz", required=True,
                   help="ΔE NPZ from sample_delta_e_PtAu.py")
    p.add_argument("--wagih-dump", default=None,
                   help="Wagih Pt_Au_20nm_GB_segregation.dump (optional reference)")
    p.add_argument("--T-list", default="500,700,900,1100",
                   help="comma-separated temperatures in K")
    p.add_argument("--Xc-min", type=float, default=1e-5)
    p.add_argument("--Xc-max", type=float, default=0.5)
    p.add_argument("--Xc-points", type=int, default=120)
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    dE_ours = load_ours(Path(args.ours_npz))
    print(f"Loaded ours: n={dE_ours.size}, mean={dE_ours.mean()*EV_TO_KJMOL:+.2f} kJ/mol, "
          f"frac(ΔE<0)={float((dE_ours<0).mean()):.3f}")
    have_wagih = args.wagih_dump is not None
    if have_wagih:
        dE_wagih = load_wagih_dump(Path(args.wagih_dump))
        print(f"Loaded Wagih: n={dE_wagih.size}, "
              f"mean={dE_wagih.mean()*EV_TO_KJMOL:+.2f} kJ/mol, "
              f"frac(ΔE<0)={float((dE_wagih<0).mean()):.3f}")

    Ts = [float(t) for t in args.T_list.split(",")]
    Xc_grid = np.logspace(np.log10(args.Xc_min), np.log10(args.Xc_max),
                          args.Xc_points)

    curves: dict[str, dict] = {"T": Ts, "X_c": Xc_grid.tolist(),
                                "ours": {}, "wagih": {}}
    for T in Ts:
        curves["ours"][f"{T:g}"] = x_gb_curve(dE_ours, T, Xc_grid).tolist()
        if have_wagih:
            curves["wagih"][f"{T:g}"] = x_gb_curve(dE_wagih, T, Xc_grid).tolist()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=140)
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd"]
    for T, c in zip(Ts, colors):
        ax.plot(Xc_grid, curves["ours"][f"{T:g}"], "-", color=c, lw=1.8,
                label=f"ours 100 Å (n={dE_ours.size})  T={T:g} K")
        if have_wagih:
            ax.plot(Xc_grid, curves["wagih"][f"{T:g}"], "--", color=c, lw=1.2,
                    alpha=0.7, label=f"Wagih (n={dE_wagih.size})  T={T:g} K")
    ax.plot([args.Xc_min, args.Xc_max], [args.Xc_min, args.Xc_max], ":",
            color="grey", lw=0.8, label=r"$X_\mathrm{GB} = X_c$ (no segregation)")
    ax.set_xscale("log")
    ax.set_xlabel(r"bulk solute fraction  $X_c$")
    ax.set_ylabel(r"GB occupancy  $X_\mathrm{GB}^\mathrm{FD}$")
    ax.set_title(r"Fermi–Dirac dilute-limit prediction — Pt(Au)")
    ax.set_ylim(0, 1.02)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7, loc="upper left", framealpha=0.95)
    fig.tight_layout()
    Path(args.out_png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_png)
    plt.close(fig)
    print(f"→ wrote {args.out_png}")

    Path(args.out_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_json).write_text(json.dumps(curves, indent=2))
    print(f"→ wrote {args.out_json}")


if __name__ == "__main__":
    main()
