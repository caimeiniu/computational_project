"""Fermi-Dirac (FD) GB-segregation predictor (Wagih 2020 eq. 2).

Given a set of per-site segregation energies {ΔE_i}_{i=1..N_GB} from
sample_delta_e.py, evaluate the dilute-limit prediction

    P_i(T, X_c) = 1 / (1 + ((1 - X_c)/X_c) · exp(ΔE_i / kT))
    X_GB^FD(T, X_c) = (1/N_GB) · Σ_i P_i(T, X_c)

on a (T, X_c) grid. Optionally overlay the same calculation done on
Wagih's 82,646 ΔE pool as a reference curve.

Sign convention: ΔE_i = E_GB^solute - E_bulk^solute (eV). Negative
means the GB site is energetically preferred → P → 1 (segregates).

Unit tests (--self-test):
    X_c → 0:    X_GB → 0
    X_c → 1:    X_GB → 1
    T → ∞:      X_GB → X_c             (no thermodynamic preference)
    T → 0:      X_GB → fraction of ΔE_i < 0
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

EV_TO_KJMOL = 96.485
KB_EV_PER_K = 8.617333262e-5  # Boltzmann constant in eV/K


def fermi_dirac(dE_eV: np.ndarray, T: float, X_c: float) -> np.ndarray:
    """Per-site occupation probability (vectorised over sites)."""
    if X_c <= 0.0:
        return np.zeros_like(dE_eV)
    if X_c >= 1.0:
        return np.ones_like(dE_eV)
    kT = KB_EV_PER_K * T  # eV; for T>0
    # exp may overflow for very positive ΔE / very small T; clip exponent
    expo = np.clip(dE_eV / kT, -700.0, 700.0)
    pref = (1.0 - X_c) / X_c
    return 1.0 / (1.0 + pref * np.exp(expo))


def x_gb(dE_eV: np.ndarray, T: float, X_c: float) -> float:
    return float(fermi_dirac(dE_eV, T, X_c).mean())


def x_gb_curve(dE_eV: np.ndarray, T: float, X_c_grid: np.ndarray) -> np.ndarray:
    return np.array([x_gb(dE_eV, T, x) for x in X_c_grid])


def load_ours(npz_path: Path) -> np.ndarray:
    return np.asarray(np.load(npz_path)["gb_delta_e"])  # eV


def load_wagih(seg_txt: Path, bulk_dat: Path) -> np.ndarray:
    e_bulk = float(open(bulk_dat).read().strip().split()[0])
    out = []
    with open(seg_txt) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                out.append(float(parts[1]) - e_bulk)
    return np.array(out)  # eV


def self_test() -> None:
    rng = np.random.default_rng(0)
    dE = rng.normal(loc=0.05, scale=0.2, size=2000)  # eV
    # X_c → 0  (with negative-tail sites, X_GB → 0 polynomially in X_c;
    # check it scales linearly: halving X_c halves X_GB)
    v1, v2 = x_gb(dE, 500.0, 1e-10), x_gb(dE, 500.0, 5e-11)
    assert v1 < 1e-3 and abs(v2 / v1 - 0.5) < 0.05, f"X_c→0 scaling: {v1}, {v2}"
    # X_c → 1
    assert x_gb(dE, 500.0, 1.0 - 1e-9) > 1.0 - 1e-3
    # T → ∞ : every site sees exp(ΔE/kT) → 1, prefactor (1-X_c)/X_c →
    #         X_GB → X_c independent of {ΔE_i}
    for X_c in [0.01, 0.1, 0.5]:
        v = x_gb(dE, 1e10, X_c)
        assert abs(v - X_c) < 1e-6, f"T→∞ failed at X_c={X_c}: {v}"
    # T → 0 : sites with ΔE<0 fully occupied; ΔE>0 empty.
    expected = float((dE < 0).mean())
    v = x_gb(dE, 1.0, 0.5)  # 1 K, 50%: dominated by ΔE sign
    assert abs(v - expected) < 1e-3, f"T→0 failed: {v} vs {expected}"
    print(f"  self-test OK  (N={dE.size}; expected fraction(ΔE<0)={expected:.3f})")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ours-npz", help="ΔE NPZ from sample_delta_e.py")
    p.add_argument("--wagih-seg", default=None,
                   help="Wagih seg_energies_Al_Mg.txt (optional reference)")
    p.add_argument("--wagih-bulk", default=None,
                   help="Wagih bulk_solute_Al_Mg.dat (optional reference)")
    p.add_argument("--T-list", default="300,500,700,900",
                   help="comma-separated temperatures in K")
    p.add_argument("--Xc-min", type=float, default=1e-5)
    p.add_argument("--Xc-max", type=float, default=0.5)
    p.add_argument("--Xc-points", type=int, default=120)
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()

    if args.self_test:
        print("Running analytic-limit self-tests...")
        self_test()

    if args.ours_npz is None:
        print("--ours-npz not provided; self-test only.")
        return

    dE_ours = load_ours(Path(args.ours_npz))
    print(f"Loaded ours: n={dE_ours.size}, mean={dE_ours.mean()*EV_TO_KJMOL:+.2f} kJ/mol, "
          f"frac(ΔE<0)={float((dE_ours<0).mean()):.3f}")
    have_wagih = args.wagih_seg and args.wagih_bulk
    if have_wagih:
        dE_wagih = load_wagih(Path(args.wagih_seg), Path(args.wagih_bulk))
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

    # Plot — log-x, linear-y. One panel per T overlay.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=140)
    colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd"]
    for T, c in zip(Ts, colors):
        ax.plot(Xc_grid, curves["ours"][f"{T:g}"], "-", color=c, lw=1.8,
                label=f"ours 200 Å (n={dE_ours.size})  T={T:g} K")
        if have_wagih:
            ax.plot(Xc_grid, curves["wagih"][f"{T:g}"], "--", color=c, lw=1.2,
                    alpha=0.7, label=f"Wagih (n={dE_wagih.size})  T={T:g} K")
    ax.plot([args.Xc_min, args.Xc_max], [args.Xc_min, args.Xc_max], ":",
            color="grey", lw=0.8, label=r"$X_\mathrm{GB} = X_c$ (no segregation)")
    ax.set_xscale("log")
    ax.set_xlabel(r"bulk solute fraction  $X_c$")
    ax.set_ylabel(r"GB occupancy  $X_\mathrm{GB}^\mathrm{FD}$")
    ax.set_title(r"Fermi–Dirac dilute-limit prediction — Al(Mg)")
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
