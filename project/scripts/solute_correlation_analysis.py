"""Mechanism analysis from existing finite-X_c HMC snapshots — three orthogonal
lines of evidence for solute-solute interactions, all cheap (no LAMMPS).

Inputs:
    - List of `(xc_label, snapshot_lmp_path)` pairs (atom_style atomic, types {1=Al, 2=Mg})
    - GB mask `.npy` (per-atom bool)
    - Reference X_c=0 spectrum `.npz` from `sample_delta_e.py` (sites' ΔE_i and 1-based ids)

Outputs:
    1. `output/g_MgMg_pair_correlation.{png,json}` — radial pair correlation
       g_MgMg(r) for GB-Mg atoms, normalized by random reference (uniform
       random GB-site placement of same N_Mg^GB count). Deviation from 1 =
       Mg-Mg correlation; range over which it deviates = interaction length.

    2. `output/site_occupation_vs_energy.{png,json}` — empirical P_i (from
       reference sites' current type in the snapshot, binary 0/1 → averaged
       per ΔE bin) overlaid with Wagih's FD prediction
       P_i^Wagih(ΔE_i; T, X_c) = 1 / (1 + ((1−X_c)/X_c) exp(ΔE_i / kT)).
       Systematic gap = Wagih breakdown localised in the spectrum.

    3. `output/site_occupation_vs_density.{png,json}` — within a fixed ΔE
       window (favourable end), P_i vs per-site n_Mg^local(r=5Å). Slope < 0
       = direct site-level repulsion: more neighbour Mg → lower occupation.

Run:
    python scripts/solute_correlation_analysis.py
        (no CLI flags by default — paths hardcoded for this project; edit the
        SNAPSHOTS / GB_MASK / REFERENCE_NPZ constants if reused on a new run)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parent.parent
SNAPSHOTS = [
    ("0.075_preseg",    REPO / "data/snapshots/hmc_T500_Xc0.075_preseg_final.lmp",           0.075),
    ("0.10_multistart", REPO / "data/snapshots/hmc_T500_Xc0.10_multistart_xgb0.3_final.lmp", 0.10),
    ("0.10_preseg",     REPO / "data/snapshots/hmc_T500_Xc0.10_preseg_final.lmp",            0.10),
    ("0.15_preseg",     REPO / "data/snapshots/hmc_T500_Xc0.15_preseg_final.lmp",            0.15),
    ("0.20_preseg",     REPO / "data/snapshots/hmc_T500_Xc0.20_preseg_final.lmp",            0.20),
    ("0.30_preseg",     REPO / "data/snapshots/hmc_T500_Xc0.30_preseg_final.lmp",            0.30),
]
GB_MASK = REPO / "data/snapshots/gb_mask_200A.npy"
REFERENCE_NPZ = Path("/cluster/scratch/cainiu/production_AlMg_200A/delta_e_results_n500_200A_tight.npz")
OUT_DIR = REPO / "output"

T_K = 500.0
KT_KJMOL = 8.314e-3 * T_K   # 4.157 kJ/mol at T=500 K

R_MAX_GR = 25.0   # Å — well past expected interaction range (~15–20 Å for elastic)
DR_GR = 0.2       # Å — bin width
R_LOCAL = 5.0     # Å — cutoff for "local Mg density" (~1st NN shell in FCC Al)
RNG_SEED = 20260429


def read_lmp(data_file: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Parse atom_style atomic data file. Returns (positions, types, box_lengths).
    Atoms section assumed `id type x y z [ix iy iz]`.
    """
    text = data_file.read_text()
    lines = text.splitlines()
    n_atoms = None
    lo = {"x": None, "y": None, "z": None}
    hi = {"x": None, "y": None, "z": None}
    atoms_start = None
    for i, line in enumerate(lines):
        s = line.split()
        if len(s) >= 2 and s[1] == "atoms":
            n_atoms = int(s[0])
        elif len(s) >= 4 and s[2] in ("xlo", "ylo", "zlo") and s[3] in ("xhi", "yhi", "zhi"):
            ax = s[2][0]
            lo[ax], hi[ax] = float(s[0]), float(s[1])
        elif line.strip() == "Atoms" or line.strip().startswith("Atoms "):
            atoms_start = i + 2
            break
    if n_atoms is None or atoms_start is None:
        raise ValueError(f"couldn't parse header in {data_file}")
    box_lengths = np.array([hi[a] - lo[a] for a in "xyz"])
    rows = []
    for raw in lines[atoms_start : atoms_start + n_atoms]:
        parts = raw.split()
        rows.append((int(parts[0]), int(parts[1]),
                     float(parts[2]), float(parts[3]), float(parts[4])))
    rows.sort(key=lambda r: r[0])
    positions = np.array([[x, y, z] for _, _, x, y, z in rows])
    types = np.array([t for _, t, _, _, _ in rows], dtype=np.int64)
    positions = positions - np.array([lo[a] for a in "xyz"])
    positions = np.mod(positions, box_lengths)
    return positions, types, box_lengths


def compute_g_r(
    positions: np.ndarray,
    box: np.ndarray,
    r_max: float,
    dr: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Radial pair correlation g(r) for the given positions, with PBC.

    g(r) = h(r) / (4π r² dr ρ N), where h(r) is the histogram of pair
    distances and ρ N gives the expected uniform-density count.

    Uses cKDTree.query_pairs for efficient PBC pair finding.
    """
    n = len(positions)
    if n < 2:
        return np.array([]), np.array([])
    tree = cKDTree(positions, boxsize=box)
    pairs = tree.query_pairs(r_max, output_type="ndarray")
    # Compute distances (with PBC) for each pair
    delta = positions[pairs[:, 0]] - positions[pairs[:, 1]]
    delta -= box * np.round(delta / box)
    dists = np.linalg.norm(delta, axis=1)
    bin_edges = np.arange(0.0, r_max + dr, dr)
    r_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    counts, _ = np.histogram(dists, bins=bin_edges)
    # Normalise: 2× counts because query_pairs gives unique pairs (i<j);
    # we want both i→j and j→i in g(r) per convention.
    counts = 2 * counts
    volume = float(np.prod(box))
    rho = n / volume
    shell_volumes = 4.0 * np.pi * r_centers**2 * dr
    g = counts / (rho * n * shell_volumes)
    return r_centers, g


def random_gb_mg_positions(
    positions: np.ndarray, gb_mask: np.ndarray, n_gb_mg: int, rng: np.random.Generator
) -> np.ndarray:
    """Pick n_gb_mg random GB sites and return their positions — i.e., a
    uniform-random Mg placement on the GB lattice, preserving GB geometry
    and X_GB count, but removing any Mg-Mg correlation."""
    gb_atom_ids = np.where(gb_mask)[0]
    chosen = rng.choice(gb_atom_ids, size=n_gb_mg, replace=False)
    return positions[chosen]


def panel_g_r(snapshots_data: list[dict]) -> dict:
    """Plot g_MgMg(r) / g_random(r) for each snapshot. Returns json-able payload."""
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    payload = {"r_axis": None, "curves": {}, "comment": (
        "g_MgMg(r) for GB-Mg atoms, divided by g_MgMg(r) of a random uniform "
        "placement on GB sites at the same N_Mg^GB count. Ratio < 1 at small r "
        "= Mg avoid each other (repulsive); > 1 = clustering. Range over which "
        "ratio deviates from 1 = interaction length scale."
    )}
    cmap = plt.cm.viridis(np.linspace(0.15, 0.85, len(snapshots_data)))
    for k, snap in enumerate(snapshots_data):
        r = snap["g_r"]["r"]
        g = snap["g_r"]["g_real"]
        g_rand = snap["g_r"]["g_random"]
        ratio = np.where(g_rand > 0, g / g_rand, np.nan)
        if payload["r_axis"] is None:
            payload["r_axis"] = r.tolist()
        payload["curves"][snap["label"]] = {
            "ratio": ratio.tolist(),
            "g_real": g.tolist(),
            "g_random": g_rand.tolist(),
            "n_gb_mg": int(snap["n_gb_mg"]),
            "x_gb": float(snap["x_gb"]),
        }
        ax.plot(r, ratio, color=cmap[k], lw=1.6,
                label=f"{snap['label']} (X_GB={snap['x_gb']:.3f})")
    ax.axhline(1.0, color="k", lw=0.7, ls="--", alpha=0.6)
    ax.set_xlim(0, R_MAX_GR)
    ax.set_xlabel("r [Å]")
    ax.set_ylabel(r"$g_\mathrm{MgMg}^{\rm HMC}(r) / g_\mathrm{MgMg}^{\rm random}(r)$")
    ax.set_title("Mg-Mg pair correlation — GB-only Mg, normalised by uniform-random GB reference")
    ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "g_MgMg_pair_correlation.png", dpi=160)
    plt.close(fig)
    return payload


def wagih_p_fd(delta_e_kjmol: np.ndarray, x_c: float) -> np.ndarray:
    """P_i^Wagih = 1 / (1 + ((1-X_c)/X_c) exp(ΔE/kT))."""
    return 1.0 / (1.0 + ((1.0 - x_c) / x_c) * np.exp(delta_e_kjmol / KT_KJMOL))


def panel_occupation_vs_energy(snapshots_data: list[dict], ref_de_kjmol: np.ndarray) -> dict:
    """For each snapshot, bin reference sites by ΔE_i, plot fraction-occupied
    vs Wagih FD curve at the snapshot's X_c."""
    fig, axes = plt.subplots(1, len(snapshots_data),
                             figsize=(3.2 * len(snapshots_data), 3.6),
                             sharey=True)
    if len(snapshots_data) == 1:
        axes = [axes]
    payload = {"comment": (
        "Empirical P_i (binary, type==2 in snapshot) averaged within ΔE bins, "
        "vs Wagih FD curve P_i^Wagih(ΔE_i; T, X_c). Systematic gap below the "
        "Wagih curve at favourable (negative) ΔE_i = strong-binding sites are "
        "occupied less than independent-site theory predicts."
    ), "snapshots": {}}
    bin_edges = np.linspace(ref_de_kjmol.min(), ref_de_kjmol.max(), 11)  # 10 bins
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    for ax, snap in zip(axes, snapshots_data):
        p_emp = snap["p_emp"]                           # (500,) binary 0/1
        x_c = snap["x_c"]
        # Per-bin mean P_i + 95% CI from binomial proportion (Wilson approx).
        bin_idx = np.digitize(ref_de_kjmol, bin_edges) - 1
        bin_idx = np.clip(bin_idx, 0, len(bin_centers) - 1)
        means, lows, highs, counts = [], [], [], []
        for b in range(len(bin_centers)):
            sel = bin_idx == b
            n = int(sel.sum())
            if n == 0:
                means.append(np.nan); lows.append(np.nan); highs.append(np.nan); counts.append(0)
                continue
            k = int(p_emp[sel].sum())
            phat = k / n
            # Normal approximation 95% CI (n=50 per bin avg, OK)
            se = np.sqrt(phat * (1 - phat) / n) if 0 < phat < 1 else 1.0 / n
            means.append(phat)
            lows.append(max(0.0, phat - 1.96 * se))
            highs.append(min(1.0, phat + 1.96 * se))
            counts.append(n)
        means = np.array(means); lows = np.array(lows); highs = np.array(highs)
        # Wagih theory line: dense ΔE grid
        de_grid = np.linspace(ref_de_kjmol.min(), ref_de_kjmol.max(), 200)
        p_theory = wagih_p_fd(de_grid, x_c)
        ax.plot(de_grid, p_theory, "k-", lw=1.5, label=f"Wagih FD (X_c={x_c:.3f})")
        ax.errorbar(bin_centers, means, yerr=[means - lows, highs - means],
                    fmt="o", color="C3", capsize=3, label="empirical (n=500 bins)")
        ax.axvline(0, color="gray", ls=":", lw=0.6, alpha=0.5)
        ax.axhline(snap["x_gb"], color="C0", ls="--", lw=0.6, alpha=0.6,
                   label=f"X_GB={snap['x_gb']:.3f}")
        ax.set_xlabel(r"$\Delta E_i$ [kJ/mol] (X_c=0 reference)")
        ax.set_title(snap["label"])
        ax.set_ylim(-0.05, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="best")
        payload["snapshots"][snap["label"]] = {
            "bin_centers_kjmol": bin_centers.tolist(),
            "p_empirical": means.tolist(),
            "p_low95": lows.tolist(),
            "p_high95": highs.tolist(),
            "n_per_bin": counts,
            "x_c": float(x_c),
            "x_gb": float(snap["x_gb"]),
        }
    axes[0].set_ylabel(r"$P_i$ (occupation probability)")
    fig.suptitle("Occupation vs ΔE_i — empirical (HMC) vs Wagih FD prediction", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "site_occupation_vs_energy.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return payload


def panel_occupation_vs_density(
    snapshots_data: list[dict], ref_de_kjmol: np.ndarray
) -> dict:
    """Within a favourable-ΔE window, plot P_i vs per-site n_Mg^local."""
    de_window = (-30.0, -5.0)   # favourable but not the deepest tail (better statistics)
    in_window = (ref_de_kjmol >= de_window[0]) & (ref_de_kjmol <= de_window[1])
    n_in = int(in_window.sum())
    fig, axes = plt.subplots(1, len(snapshots_data),
                             figsize=(3.2 * len(snapshots_data), 3.6),
                             sharey=True)
    if len(snapshots_data) == 1:
        axes = [axes]
    payload = {"comment": (
            f"P_i restricted to favourable-ΔE window {de_window} kJ/mol "
            f"({n_in}/500 sites), binned by per-site n_Mg^local(r=5Å). Slope < 0 = "
            "more neighbour Mg → lower occupation = direct site-level repulsion."
        ),
        "delta_e_window_kjmol": list(de_window),
        "n_sites_in_window": n_in,
        "r_local_angstrom": R_LOCAL,
        "snapshots": {},
    }
    for ax, snap in zip(axes, snapshots_data):
        p_emp = snap["p_emp"][in_window]
        n_mg = snap["n_mg_local"][in_window]
        if len(n_mg) == 0:
            ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center")
            continue
        # Bin by n_mg integer values
        n_mg_unique = np.arange(int(n_mg.min()), int(n_mg.max()) + 1)
        means, errs, counts = [], [], []
        for v in n_mg_unique:
            sel = n_mg == v
            n_sel = int(sel.sum())
            if n_sel == 0:
                means.append(np.nan); errs.append(np.nan); counts.append(0); continue
            ph = float(p_emp[sel].mean())
            se = np.sqrt(ph * (1 - ph) / n_sel) if 0 < ph < 1 else 1.0 / max(n_sel, 1)
            means.append(ph); errs.append(1.96 * se); counts.append(n_sel)
        means = np.array(means); errs = np.array(errs)
        ax.errorbar(n_mg_unique, means, yerr=errs, fmt="o", color="C2", capsize=3,
                    label=f"empirical ({n_in} sites)")
        # Linear fit (only points with valid means and n_sel >= 3)
        valid = ~np.isnan(means) & (np.array(counts) >= 3)
        if valid.sum() >= 2:
            coeff = np.polyfit(n_mg_unique[valid], means[valid], 1, w=1.0 / np.maximum(errs[valid], 1e-3))
            xfit = np.linspace(n_mg_unique.min(), n_mg_unique.max(), 50)
            ax.plot(xfit, np.polyval(coeff, xfit), "k--", lw=1.0,
                    label=f"slope={coeff[0]:+.3f}/Mg-neighbour")
        ax.set_xlabel(rf"$n_\mathrm{{Mg}}^\mathrm{{local}}(r\!\leq\!{R_LOCAL:.0f}\,\mathrm{{Å}})$")
        ax.set_title(snap["label"])
        ax.set_ylim(-0.05, 1.05)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, loc="best")
        payload["snapshots"][snap["label"]] = {
            "n_mg_unique": n_mg_unique.tolist(),
            "p_empirical": means.tolist(),
            "p_err95": errs.tolist(),
            "n_per_bin": counts,
            "x_c": float(snap["x_c"]),
            "x_gb": float(snap["x_gb"]),
            "slope_per_neighbour": float(coeff[0]) if valid.sum() >= 2 else None,
        }
    axes[0].set_ylabel(r"$P_i$ (within favourable ΔE window)")
    fig.suptitle(rf"Occupation vs local Mg coordination — ΔE window {de_window} kJ/mol", y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "site_occupation_vs_density.png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    return payload


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    gb_mask = np.load(GB_MASK).astype(bool)
    ref = np.load(REFERENCE_NPZ)
    ref_site_ids_1based = np.asarray(ref["gb_site_ids"], dtype=np.int64)
    ref_de_eV = np.asarray(ref["gb_delta_e"], dtype=float)
    ref_de_kjmol = ref_de_eV * 96.485

    print(f"Loaded reference: n={len(ref_site_ids_1based)} sites, "
          f"ΔE_i range [{ref_de_kjmol.min():.1f}, {ref_de_kjmol.max():.1f}] kJ/mol")

    snapshots_data = []
    for label, path, x_c in SNAPSHOTS:
        if not path.exists():
            print(f"  SKIP {label}: {path} not found", file=sys.stderr)
            continue
        print(f"--- {label} ---")
        positions, types, box = read_lmp(path)
        n_atoms = len(positions)
        if gb_mask.shape != (n_atoms,):
            raise ValueError(f"mask len {len(gb_mask)} != n_atoms={n_atoms} for {label}")
        is_mg = (types == 2)
        is_gb_mg = is_mg & gb_mask
        n_gb_mg = int(is_gb_mg.sum())
        n_gb = int(gb_mask.sum())
        x_gb = n_gb_mg / n_gb
        x_c_actual = is_mg.sum() / n_atoms
        print(f"  n_atoms={n_atoms:,}  X_c_actual={x_c_actual:.4f}  X_GB={x_gb:.4f}  N_GB_Mg={n_gb_mg:,}")

        # 1. g_MgMg(r) on GB-Mg atoms
        r, g_real = compute_g_r(positions[is_gb_mg], box, R_MAX_GR, DR_GR)
        rand_pos = random_gb_mg_positions(positions, gb_mask, n_gb_mg, rng)
        _, g_random = compute_g_r(rand_pos, box, R_MAX_GR, DR_GR)

        # 2. Empirical P_i for the 500 reference sites
        p_emp = is_mg[ref_site_ids_1based - 1].astype(int)

        # 3. Per-site local Mg density at reference sites
        mg_positions = positions[is_mg]
        mg_tree = cKDTree(mg_positions, boxsize=box)
        ref_positions = positions[ref_site_ids_1based - 1]
        # Count NEIGHBOURS within R_LOCAL — careful: reference sites that are
        # currently Mg will have 0 distance to themselves. Subtract self-count.
        neighbor_lists = mg_tree.query_ball_point(ref_positions, R_LOCAL)
        n_mg_local = np.array([len(nbrs) for nbrs in neighbor_lists]) - p_emp
        # (- p_emp removes self if site is Mg; if site is Al, no self present)

        snapshots_data.append({
            "label": label,
            "x_c": x_c,
            "x_c_actual": float(x_c_actual),
            "x_gb": x_gb,
            "n_gb_mg": n_gb_mg,
            "g_r": {"r": r, "g_real": g_real, "g_random": g_random},
            "p_emp": p_emp,
            "n_mg_local": n_mg_local,
        })

    print()
    print("=== Building 3 panels ===")
    payload_g = panel_g_r(snapshots_data)
    payload_occ_e = panel_occupation_vs_energy(snapshots_data, ref_de_kjmol)
    payload_occ_n = panel_occupation_vs_density(snapshots_data, ref_de_kjmol)

    # Persist json metadata
    out_json = {
        "g_MgMg_pair_correlation": payload_g,
        "site_occupation_vs_energy": payload_occ_e,
        "site_occupation_vs_density": payload_occ_n,
        "T_K": T_K, "kT_kjmol": KT_KJMOL,
        "R_MAX_GR": R_MAX_GR, "DR_GR": DR_GR, "R_LOCAL": R_LOCAL,
        "n_reference_sites": int(len(ref_site_ids_1based)),
    }
    (OUT_DIR / "solute_correlation_analysis.json").write_text(
        json.dumps(out_json, indent=2, default=float))
    print(f"\nWrote 3 PNGs + 1 JSON to {OUT_DIR}/")
    print("  g_MgMg_pair_correlation.png    (panel 1)")
    print("  site_occupation_vs_energy.png  (panel 2)")
    print("  site_occupation_vs_density.png (panel 3)")
    print("  solute_correlation_analysis.json (metadata for all three)")


if __name__ == "__main__":
    main()
