#!/usr/bin/env python3
"""Report figure 00 (headline) — 9-pt all-fdseed update for 2026-05-11 advisor meeting.

History anchor: the previous 00_headline_hmc_vs_wagih_T500.png (2026-05-02 evening)
showed 6 X_c points across mixed IC classes (1 fdseed-verify plateau at X_c=0.05,
5 preseg UB ▽, 1 multistart UB □ at X_c=0.10). Methodology drifted to fdseed-only
between 2026-05-07 and 2026-05-11; the headline curve is now 9 fdseed measurements
at T=500 K covering X_c ∈ [0.04, 0.40].

Plateau / UB classification (2026-05-11 evening):
  X_c=0.04  fdseed_resume    plateau  (CI95 half-width 0.00009; post-burnin imb -0.302)
  X_c=0.05  fdseed           UB       (imb -0.557, Q1->Q5 drift -0.047)
  X_c=0.06  fdseed           UB       (imb -0.622, drift -0.054)
  X_c=0.075 fdseed           UB       (imb -0.627, drift -0.065)
  X_c=0.10  fdseed           UB       (imb -0.691, drift -0.073)
  X_c=0.15  fdseed           UB       (imb -0.648, drift -0.069)
  X_c=0.20  fdseed           UB       (imb -0.678, drift -0.073)
  X_c=0.30  fdseed (timeout) UB       (imb -0.596, drift -0.055)
  X_c=0.40  fdseed           UB       (imb -0.497, drift -0.049)

Visual conventions match the previous 00_headline (per 2026-05-02 user feedback
"10个 legend, 非常的杂乱, 这个图只是要证明 break 了就行了"):
  - single reference curve: canon-FD (ours, n=500). KS p=0.89 vs Wagih's n=82,646
    Zenodo spectrum (Fig 4) → one curve = both frameworks.
  - ●  filled red circle = HMC plateau (post-burnin block-bootstrap mean,
       CI95 errorbar drawn — invisible at this scale for X_c=0.04 but present)
  - ▽  open red down-triangle = HMC upper bound (one-sided; drawn at the
       post-burnin mean per the 9-pt panel d convention, NOT at the trajectory
       end value — the mean is the conservative report value already published
       in panel (d) 9-pt)
  - legend: 3 rows only.

This is a copy + rename of canonical_fd_compare_5pt_with_multistart.py per the
``no-in-place-script-edits`` rule (feedback memory); the previous report-headline
script stays intact.

Output:
  report/figures/00_headline_hmc_vs_wagih_T500_9pt_2026-05-11.png
  output/00_headline_hmc_vs_wagih_T500_9pt_2026-05-11.json
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fermi_dirac_predict import load_ours, x_gb_canonical_curve

REPO = Path("/cluster/home/cainiu/Computational_modeling/project")
OUT = REPO / "output"
FIG_DIR = REPO / "report" / "figures"
SPECTRUM_NPZ = Path("/cluster/scratch/cainiu/production_AlMg_200A/"
                    "delta_e_results_n500_200A_tight.npz")
GB_MASK_NPY = Path("/cluster/scratch/cainiu/production_AlMg_200A/"
                   "gb_mask_200A.npy")


def load_hmc_pt(json_path: Path) -> dict:
    d = json.load(open(json_path))
    return {
        "Xc": d["X_c"],
        "X_HMC": d["x_gb"]["mean"],
        "lo": d["x_gb"]["ci95_lo"],
        "hi": d["x_gb"]["ci95_hi"],
        "n_frames": d["x_gb"]["n_frames_total"],
        "imb": (d.get("swap_decomposition", {})
                  .get("post_burnin", {}).get("imbalance_signed")),
        "source": Path(json_path).name,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plateau", nargs="*", default=[
        "hmc_T500_Xc0.04_fdseed_resume.json"],
        help="fdseed HMC JSONs deemed plateau-reads (drawn ●).")
    ap.add_argument("--upper-bound", nargs="*", default=[
        "hmc_T500_Xc0.05_fdseed.json",
        "hmc_T500_Xc0.06_fdseed.json",
        "hmc_T500_Xc0.075_fdseed.json",
        "hmc_T500_Xc0.10_fdseed.json",
        "hmc_T500_Xc0.15_fdseed.json",
        "hmc_T500_Xc0.20_fdseed.json",
        "hmc_T500_Xc0.30_fdseed.json",
        "hmc_T500_Xc0.40_fdseed.json"],
        help="fdseed HMC JSONs still in the descending phase (drawn ▽).")
    ap.add_argument("--out-prefix",
                    default="00_headline_hmc_vs_wagih_T500_9pt_2026-05-11")
    ap.add_argument("--temp", type=float, default=500.0)
    ap.add_argument("--xc-max", type=float, default=0.45,
                    help="x-axis upper bound (linear scale).")
    args = ap.parse_args()

    pl_paths = [OUT / f for f in args.plateau]
    ub_paths = [OUT / f for f in args.upper_bound]

    # ---- box geometry + spectrum ----
    gb_mask = np.load(GB_MASK_NPY).astype(bool)
    N_TOTAL = int(gb_mask.size)
    N_GB = int(gb_mask.sum())
    GB_FRAC = N_GB / N_TOTAL
    dE = load_ours(SPECTRUM_NPZ)
    print(f"ours box: N_total={N_TOTAL}, N_GB={N_GB}, GB_frac={GB_FRAC:.4f}")
    print(f"ours spectrum: n={dE.size}, mean={dE.mean()*96.485:+.2f} kJ/mol")

    # ---- single reference curve: canon-FD (ours) ----
    T = args.temp
    xc_grid = np.linspace(1e-4, args.xc_max, 400)
    x_gb_canon_ours, _ = x_gb_canonical_curve(dE, T, xc_grid, N_GB, N_TOTAL)

    # ---- HMC points ----
    pl_pts = [load_hmc_pt(p) for p in pl_paths]
    ub_pts = [load_hmc_pt(p) for p in ub_paths]

    # FD prediction at each point's X_c (for traceability + console table)
    for pt in pl_pts + ub_pts:
        pt["X_FD_canon_ours"] = float(x_gb_canonical_curve(
            dE, T, np.array([pt["Xc"]]), N_GB, N_TOTAL)[0][0])
        pt["gap_HMC_minus_canon_ours"] = pt["X_HMC"] - pt["X_FD_canon_ours"]

    # ---- JSON dump for traceability ----
    summary = {
        "T": T,
        "geometry_ours": {"N_total": N_TOTAL, "N_GB": N_GB,
                          "GB_frac": GB_FRAC,
                          "spectrum_n": int(dE.size),
                          "spectrum_mean_kjmol": float(dE.mean() * 96.485)},
        "reference_curve": {
            "X_c": xc_grid.tolist(),
            "X_GB_canon_ours": x_gb_canon_ours.tolist(),
        },
        "plateau": pl_pts,
        "upper_bound": ub_pts,
        "notes": (
            "9-pt all-fdseed headline for 2026-05-11 advisor meeting. "
            "X_c=0.04 uses the resume JSON (post-burnin block-bootstrap "
            "mean, CI95 half-width 0.00009); the 8 UB points use the "
            "post-burnin mean (not trajectory end value) so the marker "
            "co-locates with the 9-pt panel (d) reading. Wagih spectrum "
            "is intentionally omitted (KS p=0.89 ↔ ours, Fig 4) per the "
            "2026-05-02 simplification feedback."),
    }
    out_json = OUT / f"{args.out_prefix}.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"wrote {out_json}")

    # ---- plot (minimal legend per 2026-05-02 user feedback) ----
    fig, ax = plt.subplots(figsize=(7.0, 5.0))

    fd_line, = ax.plot(xc_grid, x_gb_canon_ours, "-", color="C2", lw=1.8,
                       label="Wagih FD prediction")

    # plateau: filled red ● + CI95 errorbar
    pl_container = None
    if pl_pts:
        x_pl = np.array([p["Xc"] for p in pl_pts])
        y_pl = np.array([p["X_HMC"] for p in pl_pts])
        lo_pl = np.array([p["lo"] for p in pl_pts])
        hi_pl = np.array([p["hi"] for p in pl_pts])
        pl_container = ax.errorbar(x_pl, y_pl,
                                    yerr=[y_pl - lo_pl, hi_pl - y_pl],
                                    fmt="o", color="C3", ms=8, capsize=3,
                                    elinewidth=1.0, zorder=6,
                                    label="HMC plateau")

    # upper bound: open red ▽ (down-triangle), no errorbar
    ub_handle = None
    if ub_pts:
        x_ub = np.array([p["Xc"] for p in ub_pts])
        y_ub = np.array([p["X_HMC"] for p in ub_pts])
        ub_handle, = ax.plot(x_ub, y_ub, "v", color="C3",
                              mfc="white", mew=1.4, ms=9, zorder=5,
                              label="HMC upper bound")

    ax.set_xlabel(r"total Mg fraction  $X_c$", fontsize=11)
    ax.set_ylabel(r"GB Mg fraction  $X_{\mathrm{GB}}$", fontsize=11)
    ax.set_title(
        f"Al(Mg) GB occupancy at T={T:g} K  "
        fr"($N_\mathrm{{GB}}/N_\mathrm{{tot}} = {GB_FRAC:.3f}$)",
        fontsize=12.5, pad=14)
    ax.set_xlim(0, args.xc_max)
    ax.set_ylim(0, 0.70)
    ax.grid(True, alpha=0.25, lw=0.5)
    legend_handles, legend_labels = [fd_line], ["Wagih FD prediction"]
    if pl_container is not None:
        legend_handles.append(pl_container)
        legend_labels.append("HMC plateau")
    if ub_handle is not None:
        legend_handles.append(ub_handle)
        legend_labels.append("HMC upper bound")
    ax.legend(legend_handles, legend_labels,
              loc="upper left", fontsize=10.5, framealpha=0.95)

    fig.tight_layout()
    out_png = FIG_DIR / f"{args.out_prefix}.png"
    fig.savefig(out_png, dpi=150)
    print(f"wrote {out_png}")

    # ---- console table ----
    hdr = (f"\n{'X_c':>7} {'X_HMC':>9} {'canon-O':>9} "
           f"{'gap-O':>9} {'imb':>9}  source")
    print(hdr)
    print("-" * len(hdr))
    def _row(pts, tag):
        if not pts:
            return
        print(f"({tag})")
        for pt in pts:
            print(f"{pt['Xc']:>7.3f} {pt['X_HMC']:>9.4f} "
                  f"{pt['X_FD_canon_ours']:>9.4f} "
                  f"{pt['gap_HMC_minus_canon_ours']:>+9.4f} "
                  f"{(pt['imb'] if pt['imb'] is not None else float('nan')):>+9.3f}  "
                  f"{pt['source']}")
    _row(pl_pts, "plateau ●")
    _row(ub_pts, "upper bound ▽")


if __name__ == "__main__":
    main()
