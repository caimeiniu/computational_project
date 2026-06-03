#!/usr/bin/env python3
"""Generate current Al(Mg) report figures after the 2026-06-03 job check.

Outputs:
  output/almg_report_02_prediction_hmc_vs_fd_current.{json,png}
  output/almg_report_05_ergodicity_bracket_T700_T800_current.{json,png}

This script intentionally uses repo-side snapshots and output JSONs instead of
the older /cluster/scratch paths used by legacy plotting scripts.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from fermi_dirac_predict import load_ours, x_gb_canonical_curve


REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "output"
SNAP = REPO / "data" / "snapshots"

SPECTRUM_NPZ = SNAP / "delta_e_results_n500_200A_tight.npz"
GB_MASK_NPY = SNAP / "gb_mask_200A.npy"

FIG2_PREFIX = OUT / "almg_report_02_prediction_hmc_vs_fd_current"
FIG5_PREFIX = OUT / "almg_report_05_ergodicity_bracket_T700_T800_current"


def load_run(filename: str) -> dict:
    path = OUT / filename
    with path.open() as f:
        d = json.load(f)
    d["_source_file"] = filename
    return d


def ci_half_width(d: dict) -> float:
    return 0.5 * (d["x_gb"]["ci95_hi"] - d["x_gb"]["ci95_lo"])


def postburn_quintile_drift(d: dict) -> tuple[float, float]:
    """Return last-quintile minus first-quintile drift and drift/CIHW."""
    x = np.asarray(d["series"]["x_gb"], dtype=float)
    n_burn = int(d["x_gb"]["n_frames_burnin_dropped"])
    x = x[n_burn:]
    k = max(1, int(round(len(x) * 0.2)))
    drift = float(x[-k:].mean() - x[:k].mean())
    hw = ci_half_width(d)
    return drift, drift / hw if hw else float("nan")


def run_summary(d: dict) -> dict:
    drift, drift_over_ci = postburn_quintile_drift(d)
    pb = d.get("swap_decomposition", {}).get("post_burnin", {})
    return {
        "source": d["_source_file"],
        "T": d["T"],
        "X_c": d["X_c"],
        "X_GB_mean": d["x_gb"]["mean"],
        "CI95_half_width": ci_half_width(d),
        "first_frame": d["series"]["x_gb"][0],
        "last_frame": d["series"]["x_gb"][-1],
        "postburn_Q5_minus_Q1": drift,
        "postburn_Q5_minus_Q1_over_CIHW": drift_over_ci,
        "PE_drift_eV": d["thermo_prod"]["PE_drift_eV"],
        "fwd_over_rev_postburn": pb.get("fwd_over_rev"),
        "fd_predicted": d.get("fd_predicted"),
        "gap_hmc_minus_fd": d.get("gap_hmc_minus_fd"),
        "swap_accept_rate": d["swap"]["accept_rate_overall"],
        "n_frames_total": d["x_gb"]["n_frames_total"],
        "n_frames_postburn": (
            d["x_gb"]["n_frames_total"]
            - d["x_gb"]["n_frames_burnin_dropped"]
        ),
    }


def load_geometry_and_fd(T: float, xc_max: float) -> tuple[np.ndarray, np.ndarray, dict]:
    gb_mask = np.load(GB_MASK_NPY).astype(bool)
    n_total = int(gb_mask.size)
    n_gb = int(gb_mask.sum())
    dE = load_ours(SPECTRUM_NPZ)
    xc = np.linspace(1e-4, xc_max, 500)
    xgb, xbulk = x_gb_canonical_curve(dE, T, xc, n_gb, n_total)
    meta = {
        "N_total": n_total,
        "N_GB": n_gb,
        "GB_fraction": n_gb / n_total,
        "spectrum_n": int(dE.size),
        "spectrum_mean_kjmol": float(dE.mean() * 96.485),
        "X_c": xc.tolist(),
        "X_GB_canon_fd": xgb.tolist(),
        "X_bulk_canon_fd": xbulk.tolist(),
    }
    return xc, xgb, meta


def add_points(ax, runs: list[dict], marker: str, color: str, label: str,
               *, face: str | None = None, with_ci: bool = True,
               zorder: int = 5):
    x = np.array([d["X_c"] for d in runs])
    y = np.array([d["x_gb"]["mean"] for d in runs])
    if with_ci:
        yerr = np.array([ci_half_width(d) for d in runs])
        return ax.errorbar(
            x, y, yerr=yerr, fmt=marker, color=color,
            mfc=(face if face is not None else color), mec=color, mew=1.2,
            ms=7.5, capsize=3, lw=1.0, label=label, zorder=zorder,
        )
    (handle,) = ax.plot(
        x, y, marker, color=color, mfc=(face if face is not None else color),
        mec=color, mew=1.2, ms=8.0, lw=0, label=label, zorder=zorder,
    )
    return handle


def figure_2() -> dict:
    practical = [
        load_run("hmc_T500_Xc0.01_fdseed_resume3.json"),
        load_run("hmc_T500_Xc0.04_fdseed_resume.json"),
    ]
    upper = [
        load_run("hmc_T500_Xc0.075_eq_cont.json"),
        load_run("hmc_T500_Xc0.10_fdseed_resume11.json"),
    ]
    independent_ic = [
        load_run("hmc_T500_Xc0.10_multistart_xgb0.3.json"),
    ]

    xc, fd, fd_meta = load_geometry_and_fd(T=500.0, xc_max=0.12)

    summary = {
        "purpose": "Report Fig. 2 current HMC-vs-canonical-FD prediction panel.",
        "T": 500.0,
        "fd_reference": fd_meta,
        "practical_or_caveated_points": [run_summary(d) for d in practical],
        "descending_upper_bound_points": [run_summary(d) for d in upper],
        "independent_ic_upper_bound_points": [run_summary(d) for d in independent_ic],
        "notes": (
            "T=500 Xc=0.10 uses resume11. Upper-bound markers are still "
            "drifting and should not be described as equilibrium values."
        ),
    }
    FIG2_PREFIX.with_suffix(".json").write_text(json.dumps(summary, indent=2))

    fig, ax = plt.subplots(figsize=(6.9, 4.8))
    ax.plot(xc, fd, color="#2E7D32", lw=2.0,
            label="canonical FD prediction")
    add_points(ax, practical, "o", "#1F77B4", "HMC practical anchors")
    add_points(ax, upper, "v", "#C62828", "HMC descending upper bounds",
               face="white", with_ci=False)
    add_points(ax, independent_ic, "s", "0.35", "HMC independent-IC upper bound",
               face="white", with_ci=False)

    latest = upper[-1]
    ax.annotate(
        "latest Xc=0.10",
        xy=(latest["X_c"], latest["x_gb"]["mean"]),
        xytext=(0.072, 0.165),
        arrowprops={"arrowstyle": "->", "lw": 0.8, "color": "0.25"},
        fontsize=8.5,
    )

    ax.set_xlim(0.0, 0.115)
    ax.set_ylim(0.0, 0.38)
    ax.set_xlabel(r"total Mg fraction $X_c$")
    ax.set_ylabel(r"GB Mg fraction $X_{\mathrm{GB}}$")
    ax.set_title("Al(Mg) HMC remains below independent-site FD at 500 K")
    ax.grid(True, alpha=0.25, lw=0.5)
    ax.legend(loc="upper left", fontsize=8.8, framealpha=0.95)
    fig.tight_layout()
    fig.savefig(FIG2_PREFIX.with_suffix(".png"), dpi=180)
    plt.close(fig)
    return summary


def branch(files: list[str]) -> list[dict]:
    return [load_run(f) for f in files if (OUT / f).exists()]


def plot_branch(ax, runs: list[dict], color: str, marker: str, label: str):
    x = np.arange(len(runs))
    y = np.array([d["x_gb"]["mean"] for d in runs])
    yerr = np.array([ci_half_width(d) for d in runs])
    ax.errorbar(
        x, y, yerr=yerr, marker=marker, color=color, mfc="white",
        mec=color, mew=1.2, ms=6.8, lw=1.2, capsize=2.5, label=label,
    )


def figure_5() -> dict:
    panels = [
        {
            "T": 700,
            "fd": 0.2956,
            "lower": branch([
                "hmc_T700_Xc0.10_random_targeted.json",
                "hmc_T700_Xc0.10_random_targeted_resume.json",
                "hmc_T700_Xc0.10_random_targeted_resume2.json",
            ]),
            "upper": branch([
                "hmc_T700_Xc0.10_fdseed_targeted.json",
                "hmc_T700_Xc0.10_fdseed_targeted_resume.json",
                "hmc_T700_Xc0.10_fdseed_targeted_resume2.json",
                "hmc_T700_Xc0.10_fdseed_targeted_resume3.json",
                "hmc_T700_Xc0.10_fdseed_targeted_resume4.json",
                "hmc_T700_Xc0.10_fdseed_targeted_resume5.json",
            ]),
        },
        {
            "T": 800,
            "fd": 0.2734,
            "lower": branch([
                "hmc_T800_Xc0.10_random_targeted.json",
                "hmc_T800_Xc0.10_random_targeted_resume.json",
                "hmc_T800_Xc0.10_random_targeted_resume2.json",
                "hmc_T800_Xc0.10_random_targeted_resume3.json",
            ]),
            "upper": branch([
                "hmc_T800_Xc0.10_fdseed_targeted.json",
                "hmc_T800_Xc0.10_fdseed_targeted_resume.json",
                "hmc_T800_Xc0.10_fdseed_targeted_resume2.json",
                "hmc_T800_Xc0.10_fdseed_targeted_resume3.json",
                "hmc_T800_Xc0.10_fdseed_targeted_resume4.json",
                "hmc_T800_Xc0.10_fdseed_targeted_resume5.json",
            ]),
        },
    ]

    summary = {
        "purpose": "Report Fig. 5 current two-branch ergodicity/bracket check.",
        "panels": [],
        "notes": (
            "Lower and upper branches remain separated. Latest upper branches "
            "are still descending, so these are bracket diagnostics rather "
            "than equilibrium estimates."
        ),
    }

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.8), sharey=False)
    for ax, panel in zip(axes, panels):
        lower = panel["lower"]
        upper = panel["upper"]
        plot_branch(ax, lower, "#1F77B4", "o", "random-targeted lower branch")
        plot_branch(ax, upper, "#C62828", "v", "fdseed-targeted upper branch")
        ax.axhline(panel["fd"], color="#2E7D32", ls="--", lw=1.4,
                   label="canonical FD")
        latest_gap = None
        if lower and upper:
            latest_gap = upper[-1]["x_gb"]["mean"] - lower[-1]["x_gb"]["mean"]
            y_min = min([d["x_gb"]["mean"] for d in lower + upper]) - 0.01
            y_max = panel["fd"] + 0.015
            ax.set_ylim(y_min, y_max)
            ax.text(
                0.03, 0.06,
                f"latest bracket = {latest_gap:.4f}",
                transform=ax.transAxes,
                fontsize=8.4,
                bbox={"boxstyle": "round,pad=0.25", "fc": "white",
                      "ec": "0.85", "alpha": 0.9},
            )
        ax.set_title(f"T = {panel['T']} K, Xc = 0.10")
        ax.set_xlabel("continuation block")
        ax.grid(True, alpha=0.25, lw=0.5)
        ax.set_xlim(-0.25, max(len(lower), len(upper)) - 0.75)
        summary["panels"].append({
            "T": panel["T"],
            "fd_predicted": panel["fd"],
            "latest_upper_minus_lower": latest_gap,
            "lower_branch": [run_summary(d) for d in lower],
            "upper_branch": [run_summary(d) for d in upper],
        })
    axes[0].set_ylabel(r"GB Mg fraction $X_{\mathrm{GB}}$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8.7,
               framealpha=0.95, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("700 K and 800 K two-branch checks remain open", y=1.10)
    fig.tight_layout()
    fig.savefig(FIG5_PREFIX.with_suffix(".png"), dpi=180, bbox_inches="tight")
    plt.close(fig)
    FIG5_PREFIX.with_suffix(".json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    fig2 = figure_2()
    fig5 = figure_5()
    print(f"wrote {FIG2_PREFIX.with_suffix('.json')}")
    print(f"wrote {FIG2_PREFIX.with_suffix('.png')}")
    print(f"wrote {FIG5_PREFIX.with_suffix('.json')}")
    print(f"wrote {FIG5_PREFIX.with_suffix('.png')}")
    print("latest T500 Xc=0.10:",
          fig2["descending_upper_bound_points"][-1]["X_GB_mean"])
    for panel in fig5["panels"]:
        print(f"T{panel['T']} latest bracket:",
              panel["latest_upper_minus_lower"])


if __name__ == "__main__":
    main()
