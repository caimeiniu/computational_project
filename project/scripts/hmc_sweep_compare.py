#!/usr/bin/env python3
"""Compare X_GB^HMC vs X_GB^FD across the (T=500 K, X_c) sweep.

Reads `output/hmc_T*_Xc*_xgb.json` produced by hmc_xgb_timeseries.py and
the FD prediction grid `output/fd_curves_200A_tight.json`. Emits a
two-panel plot (log-log X_GB vs X_c, plus enrichment ratio) and a CSV summary.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("/cluster/home/cainiu/Computational_modeling/project/output")

# inputs
hmc_files = [
    OUT / "hmc_T500_Xc5e-4_xgb.json",
    OUT / "hmc_T500_Xc5e-3_xgb.json",
    OUT / "hmc_T500_Xc5e-2_xgb.json",
]
fd_curves = json.load(open(OUT / "fd_curves_200A_tight.json"))
fd_xc = np.array(fd_curves["X_c"])
fd_500 = np.array(fd_curves["ours"]["500"])

# load HMC
rows = []
for p in hmc_files:
    d = json.load(open(p))
    rows.append({
        "Xc": d["X_c"],
        "X_HMC": d["x_gb"]["mean"],
        "lo": d["x_gb"]["ci95_lo"],
        "hi": d["x_gb"]["ci95_hi"],
        "fd": d["fd_predicted"],
        "accept": d["swap"]["accept_rate_overall"],
        "accepts": d["swap"]["total_accepts"],
        "first_xgb": json.load(open(p))["series"]["x_gb"][0],
        "last_xgb": json.load(open(p))["series"]["x_gb"][-1],
    })

# CSV
csv_path = OUT / "hmc_sweep_T500.csv"
with csv_path.open("w") as f:
    f.write("Xc,X_HMC,X_HMC_lo,X_HMC_hi,X_FD,gap,enrichment_HMC,enrichment_FD,"
            "accept_rate,total_accepts,xgb_first,xgb_last\n")
    for r in rows:
        f.write(f"{r['Xc']:.4e},{r['X_HMC']:.5f},{r['lo']:.5f},{r['hi']:.5f},"
                f"{r['fd']:.5f},{(r['X_HMC']-r['fd']):+.5f},"
                f"{r['X_HMC']/r['Xc']:.3f},{r['fd']/r['Xc']:.3f},"
                f"{r['accept']:.4f},{r['accepts']},"
                f"{r['first_xgb']:.5f},{r['last_xgb']:.5f}\n")
print(f"wrote {csv_path}")

# plot
xs_hmc = np.array([r["Xc"] for r in rows])
y_hmc = np.array([r["X_HMC"] for r in rows])
y_lo = np.array([r["lo"] for r in rows])
y_hi = np.array([r["hi"] for r in rows])
y_fd = np.array([r["fd"] for r in rows])

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.loglog(fd_xc, fd_500, "-", color="C0", label="X_GB^FD (ours, T=500 K)")
ax.loglog(fd_xc, fd_xc, "--", color="gray", lw=0.7, label="X_GB = X_c (no segregation)")
ax.errorbar(xs_hmc, y_hmc,
            yerr=[y_hmc - y_lo, y_hi - y_hmc],
            fmt="o", color="C3", ms=8, capsize=4,
            label="X_GB^HMC (50 ps PROD)")
ax.set_xlabel("bulk Mg fraction X_c")
ax.set_ylabel("X_GB")
ax.set_title("Al(Mg) GB occupancy — T=500 K\nHMC under-sampled: X_GB^HMC ≈ X_c")
ax.legend(loc="best", fontsize=9)
ax.grid(True, which="both", alpha=0.3)

ax = axes[1]
ax.semilogx(xs_hmc, y_hmc / xs_hmc, "o-", color="C3", ms=8,
            label="HMC enrichment X_GB^HMC / X_c")
ax.semilogx(xs_hmc, y_fd / xs_hmc, "s--", color="C0", ms=8,
            label="FD enrichment X_GB^FD / X_c")
ax.axhline(1, color="gray", lw=0.6, ls=":")
ax.set_xlabel("bulk Mg fraction X_c")
ax.set_ylabel("enrichment X_GB / X_c")
ax.set_title("HMC ≈ 1 across all X_c → no segregation reached\n"
             "FD predicts 7.4× → 159× enrichment (saturating)")
ax.set_yscale("log")
ax.legend(loc="best", fontsize=9)
ax.grid(True, which="both", alpha=0.3)

fig.tight_layout()
out_png = OUT / "hmc_vs_fd_T500_sweep.png"
fig.savefig(out_png, dpi=130)
print(f"wrote {out_png}")

print("\nsummary:")
print(f"{'X_c':>10} {'X_HMC':>10} {'X_FD':>10} {'enrich_HMC':>12} {'enrich_FD':>11} {'accepts':>10}")
for r in rows:
    print(f"{r['Xc']:>10.1e} {r['X_HMC']:>10.4f} {r['fd']:>10.4f} "
          f"{r['X_HMC']/r['Xc']:>12.3f} {r['fd']/r['Xc']:>11.3f} {r['accepts']:>10}")
