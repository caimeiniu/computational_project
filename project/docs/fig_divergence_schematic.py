#!/usr/bin/env python3
"""Schematic divergence plot — matching PPT color scheme."""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# PPT-matching palette
C_BLUE   = "#4A90D9"   # solvent blue (for Wagih curve)
C_GOLD   = "#E8A735"   # solute gold (for MC curve)
C_RED    = "#F4A0A0"    # light rose
C_GREEN  = "#27AE60"   # process green (for annotations)
C_DARK   = "#1B2A4A"   # dark navy text
C_GREY   = "#D1D5DB"   # light border
C_SHADE  = "#F0F4F8"   # very subtle shade

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 13,
    "axes.linewidth": 0,
})

fig, ax = plt.subplots(figsize=(8, 5.5))

# Smooth curves
x = np.linspace(0, 1, 300)
y_wagih = 0.62 * (1 - np.exp(-3.5 * x))
y_mc = 0.62 * (1 - np.exp(-3.5 * x)) * (1 - 0.45 * x**1.8)

# Deviation region — subtle shading
x_crit = 0.38
ax.axvspan(x_crit, 1.05, alpha=0.5, color=C_SHADE, zorder=0)

# Curves
ax.plot(x, y_wagih, "--", color=C_BLUE, lw=3, zorder=3,
        label="Wagih (independent-site)")
ax.plot(x, y_mc, "-", color=C_RED, lw=3, zorder=3,
        label="MC (with interactions)")

# Markers
xm = np.array([0.05, 0.15, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95])
yw_m = 0.62 * (1 - np.exp(-3.5 * xm))
ym_m = 0.62 * (1 - np.exp(-3.5 * xm)) * (1 - 0.45 * xm**1.8)

ax.plot(xm, yw_m, "s", color=C_BLUE, ms=9, markerfacecolor="white",
        mew=2.5, zorder=5)
ax.plot(xm, ym_m, "o", color=C_RED, ms=9, markerfacecolor="white",
        mew=2.5, zorder=5)


# ── Key finding 1: Critical concentration ──
ax.axvline(x_crit, color=C_GREEN, lw=1.8, ls=":", alpha=0.8, zorder=1)
ax.annotate(r"$X_{tot}^{crit}$",
            xy=(x_crit, 0.01),
            xytext=(x_crit + 0.03, 0.13),
            fontsize=16, color=C_GREEN, fontweight="bold", ha="center",
            arrowprops=dict(arrowstyle="-|>", color=C_GREEN, lw=2.2))


# Region labels
ax.text(0.15, 0.50, "Agreement", fontsize=13, color=C_BLUE,
        ha="center", fontweight="bold", alpha=0.5)
ax.text(0.75, 0.08, "Deviation", fontsize=13, color=C_RED,
        ha="center", fontweight="bold", alpha=0.6)

# Axis labels
ax.set_xlabel(r"Total solute concentration $X_{tot}$", fontsize=12,
              color=C_DARK, labelpad=10)
ax.set_ylabel(r"GB solute concentration $X_{GB}$", fontsize=12,
              color=C_DARK, labelpad=10)

# Ticks
ax.set_xticks([0.1, 0.9])
ax.set_xticklabels(["Low", "High"], fontsize=14, color=C_DARK)
ax.set_yticks([])
ax.set_xlim(-0.02, 1.06)
ax.set_ylim(0, 0.68)

# Spines
for spine in ax.spines.values():
    spine.set_visible(False)
ax.spines["bottom"].set_visible(True)
ax.spines["bottom"].set_color(C_GREY)
ax.spines["bottom"].set_linewidth(1.5)
ax.spines["left"].set_visible(True)
ax.spines["left"].set_color(C_GREY)
ax.spines["left"].set_linewidth(1.5)

ax.tick_params(axis="x", colors=C_DARK, length=0)
ax.tick_params(axis="y", length=0)

# Legend
legend = ax.legend(fontsize=8, loc="upper left", framealpha=0.95,
                   edgecolor=C_GREY, fancybox=True)
legend.get_frame().set_linewidth(1)

fig.savefig("project/docs/fig_divergence_schematic.png", dpi=200,
            bbox_inches="tight", facecolor="white", pad_inches=0.2)
plt.close(fig)
print("Saved -> project/docs/fig_divergence_schematic.png")
