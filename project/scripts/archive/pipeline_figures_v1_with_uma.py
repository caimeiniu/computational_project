"""
Generate three technical pipeline figures for presentation:
  Fig 1: Paper's approach (Wagih et al.) — per-site ΔE + SOAP + ML
  Fig 2: Our HMC approach — MC swap + EAM relaxation
  Fig 3: UMA integration — replace EAM with universal MLIP
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# ── Color palette (consistent across all 3 figures) ──
C_STRUCT  = '#4A90D9'   # blue — structure generation
C_CALC    = '#E8833A'   # orange — energy calculation
C_ML      = '#7B68AE'   # purple — machine learning
C_RESULT  = '#50B86C'   # green — results/output
C_UMA     = '#E05555'   # red — UMA specific
C_LIGHT   = '#F5F5F5'   # light gray background
C_ARROW   = '#333333'   # dark arrows
C_ACCENT  = '#FFD700'   # gold accent

def draw_box(ax, x, y, w, h, text, color, fontsize=11, text_color='white', alpha=0.95, bold=True):
    """Draw a rounded box with centered text."""
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle="round,pad=0.15", linewidth=2,
                         edgecolor='white', facecolor=color, alpha=alpha, zorder=3)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, text, fontsize=fontsize, ha='center', va='center',
            color=text_color, fontweight=weight, zorder=4, linespacing=1.4)

def draw_arrow(ax, x1, y1, x2, y2, color=C_ARROW, style='->', lw=2.5):
    """Draw an arrow between two points."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, lw=lw, color=color,
                                connectionstyle='arc3,rad=0'))

def draw_arrow_label(ax, x1, y1, x2, y2, label, color=C_ARROW, fontsize=9):
    """Draw arrow with a label at midpoint."""
    draw_arrow(ax, x1, y1, x2, y2, color=color)
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    ax.text(mx + 0.3, my, label, fontsize=fontsize, color=color,
            ha='left', va='center', style='italic')

def setup_ax(ax, title, xlim=(0, 12), ylim=(0, 10)):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15, color='#222222')


# ╔════════════════════════════════════════════════════════════╗
# ║  FIGURE 1: Paper's Pipeline (Wagih et al. 2020)          ║
# ╚════════════════════════════════════════════════════════════╝

fig1, ax1 = plt.subplots(1, 1, figsize=(14, 8))
setup_ax(ax1, "Wagih et al. (2020): Per-Site Segregation Energy + ML Prediction",
         xlim=(0, 14), ylim=(0, 9.5))

# Row 1: Structure generation
draw_box(ax1, 2.5, 8.2, 3.8, 1.0, "Polycrystal\n(Voronoi, 16 grains)", C_STRUCT)
draw_box(ax1, 7.0, 8.2, 3.0, 1.0, "Anneal + Relax\n(LAMMPS MD)", C_STRUCT)
draw_box(ax1, 11.0, 8.2, 2.8, 1.0, "Identify GB sites\n(a-CNA)", C_STRUCT)
draw_arrow(ax1, 4.4, 8.2, 5.5, 8.2)
draw_arrow(ax1, 8.5, 8.2, 9.6, 8.2)

# Row 2: Per-site calculation (the expensive part)
draw_arrow(ax1, 11.0, 7.7, 11.0, 6.8)
draw_box(ax1, 11.0, 6.2, 3.2, 1.0, "For EACH GB site:\nplace 1 solute → relax", C_CALC)
draw_box(ax1, 7.0, 6.2, 3.2, 1.0, "Compute ΔE_seg\n= E_GB - E_bulk", C_CALC)
draw_arrow(ax1, 9.4, 6.2, 8.6, 6.2)

# Cost annotation
ax1.text(11.0, 5.3, "~10⁵ sites × relax each\n(very expensive!)", fontsize=9,
         ha='center', color='#CC4444', style='italic')

# Row 3: ML
draw_arrow(ax1, 7.0, 5.7, 7.0, 4.8)
draw_box(ax1, 3.5, 4.2, 3.5, 1.0, "SOAP descriptors\n(1015 features, r=6 Å)", C_ML)
draw_box(ax1, 7.5, 4.2, 2.5, 1.0, "Linear\nRegression", C_ML)
draw_box(ax1, 11.0, 4.2, 3.0, 1.0, "Predict ΔE_seg\nfor any alloy", C_ML)
draw_arrow(ax1, 5.25, 4.2, 6.25, 4.2)
draw_arrow(ax1, 8.75, 4.2, 9.5, 4.2)

# Row 4: Output
draw_arrow(ax1, 11.0, 3.7, 11.0, 2.8)
draw_box(ax1, 7.0, 2.2, 3.5, 1.0, "Segregation\nEnergy Spectrum", C_RESULT)
draw_box(ax1, 11.0, 2.2, 3.0, 1.0, "Equilibrium GB\nConcentration", C_RESULT)
draw_arrow(ax1, 8.75, 2.2, 9.5, 2.2)

# Annotation: Fermi-Dirac
ax1.text(9.8, 1.8, "Fermi-Dirac\nstatistics", fontsize=8, ha='center',
         color=C_RESULT, style='italic')

# Key limitation box
ax1.text(1.0, 1.2, "Key assumptions:\n• Dilute limit (1 solute at a time)\n"
         "• No solute-solute interaction\n• Static (0 K) calculation\n"
         "• Depends on classical potential quality",
         fontsize=9, ha='left', va='center',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFF3E0', edgecolor='#E8833A', linewidth=1.5))

fig1.tight_layout()
fig1.savefig("/scratch/cainiu/UMA/Computational/docs/fig1_paper_pipeline.png",
             dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig1)
print("Fig 1 saved.")


# ╔════════════════════════════════════════════════════════════╗
# ║  FIGURE 2: Our HMC Approach (LAMMPS + EAM)               ║
# ╚════════════════════════════════════════════════════════════╝

fig2, ax2 = plt.subplots(1, 1, figsize=(14, 8))
setup_ax(ax2, "Our Approach: Hybrid Monte Carlo Segregation Simulation",
         xlim=(0, 14), ylim=(0, 9.5))

# Row 1: Structure
draw_box(ax2, 2.5, 8.2, 3.8, 1.0, "Nanocrystalline\nstructure", C_STRUCT)
draw_box(ax2, 7.0, 8.2, 3.5, 1.0, "Distribute solutes\n(random, X_tot = 5%)", C_STRUCT)
draw_arrow(ax2, 4.4, 8.2, 5.25, 8.2)

# Row 2: MC loop (the core)
draw_arrow(ax2, 7.0, 7.7, 7.0, 6.8)

# MC loop box (large rounded rect)
loop_box = FancyBboxPatch((1.5, 2.8), 11.0, 4.2,
                          boxstyle="round,pad=0.3", linewidth=3,
                          edgecolor='#E8833A', facecolor='#FFF8F0', alpha=0.5, zorder=1,
                          linestyle='--')
ax2.add_patch(loop_box)
ax2.text(7.0, 6.7, "MC Loop (repeat ~10⁴ – 10⁶ steps)", fontsize=12,
         ha='center', color='#CC5500', fontweight='bold')

# Inside MC loop
draw_box(ax2, 3.5, 5.5, 3.2, 0.9, "Random swap\nsolute ↔ solvent", C_CALC)
draw_box(ax2, 7.5, 5.5, 2.8, 0.9, "Energy minimize\n(LAMMPS + EAM)", C_CALC)
draw_box(ax2, 11.0, 5.5, 2.5, 0.9, "Compute ΔE\nof swap", C_CALC)
draw_arrow(ax2, 5.1, 5.5, 6.1, 5.5)
draw_arrow(ax2, 8.9, 5.5, 9.75, 5.5)

# Metropolis
draw_arrow(ax2, 11.0, 5.05, 11.0, 4.3)
draw_box(ax2, 7.0, 3.8, 5.5, 0.9, "Metropolis: accept if ΔE<0,\n"
         "else accept with p = exp(−ΔE/kT)", '#D4A76A', text_color='#333333')

# Loop back arrow
ax2.annotate('', xy=(3.5, 4.8), xytext=(3.5, 3.8),
             arrowprops=dict(arrowstyle='->', lw=2, color='#CC5500',
                             connectionstyle='arc3,rad=-0.5'))
ax2.text(1.8, 4.3, "next\nstep", fontsize=9, ha='center', color='#CC5500', style='italic')

# Row 3: Output
draw_arrow(ax2, 7.0, 2.8, 7.0, 2.1)
draw_box(ax2, 3.5, 1.5, 3.5, 0.9, "Equilibrium GB\nsolute distribution", C_RESULT)
draw_box(ax2, 7.5, 1.5, 3.0, 0.9, "GB concentration\nvs T, X_tot", C_RESULT)
draw_box(ax2, 11.2, 1.5, 3.0, 0.9, "Compare with\nWagih prediction", C_RESULT)
draw_arrow(ax2, 5.25, 1.5, 6.0, 1.5)
draw_arrow(ax2, 9.0, 1.5, 9.7, 1.5)

# Advantages box
ax2.text(12.5, 8.5, "Advantages:\n• Finite concentration\n• Solute-solute interactions\n"
         "• Temperature effects\n• Direct equilibrium sampling",
         fontsize=9, ha='center', va='top',
         bbox=dict(boxstyle='round,pad=0.4', facecolor='#E8F5E9', edgecolor='#50B86C', linewidth=1.5))

fig2.tight_layout()
fig2.savefig("/scratch/cainiu/UMA/Computational/docs/fig2_hmc_pipeline.png",
             dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig2)
print("Fig 2 saved.")


# ╔════════════════════════════════════════════════════════════╗
# ║  FIGURE 3: UMA MLIP Integration                          ║
# ╚════════════════════════════════════════════════════════════╝

fig3, ax3 = plt.subplots(1, 1, figsize=(14, 9))
setup_ax(ax3, "UMA MLIP Integration: Universal Potential for GB Segregation",
         xlim=(0, 14), ylim=(0, 10.5))

# Top: comparison diagram — left = classical, right = UMA
# Dividing line
ax3.plot([7, 7], [4.5, 10], '--', color='gray', lw=1.5, alpha=0.5)
ax3.text(3.5, 9.8, "Classical (EAM)", fontsize=14, ha='center',
         fontweight='bold', color=C_CALC)
ax3.text(10.5, 9.8, "UMA MLIP", fontsize=14, ha='center',
         fontweight='bold', color=C_UMA)

# Left column: Classical
draw_box(ax3, 3.5, 8.8, 4.2, 0.8, "System-specific EAM potential\n(fitted per alloy pair)", C_CALC, fontsize=10)
draw_box(ax3, 3.5, 7.5, 4.2, 0.8, "NIST repository\n(259 alloys, variable quality)", C_CALC, fontsize=10)
draw_box(ax3, 3.5, 6.2, 4.2, 0.8, "Limited to available potentials\n(no potential → no prediction)", C_CALC, fontsize=10)
draw_box(ax3, 3.5, 4.9, 4.2, 0.8, "Potential-dependent results\n(different EAMs → different ΔE)", C_CALC, fontsize=10)
draw_arrow(ax3, 3.5, 8.4, 3.5, 7.9)
draw_arrow(ax3, 3.5, 7.1, 3.5, 6.6)
draw_arrow(ax3, 3.5, 5.8, 3.5, 5.3)

# Right column: UMA
draw_box(ax3, 10.5, 8.8, 4.2, 0.8, "Universal MLIP\n(one model for all elements)", C_UMA, fontsize=10)
draw_box(ax3, 10.5, 7.5, 4.2, 0.8, "Trained on massive DFT database\n(OMat, ~100M structures)", C_UMA, fontsize=10)
draw_box(ax3, 10.5, 6.2, 4.2, 0.8, "Works for ANY binary alloy\n(no system-specific fitting)", C_UMA, fontsize=10)
draw_box(ax3, 10.5, 4.9, 4.2, 0.8, "Consistent across alloy space\n(single model, uniform accuracy)", C_UMA, fontsize=10)
draw_arrow(ax3, 10.5, 8.4, 10.5, 7.9)
draw_arrow(ax3, 10.5, 7.1, 10.5, 6.6)
draw_arrow(ax3, 10.5, 5.8, 10.5, 5.3)

# VS arrows
for y in [8.8, 7.5, 6.2, 4.9]:
    ax3.text(7.0, y, "vs", fontsize=12, ha='center', va='center',
             color='gray', fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='gray'))

# Bottom: unified pipeline
draw_arrow(ax3, 7.0, 4.2, 7.0, 3.6)
ax3.text(7.0, 3.8, "▼  Drop-in replacement  ▼", fontsize=11, ha='center',
         color=C_UMA, fontweight='bold')

# Pipeline boxes
draw_box(ax3, 2.5, 2.5, 3.5, 0.9, "Nanocrystalline\nstructure", C_STRUCT, fontsize=10)
draw_box(ax3, 6.2, 2.5, 2.8, 0.9, "MC swap\nsimulation", '#D4A76A', text_color='#333', fontsize=10)

# The key replacement
draw_box(ax3, 9.8, 2.5, 2.8, 0.9, "Energy: UMA\n(replaces EAM)", C_UMA, fontsize=10)
draw_box(ax3, 13.0, 2.5, 1.8, 0.9, "Equilibrium\nresult", C_RESULT, fontsize=10)

draw_arrow(ax3, 4.25, 2.5, 4.8, 2.5)
draw_arrow(ax3, 7.6, 2.5, 8.4, 2.5)
draw_arrow(ax3, 11.2, 2.5, 12.1, 2.5)

# Key questions box
ax3.text(7.0, 0.8, "Key questions:  ① Can UMA match EAM segregation energies?  "
         "② Does UMA enable predictions for alloys without fitted potentials?  "
         "③ Is UMA accurate enough for quantitative segregation spectra?",
         fontsize=9.5, ha='center', va='center',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFEBEE', edgecolor=C_UMA, linewidth=2))

fig3.tight_layout()
fig3.savefig("/scratch/cainiu/UMA/Computational/docs/fig3_uma_integration.png",
             dpi=200, bbox_inches='tight', facecolor='white')
plt.close(fig3)
print("Fig 3 saved.")
print("\nAll 3 figures saved to docs/")
