"""
Nature-style pipeline figures for GB segregation project.
Three figures:
  Fig 1: Wagih et al. — per-site ΔE + ML prediction framework
  Fig 2: MC verification — atomistic ground truth
  Fig 3: UMA MLIP — universal potential extension
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import matplotlib.patheffects as pe
import numpy as np

# ── Nature-style global settings ──
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.linewidth': 0.6,
    'mathtext.default': 'regular',
    'figure.dpi': 200,
})

# ── Color palette: muted, professional ──
PAL = {
    'blue':     '#3B7DD8',
    'blue_l':   '#E8F0FA',
    'orange':   '#D97B3B',
    'orange_l': '#FDF0E6',
    'green':    '#3BA55D',
    'green_l':  '#E6F5EC',
    'purple':   '#7B5EA7',
    'purple_l': '#F0EBF5',
    'red':      '#D94452',
    'red_l':    '#FDECEE',
    'gray':     '#5A6068',
    'gray_l':   '#F4F5F6',
    'gold':     '#C49A2A',
    'gold_l':   '#FBF5E3',
    'text':     '#2C2C2C',
    'subtext':  '#6B7280',
    'line':     '#B0B8C4',
}

def draw_box(ax, x, y, w, h, text, fill, edge=None, fontsize=8.5,
             text_color='white', alpha=1.0, bold=True, radius=0.15, lw=0.8):
    if edge is None:
        edge = fill
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle=f"round,pad={radius}", linewidth=lw,
                         edgecolor=edge, facecolor=fill, alpha=alpha, zorder=3)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, text, fontsize=fontsize, ha='center', va='center',
            color=text_color, fontweight=weight, zorder=4, linespacing=1.35)

def draw_box_light(ax, x, y, w, h, text, color, fontsize=8.5, bold=False):
    """Light-fill box with colored border — Nature style."""
    # Map color to light version
    light_map = {
        PAL['blue']: PAL['blue_l'], PAL['orange']: PAL['orange_l'],
        PAL['green']: PAL['green_l'], PAL['purple']: PAL['purple_l'],
        PAL['red']: PAL['red_l'], PAL['gold']: PAL['gold_l'],
        PAL['gray']: PAL['gray_l'],
    }
    fill = light_map.get(color, '#F8F8F8')
    draw_box(ax, x, y, w, h, text, fill=fill, edge=color, fontsize=fontsize,
             text_color=PAL['text'], bold=bold, lw=1.2)

def arrow(ax, x1, y1, x2, y2, color=None, lw=1.2, style='-|>'):
    if color is None:
        color = PAL['line']
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=style, lw=lw, color=color,
                                shrinkA=2, shrinkB=2))

def arrow_label(ax, x1, y1, x2, y2, label, offset=(0.15, 0.08), color=None,
                fontsize=7.5):
    if color is None:
        color = PAL['subtext']
    arrow(ax, x1, y1, x2, y2, color=PAL['line'])
    mx, my = (x1 + x2) / 2 + offset[0], (y1 + y2) / 2 + offset[1]
    ax.text(mx, my, label, fontsize=fontsize, color=color, style='italic',
            ha='left', va='center')

def panel_label(ax, label, x=0.02, y=0.97):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=14, fontweight='bold',
            color=PAL['text'], va='top', ha='left')

def setup(ax, xlim=(0, 12), ylim=(0, 10)):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect('equal')
    ax.axis('off')


# ╔════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 1: Wagih et al. — Per-site + ML framework               ║
# ╚════════════════════════════════════════════════════════════════════╝

fig1, ax = plt.subplots(figsize=(11, 7.5))
setup(ax, xlim=(-0.5, 13), ylim=(-0.8, 9.5))
panel_label(ax, 'a')

# Title
ax.text(6.25, 9.2, "Wagih et al. (2020) — Per-Site Segregation Energy + ML Prediction",
        fontsize=12, ha='center', fontweight='bold', color=PAL['text'])

# ── Row 1: Structure ──
y1 = 7.8
draw_box(ax, 1.8, y1, 2.8, 0.85, "Polycrystal\nVoronoi, 16 grains", PAL['blue'])
draw_box(ax, 5.0, y1, 2.4, 0.85, "Anneal + relax\n(LAMMPS MD)", PAL['blue'])
draw_box(ax, 8.2, y1, 2.4, 0.85, "Identify GB sites\n(a-CNA)", PAL['blue'])
arrow(ax, 3.2, y1, 3.8, y1)
arrow(ax, 6.2, y1, 7.0, y1)

# Small annotation
ax.text(1.8, 7.15, "20×20×20 nm³", fontsize=7, ha='center', color=PAL['subtext'])

# ── Row 2: Per-site computation ──
y2 = 6.1
arrow(ax, 8.2, 7.35, 8.2, 6.55)
draw_box(ax, 8.2, y2, 3.2, 0.85, "For each GB site i:\nsubstitute 1 solute → relax",
         PAL['orange'])
draw_box(ax, 4.2, y2, 2.8, 0.85, "Compute ΔE_seg(i)\nfor all ~10⁵ sites", PAL['orange'])
arrow(ax, 6.6, y2, 5.6, y2)

# Cost callout
ax.text(10.5, 5.4, "~10⁵ relaxations", fontsize=7.5, color=PAL['red'],
        ha='center', style='italic')

# Equation box for ΔE
eq_y = 4.55
ax.text(4.2, eq_y, r"$\Delta E_i^{\rm seg} = E_{\rm GB,i}^{\rm solute} - E_{\rm bulk}^{\rm solute}$",
        fontsize=11, ha='center', va='center', color=PAL['text'],
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                  edgecolor=PAL['orange'], linewidth=1.0))

# ── Row 3: ML ──
y3 = 3.3
draw_box(ax, 1.5, y3, 2.2, 0.85, "SOAP\n(1015 features)", PAL['purple'])
draw_box(ax, 4.3, y3, 2.4, 0.85, "Linear regression\nor PCA + k-means", PAL['purple'])
draw_box(ax, 7.5, y3, 2.8, 0.85, "Predict ΔE_seg(i)\nfor any alloy/site", PAL['purple'])
arrow(ax, 2.6, y3, 3.1, y3)
arrow(ax, 5.5, y3, 6.1, y3)

# ML annotation
ax.text(4.3, 2.6, "100 training points\n(accelerated model)", fontsize=7,
        ha='center', color=PAL['subtext'], style='italic')

# ── Row 4: Output ──
y4 = 1.6
arrow(ax, 7.5, 2.85, 7.5, 2.1)
draw_box(ax, 5.0, y4, 3.0, 0.85, "Segregation\nenergy spectrum", PAL['green'])
draw_box(ax, 9.0, y4, 3.0, 0.85, "Equilibrium X_GB\n(Fermi-Dirac)", PAL['green'])
arrow(ax, 6.5, y4, 7.5, y4)

# Fermi-Dirac equation
eq_fd_y = 0.35
ax.text(7.0, eq_fd_y,
        r"$P_i = \left[1 + \frac{1-X^c}{X^c}\,\exp\!\left(\frac{\Delta E_i}{k_BT}\right)\right]^{-1}$"
        r"$\qquad X^{\rm GB} = \frac{1}{N_{\rm GB}}\sum_i P_i$",
        fontsize=10, ha='center', va='center', color=PAL['text'],
        bbox=dict(boxstyle='round,pad=0.35', facecolor=PAL['green_l'],
                  edgecolor=PAL['green'], linewidth=1.0))

# Assumptions sidebar
sidebar_x = 11.5
ax.text(sidebar_x, 5.2, "Assumptions", fontsize=8, fontweight='bold',
        ha='center', color=PAL['gray'])
assumptions = ["Dilute limit", "No solute–solute\ninteraction", "Static (0 K)", "Potential-\ndependent"]
for i, txt in enumerate(assumptions):
    ay = 4.5 - i * 0.7
    ax.text(sidebar_x, ay, f"• {txt}", fontsize=7, ha='center',
            color=PAL['subtext'], linespacing=1.1)

# Dashed box around assumptions
sidebar_box = FancyBboxPatch((10.3, 2.2), 2.4, 3.5,
                              boxstyle="round,pad=0.15", linewidth=0.8,
                              edgecolor=PAL['line'], facecolor='white',
                              alpha=0.5, linestyle='--', zorder=1)
ax.add_patch(sidebar_box)

fig1.tight_layout()
fig1.savefig("/scratch/cainiu/UMA/Computational/docs/fig1_paper_pipeline.png",
             dpi=250, bbox_inches='tight', facecolor='white')
plt.close(fig1)
print("Fig 1 done.")


# ╔════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 2: MC Verification — Atomistic Ground Truth             ║
# ╚════════════════════════════════════════════════════════════════════╝

fig2, ax = plt.subplots(figsize=(11, 7.5))
setup(ax, xlim=(-0.5, 13), ylim=(-0.8, 9.5))
panel_label(ax, 'b')

ax.text(6.25, 9.2, "MC Verification — Atomistic Ground Truth for Wagih's Predictions",
        fontsize=12, ha='center', fontweight='bold', color=PAL['text'])

# ── Row 1: Setup ──
y1 = 7.8
draw_box(ax, 2.0, y1, 2.8, 0.85, "Nanocrystalline\nstructure", PAL['blue'])
draw_box(ax, 5.5, y1, 3.0, 0.85, "Distribute solutes\n(random, X_tot = 5%)", PAL['blue'])
arrow(ax, 3.4, y1, 4.0, y1)

# ── MC Loop region ──
loop = FancyBboxPatch((0.3, 2.6), 9.6, 4.5,
                       boxstyle="round,pad=0.25", linewidth=1.5,
                       edgecolor=PAL['orange'], facecolor=PAL['orange_l'],
                       alpha=0.3, linestyle=(0, (5, 3)), zorder=0)
ax.add_patch(loop)
ax.text(5.1, 6.85, "Monte Carlo loop  (~10⁴ – 10⁶ steps)", fontsize=9.5,
        ha='center', fontweight='bold', color=PAL['orange'])

arrow(ax, 5.5, 7.35, 5.1, 6.95, color=PAL['orange'])

# Inside loop
y_mc1 = 5.8
draw_box_light(ax, 2.0, y_mc1, 2.5, 0.8, "①  Random swap\nsolute ↔ solvent", PAL['orange'], bold=True)
draw_box_light(ax, 5.1, y_mc1, 2.5, 0.8, "②  Energy minimize\n(LAMMPS + EAM)", PAL['orange'], bold=True)
draw_box_light(ax, 8.2, y_mc1, 2.3, 0.8, "③  Compute ΔE\nof swap", PAL['orange'], bold=True)
arrow(ax, 3.25, y_mc1, 3.85, y_mc1, color=PAL['orange'])
arrow(ax, 6.35, y_mc1, 7.05, y_mc1, color=PAL['orange'])

# Metropolis box
y_met = 4.4
draw_box_light(ax, 5.1, y_met, 5.5, 0.8, "", PAL['gold'])
ax.text(5.1, y_met,
        r"④  Metropolis:  accept if $\Delta E < 0$;   else accept with $p = e^{-\Delta E / k_BT}$",
        fontsize=8.5, ha='center', va='center', color=PAL['text'], fontweight='bold')

arrow(ax, 8.2, 5.4, 8.2, 4.8, color=PAL['orange'])

# Loop-back arrow (curved)
ax.annotate('', xy=(2.0, 5.1), xytext=(2.0, 4.0),
            arrowprops=dict(arrowstyle='-|>', lw=1.3, color=PAL['orange'],
                            connectionstyle='arc3,rad=-0.6'))
ax.text(0.7, 4.5, "repeat", fontsize=7.5, color=PAL['orange'],
        ha='center', style='italic', rotation=90)

# ── Row 3: Convergence → Output ──
arrow(ax, 5.1, 2.6, 5.1, 2.05, color=PAL['orange'])
ax.text(5.5, 2.2, "converge", fontsize=7.5, color=PAL['subtext'], style='italic')

y_out = 1.3
draw_box(ax, 2.5, y_out, 3.0, 0.8, "Equilibrium GB\nsolute distribution", PAL['green'])
draw_box(ax, 6.5, y_out, 3.0, 0.8, "X_GB^MC  vs  X_GB^Wagih",
         PAL['green'])
arrow(ax, 4.0, y_out, 5.0, y_out)

# Comparison annotation
ax.text(6.5, 0.55,
        r"Agree $\rightarrow$ validates per-site framework"
        "\n"
        r"Disagree $\rightarrow$ solute–solute effects matter",
        fontsize=8, ha='center', color=PAL['subtext'], linespacing=1.5,
        bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                  edgecolor=PAL['line'], linewidth=0.6))

# Advantages sidebar
sidebar_x = 11.2
ax.text(sidebar_x, 7.8, "Why MC?", fontsize=8.5, fontweight='bold',
        ha='center', color=PAL['blue'])
advantages = [
    "Finite concentration",
    "Solute–solute\ninteractions",
    "Temperature-\ndependent",
    "No dilute-limit\nassumption",
    "Direct equilibrium\nsampling",
]
for i, txt in enumerate(advantages):
    ay = 7.0 - i * 0.85
    ax.text(sidebar_x, ay, f"✓ {txt}", fontsize=7, ha='center',
            color=PAL['green'], linespacing=1.1)

sidebar_box2 = FancyBboxPatch((10.0, 3.1), 2.5, 5.2,
                               boxstyle="round,pad=0.15", linewidth=0.8,
                               edgecolor=PAL['line'], facecolor='white',
                               alpha=0.5, linestyle='--', zorder=0)
ax.add_patch(sidebar_box2)

fig2.tight_layout()
fig2.savefig("/scratch/cainiu/UMA/Computational/docs/fig2_hmc_pipeline.png",
             dpi=250, bbox_inches='tight', facecolor='white')
plt.close(fig2)
print("Fig 2 done.")


# ╔════════════════════════════════════════════════════════════════════╗
# ║  FIGURE 3: UMA MLIP Integration                                 ║
# ╚════════════════════════════════════════════════════════════════════╝

fig3, ax = plt.subplots(figsize=(11, 8))
setup(ax, xlim=(-0.5, 13), ylim=(-1.2, 10))
panel_label(ax, 'c')

ax.text(6.25, 9.5, "UMA MLIP — Universal Potential as Drop-in Replacement",
        fontsize=12, ha='center', fontweight='bold', color=PAL['text'])

# ── Top: side-by-side comparison ──
# Column headers
ax.text(3.2, 8.8, "Classical (EAM)", fontsize=11, ha='center',
        fontweight='bold', color=PAL['orange'])
ax.text(9.3, 8.8, "UMA MLIP", fontsize=11, ha='center',
        fontweight='bold', color=PAL['red'])

# Divider
ax.plot([6.25, 6.25], [4.2, 8.9], '--', color=PAL['line'], lw=0.8, alpha=0.6)

# Row data
rows = [
    ("System-specific potential\n(fitted per alloy pair)",
     "Universal model\n(one model, all elements)"),
    ("NIST repository\n(259 alloys, variable quality)",
     "Trained on OMat DFT database\n(~118M structures)"),
    ("No potential → no prediction\n(limited coverage)",
     "Works for ANY binary alloy\n(no system-specific fitting)"),
    ("Different EAMs → different ΔE\n(potential-dependent)",
     "Single model → consistent\n(uniform accuracy)"),
]

for i, (left, right) in enumerate(rows):
    y = 8.0 - i * 1.1
    draw_box_light(ax, 3.2, y, 4.2, 0.8, left, PAL['orange'])
    draw_box_light(ax, 9.3, y, 4.2, 0.8, right, PAL['red'])
    # vs label
    ax.text(6.25, y, "vs", fontsize=9, ha='center', va='center',
            color=PAL['gray'], fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor=PAL['line'], linewidth=0.5))

# ── Middle: drop-in replacement arrow ──
y_drop = 3.6
ax.annotate('', xy=(6.25, y_drop), xytext=(6.25, 4.1),
            arrowprops=dict(arrowstyle='-|>', lw=2, color=PAL['red']))
ax.text(6.25, 3.95, "Drop-in replacement in MC pipeline",
        fontsize=9, ha='center', color=PAL['red'], fontweight='bold')

# ── Bottom: unified pipeline ──
y_pipe = 2.6
bw = 2.3
draw_box(ax, 1.8, y_pipe, bw, 0.75, "Nanocrystalline\nstructure", PAL['blue'], fontsize=8)
draw_box(ax, 4.5, y_pipe, bw, 0.75, "MC swap\nsimulation", PAL['gold'],
         text_color=PAL['text'], fontsize=8)

# The key box — highlighted
draw_box(ax, 7.4, y_pipe, 2.6, 0.75, "", PAL['red'])
ax.text(7.4, y_pipe, "Energy: UMA\n(replaces EAM)", fontsize=8.5,
        ha='center', va='center', color='white', fontweight='bold', zorder=4)
# Glow effect
glow = FancyBboxPatch((7.4 - 1.4, y_pipe - 0.45), 2.8, 0.9,
                       boxstyle="round,pad=0.15", linewidth=2.5,
                       edgecolor=PAL['red'], facecolor='none', alpha=0.3, zorder=2)
ax.add_patch(glow)

draw_box(ax, 10.5, y_pipe, 2.3, 0.75, "Equilibrium\nX_GB", PAL['green'], fontsize=8)

arrow(ax, 2.95, y_pipe, 3.35, y_pipe)
arrow(ax, 5.65, y_pipe, 6.1, y_pipe)
arrow(ax, 8.7, y_pipe, 9.35, y_pipe)

# ── Bottom: key questions ──
y_q = 0.9
questions = [
    r"① Can UMA match EAM segregation energies ($\Delta E_{\rm seg}$)?",
    r"② Can UMA predict alloys where NO fitted potential exists?",
    r"③ Is UMA accurate enough for quantitative segregation spectra?",
]
for i, q in enumerate(questions):
    ax.text(6.25, y_q - i * 0.55, q, fontsize=8.5, ha='center',
            color=PAL['text'])

# Question box background
q_box = FancyBboxPatch((1.0, -0.65), 10.5, 2.05,
                        boxstyle="round,pad=0.2", linewidth=1.0,
                        edgecolor=PAL['red'], facecolor=PAL['red_l'],
                        alpha=0.3, zorder=0)
ax.add_patch(q_box)
ax.text(1.3, 1.2, "Key questions", fontsize=8, fontweight='bold', color=PAL['red'])

fig3.tight_layout()
fig3.savefig("/scratch/cainiu/UMA/Computational/docs/fig3_uma_integration.png",
             dpi=250, bbox_inches='tight', facecolor='white')
plt.close(fig3)
print("Fig 3 done.")

print("\nAll 3 Nature-style figures saved to docs/")
