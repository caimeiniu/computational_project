"""
Revised Figure c: UMA replaces EAM in Wagih's per-site framework (NOT in MC).
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.linewidth': 0.6,
    'mathtext.default': 'regular',
    'figure.dpi': 200,
})

PAL = {
    'blue':     '#3B7DD8',  'blue_l':   '#E8F0FA',
    'orange':   '#D97B3B',  'orange_l': '#FDF0E6',
    'green':    '#3BA55D',  'green_l':  '#E6F5EC',
    'purple':   '#7B5EA7',  'purple_l': '#F0EBF5',
    'red':      '#D94452',  'red_l':    '#FDECEE',
    'gray':     '#5A6068',  'gray_l':   '#F4F5F6',
    'gold':     '#C49A2A',  'gold_l':   '#FBF5E3',
    'text':     '#2C2C2C',  'subtext':  '#6B7280',
    'line':     '#B0B8C4',
}

def draw_box(ax, x, y, w, h, text, fill, edge=None, fontsize=8.5,
             text_color='white', bold=True, lw=0.8, radius=0.15):
    if edge is None:
        edge = fill
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                         boxstyle=f"round,pad={radius}", linewidth=lw,
                         edgecolor=edge, facecolor=fill, zorder=3)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, text, fontsize=fontsize, ha='center', va='center',
            color=text_color, fontweight=weight, zorder=4, linespacing=1.35)

def draw_box_light(ax, x, y, w, h, text, color, fontsize=8.5, bold=False):
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


fig, ax = plt.subplots(figsize=(12, 9))
ax.set_xlim(-0.5, 13.5)
ax.set_ylim(-1.5, 11.5)
ax.set_aspect('equal')
ax.axis('off')

# Panel label
ax.text(0.02, 0.97, 'c', transform=ax.transAxes, fontsize=14,
        fontweight='bold', color=PAL['text'], va='top', ha='left')

# ── Title ──
ax.text(6.5, 11.0,
        "UMA MLIP as Drop-in Replacement in Per-Site Framework",
        fontsize=13, ha='center', fontweight='bold', color=PAL['text'])

# ═══════════════════════════════════
# TOP HALF: side-by-side comparison
# ═══════════════════════════════════

ax.text(3.2, 10.2, "Classical (EAM)", fontsize=11.5, ha='center',
        fontweight='bold', color=PAL['orange'])
ax.text(9.8, 10.2, "UMA MLIP", fontsize=11.5, ha='center',
        fontweight='bold', color=PAL['red'])

ax.plot([6.5, 6.5], [5.7, 10.3], '--', color=PAL['line'], lw=0.8, alpha=0.6)

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
    y = 9.3 - i * 1.1
    draw_box_light(ax, 3.2, y, 4.5, 0.8, left, PAL['orange'])
    draw_box_light(ax, 9.8, y, 4.5, 0.8, right, PAL['red'])
    ax.text(6.5, y, "vs", fontsize=9, ha='center', va='center',
            color=PAL['gray'], fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.1', facecolor='white',
                      edgecolor=PAL['line'], linewidth=0.5))

# ═══════════════════════════════════
# BOTTOM HALF: per-site pipeline
# ═══════════════════════════════════

# Separator
ax.text(6.5, 5.35, "▼  Replace EAM with UMA in Wagih's per-site framework  ▼",
        fontsize=9.5, ha='center', color=PAL['red'], fontweight='bold')

# ── Pipeline (mirrors Fig a structure) ──
# Row 1: structure
yr1 = 4.3
draw_box(ax, 2.0, yr1, 2.6, 0.75, "Polycrystal\n+ anneal", PAL['blue'], fontsize=8)
draw_box(ax, 5.0, yr1, 2.2, 0.75, "Identify GB\nsites (a-CNA)", PAL['blue'], fontsize=8)
arrow(ax, 3.3, yr1, 3.9, yr1)

# Row 2: per-site — this is where UMA comes in
yr2 = 3.0
arrow(ax, 5.0, 3.9, 5.0, 3.4)

# EAM version (left, faded out)
draw_box(ax, 2.5, yr2, 3.2, 0.75, "Per-site relax\n(LAMMPS + EAM)", '#C4C4C4',
         fontsize=8, text_color='#999999')

# UMA version (right, highlighted)
draw_box(ax, 7.8, yr2, 3.2, 0.75, "Per-site relax\n(ASE + UMA)", PAL['red'], fontsize=8.5)
# Glow
glow = FancyBboxPatch((7.8 - 1.7, yr2 - 0.45), 3.4, 0.9,
                       boxstyle="round,pad=0.15", linewidth=2.5,
                       edgecolor=PAL['red'], facecolor='none', alpha=0.25, zorder=2)
ax.add_patch(glow)

# Arrow showing replacement
ax.annotate('', xy=(6.1, yr2), xytext=(4.2, yr2),
            arrowprops=dict(arrowstyle='-|>', lw=2, color=PAL['red'],
                            connectionstyle='arc3,rad=0'))
ax.text(5.15, yr2 + 0.35, "replace", fontsize=8, ha='center',
        color=PAL['red'], fontweight='bold')

# Row 3: compute ΔE
yr3 = 1.6
arrow(ax, 7.8, 2.6, 7.8, 2.05)
draw_box(ax, 5.0, yr3, 3.0, 0.75, "Compute ΔE_seg(i)\nfor all GB sites", PAL['green'],
         fontsize=8)
draw_box(ax, 9.5, yr3, 3.2, 0.75, "Segregation spectrum\n+ Fermi-Dirac → X_GB", PAL['green'],
         fontsize=8)
arrow(ax, 6.5, yr3, 7.9, yr3)

# Equation
ax.text(6.5, 0.45,
        r"$\Delta E_i^{\rm seg} = E_{\rm GB,i}^{\rm UMA} - E_{\rm bulk}^{\rm UMA}$"
        r"$\qquad\longrightarrow\qquad$"
        r"$P_i = \left[1 + \frac{1-X^c}{X^c}\exp\!\left(\frac{\Delta E_i}{k_BT}\right)\right]^{-1}$",
        fontsize=10.5, ha='center', va='center', color=PAL['text'],
        bbox=dict(boxstyle='round,pad=0.35', facecolor=PAL['green_l'],
                  edgecolor=PAL['green'], linewidth=1.0))

# ── Key questions ──
q_box = FancyBboxPatch((0.5, -1.4), 12.0, 1.3,
                        boxstyle="round,pad=0.2", linewidth=1.0,
                        edgecolor=PAL['red'], facecolor=PAL['red_l'],
                        alpha=0.3, zorder=0)
ax.add_patch(q_box)
ax.text(1.0, -0.35, "Key questions", fontsize=8.5, fontweight='bold', color=PAL['red'])

questions = [
    r"① Same structure, same sites — does UMA ΔE_seg match EAM ΔE_seg?",
    r"② Can UMA predict segregation for alloys with NO available EAM potential?",
    r"③ Is UMA fast enough to skip the ML acceleration step entirely?",
]
for i, q in enumerate(questions):
    ax.text(6.5, -0.55 - i * 0.4, q, fontsize=8.5, ha='center', color=PAL['text'])

fig.tight_layout()
fig.savefig("/scratch/cainiu/UMA/Computational/docs/fig3_uma_integration.png",
            dpi=250, bbox_inches='tight', facecolor='white')
plt.close(fig)
print("Fig c (revised) saved.")
