"""
Demo: 2D illustration of the polycrystal generation + GB segregation workflow.
Generates a figure showing:
  (a) Voronoi tessellation with colored grains
  (b) a-CNA: bulk vs GB atom identification
  (c) Solute placement at GB site + relaxation → ΔE_seg
  (d) MC swap procedure
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.spatial import Voronoi, voronoi_plot_2d

# ── reproducible random seed ──
np.random.seed(42)

# ── parameters ──
box_size = 20.0     # nm
n_grains = 16
lattice_a = 0.405   # nm (Al)
n_atoms_per_row = int(box_size / lattice_a)

fig, axes = plt.subplots(2, 2, figsize=(14, 13))

# ═══════════════════════════════════════════
# Panel (a): Voronoi tessellation — grains
# ═══════════════════════════════════════════
ax = axes[0, 0]
ax.set_title("(a) Voronoi Polycrystal Generation", fontsize=13, fontweight='bold')

# Random grain centers
grain_centers = np.random.uniform(2, box_size - 2, size=(n_grains, 2))

# Mirror points for periodic Voronoi
mirrors = []
for dx in [-box_size, 0, box_size]:
    for dy in [-box_size, 0, box_size]:
        if dx == 0 and dy == 0:
            continue
        mirrors.append(grain_centers + np.array([dx, dy]))
all_points = np.vstack([grain_centers] + mirrors)

vor = Voronoi(all_points)

# Color each grain
colors = plt.cm.Set3(np.linspace(0, 1, n_grains))
for i, region_idx in enumerate(vor.point_region[:n_grains]):
    region = vor.regions[region_idx]
    if -1 in region or len(region) == 0:
        continue
    polygon = [vor.vertices[v] for v in region]
    polygon = np.array(polygon)
    ax.fill(*polygon.T, color=colors[i % n_grains], alpha=0.6, edgecolor='k', linewidth=1.5)

# Plot grain centers
ax.plot(grain_centers[:, 0], grain_centers[:, 1], 'k^', markersize=6, label='Grain centers')

ax.set_xlim(0, box_size)
ax.set_ylim(0, box_size)
ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_aspect('equal')
ax.legend(loc='upper right', fontsize=9)
ax.text(0.5, -0.5, f"{n_grains} grains, {box_size}×{box_size} nm² box",
        fontsize=10, ha='left', color='gray')

# ═══════════════════════════════════════════
# Panel (b): Atom-level view — bulk vs GB
# ═══════════════════════════════════════════
ax = axes[0, 1]
ax.set_title("(b) GB Site Identification (a-CNA)", fontsize=13, fontweight='bold')

# Generate a small 2D "lattice" to illustrate
# Use 2 grains with different orientations for clarity
region_size = 6.0  # nm
spacing = 0.35     # nm

# Grain 1: no rotation (left half)
x1 = np.arange(0, region_size / 2, spacing)
y1 = np.arange(0, region_size, spacing)
xx1, yy1 = np.meshgrid(x1, y1)
atoms1_x = xx1.ravel()
atoms1_y = yy1.ravel()

# Grain 2: rotated 25° (right half)
x2 = np.arange(region_size / 2, region_size, spacing)
y2 = np.arange(-1, region_size + 1, spacing)
xx2, yy2 = np.meshgrid(x2, y2)
angle = np.radians(25)
cx, cy = region_size * 0.75, region_size / 2
rot_x = cx + (xx2.ravel() - cx) * np.cos(angle) - (yy2.ravel() - cy) * np.sin(angle)
rot_y = cy + (xx2.ravel() - cx) * np.sin(angle) + (yy2.ravel() - cy) * np.cos(angle)
mask2 = (rot_x > region_size / 2 - 0.1) & (rot_x < region_size + 0.1) & \
        (rot_y > -0.1) & (rot_y < region_size + 0.1)
atoms2_x = rot_x[mask2]
atoms2_y = rot_y[mask2]

# Classify: atoms near the boundary (x ≈ region_size/2) as GB
all_x = np.concatenate([atoms1_x, atoms2_x])
all_y = np.concatenate([atoms1_y, atoms2_y])

gb_width = 0.6  # nm
is_gb = np.abs(all_x - region_size / 2) < gb_width

# Plot bulk atoms
ax.scatter(all_x[~is_gb], all_y[~is_gb], s=20, c='steelblue', alpha=0.7,
           label='Bulk (FCC)', zorder=2)
# Plot GB atoms
ax.scatter(all_x[is_gb], all_y[is_gb], s=25, c='red', alpha=0.9,
           label='GB sites', zorder=3)

# GB region shading
ax.axvspan(region_size / 2 - gb_width, region_size / 2 + gb_width,
           color='red', alpha=0.08, zorder=1)

# Labels
ax.text(region_size * 0.2, region_size * 0.95, "Grain 1\n(θ=0°)",
        fontsize=11, ha='center', fontweight='bold', color='steelblue')
ax.text(region_size * 0.78, region_size * 0.95, "Grain 2\n(θ=25°)",
        fontsize=11, ha='center', fontweight='bold', color='steelblue')

ax.set_xlim(-0.2, region_size + 0.2)
ax.set_ylim(-0.2, region_size + 0.2)
ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_aspect('equal')
ax.legend(loc='lower right', fontsize=9)

n_gb = np.sum(is_gb)
n_total = len(all_x)
ax.text(0.1, -0.5, f"GB fraction: {n_gb}/{n_total} = {n_gb/n_total:.1%}",
        fontsize=10, ha='left', color='gray')

# ═══════════════════════════════════════════
# Panel (c): Segregation energy calculation
# ═══════════════════════════════════════════
ax = axes[1, 0]
ax.set_title("(c) Segregation Energy ΔE_seg", fontsize=13, fontweight='bold')
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_aspect('equal')
ax.axis('off')

# Left: solute at GB
# Draw schematic atoms
def draw_atom_cluster(ax, cx, cy, label, highlight_idx=None, highlight_color='orange'):
    """Draw a small cluster of atoms around (cx, cy)."""
    positions = []
    for i in range(-2, 3):
        for j in range(-2, 3):
            x = cx + i * 0.45
            y = cy + j * 0.45
            positions.append((x, y))

    for idx, (x, y) in enumerate(positions):
        if highlight_idx is not None and idx == highlight_idx:
            ax.plot(x, y, 'o', color=highlight_color, markersize=12, zorder=5,
                    markeredgecolor='k', markeredgewidth=1.5)
        else:
            ax.plot(x, y, 'o', color='steelblue', markersize=8, zorder=3,
                    markeredgecolor='k', markeredgewidth=0.5, alpha=0.6)
    ax.text(cx, cy - 1.6, label, fontsize=10, ha='center', fontweight='bold')

# GB site scenario
draw_atom_cluster(ax, 2.5, 6.5, "Solute at GB site", highlight_idx=12, highlight_color='orange')
ax.text(2.5, 4.3, r"$E_{\rm GB}^{\rm solute}$", fontsize=14, ha='center', color='darkorange')

# Bulk site scenario
draw_atom_cluster(ax, 7.5, 6.5, "Solute at bulk site", highlight_idx=12, highlight_color='limegreen')
ax.text(7.5, 4.3, r"$E_{\rm bulk}^{\rm solute}$", fontsize=14, ha='center', color='green')

# Arrow and equation
ax.annotate('', xy=(6.0, 3.2), xytext=(4.0, 3.2),
            arrowprops=dict(arrowstyle='->', lw=2, color='black'))
ax.text(5.0, 3.6, "minus", fontsize=10, ha='center')

# Equation
ax.text(5.0, 2.0, r"$\Delta E_{\rm seg} = E_{\rm GB}^{\rm solute} - E_{\rm bulk}^{\rm solute}$",
        fontsize=15, ha='center', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', edgecolor='gray'))
ax.text(5.0, 0.8, r"$\Delta E < 0$: segregation favorable" + "\n"
        r"$\Delta E > 0$: segregation unfavorable",
        fontsize=10, ha='center', color='gray')

# Legend entries
ax.plot([], [], 'o', color='orange', markersize=10, markeredgecolor='k', label='Solute atom')
ax.plot([], [], 'o', color='steelblue', markersize=8, markeredgecolor='k', label='Solvent atom')
ax.legend(loc='upper center', fontsize=9, ncol=2)

# ═══════════════════════════════════════════
# Panel (d): MC swap procedure
# ═══════════════════════════════════════════
ax = axes[1, 1]
ax.set_title("(d) Monte Carlo Swap Simulation", fontsize=13, fontweight='bold')
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.set_aspect('equal')
ax.axis('off')

# Step 1: initial config (top)
y_top = 7.8
ax.text(2.0, y_top + 0.5, "Before swap", fontsize=11, ha='center', fontweight='bold')
# GB region
for i in range(7):
    x = 0.5 + i * 0.5
    color = 'orange' if i == 3 else 'steelblue'
    ax.plot(x, y_top, 'o', color=color, markersize=14, markeredgecolor='k', markeredgewidth=1)
ax.text(0.5, y_top - 0.6, "← GB →", fontsize=9, ha='left', color='red')
# Bulk region
for i in range(7):
    x = 5.5 + i * 0.5
    color = 'steelblue'
    ax.plot(x, y_top, 'o', color=color, markersize=14, markeredgecolor='k', markeredgewidth=1)
ax.text(5.5, y_top - 0.6, "← Bulk →", fontsize=9, ha='left', color='blue')

# Arrow down: swap
ax.annotate('', xy=(5.0, 6.2), xytext=(5.0, 7.0),
            arrowprops=dict(arrowstyle='->', lw=2, color='red'))
ax.text(5.5, 6.5, "Random swap\nsolute ↔ solvent", fontsize=9, ha='left', color='red')

# Step 2: after swap (middle)
y_mid = 5.5
ax.text(2.0, y_mid + 0.5, "After swap", fontsize=11, ha='center', fontweight='bold')
for i in range(7):
    x = 0.5 + i * 0.5
    color = 'steelblue'
    ax.plot(x, y_mid, 'o', color=color, markersize=14, markeredgecolor='k', markeredgewidth=1)
for i in range(7):
    x = 5.5 + i * 0.5
    color = 'orange' if i == 2 else 'steelblue'
    ax.plot(x, y_mid, 'o', color=color, markersize=14, markeredgecolor='k', markeredgewidth=1)

# Arrow down: relax + accept/reject
ax.annotate('', xy=(5.0, 3.8), xytext=(5.0, 4.7),
            arrowprops=dict(arrowstyle='->', lw=2, color='darkgreen'))
ax.text(5.5, 4.2, "Energy minimize\n→ compute ΔE", fontsize=9, ha='left', color='darkgreen')

# Metropolis criterion box
metro_text = (
    "Metropolis criterion:\n"
    r"if $\Delta E < 0$: accept" + "\n"
    r"if $\Delta E > 0$: accept with $p = e^{-\Delta E / k_B T}$"
)
ax.text(5.0, 2.5, metro_text, fontsize=11, ha='center',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='lightyellow', edgecolor='gray'))

# Bottom: repeat arrow
ax.annotate('', xy=(1.5, 1.0), xytext=(8.5, 1.0),
            arrowprops=dict(arrowstyle='<->', lw=2, color='purple'))
ax.text(5.0, 0.5, "Repeat ~10⁴–10⁶ steps → equilibrium segregation",
        fontsize=10, ha='center', color='purple', fontweight='bold')

plt.tight_layout()
plt.savefig("/scratch/cainiu/UMA/Computational/output/method_overview.png", dpi=150, bbox_inches='tight')
plt.close()
print("Figure saved to output/method_overview.png")
