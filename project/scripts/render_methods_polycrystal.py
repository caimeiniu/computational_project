#!/usr/bin/env python3
"""Render the 16-grain Al polycrystal for the Methods figure (right panel
of Fig.~\\ref{fig:method}) in the GB-segregation report.

Pipeline:
    1. Load the annealed pure-Al polycrystal LAMMPS data file
    2. Run adaptive common-neighbour analysis (a-CNA) → classify each
       atom as FCC (bulk) or non-FCC (GB)
    3. Recolour: bulk atoms in light grey, GB atoms in warm red
    4. Render perspective view to PNG

INPUT:
    /cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g_annealed.lmp

OUTPUT:
    report/figures/polycrystal_geometry.png
        (1600 × 1200 px, white background; sized for ~6 cm width in a
        double-column figure)

USAGE:
    OVITO Pro (recommended, supports vector export):
        python3 scripts/render_methods_polycrystal.py
    OVITO Basic (free, raster only):
        # same command — OVITO Basic also exposes the Python API for
        # the modifiers used here.
    On Euler:
        module load ovito  # if available; otherwise pip install ovito
        python3 scripts/render_methods_polycrystal.py

CAMERA TUNING:
    The default camera looks at the 200 Å cube from a corner-ish elevation
    so the GB network is visible through the cube interior. If the result
    is unsatisfying, edit `vp.camera_pos` and `vp.camera_dir` near the end
    of `main()` and re-run. The camera_pos = (350, -300, 270) is in Å and
    points the camera roughly towards the box centre at (100, 100, 100).
"""
from __future__ import annotations
from pathlib import Path
import sys

import numpy as np

try:
    from ovito.io import import_file
    from ovito.modifiers import CommonNeighborAnalysisModifier
    from ovito.vis import Viewport, OpenGLRenderer
except ImportError as e:
    sys.exit(
        f"OVITO Python module not available: {e}\n"
        "Install with one of:\n"
        "  pip install ovito\n"
        "  conda install -c https://conda.ovito.org ovito-pro\n"
        "Or, on Euler, try `module spider ovito` and load the matching module.\n"
    )


REPO = Path("/cluster/home/cainiu/Computational_modeling/project")
INPUT_LMP = Path(
    "/cluster/scratch/cainiu/production_AlMg_200A/"
    "poly_Al_200A_16g_annealed.lmp"
)
OUTPUT_PNG = REPO / "report" / "figures" / "polycrystal_geometry.png"


def main() -> None:
    if not INPUT_LMP.exists():
        sys.exit(f"input not found: {INPUT_LMP}")
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)

    # ---- 1. Load the annealed polycrystal -----------------------------
    pipeline = import_file(str(INPUT_LMP))

    # ---- 2. Adaptive CNA: tag atoms by local crystal structure --------
    # FCC = bulk Al; everything else (HCP / BCC / ICO / OTHER) = GB.
    pipeline.modifiers.append(
        CommonNeighborAnalysisModifier(
            mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff,
        )
    )

    # Compute one frame so we can read the Structure Type property and
    # set per-atom colours manually (more control than ColorCodingModifier).
    data = pipeline.compute()
    structure = np.asarray(data.particles["Structure Type"][...])
    n_atoms = data.particles.count
    fcc_type = CommonNeighborAnalysisModifier.Type.FCC
    n_fcc = int((structure == fcc_type).sum())
    n_gb = n_atoms - n_fcc
    print(f"a-CNA: {n_fcc} FCC bulk + {n_gb} non-FCC GB out of "
          f"{n_atoms} atoms ({n_gb/n_atoms*100:.2f} % GB)")

    # ---- 3. Recolour atoms ---------------------------------------------
    colors = np.empty((n_atoms, 3), dtype=np.float32)
    fcc_mask = (structure == fcc_type)
    colors[fcc_mask]  = (0.82, 0.82, 0.82)   # bulk Al   — light grey
    colors[~fcc_mask] = (0.85, 0.30, 0.20)   # GB atoms  — warm red

    # Inject as the per-particle Color property; OVITO uses this to
    # override the default per-type colour map at render time.
    data.particles_.create_property("Color", data=colors)

    # Render-time visualisation tweaks: slightly smaller particles so the
    # GB network is not obscured by overlapping bulk spheres.
    pipeline.source.data.particles_.vis.radius = 1.05

    pipeline.add_to_scene()

    # ---- 4. Camera + render --------------------------------------------
    # Box is 200 Å on a side, centred on (100, 100, 100). Camera at
    # (350, -300, 270) Å, pointing back towards the box → roughly
    # corner-ish elevation, ~35° field of view.
    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    vp.camera_pos = (350.0, -300.0, 270.0)
    vp.camera_dir = (-0.65, 0.55, -0.52)
    vp.fov = float(np.radians(35.0))

    vp.render_image(
        filename=str(OUTPUT_PNG),
        size=(1600, 1200),
        background=(1.0, 1.0, 1.0),
        renderer=OpenGLRenderer(),
    )
    print(f"wrote {OUTPUT_PNG}")
    print("Upload this PNG to Overleaf next to main.tex; "
          "methods_draft.tex's \\includegraphics will resolve.")


if __name__ == "__main__":
    main()
