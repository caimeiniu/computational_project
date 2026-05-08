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
import warnings

# Suppress the "OVITO PyPI in conda env" warning that fires on every import
# when ovito was pip-installed inside a miniconda environment. Cosmetic only.
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI.*")

import numpy as np

try:
    from ovito.io import import_file
    from ovito.modifiers import (
        CommonNeighborAnalysisModifier,
        PythonScriptModifier,
    )
    from ovito.vis import Viewport, TachyonRenderer
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

    # ---- 3. Per-atom colour + transparency via PythonScriptModifier ----
    # We need this in the modifier chain (not on a one-shot data.compute())
    # because viewport.render_image() re-runs the pipeline; properties set
    # on a stale frame would be discarded.
    def colour_and_fade(frame, data):
        structure = np.asarray(data.particles["Structure Type"][...])
        n = data.particles.count
        fcc_mask = (structure == CommonNeighborAnalysisModifier.Type.FCC)
        n_gb = n - int(fcc_mask.sum())
        colors = np.empty((n, 3), dtype=np.float32)
        transp = np.empty(n, dtype=np.float32)
        colors[fcc_mask]  = (0.85, 0.85, 0.85)   # bulk grey
        colors[~fcc_mask] = (0.27, 0.46, 0.71)   # GB steel blue
        transp[fcc_mask]  = 0.82                  # bulk almost see-through
        transp[~fcc_mask] = 0.00                  # GB fully opaque
        data.particles_.create_property("Color", data=colors)
        data.particles_.create_property("Transparency", data=transp)
        # Print only on the first frame call for a clean log.
        if frame == 0:
            print(f"a-CNA: {n - n_gb} FCC bulk + {n_gb} non-FCC GB out of "
                  f"{n} atoms ({n_gb/n*100:.2f} % GB)")

    pipeline.modifiers.append(
        PythonScriptModifier(function=colour_and_fade)
    )

    # Render-time particle size: slightly smaller than default so bulk
    # spheres don't fully occlude the GB network even at 0.82 transparency.
    pipeline.source.data.particles_.vis.radius = 0.95

    pipeline.add_to_scene()

    # ---- 4. Camera + render --------------------------------------------
    # Box is 200 Å on a side, centred on (100, 100, 100). We use a
    # corner-ish 3/4 view, then zoom_all() auto-fits the cube into the
    # viewport. The exact camera_pos / camera_dir act as direction hints;
    # zoom_all overrides distance so the cube fills the frame.
    vp = Viewport()
    vp.type = Viewport.Type.Perspective
    vp.camera_pos = (450.0, -200.0, 350.0)
    vp.camera_dir = (-0.62, 0.55, -0.56)
    vp.fov = float(np.radians(35.0))
    vp.zoom_all()

    # TachyonRenderer is CPU raytraced — slower than OpenGL (~30 s vs ~5 s)
    # but works headlessly (no X / OpenGL context required), which matters
    # on Euler login / compute nodes.
    vp.render_image(
        filename=str(OUTPUT_PNG),
        size=(1600, 1200),
        background=(1.0, 1.0, 1.0),
        renderer=TachyonRenderer(),
    )
    print(f"wrote {OUTPUT_PNG}")
    print("Upload this PNG to Overleaf next to main.tex; "
          "methods_draft.tex's \\includegraphics will resolve.")


if __name__ == "__main__":
    main()
