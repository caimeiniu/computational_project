#!/usr/bin/env python3
"""Render panel (f): OVITO snapshot of a converged HMC final config showing
GB segregation (Mg clustered at GB planes inside a polycrystalline Al
matrix).

Usage:
    python3 scripts/render_panel_f.py \
        --lmp project/data/snapshots/hmc_T500_Xc0.20_preseg_final.lmp \
        --out project/output/panel_f_xc0.20_preseg.png

Requires the `ovito` Python package (pip install ovito --user).

Renders headlessly via TachyonRenderer (software ray-tracer; OpenGL
requires an X server, not available on Euler compute / login nodes).

Default visual scheme:
  - Al (type 1): small radius, light gray — provides matrix context
    without obscuring GB structure
  - Mg (type 2): larger radius, orange-red — segregated atoms at GBs
    pop visually as 2D sheets inside the 3D box
  - White background, perspective camera, auto-fit zoom
"""
from __future__ import annotations
import argparse
from pathlib import Path
import warnings

# Silence the PyPI-vs-Anaconda warning that the OVITO pip package emits at
# import time on managed Python installs.
warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")

from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lmp", type=Path, required=True,
                    help="LAMMPS data file (atom_style atomic, types 1=Al "
                         "2=Mg).")
    ap.add_argument("--out", type=Path, required=True,
                    help="Output PNG path.")
    ap.add_argument("--size", type=int, nargs=2, default=(1600, 1200),
                    metavar=("W", "H"))
    ap.add_argument("--al-radius", type=float, default=0.6,
                    help="Sphere radius (Å) for Al atoms; smaller = less "
                         "matrix occlusion of the GB pattern. (default 0.6)")
    ap.add_argument("--mg-radius", type=float, default=1.3,
                    help="Sphere radius (Å) for Mg atoms; larger = more "
                         "GB structure visibility. (default 1.3)")
    ap.add_argument("--hide-al", action="store_true",
                    help="Render only Mg atoms; emphasizes the GB-Mg "
                         "structure but loses the polycrystal matrix.")
    ap.add_argument("--view", choices=["perspective", "front", "side", "top"],
                    default="perspective",
                    help="Camera orientation. (default perspective)")
    ap.add_argument("--bg", type=float, nargs=3, default=(1.0, 1.0, 1.0),
                    metavar=("R", "G", "B"),
                    help="Background colour. (default white)")
    args = ap.parse_args()

    pipeline = import_file(str(args.lmp))

    # Configure type appearance at the SOURCE level (not on a computed
    # copy) — modifications to `pipeline.source.data` propagate through
    # the pipeline and are honoured by the renderer. OVITO uses CoW data
    # ownership; the trailing-`_` accessors (e.g. `particles_`,
    # `particle_types_`, `type_by_id_`) return a mutable handle. Setting
    # colors on the result of `pipeline.compute()` would have no effect
    # on the rendered scene.
    src_types = pipeline.source.data.particles_.particle_types_

    def _set_type(type_id: int, name: str, color, radius: float) -> None:
        t = src_types.type_by_id_(type_id)   # mutable type handle
        t.name = name
        t.color = color
        t.radius = radius

    _set_type(1, "Al", (0.70, 0.70, 0.74),
              0.0 if args.hide_al else args.al_radius)
    _set_type(2, "Mg", (0.95, 0.42, 0.20), args.mg_radius)

    pipeline.add_to_scene()

    vp_type = {
        "perspective": Viewport.Type.Perspective,
        "front":       Viewport.Type.Front,
        "side":        Viewport.Type.Right,
        "top":         Viewport.Type.Top,
    }[args.view]
    vp = Viewport(type=vp_type)
    vp.zoom_all()

    out_path = args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    vp.render_image(
        filename=str(out_path),
        size=tuple(args.size),
        background=tuple(args.bg),
        renderer=TachyonRenderer(),
    )
    print(f"wrote {out_path}  ({args.size[0]}×{args.size[1]} px)")


if __name__ == "__main__":
    main()
