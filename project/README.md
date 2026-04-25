# GB Segregation — HMC probe of the dilute-limit assumption

Grain boundary (GB) solute segregation study using Hybrid Monte Carlo to test where
Wagih's independent-site (non-interacting solute) framework breaks down as a function
of solute concentration and temperature.

## Reference

Wagih, Larsen & Schuh, "Learning grain boundary segregation energy spectra in polycrystals",
*Nature Communications* 11:6376 (2020). https://doi.org/10.1038/s41467-020-20083-6

## Project structure

```
Computational/
├── scripts/     # Python scripts (polycrystal generation, HMC, analysis)
├── data/        # Input data (potential files, structure files)
├── output/      # Simulation results (gitignored, large files)
├── docs/        # Notes, references
├── CHANGELOG.md # Development log
└── README.md    # This file
```

## Environment

```bash
conda activate gb-seg
```

## Workflow

1. Generate polycrystalline structure (Voronoi tessellation + thermal annealing)
2. Identify GB sites via structural analysis (a-CNA)
3. Compute per-site segregation energies ΔE_i (single solute, 0 K, CG relax) → spectrum
4. Hybrid Monte Carlo on a `(T, X_c)` grid to sample equilibrium occupation
5. Compare `X_GB(T, X_c)` and site-resolved `P_i` from HMC vs Fermi-Dirac prediction
   built from the 0 K ΔE spectrum → locate the dilute-limit breakdown concentration
6. Solute-solute `g(r)` at GB + local-density vs `P_i` correlation to diagnose the
   physical origin of the breakdown

## 3D Cu-Ni prototype

The 2D columnar Cu-Ni advisor example is preserved under `data/examples/nc_swap_CuNi/`.
The current 3D Cu-Ni prototype lives under `data/cuni_3d/` and is driven by:

```bash
python3 scripts/generate_cuni_3d_polycrystal.py --box-nm 10 --grains 8
cd data/cuni_3d/lammps && lmp -in in_relax_cuni_3d.lammps && cd ../../..
python3 scripts/select_cuni_deltae_sites.py --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps --samples 500
python3 scripts/write_cuni_deltae_lammps_jobs.py --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps
```

Run the generated LAMMPS jobs, then analyze:

```bash
python3 scripts/analyze_cuni_deltae.py
```

The analyzer reports ΔE in kJ/mol and overlays the sampled 3D Cu-Ni spectrum with the
approximate Wagih/Larsen/Schuh Ni(Cu) reference noted in `docs/paper_notes.md`.
