# GB Segregation with UMA MLIP

Grain boundary (GB) solute segregation study using Hybrid Monte Carlo and UMA universal machine learning interatomic potential.

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
2. Identify GB sites via structural analysis (CNA/PTM)
3. Compute segregation energies (substitutional solute at each GB site)
4. Hybrid Monte Carlo to sample equilibrium segregation state
5. Compare: classical potential vs UMA MLIP
