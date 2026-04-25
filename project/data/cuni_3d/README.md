# 3D Cu-Ni ΔE Workflow

This folder is for the 3D Cu-Ni prototype: generate a true 3D Voronoi FCC Cu
polycrystal, sample 500 grain-boundary sites, calculate single-solute Cu->Ni
segregation energies, and compare the spectrum to Wagih/Larsen/Schuh Ni(Cu).

## Local setup

The scripts only need Python, NumPy, SciPy, and Matplotlib. The actual energy
relaxations require LAMMPS with `eam/alloy`.

## Steps

From the repo root:

### One-command Euler run

On Euler, from the repo root:

```bash
sbatch -n 16 --time=24:00:00 --output=cuni3d_%j.out scripts/euler_cuni3d_pipeline.sbatch
```

Check progress:

```bash
squeue -u $USER
tail -f cuni3d_<jobid>.out
```

The batch script runs in `/cluster/scratch/$USER/cuni3d_<jobid>` and copies compact
results back to `data/cuni_3d/euler_results_<jobid>/`.

### Manual/local steps

```bash
python3 scripts/generate_cuni_3d_polycrystal.py --box-nm 10 --grains 8
```

Anneal and minimize the raw Voronoi structure on Euler or another machine with LAMMPS:

```bash
cd data/cuni_3d/lammps
lmp -in in_relax_cuni_3d.lammps
cd ../../..
```

Then sample 500 sites from the relaxed structure and write the ΔE jobs:

```bash
python3 scripts/select_cuni_deltae_sites.py --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps --samples 500
python3 scripts/write_cuni_deltae_lammps_jobs.py --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps
```

Run the generated LAMMPS jobs:

```bash
cd data/cuni_3d/lammps/deltae_jobs
./run_all_deltae.sh
```

Collect and compare:

```bash
python3 scripts/analyze_cuni_deltae.py
```

Outputs:

- `cuni_3d_polycrystal.lammps`: pure Cu 3D Voronoi structure.
- `lammps/cuni_3d_relaxed.lammps`: annealed/minimized pure Cu 3D structure.
- `cuni_3d_deltae_sites.csv`: one bulk reference plus 500 sampled GB sites.
- `lammps/deltae_jobs/`: one minimization input per sampled site.
- `cuni_3d_deltae_results.csv`: relaxed ΔE values in kJ/mol.
- `cuni_3d_deltae_vs_wagih.png`: histogram compared with the approximate Wagih Ni(Cu)
  spectrum from the project notes.

## Notes

Wagih's public Zenodo dataset is large, about 4 GB, so this workflow uses the
project-note Ni(Cu) reference values by default: μ ≈ -2 kJ/mol, σ ≈ 8 kJ/mol,
range ≈ [-30, 30] kJ/mol. If we later download the full database, replace the
approximate curve with the site-resolved Ni(Cu) distribution.
