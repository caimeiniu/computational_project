# GB Segregation — HMC probe of the dilute-limit assumption

Grain boundary (GB) solute segregation study using Hybrid Monte Carlo to
test where Wagih's independent-site (non-interacting solute) framework
breaks down as a function of solute concentration and temperature.

## Reference

Wagih, Larsen & Schuh, "Learning grain boundary segregation energy
spectra in polycrystals", *Nature Communications* 11:6376 (2020).
https://doi.org/10.1038/s41467-020-20083-6

## Project structure

```
project/
├── scripts/      Python: generation, GB ID, ΔE sampling, fits, FD predictor
├── data/
│   ├── decks/    LAMMPS input decks + SLURM submit shells
│   └── potentials/  EAM/fs files
├── output/       fits + figures (PNG/JSON kept; large dumps gitignored)
├── docs/         notes, paper PDF
├── CHANGELOG.md  development log (newest first)
└── README.md     this file
```

## Environment

Tested on ETHZ Euler (SLURM); replace `module load …` and the
`/cluster/scratch/$USER` paths with your cluster's equivalents on other
systems.

```bash
# LAMMPS 20240829 with MANYBODY + MC + VORONOI + KIM packages.
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4
# Python ≥ 3.7 (uses `from __future__ import annotations`).
# Required Python packages: numpy, scipy, matplotlib (any reasonable versions).
# Install however you prefer: conda env, venv, pip --user, ...
```

Quick environment sanity check (≈30 s, no SLURM allocation needed):

```bash
python scripts/fermi_dirac_predict.py --self-test --out-png /tmp/x.png --out-json /tmp/x.json
# Should print "self-test OK"; produces no real output, only the analytic-limit checks.
```

Design rationale for every choice (CG tolerances, anneal protocol, n=500
sample size, FD-before-HMC ordering, etc.) is in `CHANGELOG.md` (newest
entries first).

## Pipeline overview

```
generate_polycrystal.py  →  anneal_AlMg.lammps  →  gb_identify.py  →  sample_delta_e.py
       Voronoi              Wagih-style anneal       a-CNA mask         per-site ΔE_i [NPZ]
                                                                              │
                                                                              ▼
                                                                  fit_delta_e_spectrum.py
                                                                  fermi_dirac_predict.py
                                                                  (downstream:  HMC scan)
```

## Validated scripts (ready for any FCC/BCC/HCP binary alloy)

The following are validated end-to-end on Al(Mg) Mendelev 2009 against
the Wagih 2020 dataset (KS p=0.892 vs Wagih's 82,646-site Zenodo pool;
6-statistic bootstrap CI all inside; paired per-site PE residual −0.035 meV
on Wagih's own structure). Cu(Ni) requires only changing the potential
file and a few CLI flags — see *Cu(Ni) quickstart* below.

| Script / deck | Role | Generic over alloy? |
|---------------|------|---------------------|
| `scripts/generate_polycrystal.py` | Voronoi 3D polycrystal, FCC/BCC/HCP | yes — `--structure --lattice-a --grains --box` |
| `scripts/gb_identify.py` | a-CNA → per-atom GB bool mask | yes — `--parent {fcc,bcc,hcp} --lattice-a` |
| `scripts/sample_delta_e.py` | per-site ΔE_seg via single-substitution + tight CG | yes — `--elements "Cu Ni" --masses "63.546 58.6934"` |
| `scripts/fit_delta_e_spectrum.py` | skew-normal fit of ΔE spectrum | yes — pure NPZ → fit |
| `scripts/fermi_dirac_predict.py` | X_GB^FD(T, X_c) curves from ΔE NPZ | yes — pure NPZ → curves |
| `data/decks/anneal_AlMg.lammps` | Wagih-style anneal (CG → ramp → NPT hold → cool → CG) | yes — `-var EL1 -var EL2 -var MASS1 -var MASS2 -var T_HOLD` |

## Validation-only (Al(Mg)-specific, not relevant to Cu(Ni) unless they have a Wagih-style reference dataset)

| Script | Purpose | Why Al(Mg)-only |
|--------|---------|-----------------|
| `scripts/compare_vs_wagih.py` | overlay our spectrum vs Wagih Zenodo pool, KS test | needs Wagih's `seg_energies_Al_Mg.txt` + `bulk_solute_Al_Mg.dat` |
| `scripts/bootstrap_vs_wagih.py` | 6-statistic bootstrap CI of n=500 vs 82,646 pool | same |
| `scripts/paired_pipeline_residual.py` | per-site PE residual on Wagih's own dump | needs Wagih's structure + per-site E_GB |
| `scripts/wagih_dump_to_data.py` | convert Wagih's dump to LAMMPS data file for re-substitution | one-time helper |

## Not yet validated (do not depend on)

| Item | Status |
|------|--------|
| `data/decks/hmc_AlMg.lammps` + `submit_hmc_dryrun.sh` | dry-run pending; merge after swap acceptance rate (target 5–30 %), PE plateau, and X_GB(t) post-pipeline are verified |
| HMC post-processing (`hmc_xgb_timeseries.py`) | not written |

## Cu(Ni) quickstart

### Conventions

- **Atom types**: type 1 = solvent (the polycrystal host, here Cu);
  type 2 = solute (substituted in single-site, here Ni). All scripts
  assume this; `sample_delta_e.py --solute-type` defaults to 2.
- **Element order in EAM**: `--elements "Cu Ni"` and `pair_coeff`
  ordering must match the EAM/fs file header line 4. Verify before
  running anything:

  ```bash
  head -5 your-Cu-Ni.eam.fs
  # Look for the line that starts with `2  <El1>  <El2>` — that's the
  # order. If it says `2  Ni  Cu`, use `--elements "Ni Cu"` everywhere
  # (and remember type 1 is then Ni in your data file).
  ```

- **Cu-Ni EAM/fs source**: NIST Interatomic Potentials Repository
  (https://www.ctcms.nist.gov/potentials/) hosts several Cu-Ni alloy
  potentials (Foiles 1986, Onat & Durukanoglu 2014, …). Pick one and
  drop into `data/potentials/`.

### Workflow

```bash
# ----- paths (edit for your cluster / clone location) -----
ALLOY=CuNi
SCRATCH=/cluster/scratch/$USER/$ALLOY                                 # Euler example
PROJECT=/cluster/home/$USER/Computational_modeling/project            # your clone
POT=$PROJECT/data/potentials/Cu-Ni.eam.fs                             # provide yourself
mkdir -p "$SCRATCH" && cd "$SCRATCH"

# ----- Step 1: generate 20³ nm³, 16-grain FCC polycrystal -----
# Cu lattice 3.615 Å (room T). Cheap (~1 min on login node).
python "$PROJECT/scripts/generate_polycrystal.py" \
  --structure fcc --box 200 --grains 16 --lattice-a 3.615 \
  --structure-seed 1 --types 2 \
  --out poly_Cu_200A_16g.lmp

# ----- Step 2: anneal (Wagih-style) — submit via SLURM -----
# Copy the AlMg submit shell as a starting template, then edit the variables
# (RUN_DIR, DATAFILE, OUTSTUB, POTENTIAL, and add EL1/EL2/MASS1/MASS2/T_HOLD).
# T_HOLD = 0.4·T_melt(Cu=1357 K) ≈ 540 K.
cp "$PROJECT/data/decks/submit_anneal_200A.sh" submit_anneal_Cu_200A.sh
# … edit the file:
#   RUN_DIR=$SCRATCH
#   DATAFILE=poly_Cu_200A_16g.lmp
#   OUTSTUB=poly_Cu_200A_16g
#   POTENTIAL=$POT
#   T_HOLD=540.0
# … and add to the `srun -n ... lmp ...` invocation:
#   -var EL1 Cu -var EL2 Ni -var MASS1 63.546 -var MASS2 58.6934
sbatch submit_anneal_Cu_200A.sh   # ≈ 5–6 h on 32 cores, 8 h budget

# ----- Step 3: GB identification (login node, ~1 min) -----
python "$PROJECT/scripts/gb_identify.py" \
  poly_Cu_200A_16g_annealed.lmp \
  --parent fcc --lattice-a 3.615 \
  --out-mask gb_mask_Cu_200A.npy \
  --out-report gb_info_Cu_200A.json

# ----- Step 4: per-site ΔE_seg (n=500 GB sites, tight CG) — submit via SLURM -----
# DO NOT run sample_delta_e.py directly on the login node — it spawns
# `mpirun lmp` per site and will be killed. Use sbatch.
# Copy the AlMg submit shell as template:
cp "$PROJECT/data/decks/submit_delta_e_200A.sh" submit_delta_e_Cu_200A.sh
# … edit the file:
#   RUN_DIR=$SCRATCH
#   ANNEALED=$RUN_DIR/poly_Cu_200A_16g_annealed.lmp
#   GB_MASK=$RUN_DIR/gb_mask_Cu_200A.npy
#   POTENTIAL=$POT
# … and append to the `python "$DRIVER" \` block:
#   --elements "Cu Ni" --masses "63.546 58.6934"
# Optionally drop the conda activate line if you're not using a conda env.
sbatch submit_delta_e_Cu_200A.sh   # ≈ 45 min on 32 cores, default 8 h budget

# ----- Step 5: skew-normal fit + Fermi-Dirac predictions (login node) -----
python "$PROJECT/scripts/fit_delta_e_spectrum.py" \
  --npz delta_e_results_Cu_200A.npz \
  --out-png delta_e_spectrum_Cu_200A.png \
  --out-json delta_e_fit_Cu_200A.json

python "$PROJECT/scripts/fermi_dirac_predict.py" \
  --ours-npz delta_e_results_Cu_200A.npz \
  --T-list 500,700,900,1100 \
  --out-png fd_curves_Cu_200A.png \
  --out-json fd_curves_Cu_200A.json
```

### Notes for Cu(Ni)

- **Lattice parameters at room T**: Cu **3.615 Å**, Ni **3.524 Å**.
  Use the solvent (host = Cu) lattice for the polycrystal: **3.615 Å**.
- **T_hold ranges** (0.3–0.5·T_melt(Cu)): **407–679 K**. We use 0.4·T_melt
  on Al(Mg); 540 K is the analogous choice for Cu.
- **a-CNA cutoff**: default `0.854 a` for FCC is correct for both Cu and
  Ni; no manual `--cutoff` needed.
- The `compare_vs_wagih.py` / `bootstrap_vs_wagih.py` /
  `paired_pipeline_residual.py` scripts apply only if you have a
  reference per-site Cu(Ni) dataset to validate against. For novel
  Cu(Ni) predictions, the pipeline ends at `fermi_dirac_predict.py`.
- **Expected timing on 32 cores** (Euler):
  | step | wall time |
  |------|-----------|
  | generate_polycrystal | ~1 min, login node |
  | anneal (Wagih protocol) | ~5–6 h on 32 cores |
  | gb_identify | ~1 min, login node |
  | sample_delta_e (n=500, tight CG) | ~45 min on 32 cores |
  | fit + Fermi-Dirac | seconds, login node |

## What this project is doing on Al(Mg) (dilute-limit study)

The Cu(Ni) quickstart above ends at **step 5** (the Fermi-Dirac
prediction). On Al(Mg) we are continuing past that point to test where
the dilute-limit assumption breaks down:

6. Hybrid Monte Carlo on a `(T, X_c)` grid sampling the equilibrium
   GB occupation directly (LAMMPS `fix atom/swap`).
7. Compare `X_GB(T, X_c)` from HMC vs the Fermi-Dirac prediction →
   locate the breakdown concentration where independent-site assumption fails.
8. Solute-solute g(r) at GB + local-density vs `P_i` correlation to
   diagnose the physical origin of the breakdown.

Steps 6–8 will be merged into this repo once validated on Al(Mg) (HMC
deck currently in dry-run; see *Not yet validated* table above).
