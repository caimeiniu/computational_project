#!/bin/bash
# Submit the Wagih-style Al polycrystal anneal to SLURM on Euler — production
# 20³ nm³ / 16 grains / ~476k atoms box.
#
# Usage (from anywhere):
#   sbatch /cluster/home/cainiu/Computational_modeling/project/data/decks/submit_anneal_200A.sh
#
# Timing basis (job 64663404): 1.74 ns/day on 32 cores (49.7 s / 1000 steps).
# Full anneal is 5 ps ramp + 250 ps hold + 124 ps cool = 379 ps ≈ 5.2 h MD;
# requesting 8 h to absorb CG + startup + any long-tail stragglers.

#SBATCH --job-name=anneal_AlMg_200A
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/production_AlMg_200A/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/production_AlMg_200A/%x-%j.err

set -euo pipefail

# ----- environment -----
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

# ----- job parameters -----
RUN_DIR=/cluster/scratch/cainiu/production_AlMg_200A
DECK=/cluster/home/cainiu/Computational_modeling/project/data/decks/anneal_AlMg.lammps
DATAFILE=poly_Al_200A_16g.lmp
OUTSTUB=poly_Al_200A_16g
T_HOLD=373.0        # K  (0.4 T_melt for Al ≈ 933 K) — matches prototype
T_HOLD_PS=250       # ps
COOL_RATE=3.0       # K/ps
VELOCITY_SEED=20240101  # same as prototype, so only box size differs

# ----- bookkeeping -----
cd "$RUN_DIR"
echo "=========================================================="
echo "Job ID       : $SLURM_JOB_ID"
echo "Job name     : $SLURM_JOB_NAME"
echo "Node list    : $SLURM_JOB_NODELIST"
echo "MPI tasks    : $SLURM_NTASKS"
echo "Run dir      : $(pwd)"
echo "Deck         : $DECK"
echo "Datafile     : $DATAFILE (~476k atoms, 20³ nm³, 16 grains)"
echo "Potential    : Al-Mg.eam.fs (Mendelev 2009)"
echo "T_hold       : $T_HOLD K (hold $T_HOLD_PS ps, cool $COOL_RATE K/ps)"
echo "Started      : $(date)"
echo "=========================================================="

# ----- run -----
srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var datafile      "$DATAFILE" \
  -var outstub       "$OUTSTUB" \
  -var T_HOLD        "$T_HOLD" \
  -var T_HOLD_PS     "$T_HOLD_PS" \
  -var COOL_RATE     "$COOL_RATE" \
  -var VELOCITY_SEED "$VELOCITY_SEED"

echo "=========================================================="
echo "Finished     : $(date)"
echo "=========================================================="
