#!/bin/bash
# Submit the Wagih-style Al polycrystal anneal to SLURM on Euler.
#
# Usage (from anywhere):
#   sbatch /cluster/home/cainiu/Computational_modeling/project/data/decks/submit_anneal.sh
#
# The script cd's into the scratch run directory, so slurm-*.out lands there.
# Override parameters via sbatch --export=ALL,KEY=VALUE or by editing below.

#SBATCH --job-name=anneal_AlMg_100A
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/prototype_AlMg_100A/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/prototype_AlMg_100A/%x-%j.err

set -euo pipefail

# ----- environment -----
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

# ----- job parameters (edit here to retune) -----
RUN_DIR=/cluster/scratch/cainiu/prototype_AlMg_100A
DECK=/cluster/home/cainiu/Computational_modeling/project/data/decks/anneal_AlMg.lammps
DATAFILE=poly_Al_100A_8g.lmp
OUTSTUB=poly_Al_100A_8g
T_HOLD=373.0        # K  (0.4 T_melt for Al ≈ 933 K)
T_HOLD_PS=250       # ps
COOL_RATE=3.0       # K/ps
VELOCITY_SEED=20240101

# ----- bookkeeping -----
cd "$RUN_DIR"
echo "=========================================================="
echo "Job ID       : $SLURM_JOB_ID"
echo "Job name     : $SLURM_JOB_NAME"
echo "Node list    : $SLURM_JOB_NODELIST"
echo "MPI tasks    : $SLURM_NTASKS"
echo "Run dir      : $(pwd)"
echo "Deck         : $DECK"
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
