#!/bin/bash
# Submit the Wagih-style Pt polycrystal anneal to SLURM on Euler — prototype
# 10³ nm³ / 8 grains / ~64k atoms box (mirrors prototype_AlMg_100A).
#
# Usage (from anywhere):
#   sbatch "$HOME/computational_project/project/PtAu/data/decks/submit_anneal_PtAu_100A.sh"
#
# Timing estimate: ~1 h on 16 cores for 549 ps MD at ~64k atoms (NPT cool to
# 0 K is longer than Al-Mg because T_hold is higher → more cool ps); 6 h budget.

#SBATCH --job-name=anneal_PtAu_100A
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/%u/prototype_PtAu_100A/%x-%j.out
#SBATCH --error=/cluster/scratch/%u/prototype_PtAu_100A/%x-%j.err

set -euo pipefail

# ----- environment -----
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

# ----- job parameters -----
PROJECT="${PROJECT:-$HOME/computational_project/project}"
RUN_DIR="${RUN_DIR:-/cluster/scratch/$USER/prototype_PtAu_100A}"
DECK="$PROJECT/PtAu/data/decks/anneal_PtAu.lammps"
POTENTIAL="$PROJECT/PtAu/data/potentials/PtAu.eam.alloy"
DATAFILE=poly_Pt_100A_8g.lmp
OUTSTUB=poly_Pt_100A_8g
T_HOLD=816.0        # K  (0.4 T_melt for Pt ≈ 2041 K) — middle Wagih range
T_HOLD_PS=250       # ps
COOL_RATE=3.0       # K/ps  -> 816/3 = 272 ps cool segment
VELOCITY_SEED=20240101  # same as Al-Mg prototype so only alloy/lattice/structure vary

# ----- bookkeeping -----
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
echo "=========================================================="
echo "Job ID       : $SLURM_JOB_ID"
echo "Job name     : $SLURM_JOB_NAME"
echo "Node list    : $SLURM_JOB_NODELIST"
echo "MPI tasks    : $SLURM_NTASKS"
echo "Run dir      : $(pwd)"
echo "Deck         : $DECK"
echo "Datafile     : $DATAFILE (~64k atoms, 10³ nm³, 8 grains, Pt host)"
echo "Potential    : PtAu.eam.alloy (O'Brien 2017)"
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
  -var VELOCITY_SEED "$VELOCITY_SEED" \
  -var POTFILE       "$POTENTIAL" \
  -var EL1           Pt \
  -var EL2           Au \
  -var MASS1         195.0900 \
  -var MASS2         196.9665

echo "=========================================================="
echo "Finished     : $(date)"
echo "=========================================================="
