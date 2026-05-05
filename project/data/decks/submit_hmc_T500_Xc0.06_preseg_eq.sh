#!/bin/bash
# Equilibration sweep job (2026-05-05): T=500 K, X_c=0.06.
# Re-runs from canonical preseg IC (X_GB(0)=0.3206) using v2 deck which
# writes a checkpoint every 5 ps. Replaces the 188-ps 16-rank TIMEOUT
# run (job 65259659) with a checkpoint-aware trajectory.
#
# 16-rank for parallel throughput (33 % public QOS quota). 24 h budget;
# ~10 ps/h on 16 ranks yields ~240 ps PROD before TIMEOUT, then
# resumable via hmc_AlMg_resume.lammps if 300 ps not reached.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.06_preseg_eq
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/hmc_AlMg/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/hmc_AlMg/%x-%j.err

set -euo pipefail

module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

REPO=/cluster/home/cainiu/Computational_modeling/project
RUN_DIR=/cluster/scratch/cainiu/hmc_AlMg
DECK=$REPO/data/decks/hmc_AlMg_v2.lammps
PRESEG=/cluster/scratch/cainiu/production_AlMg_200A/poly_AlMg_200A_preseg_Xc0.06.lmp

if [ ! -f "$PRESEG" ]; then
  echo "ERROR: pre-segregated data file not found: $PRESEG" >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$PRESEG" $(basename "$PRESEG")

T=500.0
XC=0.06
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
RESTART_PS=5
OUTSTUB="hmc_T${T%.*}_Xc0.06_preseg_eq"
FD_PRED=0.2611

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : equilibration sweep, 16-rank, v2 deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Initial Mg     : PRE-SEGREGATED (canonical preseg IC)"
echo "canon-FD pred  : $FD_PRED (ours)"
echo "Started        : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var annealed       $(basename "$PRESEG") \
  -var outstub        "$OUTSTUB" \
  -var T              "$T" \
  -var XC             "$XC" \
  -var EQUIL_PS       "$EQUIL_PS" \
  -var PROD_PS        "$PROD_PS" \
  -var RESTART_PS     "$RESTART_PS" \
  -var SWAPS_PER_CALL "$SWAPS_PER_CALL" \
  -var SKIP_PLACE     1

echo "=========================================================="
echo "Sim done at $(date); auto-post-processing..."
echo "=========================================================="

python3 $REPO/scripts/hmc_xgb_timeseries.py \
  --stub $RUN_DIR/$OUTSTUB \
  --gb-mask $REPO/data/snapshots/gb_mask_200A.npy \
  --xc $XC \
  --temp $T \
  --fd-pred $FD_PRED \
  --out-prefix $REPO/output/$OUTSTUB \
  || echo "POST-PROCESS FAILED but sim succeeded; rerun manually"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
