#!/bin/bash
# Random-IC HMC run (2026-05-13): T=800 K, X_c=0.10.
#
# Bracket-from-below companion to submit_hmc_T800_Xc0.10_fdseed.sh
# (the bracket-from-above half). Together: two-IC convergence test
# at T=800 K. Advisor's suggestion to validate HMC ergodicity at a T
# where swap acceptance is high enough to mix in 24 h.
#
# Initial condition: poly_AlMg_200A_random_Xc0.1.lmp from
# pre_segregate.py with --xgb-init 0.10 (= X_c, so each atom GB or
# bulk has uniform prob X_c of being Mg). X_GB(0) = 0.1000 = X_bulk(0).
# This is BELOW the FD prediction X_GB^FD(800 K, 0.10) = 0.2734, so
# HMC will ASCEND.
#
# Convergence outcomes:
#   ascend → X_GB^∞ ≈ 0.2734 AND fdseed companion stays at 0.2734
#     → HMC ergodic at T=800 K, FD recovered, no breakdown at high T
#   ascend → X_GB^∞ < 0.2734 AND fdseed companion descends to same
#     → HMC ergodic, breakdown is real and persists at high T
#   ascend → X_GB^∞ ≠ fdseed end value
#     → HMC not yet ergodic at T=800 K — schedule resume

#SBATCH --job-name=hmc_AlMg_T800_Xc0.10_random
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
SNAPSHOT=$RUN_DIR/poly_AlMg_200A_random_Xc0.1.lmp

if [ ! -f "$SNAPSHOT" ]; then
  echo "ERROR: snapshot not found: $SNAPSHOT" >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

T=800.0
XC=0.10
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
RESTART_PS=5
OUTSTUB="hmc_T800_Xc0.10_random"
FD_PRED=0.2734

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : Random/uniform IC, 16-rank, v2 deck"
echo "T              : $T K  (0.86 T_melt(Al=933 K))"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Initial X_GB   : 0.1000 (= X_c, uniform Mg)"
echo "canon-FD pred  : $FD_PRED (ours; HMC ascends from below)"
echo "Companion      : submit_hmc_T800_Xc0.10_fdseed.sh (X_GB(0)=0.2734)"
echo "Started        : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var annealed       $(basename "$SNAPSHOT") \
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
