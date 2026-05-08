#!/bin/bash
# FD-seeded HMC run (2026-05-08 evening): T=500 K, X_c=0.075.
# Headline-curve methodological cleanup batch (CHANGELOG 2026-05-08 evening).
#
# Initial condition: poly_AlMg_200A_fdseed_T500K_Xc0.075.lmp from
# build_fdseed_inits.py, X_GB(0) = X_GB^FD = 0.3007 at random GB sites
# (N_Mg(GB) = 26773) and the remaining Mg at random bulk sites
# (N_Mg(bulk) = 8906). Mass conservation: 26773 + 8906 = 35679 = X_c * N_total.
#
# Purpose: replace the existing eq_cont point at X_c=0.075 (preseg
# init + 300 ps continuation; post-burnin imbalance -0.548, semi-
# upper-bound) with a plateau-pinned fdseed value so the 8-point
# headline curve is methodologically uniform (all fdseed). The
# previous trajectory was the closest-to-plateau of the dilute trio
# but still descending, and methodologically diverged from the
# fdseed runs at X_c >= 0.10.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.075_fdseed
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
SNAPSHOT=$RUN_DIR/poly_AlMg_200A_fdseed_T500K_Xc0.075.lmp

if [ ! -f "$SNAPSHOT" ]; then
  echo "ERROR: snapshot not found: $SNAPSHOT" >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

T=500.0
XC=0.075
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.075_fdseed"
FD_PRED=0.3007

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : Fermi-Dirac-seeded init, 16-rank, v2 deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Initial X_GB   : X_GB^FD = $FD_PRED (random GB+bulk placement)"
echo "canon-FD pred  : $FD_PRED (ours)"
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
