#!/bin/bash
# FD-seeded HMC run (2026-05-10): T=500 K, X_c=0.40.
# Saturation-edge probe to close the broken-band right edge in panel (d).
#
# Initial condition: poly_AlMg_200A_fdseed_T500K_Xc0.4.lmp from
# build_fdseed_inits.py, which places X_GB(0) = X_GB^FD = 0.5848 at
# random GB sites and the remaining Mg at random bulk sites
# (X_bulk(0) = 0.3574). Total Mg inventory N_Mg = 190 286 = 0.4 · N_total.
#
# Why X_c = 0.40: at T = 500 K, the broken band's right edge has not been
# pinned. X_c = 0.30 sits at gap-O = -0.046 (smallest of the 8-pt panel),
# trending toward closure. X_c = 0.40 lets us bracket the closure point
# above 0.30 (predicted: gap-O \approx -0.02 to 0 if Wagih's
# dilute-limit assumption recovers as Mg saturates the favourable GB
# sites). Geometric ceiling at X_c = 0.40 is X_c · N_total / N_GB =
# 2.137 (capped at 1.0), so FD-seed at 0.585 has plenty of headroom
# and the run is not ceiling-locked.
#
# Why this initial condition: any net descent of X_GB(t) directly
# evidences X_GB^infinity < X_GB^FD (the central project claim); does
# not require reaching equilibrium. At X_c = 0.40 the trajectory may
# also ASCEND if the true equilibrium > FD prediction (only physically
# possible if site-site repulsion is now dominated by saturation /
# entropy mixing - watch for sign of d X_GB / d t).

#SBATCH --job-name=hmc_AlMg_T500_Xc0.40_fdseed
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
SNAPSHOT=$RUN_DIR/poly_AlMg_200A_fdseed_T500K_Xc0.4.lmp

if [ ! -f "$SNAPSHOT" ]; then
  echo "ERROR: snapshot not found: $SNAPSHOT" >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

T=500.0
XC=0.40
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.40_fdseed"
FD_PRED=0.5848

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
