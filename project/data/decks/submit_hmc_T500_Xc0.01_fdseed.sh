#!/bin/bash
# FD-seeded HMC run (2026-05-11): T=500 K, X_c=0.01 — far-dilute break-point probe.
#
# Initial condition: poly_AlMg_200A_fdseed_T500K_Xc0.01.lmp from
# build_fdseed_inits.py (regenerated 2026-05-11 with deterministic
# seed=20260505), X_GB(0) = X_GB^FD = 0.0526 at random GB sites
# (N_Mg(GB) = 4,683) and the remaining Mg at random bulk sites
# (N_Mg(bulk) = 74). Total Mg only 4,757 atoms (vs 19,029 at X_c=0.04).
#
# Caveats unique to this far-dilute point:
#   - Closed-box ceiling at X_c=0.01 is 0.0534, so the FD seed sits
#     0.0008 below ceiling (far below the 0.02 "viable threshold" used
#     at X_c=0.04). However the descent room from FD to 0 is 0.0526,
#     so any breakdown that pushes Mg into bulk is observable in
#     principle — this is NOT a mass-conservation-pinned point.
#   - Few-atom statistics: 4,683 GB Mg vs 89,042 GB sites; per-frame
#     X_GB single-bin noise sqrt(p(1-p)/N_GB) ≈ 0.0007. Block-bootstrap
#     over 240 frames should still give CI95 half-width ~0.0001
#     (same order as X_c=0.04 plateau).
#   - Swap rate will be the lowest of any HMC point we've run
#     (≈ 2 %, scaling with sqrt(X_GB^FD)); EXPECT this to need a
#     resume even more than X_c=0.04 did.
#
# Purpose: directly test Wagih's dilute-limit recovery claim. If the
# site-independent FD becomes exact as X_c → 0 (Mg-Mg pair density
# ∝ (X_GB^FD)² → 0), gap_O at X_c=0.01 should be ≤ 0.005 (within CI95
# of zero). Pair-density scaling from X_c=0.04 (gap -0.060) predicts
# gap_O at X_c=0.01 ≈ -0.005. A nonzero gap at X_c=0.01 would imply
# the breakdown has a floor that does NOT vanish in the dilute limit,
# which is a bigger physics story than the current V-band finding.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.01_fdseed
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
SNAPSHOT=$RUN_DIR/poly_AlMg_200A_fdseed_T500K_Xc0.01.lmp

if [ ! -f "$SNAPSHOT" ]; then
  echo "ERROR: snapshot not found: $SNAPSHOT" >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

T=500.0
XC=0.01
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.01_fdseed"
FD_PRED=0.0526

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : Fermi-Dirac-seeded init, 16-rank, v2 deck; far-dilute break probe"
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
