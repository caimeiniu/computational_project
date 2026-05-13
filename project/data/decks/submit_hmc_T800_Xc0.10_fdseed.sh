#!/bin/bash
# FD-seeded HMC run (2026-05-13): T=800 K, X_c=0.10.
#
# Half of a two-IC convergence test at T=800 K (the other half is the
# bracket-from-below companion submit_hmc_T800_Xc0.10_random.sh).
# Advisor's suggestion: at higher T, HMC should ergodicize fast enough
# that random IC (ascends from X_GB=X_c=0.10) and FD-seeded IC
# (descends or stays at X_GB=X_GB^FD=0.2734) converge to the same
# equilibrium. If they do → HMC sampler is correct → low-T (T=500 K)
# breakdown findings are real physics, not sampler failure. If they
# don't → fundamental sampling problem.
#
# Why T=800 K: 0.86 × T_melt(Al=933 K) — well above existing T=700 K
# data, comfortably below GB-premelting risk. kT/kT(500) = 1.6 → swap
# rate substantially higher than the T=500 K runs.
#
# Why X_c=0.10: central pivot of T-axis subsweep, plateau ● already
# established at both T=500 K and T=700 K with fdseed IC; this is the
# T=800 K third point completing the trio.
#
# Initial condition: poly_AlMg_200A_fdseed_T800K_Xc0.1.lmp from
# build_fdseed_inits.py, X_GB(0) = X_GB^FD(800 K, 0.10) = 0.2734
# (vs T=500 K: 0.3519, 22% weaker as expected).

#SBATCH --job-name=hmc_AlMg_T800_Xc0.10_fdseed
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
SNAPSHOT=$RUN_DIR/poly_AlMg_200A_fdseed_T800K_Xc0.1.lmp

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
OUTSTUB="hmc_T800_Xc0.10_fdseed"
FD_PRED=0.2734

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : Fermi-Dirac-seeded init, 16-rank, v2 deck"
echo "T              : $T K  (0.86 T_melt(Al=933 K))"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Initial X_GB   : X_GB^FD = $FD_PRED (random GB+bulk placement)"
echo "canon-FD pred  : $FD_PRED (ours)"
echo "Companion      : submit_hmc_T800_Xc0.10_random.sh (X_GB(0)=0.10)"
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
