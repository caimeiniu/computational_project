#!/bin/bash
# 32-rank rerun of X_c=0.06 preseg with embedded auto-post-process.
# Supersedes the 16-rank job 65259659 (TIMEOUT at 188 ps PROD on 2026-05-05).
#
# Reason for rerun: job 65259659 ended descending at X_GB=0.2284
# (last-10 mean 0.2295) — already below canon-FD ours = 0.2611, but no
# `_final.lmp` written and trajectory not yet plateaued. A clean 32-rank
# run reaches 300 ps in ~10 h (mirrors X_c=0.075 wallclock 10h20m), giving
# equilibrated breakdown evidence at the threshold + a `_final.lmp` for
# mechanism analysis (master-figure panel f candidate at the threshold).
#
# Fairshare healthy (LevelFS 0.55-0.87 across partitions, was 0.089 on
# 2026-05-01); queue empty as of 2026-05-05 morning, so 32-rank should
# start without contention.
#
# Pre-seg data file built by scripts/pre_segregate.py with seed 20260430.
# X_GB(0) = 0.3206 (all 28,543 Mg in GB; X_c=0.06 < GB_frac=0.1872 means
# GB capacity is sufficient, no bulk overflow).
# `-var SKIP_PLACE 1` skips the deck's random Mg placement.
#
# 24 h budget mirrors X_c=0.075 (which used 10h20m of 24h).

#SBATCH --job-name=hmc_AlMg_T500_Xc0.06_preseg
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/hmc_AlMg/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/hmc_AlMg/%x-%j.err

set -euo pipefail

module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

REPO=/cluster/home/cainiu/Computational_modeling/project
RUN_DIR=/cluster/scratch/cainiu/hmc_AlMg
DECK=$REPO/data/decks/hmc_AlMg.lammps
PRESEG=/cluster/scratch/cainiu/production_AlMg_200A/poly_AlMg_200A_preseg_Xc0.06.lmp

if [ ! -f "$PRESEG" ]; then
  echo "ERROR: pre-segregated data file not found: $PRESEG" >&2
  echo "Build it with: python scripts/pre_segregate.py --xc 0.06 ..." >&2
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
OUTSTUB="hmc_T${T%.*}_Xc0.06_preseg"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "MPI ranks      : $SLURM_NTASKS (32-rank, autopp variant)"
echo "Initial Mg     : PRE-SEGREGATED (X_GB(0) = 0.3206, all in GB)"
echo "canon-FD pred  : X_GB ≈ 0.2611 (ours)  /  0.2727 (Wagih)"
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
  --fd-pred 0.2611 \
  --out-prefix $REPO/output/$OUTSTUB \
  || echo "POST-PROCESS FAILED but sim succeeded; rerun manually"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
