#!/bin/bash
# Threshold binary-search continuation (CHANGELOG 2026-04-30 morning):
# X_c=0.075 preseg UB descended to 0.2543 < canon-FD 0.3007 → breakdown
# at 0.075. Threshold X_c* now ∈ (0.05, 0.075]. X_c=0.06 = next midpoint.
#
# Pre-seg data file built by scripts/pre_segregate.py with seed 20260430.
# X_GB(0) = 0.3206 (all 28,543 Mg in GB; X_c=0.06 < GB_frac=0.1872 means
# GB capacity is sufficient, no bulk overflow).
# `-var SKIP_PLACE 1` skips the deck's random Mg placement.
#
# 24 h budget mirrors X_c=0.075 (which used 10h 20min of 24h).

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

RUN_DIR=/cluster/scratch/cainiu/hmc_AlMg
DECK=/cluster/home/cainiu/Computational_modeling/project/data/decks/hmc_AlMg.lammps
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
echo "Finished       : $(date)"
echo "=========================================================="
