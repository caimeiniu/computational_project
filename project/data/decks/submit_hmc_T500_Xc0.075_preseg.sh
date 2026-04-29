#!/bin/bash
# Bonus stretch run (CHANGELOG 2026-04-29 evening): X_c=0.075 preseg HMC
# at T=500 K to probe the breakdown threshold X_c* in (0.05, 0.10).
# X_c=0.05 has Wagih agreement (0.238 vs FD 0.228), X_c=0.10 multistart
# UB 0.228 ≪ FD 0.352. 0.075 is the binary-search midpoint.
#
# Pre-seg data file built by scripts/pre_segregate.py with seed 20260429.
# X_GB(0) = 0.4007 (all 35,679 Mg in GB; X_c=0.075 < GB_frac=0.1872 means
# GB capacity is sufficient, no bulk overflow).
# `-var SKIP_PLACE 1` skips the deck's random Mg placement.
#
# 24 h budget — generous buffer relative to past ~5–9 h preseg runs;
# extra room because behavior near the threshold is unknown.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.075_preseg
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
PRESEG=/cluster/scratch/cainiu/production_AlMg_200A/poly_AlMg_200A_preseg_Xc0.075.lmp

if [ ! -f "$PRESEG" ]; then
  echo "ERROR: pre-segregated data file not found: $PRESEG" >&2
  echo "Build it with: python scripts/pre_segregate.py --xc 0.075 ..." >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$PRESEG" $(basename "$PRESEG")

T=500.0
XC=0.075
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
OUTSTUB="hmc_T${T%.*}_Xc0.075_preseg"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "Initial Mg     : PRE-SEGREGATED (X_GB(0) = 0.4007, all in GB)"
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
