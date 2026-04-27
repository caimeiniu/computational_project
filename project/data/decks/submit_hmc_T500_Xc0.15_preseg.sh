#!/bin/bash
# Preseg restart at T=500 K, X_c=0.15 (random-IC sweep failed kinetically
# even at X_c >= 0.10 — see CHANGELOG 2026-04-27). Starts above canon-FD
# (X_GB(0) at ceiling or =1 if X_c > GB_frac) and descends.
#
# Pre-seg data file built by scripts/pre_segregate.py with seed 20260427.
# `-var SKIP_PLACE 1` skips the deck's random Mg placement (data file
# already has type-2 atoms).

#SBATCH --job-name=hmc_AlMg_T500_Xc0.15_preseg
#SBATCH --time=12:00:00
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
PRESEG=/cluster/scratch/cainiu/production_AlMg_200A/poly_AlMg_200A_preseg_Xc0.15.lmp

if [ ! -f "$PRESEG" ]; then
  echo "ERROR: pre-segregated data file not found: $PRESEG" >&2
  echo "Build it with: python scripts/pre_segregate.py --xc 0.15 ..." >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$PRESEG" $(basename "$PRESEG")

T=500.0
XC=0.15
EQUIL_PS=10
PROD_PS=200
SWAPS_PER_CALL=100
OUTSTUB="hmc_T${T%.*}_Xc0.15_preseg"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "Initial Mg     : PRE-SEGREGATED (see preseg lmp header for X_GB(0)/X_bulk(0))"
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
