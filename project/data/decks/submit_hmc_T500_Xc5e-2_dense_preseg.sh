#!/bin/bash
# A2: equilibration verification at T=500 K, X_c=5e-2, PRE-SEGREGATED
# initial state — all 23 786 Mg atoms placed at randomly-chosen GB sites
# so X_GB(0) ≈ 0.2671 (= closed-box ceiling, well above canon-FD's 0.228).
#
# Pre-seg data file built once by scripts/pre_segregate.py; this submit
# expects it at scratch/production_AlMg_200A/poly_AlMg_200A_preseg_Xc5e-2.lmp.
#
# Same dense-sampling params as the random-IC twin (SWAPS=100, PROD=200).
# `-var SKIP_PLACE 1` tells the deck not to do its own random Mg placement;
# the data file already has type-2 atoms.

#SBATCH --job-name=hmc_AlMg_T500_Xc5e-2_dense_preseg
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/hmc_AlMg/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/hmc_AlMg/%x-%j.err

set -euo pipefail

module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

RUN_DIR=/cluster/scratch/cainiu/hmc_AlMg
DECK=/cluster/home/cainiu/Computational_modeling/project/data/decks/hmc_AlMg.lammps
PRESEG=/cluster/scratch/cainiu/production_AlMg_200A/poly_AlMg_200A_preseg_Xc5e-2.lmp

if [ ! -f "$PRESEG" ]; then
  echo "ERROR: pre-segregated data file not found: $PRESEG" >&2
  echo "Build it with: python scripts/pre_segregate.py --xc 0.05 ..." >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$PRESEG" poly_AlMg_200A_preseg_Xc5e-2.lmp

T=500.0
XC=0.05
EQUIL_PS=10
PROD_PS=200
SWAPS_PER_CALL=100
OUTSTUB="hmc_T${T%.*}_Xc${XC}_dense_preseg"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC  (canon-FD predicts X_GB ≈ 0.228)"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "Initial Mg     : PRE-SEGREGATED  X_GB(0) = 23786/89042 ≈ 0.267 (ceiling)"
echo "Started        : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var annealed       poly_AlMg_200A_preseg_Xc5e-2.lmp \
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
