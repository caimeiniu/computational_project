#!/bin/bash
# B-sweep point at T=500 K, X_c=0.30. Production-grade dense sampling
# (SWAPS_PER_CALL=100, PROD_PS=200). Random IC. This X_c lies in the
# regime where canonical FD and grand-canonical FD differ meaningfully
# AND solute-solute interactions are physically expected — the actual
# dilute-limit-breakdown probe.

#SBATCH --job-name=hmc_AlMg_T500_Xc3e-1
#SBATCH --time=08:00:00
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
SRC=/cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g_annealed.lmp

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$SRC" poly_Al_200A_16g_annealed.lmp

T=500.0
XC=0.30
EQUIL_PS=10
PROD_PS=200
SWAPS_PER_CALL=100
OUTSTUB="hmc_T${T%.*}_Xc${XC}"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC  (B-sweep production point)"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "Initial Mg     : RANDOM"
echo "Started        : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var annealed       poly_Al_200A_16g_annealed.lmp \
  -var outstub        "$OUTSTUB" \
  -var T              "$T" \
  -var XC             "$XC" \
  -var EQUIL_PS       "$EQUIL_PS" \
  -var PROD_PS        "$PROD_PS" \
  -var SWAPS_PER_CALL "$SWAPS_PER_CALL"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
