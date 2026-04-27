#!/bin/bash
# HMC point at T=500 K, X_c=5e-2 (well above FD knee, expect X_GB^FD ≈ 0.43).
# Most informative single point for breakdown detection — if X_GB^HMC
# disagrees with FD here, dilute-limit assumption has broken down.

#SBATCH --job-name=hmc_AlMg_T500_Xc5e-2
#SBATCH --time=02:00:00
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
SRC=/cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g_annealed.lmp

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$SRC" poly_Al_200A_16g_annealed.lmp

T=500.0
XC=0.05            # 5e-2, well above FD knee at 500 K (FD predicts X_GB ≈ 0.43)
EQUIL_PS=10
PROD_PS=50
OUTSTUB="hmc_T${T%.*}_Xc${XC}"

echo "=========================================================="
echo "Job ID    : $SLURM_JOB_ID"
echo "Node list : $SLURM_JOB_NODELIST"
echo "MPI tasks : $SLURM_NTASKS"
echo "T         : $T K"
echo "X_c       : $XC  (FD prediction X_GB^FD ≈ 0.43)"
echo "EQUIL_PS  : $EQUIL_PS  ;  PROD_PS : $PROD_PS"
echo "Started   : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var annealed       poly_Al_200A_16g_annealed.lmp \
  -var outstub        "$OUTSTUB" \
  -var T              "$T" \
  -var XC             "$XC" \
  -var EQUIL_PS       "$EQUIL_PS" \
  -var PROD_PS        "$PROD_PS"

echo "=========================================================="
echo "Finished  : $(date)"
echo "=========================================================="
