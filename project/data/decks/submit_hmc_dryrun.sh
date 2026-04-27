#!/bin/bash
# HMC dry-run on the production 200A Al(Mg) box.
# Single (T, X_c) point chosen near the FD knee for T=500 K → verifies:
#   (a) `fix atom/swap` runs to completion;
#   (b) Metropolis acceptance rate lands in 5–30 % regime;
#   (c) PE plateaus before measurement window;
#   (d) X_GB(t) post-process pipeline (gb_mask cross + block average) works.
#
# Outputs land in scratch dir; deck path is absolute so we can sbatch from
# anywhere.

#SBATCH --job-name=hmc_dry_AlMg_T500_Xc5e-3
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

# Source structure: the production 200A annealed pure-Al box.
SRC=/cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g_annealed.lmp

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

# Stage the annealed structure (symlink — read-only, saves disk).
ln -sf "$SRC" poly_Al_200A_16g_annealed.lmp

T=500.0          # K
XC=0.005         # 5e-3, near 500 K FD knee (X_GB^FD ≈ 0.21)
EQUIL_PS=10
PROD_PS=50
OUTSTUB="hmc_dry_T${T%.*}_Xc${XC}"

echo "=========================================================="
echo "Job ID    : $SLURM_JOB_ID"
echo "Node list : $SLURM_JOB_NODELIST"
echo "MPI tasks : $SLURM_NTASKS"
echo "Run dir   : $(pwd)"
echo "Deck      : $DECK"
echo "Source    : $SRC"
echo "T         : $T K"
echo "X_c       : $XC"
echo "EQUIL_PS  : $EQUIL_PS"
echo "PROD_PS   : $PROD_PS"
echo "outstub   : $OUTSTUB"
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
