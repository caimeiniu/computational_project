#!/bin/bash
# X_c=0.10 multi-start preseg at T=500 K. Initial condition X_GB(0)=0.300
# (BELOW canon-FD prediction 0.357), built by:
#   pre_segregate.py --xc 0.10 --xgb-init 0.30 --seed 20260428
#
# Companion to the descending-IC run hmc_T500_Xc0.10_preseg (X_GB(0)=0.534
# → 0.375 after 12 h, fwd/rev=0.069 still descending). Together the two
# trajectories form a two-sided bracket: ascending from 0.30 + descending
# from 0.375 → meeting point ≈ true equilibrium.
#
# Outcomes:
#   meet ≈ 0.36 (= canon-FD)  →  Wagih works at X_c=0.10; breakdown > 0.10
#   meet < 0.36 (e.g. 0.32)   →  Wagih over-predicts at X_c=0.10; breakdown ≤ 0.10
#   not meet (ascend → 0.34, descend → 0.37 untouched) → bracket [0.34, 0.37]
#
# `-var SKIP_PLACE 1` skips the deck's random Mg placement (data file
# already has type-2 atoms).

#SBATCH --job-name=hmc_AlMg_T500_Xc0.10_multistart_xgb0.3
#SBATCH --time=14:00:00
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
PRESEG=/cluster/scratch/cainiu/production_AlMg_200A/poly_AlMg_200A_multistart_Xc0.10_xgb0.30.lmp

if [ ! -f "$PRESEG" ]; then
  echo "ERROR: multi-start data file not found: $PRESEG" >&2
  echo "Build it with:" >&2
  echo "  python scripts/pre_segregate.py --xc 0.10 --xgb-init 0.30 \\" >&2
  echo "    --in-data /cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g_annealed.lmp \\" >&2
  echo "    --gb-mask /cluster/scratch/cainiu/production_AlMg_200A/gb_mask_200A.npy \\" >&2
  echo "    --seed 20260428 --out-data $PRESEG" >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$PRESEG" $(basename "$PRESEG")

T=500.0
XC=0.10
EQUIL_PS=10
# 300 ps at observed 25 ps/h ≈ 12 h compute; --time=14h gives 2 h headroom.
# Comparison with descending preseg (PROD_PS=200) is fair within the
# overlapping 200 ps window; the extra 100 ps lets fwd/rev → 1 after the
# trajectory has plateaued (so we can detect "reached equilibrium" rather
# than "still ascending at end of run").
PROD_PS=300
SWAPS_PER_CALL=100
OUTSTUB="hmc_T${T%.*}_Xc0.10_multistart_xgb0.3"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "Initial Mg     : MULTI-START X_GB(0)=0.300 (see preseg lmp header)"
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
