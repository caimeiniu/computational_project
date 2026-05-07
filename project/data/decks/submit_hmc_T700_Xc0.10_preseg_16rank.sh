#!/bin/bash
# 16-rank variant of submit_hmc_T700_Xc0.10_preseg.sh (CHANGELOG 2026-05-01
# evening): original 32-rank job 65227656 cancelled after >30 h pending
# (CVE-2026-31431 cluster reboot wave + fairshare debt LevelFS=0.089).
# 16 ranks = 33% of public QOS cap (vs 67% at 32 ranks) → two such jobs
# can run concurrently; expected wallclock ~17 h, 20 h budget.
#
# T-axis exploration (CHANGELOG 2026-04-30 afternoon): T=500 K story is
# done — breakdown threshold X_c* ≤ 0.075. Now probe whether the
# breakdown persists at higher T.
#
# Pre-seg data file = same poly_AlMg_200A_preseg_Xc0.10.lmp used at T=500
# (preseg IC is purely positions+types; T enters only via the LAMMPS deck
# `-var T 700`). X_GB(0) = 0.40.
#
# canon-FD at T=700, X_c=0.10 = 0.2956 (vs 0.3519 at T=500). Trajectory
# starts at 0.40 → descends regardless of breakdown vs Wagih (Wagih's
# prediction is already below preseg IC).
#
# Auto-post-processing embedded at end of script: hmc_xgb_timeseries.py
# runs in the same allocation, writes output/<stub>.{json,png}.

#SBATCH --job-name=hmc_AlMg_T700_Xc0.10_preseg_16r
#SBATCH --time=20:00:00
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
DECK=$REPO/data/decks/hmc_AlMg.lammps
PRESEG=/cluster/scratch/cainiu/production_AlMg_200A/poly_AlMg_200A_preseg_Xc0.10.lmp

if [ ! -f "$PRESEG" ]; then
  echo "ERROR: pre-segregated data file not found: $PRESEG" >&2
  exit 1
fi

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
ln -sf "$PRESEG" $(basename "$PRESEG")

T=700.0
XC=0.10
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
OUTSTUB="hmc_T${T%.*}_Xc0.10_preseg"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS"
echo "PROD_PS        : $PROD_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "MPI ranks      : $SLURM_NTASKS (16-rank variant)"
echo "Initial Mg     : PRE-SEGREGATED (X_GB(0) = 0.40, all in GB)"
echo "canon-FD pred  : X_GB ≈ 0.2956 (ours) at T=700, X_c=0.10"
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
  --fd-pred 0.2956 \
  --out-prefix $REPO/output/$OUTSTUB \
  || echo "POST-PROCESS FAILED but sim succeeded; rerun manually"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
