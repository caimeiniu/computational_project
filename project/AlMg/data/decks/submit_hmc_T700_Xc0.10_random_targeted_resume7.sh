#!/bin/bash
# Resume7 of the T=700 K, X_c=0.10 random-targeted LB chain.
# resume6 completed cleanly at X_GB=0.1988; bracket still open vs UB resume10:
# Q5-Q1 = -0.00071 (~6.8 CI half-widths), fwd/rev = 0.995.

#SBATCH --job-name=hmc_AlMg_T700_Xc0.10_random_targeted_resume7
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/home/cainiu/Computational_modeling/project/AlMg/results/%x-%j.out
#SBATCH --error=/cluster/home/cainiu/Computational_modeling/project/AlMg/results/%x-%j.err

set -euo pipefail

module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

REPO=/cluster/home/cainiu/Computational_modeling/project
RUN_DIR=/cluster/fs/scratchnv/06/cainiu/hmc_AlMg
DECK=$REPO/data/decks/hmc_AlMg_resume_targeted.lammps
RSTFILE=$RUN_DIR/hmc_T700_Xc0.10_random_targeted_resume6.rst2

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  ls -ld /cluster/fs/scratchnv/06/cainiu "$RUN_DIR" >&2 || true
  exit 1
fi

mkdir -p "$REPO/AlMg/results"
cd "$RUN_DIR"

T=700.0
XC=0.10
EXTRA_PS=300
SWAPS_PER_CALL=50
RESTART_PS=5
OUTSTUB="hmc_T700_Xc0.10_random_targeted_resume7"
FD_PRED=0.2956

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : targeted RESUME7 of T=700 K X_c=0.10 random/LB"
echo "Deck           : $DECK"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL  (per fix; total = 2x = 100 / 100 steps)"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "Previous X_GB  : resume6 = 0.1988; bracket still open"
echo "UB partner     : hmc_T700_Xc0.10_fdseed_targeted_resume11"
echo "canon-FD pred  : $FD_PRED"
echo "Started        : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var rstfile        "$(basename "$RSTFILE")" \
  -var outstub        "$OUTSTUB" \
  -var T              "$T" \
  -var XC             "$XC" \
  -var EXTRA_PS       "$EXTRA_PS" \
  -var RESTART_PS     "$RESTART_PS" \
  -var SWAPS_PER_CALL "$SWAPS_PER_CALL"

echo "=========================================================="
echo "Sim done at $(date); auto-post-processing (targeted variant)..."
echo "=========================================================="

python3 "$REPO/scripts/hmc_xgb_timeseries_targeted.py" \
  --stub "$RUN_DIR/$OUTSTUB" \
  --gb-mask "$REPO/AlMg/data/snapshots/gb_mask_200A.npy" \
  --xc "$XC" \
  --temp "$T" \
  --fd-pred "$FD_PRED" \
  --out-prefix "$REPO/AlMg/results/$OUTSTUB" \
  || echo "POST-PROCESS FAILED but sim succeeded; rerun manually"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
