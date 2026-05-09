#!/bin/bash
# Resume of T=700 K, X_c=0.10 fdseed run (Task C, 2026-05-09).
# Original run (job 65732280, submit_hmc_T700_Xc0.10_fdseed.sh) hit the
# 24h SLURM wall at step 272100 / 310000 (~88%). X_GB(t) was still
# monotonically descending from FD-seed (0.295 at t=0 -> 0.219 at last
# 50 frames), fwd/rev imbalance = -0.568 (REV >> FWD), so the reported
# mean 0.2369 oversamples the relaxation phase. Need to drive the
# system to steady state before reading off the T-axis upper point.
#
# Resume picks rst1 (08:10, newer than rst2 at 07:43; timeout was 08:12,
# so we lose ~2 min of compute). Uses hmc_AlMg_resume.lammps which
# read_restart's positions/velocities/box and re-attaches NVT + atom/swap
# without re-doing CG or EQUIL.
#
# Note: f_hmc[1]/f_hmc[2] counters reset on fix re-attach. Post-process
# below treats the resume window as its own production block; combine
# with the original run's counts manually if a global accept rate is
# wanted.

#SBATCH --job-name=hmc_AlMg_T700_Xc0.10_fdseed_resume
#SBATCH --time=24:00:00
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
DECK=$REPO/data/decks/hmc_AlMg_resume.lammps
RSTFILE=$RUN_DIR/hmc_T700_Xc0.10_fdseed.rst1

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  exit 1
fi

cd "$RUN_DIR"

T=700.0
XC=0.10
EXTRA_PS=300
RESTART_PS=5
OUTSTUB="hmc_T700_Xc0.10_fdseed_resume"
FD_PRED=0.2956

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : RESUME of fdseed Task C; resume deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "X_GB at resume : ~0.219 (last 50 frames of original run)"
echo "canon-FD pred  : $FD_PRED (ours)"
echo "Started        : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var rstfile    $(basename "$RSTFILE") \
  -var outstub    "$OUTSTUB" \
  -var T          "$T" \
  -var XC         "$XC" \
  -var EXTRA_PS   "$EXTRA_PS" \
  -var RESTART_PS "$RESTART_PS"

echo "=========================================================="
echo "Sim done at $(date); auto-post-processing..."
echo "=========================================================="

python3 $REPO/scripts/hmc_xgb_timeseries.py \
  --stub $RUN_DIR/$OUTSTUB \
  --gb-mask $REPO/data/snapshots/gb_mask_200A.npy \
  --xc $XC \
  --temp $T \
  --fd-pred $FD_PRED \
  --out-prefix $REPO/output/$OUTSTUB \
  || echo "POST-PROCESS FAILED but sim succeeded; rerun manually"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
