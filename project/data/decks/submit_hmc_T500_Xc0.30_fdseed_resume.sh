#!/bin/bash
# Resume of T=500 K, X_c=0.30 fdseed run (2026-05-07 saturation interior).
# Original run (job 65485900, submit_hmc_T500_Xc0.30_fdseed.sh) was a
# TIMEOUT at step 265100 / 310000 (~85%) -- written 2026-05-08 07:37.
# Post-process from the timeout-residue dump (CHANGELOG 2026-05-08 early
# afternoon, panel (d) 7-pt) gave X_GB^HMC = 0.4874 [0.4830, 0.4918] with
# fwd/rev imbalance -0.596 (REV > FWD) and Q1->Q5 drift -0.055 across the
# production window (0.524 -> 0.470), so the reported mean is a STRICT
# UPPER BOUND on the true equilibrium X_GB at X_c=0.30, T=500 K.
#
# Why resume now: X_c=0.30 is the second-shallowest of the 9-pt panel (d)
# all-fdseed (CHANGELOG 2026-05-11 evening), gap_O = -0.0462 - just above
# X_c=0.40's -0.0406. Together these two form the saturation arm of the
# V-shaped breakdown band; both are currently upper bounds. Pairing the
# X_c=0.40 resume submitted in this batch with the X_c=0.30 resume gives
# us two consecutive plateau reads on the saturation arm, replacing the
# upper-bound caveat on the right half of panel (d).
#
# Starting imb -0.596 is comparable to T=700 X_c=0.10's -0.568 (which
# plateaued in 12.5 h with imb -> -0.453, CI95 hw 0.0003, CHANGELOG
# 2026-05-10 afternoon). Swap accept at X_c=0.30 was 8.2% in the original
# window, between dilute (~5%) and saturation (~11%) -> 24 h resume budget
# should be sufficient.
#
# Resume picks rst1 (2026-05-08 08:05, NEWER than rst2 at 07:37; the
# original log's last thermo line was step 265100 ~ 08:11 SLURM-wall
# hit time, so rst1 is the more recent of the alternating checkpoint
# pair and the right resume seed). If rst1 turns out to be corrupted,
# fall back to rst2.
#
# Uses hmc_AlMg_resume.lammps which read_restart's positions/velocities/
# box and re-attaches NVT + atom/swap without re-doing CG or EQUIL.
#
# Note: f_hmc[1]/f_hmc[2] counters reset on fix re-attach. Post-process
# below treats the resume window as its own production block; the
# original window's stats stay in hmc_T500_Xc0.30_fdseed.json.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.30_fdseed_resume
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
RSTFILE=$RUN_DIR/hmc_T500_Xc0.30_fdseed.rst1

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  exit 1
fi

cd "$RUN_DIR"

T=500.0
XC=0.30
EXTRA_PS=300
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.30_fdseed_resume"
FD_PRED=0.5337

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : RESUME of fdseed saturation interior; resume deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "X_GB at resume : ~0.465 (last 10 frames of original timeout, upper bound)"
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
