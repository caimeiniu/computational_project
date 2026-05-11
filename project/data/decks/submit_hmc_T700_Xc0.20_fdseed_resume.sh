#!/bin/bash
# Resume of T=700 K, X_c=0.20 fdseed run (2026-05-10 saturation-edge T-axis).
# Original run (job 66027275, submit_hmc_T700_Xc0.20_fdseed.sh) COMPLETED
# at the 24h SLURM wall with auto-postprocess output
# hmc_T700_Xc0.20_fdseed.{json,png} written 2026-05-11 15:14. However,
# post-burn-in fwd/rev imbalance was -0.570 (REV > FWD) with Q1->Q5
# drift -0.071 across the production window (0.405 -> 0.335), so the
# reported mean X_GB^HMC = 0.3563 [0.3512, 0.3613] oversamples the
# relaxation phase: it is a STRICT UPPER BOUND on the true equilibrium
# X_GB at X_c=0.20, T=700 K.
#
# Why resume now: this point is the high-T anchor of the saturation-
# edge T-axis (panel_d_T_axis_X_c0.20_2pt produced 2026-05-11 evening;
# both T=500 and T=700 are currently upper bounds, gap_O = -0.066 and
# -0.064 respectively). A T=700 resume here is the direct analogue of
# the T=700 X_c=0.10 resume (CHANGELOG 2026-05-10 afternoon) that
# plateaued cleanly in 12.5 h with imb -0.568 -> -0.453 and CI95
# half-width 0.0003. The X_c=0.20 starting imb is slightly less
# imbalanced (-0.570 vs -0.568) and swap accept is slightly higher
# (11.2% vs 9.2%) -> plateau probability is at least as high. If this
# resume plateaus, the saturation-edge T-axis has a confirmed plateau
# at T=700, mirroring the breakdown-apex T-axis at X_c=0.10.
#
# Resume picks rst2 (15:12, newer than rst1 at 14:50); the original
# write_data wrote _final.lmp at 15:14 (2 min after rst2). Uses
# hmc_AlMg_resume.lammps which read_restart's positions/velocities/box
# and re-attaches NVT + atom/swap without re-doing CG or EQUIL.
#
# Note: f_hmc[1]/f_hmc[2] counters reset on fix re-attach. Post-process
# below treats the resume window as its own production block; the
# original window's stats stay in hmc_T700_Xc0.20_fdseed.json.

#SBATCH --job-name=hmc_AlMg_T700_Xc0.20_fdseed_resume
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
RSTFILE=$RUN_DIR/hmc_T700_Xc0.20_fdseed.rst2

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  exit 1
fi

cd "$RUN_DIR"

T=700.0
XC=0.20
EXTRA_PS=300
RESTART_PS=5
OUTSTUB="hmc_T700_Xc0.20_fdseed_resume"
FD_PRED=0.4198

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : RESUME of fdseed saturation-edge T-axis high-T; resume deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "X_GB at resume : ~0.330 (last 10 frames of original run, upper bound)"
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
