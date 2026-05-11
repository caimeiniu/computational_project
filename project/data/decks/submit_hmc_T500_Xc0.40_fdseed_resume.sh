#!/bin/bash
# Resume of T=500 K, X_c=0.40 fdseed run (2026-05-09 saturation-edge).
# Original run (job 66002267, submit_hmc_T500_Xc0.40_fdseed.sh) COMPLETED
# at the 24h SLURM wall with auto-postprocess output
# hmc_T500_Xc0.40_fdseed.{json,png} written 2026-05-11 09:39. However,
# post-burn-in fwd/rev imbalance was -0.497 (REV > FWD) with Q1->Q5
# drift -0.049 across the production window (0.577 -> 0.528), so the
# reported mean X_GB^HMC = 0.5443 [0.5405, 0.5480] oversamples the
# relaxation phase: it is a STRICT UPPER BOUND on the true equilibrium
# X_GB at X_c=0.40, T=500 K.
#
# Why resume now: X_c=0.40 is the saturation-edge anchor of the
# headline 9-pt panel (d) all-fdseed (CHANGELOG 2026-05-11 evening).
# Its gap_O = -0.0406 is the shallowest of the nine red circles; if a
# resume confirms the gap survives below FD by more than CI95 at the
# saturation edge, the V-shape of the breakdown band has a confirmed
# saturation arm. If the resume drives X_GB further down (likely,
# given imb -0.497), the gap deepens and the band-edge claim
# strengthens.
#
# Recipe is the direct analogue of the two prior fdseed resumes that
# plateaued:
#   - T=700 X_c=0.10: imb -0.568 -> -0.453, CI95 hw 0.0044 -> 0.0003
#     (CHANGELOG 2026-05-10 afternoon)
#   - T=500 X_c=0.04: imb -0.469 -> -0.302, CI95 hw 0.0025 -> 0.00009
#     (CHANGELOG 2026-05-11 morning post-job; job 66066111)
# X_c=0.40 starts from imb -0.497, between the two precedents.
#
# Resume picks rst2 (09:36, newer than rst1 at 09:13); the original
# write_data wrote _final.lmp at 09:38 (2 min after rst2). Uses
# hmc_AlMg_resume.lammps which read_restart's positions/velocities/box
# and re-attaches NVT + atom/swap without re-doing CG or EQUIL.
#
# Note: f_hmc[1]/f_hmc[2] counters reset on fix re-attach. Post-process
# below treats the resume window as its own production block; the
# original window's stats stay in hmc_T500_Xc0.40_fdseed.json.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.40_fdseed_resume
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
RSTFILE=$RUN_DIR/hmc_T500_Xc0.40_fdseed.rst2

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  exit 1
fi

cd "$RUN_DIR"

T=500.0
XC=0.40
EXTRA_PS=300
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.40_fdseed_resume"
FD_PRED=0.5848

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : RESUME of fdseed saturation-edge; resume deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "X_GB at resume : ~0.524 (last 10 frames of original run, upper bound)"
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
