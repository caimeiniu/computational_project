#!/bin/bash
# Resume of T=500 K, X_c=0.04 fdseed run (Task A, 2026-05-08).
# Original run (job 65732267, submit_hmc_T500_Xc0.04_fdseed.sh) COMPLETED
# naturally at step 310 000 (PROD_PS=300 ps). However post-burn-in fwd/rev
# imbalance was -0.469 (REV > FWD), so X_GB(t) was still net-descending
# at the end of the original window. The reported mean
# X_GB^HMC = 0.1434 [0.1409, 0.1459] is a strict upper bound on the
# true equilibrium X_GB^infinity at X_c=0.04, T=500 K.
#
# Why resume now: X_c=0.04/0.05/0.06 are the three dilute kinetic-caveat
# points in the publication-ready 8-pt panel (d) all-fdseed (CHANGELOG
# 2026-05-10 morning). X_c=0.04 has the LEAST-reverse imbalance of the
# trio (-0.469 vs -0.557 / -0.622 for X_c=0.05/0.06), making it the most
# likely to plateau in a 24-h resume. Same starting condition as the
# T=700 X_c=0.10 resume (imbalance -0.568 -> plateaued to -0.453 with
# CI95 half-width 0.0003 within 12.5 h, CHANGELOG 2026-05-10 afternoon).
# If this resume plateaus, the dilute-edge of the headline panel (d)
# becomes a plateau read instead of an upper bound, simplifying the
# figure caption.
#
# Resume picks rst2 (16:09, newer than rst1 at 15:51); the original
# write_data wrote _final.lmp at 16:11 (2 min after rst2). Uses
# hmc_AlMg_resume.lammps which read_restart's positions/velocities/box
# and re-attaches NVT + atom/swap without re-doing CG or EQUIL.
#
# Note: f_hmc[1]/f_hmc[2] counters reset on fix re-attach. Post-process
# below treats the resume window as its own production block; the
# original window's stats stay in hmc_T500_Xc0.04_fdseed.json.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.04_fdseed_resume
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
RSTFILE=$RUN_DIR/hmc_T500_Xc0.04_fdseed.rst2

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  exit 1
fi

cd "$RUN_DIR"

T=500.0
XC=0.04
EXTRA_PS=300
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.04_fdseed_resume"
FD_PRED=0.1912

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : RESUME of fdseed Task A; resume deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "X_GB at resume : ~0.143 (post-burnin mean of original run, upper bound)"
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
