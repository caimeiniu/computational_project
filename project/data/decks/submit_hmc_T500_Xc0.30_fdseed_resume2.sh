#!/bin/bash
# Second resume of T=500 K, X_c=0.30 fdseed run (saturation interior).
# Lineage:
#   - original: job 65485900 (submit_hmc_T500_Xc0.30_fdseed.sh)
#       TIMEOUT at step 265100 / 310000; X_GB^HMC = 0.4874 [0.4830, 0.4918]
#       (UB; imb=-0.596, drift=-0.055 across PROD).
#   - resume1: job 66167089 (submit_hmc_T500_Xc0.30_fdseed_resume.sh)
#       COMPLETED 2026-05-12 04:08, 11:15:28 wall, EXTRA_PS=300.
#       Post-process (output/hmc_T500_Xc0.30_fdseed_resume.json):
#         x_gb mean = 0.4613 [0.4609, 0.4616]  (Q1..Q5 = [0.4643,0.4630,
#         0.4619,0.4607,0.4672], drift +0.003)
#         fwd=136, rev=551, net=-415, imb_signed=-0.604, fwd/rev=0.247
#         accept_rate_overall = 9.95 %
#         fd_predicted = 0.5337, gap = -0.0724
#       Verdict: descent room recovered (UB lowered by 0.026 from
#       resume1) but |imb|=0.604 still WAY above the 0.3 plateau
#       threshold and dwarfs the 0.46 imb that T=700 X_c=0.10 plateaued
#       at. Q5 actually rebounded above Q4 (0.467 vs 0.461) -> the
#       trajectory bounces around the descending mean rather than
#       sitting on a plateau. Conclusion: this point still needs
#       another 24 h block of descent before reading.
#
# Why submit now: X_c=0.30 is the saturation interior anchor of panel
# (d) - one of the two right-arm points (with X_c=0.40) that need
# plateau classification to convert the headline figure from "9 UB
# triangles + 1 plateau" toward "saturation arm closes". The other
# two overnight resumes (X_c=0.40, T=700 X_c=0.20) had imb_signed
# -0.42 and -0.46 respectively - both LESS deep than X_c=0.30's -0.60,
# so X_c=0.30 has the most descent room remaining and is the best
# single use of the 16-CPU slot today (cf. quota: 32/48 already in
# use by the dilute-arm running pair X_c=0.01/0.03).
#
# Resume picks rst2 (2026-05-12 04:08, NEWER than rst1 at 03:56).
# The resume1 _final.lmp also exists (also 04:08) but read_restart
# preserves velocities + box state more faithfully than read_data.
#
# Uses hmc_AlMg_resume.lammps; same fixes/dump cadence as resume1.
# f_hmc[1]/f_hmc[2] reset on re-attach as before; post-process treats
# the resume2 window as its own production block, written to
# output/hmc_T500_Xc0.30_fdseed_resume2.{json,png}. Resume1's window
# stats stay in hmc_T500_Xc0.30_fdseed_resume.json.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.30_fdseed_resume2
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
RSTFILE=$RUN_DIR/hmc_T500_Xc0.30_fdseed_resume.rst2

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  exit 1
fi

cd "$RUN_DIR"

T=500.0
XC=0.30
EXTRA_PS=300
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.30_fdseed_resume2"
FD_PRED=0.5337

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : RESUME2 of fdseed saturation interior; resume deck"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "X_GB at resume : ~0.461 (resume1 post-burnin mean; UB, imb=-0.604)"
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
