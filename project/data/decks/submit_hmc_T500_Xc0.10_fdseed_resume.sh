#!/bin/bash
# Resume of T=500 K, X_c=0.10 fdseed run (headline panel-d central point).
# Original run (job ID see hmc_T500_Xc0.10_fdseed-*.out, COMPLETED 2026-05-07
# 09:40 CEST) reached PROD_PS=300 ps but X_GB(t) was still descending at the
# end: start (fdseed) 0.351 → 50% 0.283 → 80% 0.263 → end 0.253. The
# reported mean X_GB^HMC = 0.279 (CI95 [..., ...]) is a strict upper bound on
# the true equilibrium X_GB^∞ at (T=500 K, X_c=0.10).
#
# Why resume now (2026-05-13 afternoon, post-Pt(Au) handoff):
#
# 1. PURPOSE A — headline ▽ → ● conversion. X_c=0.10 is the central point
#    of the 9-pt panel (d) at T=500 K and the point most often quoted as
#    "broken band maximum". Currently it is drawn as ▽ (UB) per the 2026-05-11
#    lockdown classification (drift -0.073/window, imbalance -0.691). The
#    X_c=0.04 resume (job 65732267 → resume) successfully went 0.143 → 0.131
#    plateau, demonstrating the resume recipe. If this X_c=0.10 resume also
#    plateaus, the headline central point becomes ●, simplifying the caption
#    and tightening the broken-band claim numerically.
#
# 2. PURPOSE B — supply a converged snapshot for the three mechanism figures
#    (Mg-Mg clustering 01, occupation P_i 02, Mg-Mg repulsion 03). Until now
#    these figures sourced a snapshot of X_c=0.075 which is NOT in plateau.
#    Replacing the source with X_c=0.10 fdseed_resume _final.lmp (whose run
#    will have ~ extra 300 ps PROD beyond the original 300 ps) puts the
#    mechanism figures at a properly equilibrated state in the middle of the
#    broken band, which is exactly where Mg-Mg interaction effects are most
#    visible.
#
# Outcome branches:
# - plateau (drift in last-1/3 < ~0.005): convert ▽ → ●; use _final.lmp for
#   mechanism figs.
# - still drifting (mean continues descending): tighter ▽ at lower mean;
#   _final.lmp still a better-than-X_c=0.075 mechanism-fig source (closer to
#   equilibrium).
#
# Resume picks rst2 (2026-05-07 09:37, newer than rst1 at 09:14); original
# write_data wrote _final.lmp at 09:40 (3 min after rst2). Uses
# hmc_AlMg_resume.lammps which read_restart's positions/velocities/box and
# re-attaches NVT + atom/swap without re-doing CG or EQUIL.
#
# Note: f_hmc[1]/f_hmc[2] counters reset on fix re-attach. Post-process
# below treats the resume window as its own production block; the original
# window's stats stay in hmc_T500_Xc0.10_fdseed.json.

#SBATCH --job-name=hmc_AlMg_T500_Xc0.10_fdseed_resume
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
RSTFILE=$RUN_DIR/hmc_T500_Xc0.10_fdseed.rst2

if [ ! -f "$RSTFILE" ]; then
  echo "ERROR: restart file not found: $RSTFILE" >&2
  exit 1
fi

cd "$RUN_DIR"

T=500.0
XC=0.10
EXTRA_PS=300
RESTART_PS=5
OUTSTUB="hmc_T500_Xc0.10_fdseed_resume"
FD_PRED=0.3519

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Variant        : RESUME of fdseed central-headline point"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EXTRA_PS       : $EXTRA_PS"
echo "RESTART_PS     : $RESTART_PS"
echo "MPI ranks      : $SLURM_NTASKS"
echo "Resume from    : $(basename "$RSTFILE") ($(stat -c %y "$RSTFILE"))"
echo "X_GB at resume : ~0.253 (end of original PROD; reported mean 0.279 was UB)"
echo "canon-FD pred  : $FD_PRED (ours, X_c=0.10, T=500 K)"
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
