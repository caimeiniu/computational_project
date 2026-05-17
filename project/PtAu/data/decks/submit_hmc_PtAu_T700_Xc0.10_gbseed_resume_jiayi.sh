#!/bin/bash
# Resume the Pt(Au) GB-seeded reverse initial-condition check at
# T=700 K, X_c=0.10. This continues only the over-segregated-start run,
# whose first 600 ps remained above the random-start plateau.

#SBATCH --job-name=hmc_PtAu_Xc0.10_gbseed_resume
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/jiayfu/prototype_PtAu_100A/%x-%j.out
#SBATCH --error=/cluster/scratch/jiayfu/prototype_PtAu_100A/%x-%j.err

set -euo pipefail

module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4 python/3.11.6

PROJECT=/cluster/home/jiayfu/computational_project/project
RUN_DIR=/cluster/scratch/jiayfu/prototype_PtAu_100A
DECK=$PROJECT/PtAu/data/decks/hmc_PtAu_resume.lammps
GB_MASK=$RUN_DIR/gb_mask_PtAu_100A.npy
POTENTIAL=$PROJECT/PtAu/data/potentials/PtAu.eam.alloy

T=700.0
XC=0.10
PROD_PS=${PROD_PS:-600}
SWAPS_PER_CALL=${SWAPS_PER_CALL:-100}
RESTART_PS=${RESTART_PS:-5}
OLD_STUB=hmc_PtAu_T700_Xc0.10_gbseed
OUTSTUB=hmc_PtAu_T700_Xc0.10_gbseed_resume

if [[ ! -d "$RUN_DIR" ]]; then
  echo "ERROR: run directory missing: $RUN_DIR" >&2
  exit 1
fi

if [[ ! -f "$GB_MASK" ]]; then
  echo "ERROR: GB mask missing: $GB_MASK" >&2
  exit 1
fi

RESTART=$(ls -t "$RUN_DIR/${OLD_STUB}".rst1 "$RUN_DIR/${OLD_STUB}".rst2 2>/dev/null | head -n 1 || true)
if [[ -z "$RESTART" ]]; then
  echo "ERROR: no restart found for ${OLD_STUB}.rst1/.rst2 in $RUN_DIR" >&2
  exit 1
fi

cd "$RUN_DIR"

echo "=========================================================="
echo "Job ID        : $SLURM_JOB_ID"
echo "Variant       : Pt(Au) GB-seeded Xc=0.10 resume"
echo "T             : $T K"
echo "X_c           : $XC"
echo "Extra prod ps : $PROD_PS"
echo "Run dir       : $RUN_DIR"
echo "Deck          : $DECK"
echo "Restart       : $RESTART"
echo "GB mask       : $GB_MASK"
echo "Started       : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var restart "$RESTART" \
  -var outstub "$OUTSTUB" \
  -var T "$T" \
  -var PROD_PS "$PROD_PS" \
  -var RESTART_PS "$RESTART_PS" \
  -var SWAPS_PER_CALL "$SWAPS_PER_CALL" \
  -var POTFILE "$POTENTIAL" \
  -var SWAP_SEED 20260631

echo "=========================================================="
echo "Resume sim done at $(date); auto-post-processing..."
echo "=========================================================="

mkdir -p "$PROJECT/PtAu/output"
python3 "$PROJECT/scripts/hmc_xgb_timeseries.py" \
  --stub "$RUN_DIR/$OLD_STUB" \
  --stub "$RUN_DIR/$OUTSTUB" \
  --gb-mask "$GB_MASK" \
  --xc "$XC" \
  --temp "$T" \
  --fd-pred nan \
  --out-prefix "$PROJECT/PtAu/output/$OUTSTUB" \
  || echo "POST-PROCESS FAILED but sim succeeded; rerun manually"

echo "=========================================================="
echo "Finished      : $(date)"
echo "=========================================================="
