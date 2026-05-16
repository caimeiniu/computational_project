#!/bin/bash
# Pt(Au) 700 K HMC concentration bracket, following the Al(Mg) FD-first
# sweep idea but centered below the already-converged Xc=0.10 Pt(Au) point.
#
# Run as:
#   sbatch submit_hmc_PtAu_T700_bracket_jiayi.sh
#
# Optional overrides:
#   PROD_PS=300 sbatch ...
#   XC_LIST="0.005 0.01 0.03 0.05 0.07" sbatch --array=0-4%3 ...

#SBATCH --job-name=hmc_PtAu_T700_bracket
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --array=0-4%3
#SBATCH --output=/cluster/scratch/jiayfu/prototype_PtAu_100A/%x-%A_%a.out
#SBATCH --error=/cluster/scratch/jiayfu/prototype_PtAu_100A/%x-%A_%a.err

set -euo pipefail

module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4 python/3.11.6

PROJECT=/cluster/home/jiayfu/computational_project/project
RUN_DIR=/cluster/scratch/jiayfu/prototype_PtAu_100A
DECK=$PROJECT/PtAu/data/decks/hmc_PtAu.lammps
SNAPSHOT=$RUN_DIR/poly_Pt_100A_8g_annealed.lmp
GB_MASK=$RUN_DIR/gb_mask_PtAu_100A.npy
POTENTIAL=$PROJECT/PtAu/data/potentials/PtAu.eam.alloy

T=700.0
XC_LIST=${XC_LIST:-"0.005 0.01 0.03 0.05 0.07"}
read -r -a XCS <<< "$XC_LIST"

if [[ ${SLURM_ARRAY_TASK_ID:-0} -ge ${#XCS[@]} ]]; then
  echo "ERROR: SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:-unset} outside XC_LIST length ${#XCS[@]}" >&2
  exit 1
fi

XC=${XCS[$SLURM_ARRAY_TASK_ID]}
XC_TAG=${XC/./p}

EQUIL_PS=${EQUIL_PS:-10}
PROD_PS=${PROD_PS:-600}
SWAPS_PER_CALL=${SWAPS_PER_CALL:-100}
RESTART_PS=${RESTART_PS:-5}
OUTSTUB=hmc_PtAu_T700_Xc${XC}_random

if [[ ! -f "$SNAPSHOT" ]]; then
  echo "ERROR: annealed structure missing: $SNAPSHOT" >&2
  exit 1
fi

if [[ ! -f "$GB_MASK" ]]; then
  echo "ERROR: GB mask missing: $GB_MASK" >&2
  exit 1
fi

cd "$RUN_DIR"

echo "=========================================================="
echo "Job ID        : $SLURM_JOB_ID"
echo "Array task    : ${SLURM_ARRAY_TASK_ID:-0}"
echo "Variant       : Pt(Au) 700 K concentration bracket"
echo "T             : $T K"
echo "X_total target: $XC"
echo "Run dir       : $RUN_DIR"
echo "Deck          : $DECK"
echo "Snapshot      : $SNAPSHOT"
echo "GB mask       : $GB_MASK"
echo "Outstub       : $OUTSTUB"
echo "Equil ps      : $EQUIL_PS"
echo "Prod ps       : $PROD_PS"
echo "Started       : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var annealed "$SNAPSHOT" \
  -var outstub "$OUTSTUB" \
  -var T "$T" \
  -var XC "$XC" \
  -var EQUIL_PS "$EQUIL_PS" \
  -var PROD_PS "$PROD_PS" \
  -var RESTART_PS "$RESTART_PS" \
  -var SWAPS_PER_CALL "$SWAPS_PER_CALL" \
  -var POTFILE "$POTENTIAL"

echo "=========================================================="
echo "Sim done at $(date); auto-post-processing..."
echo "=========================================================="

mkdir -p "$PROJECT/PtAu/output"
python3 "$PROJECT/scripts/hmc_xgb_timeseries.py" \
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
