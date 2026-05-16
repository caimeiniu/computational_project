#!/bin/bash
# First Pt(Au) HMC point: T=700 K, X_c=0.10, random/uniform Au initial condition.
#
# FD prediction from PtAu/output/fd_curves_PtAu_100A.json is approximately
# X_GB^FD = 0.16 for our 100 A Pt(Au) spectrum, so this random IC starts below
# FD at X_GB(0) ~= X_c = 0.10 and should climb if HMC equilibrates.

#SBATCH --job-name=hmc_PtAu_T700_Xc0.10_random
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
DECK=$PROJECT/PtAu/data/decks/hmc_PtAu.lammps
SNAPSHOT=$RUN_DIR/poly_Pt_100A_8g_annealed.lmp
GB_MASK=$RUN_DIR/gb_mask_PtAu_100A.npy
POTENTIAL=$PROJECT/PtAu/data/potentials/PtAu.eam.alloy

if [[ ! -f "$SNAPSHOT" ]]; then
  echo "ERROR: annealed structure missing: $SNAPSHOT" >&2
  exit 1
fi

if [[ ! -f "$GB_MASK" ]]; then
  echo "ERROR: GB mask missing: $GB_MASK" >&2
  echo "Run submit_delta_e_PtAu_100A_jiayi.sh first, or run gb_identify.py." >&2
  exit 1
fi

cd "$RUN_DIR"

T=700.0
XC=0.10
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
RESTART_PS=5
OUTSTUB=hmc_PtAu_T700_Xc0.10_random
FD_PRED=0.1605

echo "=========================================================="
echo "Job ID        : $SLURM_JOB_ID"
echo "Variant       : Pt(Au) random/uniform IC"
echo "T             : $T K"
echo "X_c           : $XC"
echo "FD pred       : $FD_PRED"
echo "Run dir       : $RUN_DIR"
echo "Deck          : $DECK"
echo "Snapshot      : $SNAPSHOT"
echo "GB mask       : $GB_MASK"
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
  --fd-pred "$FD_PRED" \
  --out-prefix "$PROJECT/PtAu/output/$OUTSTUB" \
  || echo "POST-PROCESS FAILED but sim succeeded; rerun manually"

echo "=========================================================="
echo "Finished      : $(date)"
echo "=========================================================="
