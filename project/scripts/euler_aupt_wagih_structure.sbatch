#!/usr/bin/env bash
#SBATCH --job-name=aupt_wagih_structure
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --output=aupt_wagih_structure_%j.out

set -euo pipefail

module load stack/2024-06
module load python/3.11.6
module load openmpi/4.1.6
module load lammps/20230802

PROJECT_DIR="$(pwd)"
LMP_CMD="mpirun -np ${SLURM_NTASKS} lmp"
SOLVENT="${SOLVENT:-Au}"
SOLUTE="${SOLUTE:-Pt}"
WORKFLOW_NAME="${WORKFLOW_NAME:-aupt_wagih_structure}"
DATA_DIR="data/${WORKFLOW_NAME}"
RUN_DIR="/cluster/scratch/${USER}/${WORKFLOW_NAME}_${SLURM_JOB_ID}"
POTENTIAL="${POTENTIAL:-data/potentials/PtAu.eam.alloy}"
WAGIH_DUMP="${WAGIH_DUMP:?Set WAGIH_DUMP to the Au_Pt_20nm_GB_segregation.dump path}"
SAMPLES="${SAMPLES:-1000}"
BULK_REFERENCES="${BULK_REFERENCES:-50}"

echo "Project: ${PROJECT_DIR}"
echo "Scratch run directory: ${RUN_DIR}"
echo "LAMMPS command: ${LMP_CMD}"
echo "Alloy: ${SOLVENT}(${SOLUTE})"
echo "Wagih dump: ${WAGIH_DUMP}"
echo "Samples/bulk refs: ${SAMPLES} / ${BULK_REFERENCES}"
echo "Potential: ${POTENTIAL}"

python3 - <<'PY'
import numpy
import scipy
import matplotlib
print(f"NumPy: {numpy.__version__}")
print(f"SciPy: {scipy.__version__}")
print(f"Matplotlib: {matplotlib.__version__}")
PY

mkdir -p "${RUN_DIR}"
rsync -a \
  --exclude ".git" \
  --exclude "__pycache__" \
  --exclude "archive" \
  "${PROJECT_DIR}/" "${RUN_DIR}/"

cd "${RUN_DIR}"
mkdir -p "${DATA_DIR}/lammps"

echo "Step 1/4: convert Wagih dump structure to LAMMPS data"
python3 scripts/convert_wagih_dump_to_lammps_data.py \
  --dump "${WAGIH_DUMP}" \
  --output "${DATA_DIR}/lammps/${WORKFLOW_NAME}_structure.lammps" \
  --solvent "${SOLVENT}" \
  --solute "${SOLUTE}"

echo "Step 2/4: sample sites from Wagih nonzero segregation-energy atoms"
python3 scripts/select_wagih_dump_deltae_sites.py \
  --dump "${WAGIH_DUMP}" \
  --samples "${SAMPLES}" \
  --bulk-references "${BULK_REFERENCES}" \
  --output "${DATA_DIR}/${WORKFLOW_NAME}_deltae_sites.csv"

echo "Step 3/4: write and run substitution Delta E jobs"
python3 scripts/write_cuni_deltae_lammps_jobs.py \
  --sites "${DATA_DIR}/${WORKFLOW_NAME}_deltae_sites.csv" \
  --data "${DATA_DIR}/lammps/${WORKFLOW_NAME}_structure.lammps" \
  --potential "${POTENTIAL}" \
  --jobs-dir "${DATA_DIR}/lammps/deltae_jobs" \
  --lammps-command "${LMP_CMD}" \
  --solvent "${SOLVENT}" \
  --solute "${SOLUTE}"

cd "${DATA_DIR}/lammps/deltae_jobs"
./run_all_deltae.sh
cd "${RUN_DIR}"

echo "Step 4/4: analyze Delta E spectrum"
python3 scripts/analyze_cuni_deltae.py \
  --sites "${DATA_DIR}/${WORKFLOW_NAME}_deltae_sites.csv" \
  --jobs-dir "${DATA_DIR}/lammps/deltae_jobs" \
  --output "${DATA_DIR}/${WORKFLOW_NAME}_deltae_results.csv" \
  --figure "${DATA_DIR}/${WORKFLOW_NAME}_deltae_vs_wagih.png" \
  --summary "${DATA_DIR}/${WORKFLOW_NAME}_deltae_summary.txt"

echo "Copying compact results back to project directory"
mkdir -p "${PROJECT_DIR}/${DATA_DIR}/euler_results_${SLURM_JOB_ID}"
rsync -a \
  "${DATA_DIR}/${WORKFLOW_NAME}_deltae_sites.csv" \
  "${DATA_DIR}/${WORKFLOW_NAME}_deltae_results.csv" \
  "${DATA_DIR}/${WORKFLOW_NAME}_deltae_summary.txt" \
  "${DATA_DIR}/${WORKFLOW_NAME}_deltae_vs_wagih.png" \
  "${PROJECT_DIR}/${DATA_DIR}/euler_results_${SLURM_JOB_ID}/"

echo "Done."
echo "Full scratch run remains at: ${RUN_DIR}"
echo "Compact results copied to: ${PROJECT_DIR}/${DATA_DIR}/euler_results_${SLURM_JOB_ID}"
