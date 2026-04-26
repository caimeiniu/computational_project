#!/usr/bin/env bash
#SBATCH --job-name=cuni3d_resume
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --output=cuni3d_resume_%j.out

set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "Usage: sbatch scripts/euler_cuni3d_resume_after_relax.sbatch /cluster/scratch/USER/cuni3d_OLDJOBID"
  exit 2
fi

module load stack/2024-06
module load python/3.11.6
module load openmpi/4.1.6
module load lammps/20230802

PROJECT_DIR="$(pwd)"
RUN_DIR="$1"
LMP_CMD="mpirun -np ${SLURM_NTASKS} lmp"

echo "Project: ${PROJECT_DIR}"
echo "Resuming scratch run directory: ${RUN_DIR}"
echo "LAMMPS command: ${LMP_CMD}"

test -f "${RUN_DIR}/data/cuni_3d/lammps/cuni_3d_relaxed.lammps"

# Copy fixed Python scripts into the old scratch run without touching LAMMPS outputs.
rsync -a "${PROJECT_DIR}/scripts/" "${RUN_DIR}/scripts/"

cd "${RUN_DIR}"

echo "Step 3/5: select 500 GB sites from relaxed structure"
python3 scripts/select_cuni_deltae_sites.py \
  --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps \
  --samples 500

echo "Step 4/5: write and run Cu->Ni Delta E minimization jobs"
python3 scripts/write_cuni_deltae_lammps_jobs.py \
  --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps \
  --lammps-command "${LMP_CMD}"

cd data/cuni_3d/lammps/deltae_jobs
./run_all_deltae.sh
cd "${RUN_DIR}"

echo "Step 5/5: analyze Delta E spectrum and compare with Wagih Ni(Cu)"
python3 scripts/analyze_cuni_deltae.py

echo "Copying compact results back to project directory"
mkdir -p "${PROJECT_DIR}/data/cuni_3d/euler_results_resume_${SLURM_JOB_ID}"
rsync -a \
  data/cuni_3d/cuni_3d_deltae_sites.csv \
  data/cuni_3d/cuni_3d_deltae_results.csv \
  data/cuni_3d/cuni_3d_deltae_summary.txt \
  data/cuni_3d/cuni_3d_deltae_vs_wagih.png \
  data/cuni_3d/lammps/log_relax_cuni_3d.lammps \
  "${PROJECT_DIR}/data/cuni_3d/euler_results_resume_${SLURM_JOB_ID}/"

echo "Done."
echo "Full scratch run remains at: ${RUN_DIR}"
echo "Compact results copied to: ${PROJECT_DIR}/data/cuni_3d/euler_results_resume_${SLURM_JOB_ID}"
