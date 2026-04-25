#!/usr/bin/env bash
#SBATCH --job-name=cuni3d_deltae
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --output=cuni3d_%j.out

set -euo pipefail

module load stack/2024-06
module load openmpi/4.1.6
module load lammps/20230802

PROJECT_DIR="$(pwd)"
RUN_DIR="/cluster/scratch/${USER}/cuni3d_${SLURM_JOB_ID}"
LMP_CMD="mpirun -np ${SLURM_NTASKS} lmp"

echo "Project: ${PROJECT_DIR}"
echo "Scratch run directory: ${RUN_DIR}"
echo "LAMMPS command: ${LMP_CMD}"

mkdir -p "${RUN_DIR}"

# Keep the run self-contained in scratch, but avoid copying bulky/generated outputs
# from previous attempts.
rsync -a \
  --exclude ".git" \
  --exclude "__pycache__" \
  --exclude "data/cuni_3d/lammps/deltae_jobs/log_*.lammps" \
  --exclude "data/cuni_3d/lammps/deltae_jobs/relaxed_*.lammps" \
  --exclude "data/cuni_3d/lammps/log_*.lammps" \
  --exclude "data/cuni_3d/lammps/cuni_3d_*.lammps" \
  "${PROJECT_DIR}/" "${RUN_DIR}/"

cd "${RUN_DIR}"

echo "Step 1/5: generate 3D Cu polycrystal"
python3 scripts/generate_cuni_3d_polycrystal.py --box-nm 10 --grains 8

echo "Step 2/5: anneal/minimize 3D polycrystal"
cd data/cuni_3d/lammps
${LMP_CMD} -in in_relax_cuni_3d.lammps
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
mkdir -p "${PROJECT_DIR}/data/cuni_3d/euler_results_${SLURM_JOB_ID}"
rsync -a \
  data/cuni_3d/cuni_3d_deltae_sites.csv \
  data/cuni_3d/cuni_3d_deltae_results.csv \
  data/cuni_3d/cuni_3d_deltae_summary.txt \
  data/cuni_3d/cuni_3d_deltae_vs_wagih.png \
  data/cuni_3d/lammps/log_relax_cuni_3d.lammps \
  "${PROJECT_DIR}/data/cuni_3d/euler_results_${SLURM_JOB_ID}/"

echo "Done."
echo "Full scratch run remains at: ${RUN_DIR}"
echo "Compact results copied to: ${PROJECT_DIR}/data/cuni_3d/euler_results_${SLURM_JOB_ID}"
