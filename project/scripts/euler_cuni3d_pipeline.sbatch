#!/usr/bin/env bash
#SBATCH --job-name=cuni3d_deltae
#SBATCH --ntasks=16
#SBATCH --time=24:00:00
#SBATCH --output=cuni3d_%j.out

set -euo pipefail

module load stack/2024-06
module load python/3.11.6
module load openmpi/4.1.6
module load lammps/20230802

PROJECT_DIR="$(pwd)"
LMP_CMD="mpirun -np ${SLURM_NTASKS} lmp"
SOLVENT="${SOLVENT:-Cu}"
SOLUTE="${SOLUTE:-Ni}"
WORKFLOW_NAME="${WORKFLOW_NAME:-$(printf '%s%s_3d' "${SOLVENT}" "${SOLUTE}" | tr '[:upper:]' '[:lower:]')}"
DATA_DIR="data/${WORKFLOW_NAME}"
RUN_DIR="/cluster/scratch/${USER}/${WORKFLOW_NAME}_${SLURM_JOB_ID}"
POTENTIAL="${POTENTIAL:-data/examples/nc_swap_CuNi/Cu_Ni_Fischer_2018.eam.alloy}"
BOX_NM="${BOX_NM:-10}"
GRAINS="${GRAINS:-8}"
SAMPLES="${SAMPLES:-500}"
BULK_REFERENCES="${BULK_REFERENCES:-20}"
case "${SOLVENT}" in
  Au)
    LATTICE="${LATTICE:-4.0782}"
    SOLVENT_MASS="196.96657"
    ANNEAL_TEMP="${ANNEAL_TEMP:-650.0}"
    ;;
  Cu)
    LATTICE="${LATTICE:-3.615}"
    SOLVENT_MASS="63.546"
    ANNEAL_TEMP="${ANNEAL_TEMP:-650.0}"
    ;;
  Ni)
    LATTICE="${LATTICE:-3.520}"
    SOLVENT_MASS="58.693"
    ANNEAL_TEMP="${ANNEAL_TEMP:-800.0}"
    ;;
  Pt)
    LATTICE="${LATTICE:-3.9231}"
    SOLVENT_MASS="195.084"
    ANNEAL_TEMP="${ANNEAL_TEMP:-700.0}"
    ;;
  *)
    echo "Unsupported SOLVENT=${SOLVENT}; expected Au, Cu, Ni, or Pt"
    exit 2
    ;;
esac
case "${SOLUTE}" in
  Au) SOLUTE_MASS="196.96657" ;;
  Cu) SOLUTE_MASS="63.546" ;;
  Ni) SOLUTE_MASS="58.693" ;;
  Pt) SOLUTE_MASS="195.084" ;;
  *)
    echo "Unsupported SOLUTE=${SOLUTE}; expected Au, Cu, Ni, or Pt"
    exit 2
    ;;
esac
if [ "${SOLVENT}" = "${SOLUTE}" ]; then
  echo "SOLVENT and SOLUTE must be different"
  exit 2
fi

echo "Project: ${PROJECT_DIR}"
echo "Scratch run directory: ${RUN_DIR}"
echo "Data directory: ${DATA_DIR}"
echo "LAMMPS command: ${LMP_CMD}"
echo "Alloy: ${SOLVENT}(${SOLUTE})"
echo "Box/grains/samples/bulk refs/lattice: ${BOX_NM} nm / ${GRAINS} / ${SAMPLES} / ${BULK_REFERENCES} / ${LATTICE} A"
echo "Potential: ${POTENTIAL}"
echo "Python: $(which python3)"
python3 - <<'PY'
import numpy
import scipy
import matplotlib
print(f"NumPy: {numpy.__version__}")
print(f"SciPy: {scipy.__version__}")
print(f"Matplotlib: {matplotlib.__version__}")
PY

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
mkdir -p "${DATA_DIR}/lammps"

echo "Step 1/5: generate 3D ${SOLVENT} polycrystal"
python3 scripts/generate_cuni_3d_polycrystal.py \
  --box-nm "${BOX_NM}" \
  --grains "${GRAINS}" \
  --solvent "${SOLVENT}" \
  --solute "${SOLUTE}" \
  --lattice "${LATTICE}" \
  --output "${DATA_DIR}/${WORKFLOW_NAME}_polycrystal.lammps" \
  --metadata "${DATA_DIR}/${WORKFLOW_NAME}_polycrystal_sites.csv"

cat > "${DATA_DIR}/lammps/in_relax_${WORKFLOW_NAME}.lammps" <<EOF
log log_relax_${WORKFLOW_NAME}.lammps

units metal
dimension 3
boundary p p p
atom_style atomic

read_data ../${WORKFLOW_NAME}_polycrystal.lammps

mass 1 ${SOLVENT_MASS}
mass 2 ${SOLUTE_MASS}

pair_style eam/alloy
pair_coeff * * ../../../${POTENTIAL} ${SOLVENT} ${SOLUTE}

neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes

timestep 0.001
thermo 1000
thermo_style custom step time temp pe ke press vol

min_style cg
minimize 1.0e-10 1.0e-10 10000 100000
write_data ${WORKFLOW_NAME}_initial_minimized.lammps

variable anneal_temp equal ${ANNEAL_TEMP}
variable final_temp equal 1.0
variable anneal_steps equal 250000
variable cool_steps equal ceil((\${anneal_temp}-\${final_temp})/3.0/0.001)

velocity all create \${anneal_temp} 20260425 mom yes rot yes dist gaussian

fix anneal all npt temp \${anneal_temp} \${anneal_temp} 0.1 iso 0.0 0.0 1.0
run \${anneal_steps}
unfix anneal

fix cool all npt temp \${anneal_temp} \${final_temp} 0.1 iso 0.0 0.0 1.0
run \${cool_steps}
unfix cool

minimize 1.0e-12 1.0e-12 10000 100000
write_data ${WORKFLOW_NAME}_relaxed.lammps
EOF

echo "Step 2/5: anneal/minimize 3D polycrystal"
cd "${DATA_DIR}/lammps"
${LMP_CMD} -in "in_relax_${WORKFLOW_NAME}.lammps"
cd "${RUN_DIR}"

echo "Step 3/5: select 500 GB sites from relaxed structure"
python3 scripts/select_cuni_deltae_sites.py \
  --data "${DATA_DIR}/lammps/${WORKFLOW_NAME}_relaxed.lammps" \
  --samples "${SAMPLES}" \
  --bulk-references "${BULK_REFERENCES}" \
  --lattice "${LATTICE}" \
  --output "${DATA_DIR}/${WORKFLOW_NAME}_deltae_sites.csv"

echo "Step 4/5: write and run Cu->Ni Delta E minimization jobs"
python3 scripts/write_cuni_deltae_lammps_jobs.py \
  --sites "${DATA_DIR}/${WORKFLOW_NAME}_deltae_sites.csv" \
  --data "${DATA_DIR}/lammps/${WORKFLOW_NAME}_relaxed.lammps" \
  --potential "${POTENTIAL}" \
  --jobs-dir "${DATA_DIR}/lammps/deltae_jobs" \
  --lammps-command "${LMP_CMD}" \
  --solvent "${SOLVENT}" \
  --solute "${SOLUTE}"

cd "${DATA_DIR}/lammps/deltae_jobs"
./run_all_deltae.sh
cd "${RUN_DIR}"

echo "Step 5/5: analyze Delta E spectrum and compare with Wagih Ni(Cu)"
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
  "${DATA_DIR}/lammps/log_relax_${WORKFLOW_NAME}.lammps" \
  "${PROJECT_DIR}/${DATA_DIR}/euler_results_${SLURM_JOB_ID}/"

echo "Done."
echo "Full scratch run remains at: ${RUN_DIR}"
echo "Compact results copied to: ${PROJECT_DIR}/${DATA_DIR}/euler_results_${SLURM_JOB_ID}"
