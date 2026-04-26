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
RUN_DIR="/cluster/scratch/${USER}/cuni3d_${SLURM_JOB_ID}"
LMP_CMD="mpirun -np ${SLURM_NTASKS} lmp"
SOLVENT="${SOLVENT:-Cu}"
SOLUTE="${SOLUTE:-Ni}"
BOX_NM="${BOX_NM:-10}"
GRAINS="${GRAINS:-8}"
SAMPLES="${SAMPLES:-500}"
case "${SOLVENT}" in
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
  *)
    echo "Unsupported SOLVENT=${SOLVENT}; expected Cu or Ni"
    exit 2
    ;;
esac
case "${SOLUTE}" in
  Cu) SOLUTE_MASS="63.546" ;;
  Ni) SOLUTE_MASS="58.693" ;;
  *)
    echo "Unsupported SOLUTE=${SOLUTE}; expected Cu or Ni"
    exit 2
    ;;
esac
if [ "${SOLVENT}" = "${SOLUTE}" ]; then
  echo "SOLVENT and SOLUTE must be different"
  exit 2
fi

echo "Project: ${PROJECT_DIR}"
echo "Scratch run directory: ${RUN_DIR}"
echo "LAMMPS command: ${LMP_CMD}"
echo "Alloy: ${SOLVENT}(${SOLUTE})"
echo "Box/grains/samples/lattice: ${BOX_NM} nm / ${GRAINS} / ${SAMPLES} / ${LATTICE} A"
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

echo "Step 1/5: generate 3D ${SOLVENT} polycrystal"
python3 scripts/generate_cuni_3d_polycrystal.py \
  --box-nm "${BOX_NM}" \
  --grains "${GRAINS}" \
  --solvent "${SOLVENT}" \
  --solute "${SOLUTE}" \
  --lattice "${LATTICE}"

cat > data/cuni_3d/lammps/in_relax_cuni_3d.lammps <<EOF
log log_relax_cuni_3d.lammps

units metal
dimension 3
boundary p p p
atom_style atomic

read_data ../cuni_3d_polycrystal.lammps

mass 1 ${SOLVENT_MASS}
mass 2 ${SOLUTE_MASS}

pair_style eam/alloy
pair_coeff * * ../../examples/nc_swap_CuNi/Cu_Ni_Fischer_2018.eam.alloy ${SOLVENT} ${SOLUTE}

neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes

timestep 0.001
thermo 1000
thermo_style custom step time temp pe ke press vol

min_style cg
minimize 1.0e-10 1.0e-10 10000 100000
write_data cuni_3d_initial_minimized.lammps

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
write_data cuni_3d_relaxed.lammps
EOF

echo "Step 2/5: anneal/minimize 3D polycrystal"
cd data/cuni_3d/lammps
${LMP_CMD} -in in_relax_cuni_3d.lammps
cd "${RUN_DIR}"

echo "Step 3/5: select 500 GB sites from relaxed structure"
python3 scripts/select_cuni_deltae_sites.py \
  --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps \
  --samples "${SAMPLES}" \
  --lattice "${LATTICE}"

echo "Step 4/5: write and run Cu->Ni Delta E minimization jobs"
python3 scripts/write_cuni_deltae_lammps_jobs.py \
  --data data/cuni_3d/lammps/cuni_3d_relaxed.lammps \
  --lammps-command "${LMP_CMD}" \
  --solvent "${SOLVENT}" \
  --solute "${SOLUTE}"

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
