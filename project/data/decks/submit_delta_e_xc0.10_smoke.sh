#!/bin/bash
# Smoke test for sample_delta_e_finite_xc.py — 50 GB sites + 5 bulk refs on
# the X_c=0.10 preseg final config. Validates: LAMMPS deck builds, type
# parser returns correct types, single-Mg insertion CG-relaxes cleanly,
# spectrum has reasonable shape.
#
# Substrate carries pre-existing Mg (X_c=0.10, X_gb=0.375) — this is the
# whole point: testing whether ΔE_eff(i; X_c=0.10) shifts vs the X_c=0
# reference spectrum (n=500, mean ≈ −6.9 kJ/mol).
#
# Timing basis: 200A pure-Al ran 36 s/site on 16 cores → 55 sites ≈ 33 min;
# 1 h budget for variance + queue + finite-X_c overhead. Production sweep
# (3 X_c × 500 sites) will follow once smoke is green.
#
# Usage:
#   sbatch project/data/decks/submit_delta_e_xc0.10_smoke.sh

#SBATCH --job-name=delta_e_xc0.10_smoke
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/delta_e_xc0.10_smoke/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/delta_e_xc0.10_smoke/%x-%j.err

set -euo pipefail

# ----- environment -----
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4
source /cluster/home/cainiu/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# ----- job parameters -----
PROJECT=/cluster/home/cainiu/Computational_modeling/project
RUN_DIR=/cluster/scratch/cainiu/delta_e_xc0.10_smoke
mkdir -p "$RUN_DIR"

BACKGROUND=$PROJECT/data/snapshots/hmc_T500_Xc0.10_preseg_final.lmp
GB_MASK=$PROJECT/data/snapshots/gb_mask_200A.npy
POTENTIAL=$PROJECT/data/potentials/Al-Mg.eam.fs
DRIVER=$PROJECT/scripts/sample_delta_e_finite_xc.py

N_GB=50
N_BULK=5
SEED=42                # match the n=500 reference seed for cross-comparison

# ----- bookkeeping -----
cd "$RUN_DIR"
echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Node           : $SLURM_JOB_NODELIST"
echo "MPI tasks      : $SLURM_NTASKS"
echo "Run dir        : $(pwd)"
echo "Background     : $BACKGROUND  (~476k atoms, X_c=0.10 preseg final)"
echo "GB mask        : $GB_MASK"
echo "Driver         : $DRIVER"
echo "N_GB, N_BULK   : $N_GB, $N_BULK   (seed=$SEED)"
echo "Started        : $(date)"
echo "=========================================================="

# ----- run -----
python "$DRIVER" \
    --background-data "$BACKGROUND" \
    --gb-mask "$GB_MASK" \
    --potential "$POTENTIAL" \
    --n-gb "$N_GB" --n-bulk "$N_BULK" --seed "$SEED" \
    --xc-label 0.10 \
    --lmp lmp --mpi-ranks "$SLURM_NTASKS" --mpi-cmd mpirun \
    --work-dir "$RUN_DIR/delta_e_run" \
    --out-npz  "$RUN_DIR/delta_e_results_xc0.10_smoke.npz" \
    --out-json "$RUN_DIR/delta_e_meta_xc0.10_smoke.json"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
