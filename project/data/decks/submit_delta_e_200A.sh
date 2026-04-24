#!/bin/bash
# Submit Phase 3 per-site ΔE_seg sampling on the 20³ nm³ production polycrystal.
#
# Assumes anneal + gb_identify have already run and produced:
#   <RUN_DIR>/poly_Al_200A_16g_annealed.lmp
#   <RUN_DIR>/gb_mask_200A.npy
#
# Usage:
#   sbatch /cluster/home/cainiu/Computational_modeling/project/data/decks/submit_delta_e_200A.sh
#
# Timing basis: prototype (59k atoms) was 3.1 s/site × 16 cores × 510 sites =
# 26.6 min. Production is ~8× atoms; CG cost scales ~linearly → estimate
# ~25 s/site × 510 sites ≈ 3.5 h. Requesting 6 h to absorb any long-tail CG.

#SBATCH --job-name=delta_e_AlMg_200A_n500
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/production_AlMg_200A/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/production_AlMg_200A/%x-%j.err

set -euo pipefail

# ----- environment -----
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4
source /cluster/home/cainiu/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# ----- job parameters -----
RUN_DIR=/cluster/scratch/cainiu/production_AlMg_200A
ANNEALED=$RUN_DIR/poly_Al_200A_16g_annealed.lmp
GB_MASK=$RUN_DIR/gb_mask_200A.npy
POTENTIAL=/cluster/home/cainiu/Computational_modeling/project/data/potentials/Al-Mg.eam.fs
DRIVER=/cluster/home/cainiu/Computational_modeling/project/scripts/sample_delta_e.py

N_GB=500
N_BULK=10
SEED=42                # same as prototype for reproducibility across scans
ETOL=1.0e-8
FTOL=1.0e-10

# ----- bookkeeping -----
cd "$RUN_DIR"
echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Node           : $SLURM_JOB_NODELIST"
echo "MPI tasks      : $SLURM_NTASKS"
echo "Run dir        : $(pwd)"
echo "Annealed file  : $ANNEALED (~476k atoms, 20³ nm³, 16 grains)"
echo "GB mask        : $GB_MASK"
echo "N_GB, N_BULK   : $N_GB, $N_BULK   (seed=$SEED)"
echo "Started        : $(date)"
echo "=========================================================="

# ----- run -----
python "$DRIVER" \
    --annealed "$ANNEALED" \
    --gb-mask "$GB_MASK" \
    --potential "$POTENTIAL" \
    --n-gb "$N_GB" --n-bulk "$N_BULK" --seed "$SEED" \
    --etol "$ETOL" --ftol "$FTOL" \
    --lmp lmp --mpi-ranks "$SLURM_NTASKS" --mpi-cmd mpirun \
    --work-dir "$RUN_DIR/delta_e_run_n500_200A" \
    --out-npz  "$RUN_DIR/delta_e_results_n500_200A.npz" \
    --out-json "$RUN_DIR/delta_e_meta_n500_200A.json"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
