#!/bin/bash
# Submit Phase 3 per-site ΔE_seg sampling on Euler.
#
# Assumes anneal + gb_identify have already run and produced:
#   <RUN_DIR>/poly_Al_100A_8g_annealed.lmp
#   <RUN_DIR>/gb_mask.npy
#
# Usage:
#   sbatch /cluster/home/cainiu/Computational_modeling/project/data/decks/submit_delta_e.sh
#
# To rescale, edit N_GB / N_BULK / RUN_DIR / SEED below.

#SBATCH --job-name=delta_e_AlMg_500
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/prototype_AlMg_100A/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/prototype_AlMg_100A/%x-%j.err

set -euo pipefail

# ----- environment -----
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4
# Python env with numpy + scipy for the driver script.
source /cluster/home/cainiu/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# ----- job parameters -----
RUN_DIR=/cluster/scratch/cainiu/prototype_AlMg_100A
ANNEALED=$RUN_DIR/poly_Al_100A_8g_annealed.lmp
GB_MASK=$RUN_DIR/gb_mask.npy
POTENTIAL=/cluster/home/cainiu/Computational_modeling/project/data/potentials/Al-Mg.eam.fs
DRIVER=/cluster/home/cainiu/Computational_modeling/project/scripts/sample_delta_e.py

N_GB=500
N_BULK=10
SEED=42
# CG tolerances: tighter than anneal (ΔE accuracy matters at sub-meV level)
ETOL=1.0e-8
FTOL=1.0e-10

# ----- bookkeeping -----
cd "$RUN_DIR"
echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Node           : $SLURM_JOB_NODELIST"
echo "MPI tasks      : $SLURM_NTASKS"
echo "Run dir        : $(pwd)"
echo "Annealed file  : $ANNEALED"
echo "GB mask        : $GB_MASK"
echo "N_GB, N_BULK   : $N_GB, $N_BULK   (seed=$SEED)"
echo "Started        : $(date)"
echo "=========================================================="

# ----- sanity: did gb_identify run already? -----
if [[ ! -f "$GB_MASK" ]]; then
    echo "gb_mask.npy missing — running gb_identify.py first" >&2
    python /cluster/home/cainiu/Computational_modeling/project/scripts/gb_identify.py \
        "$ANNEALED" \
        --lattice-a 4.05 \
        --out-mask "$GB_MASK" \
        --out-report "$RUN_DIR/gb_info.json" \
        --out-dump "$RUN_DIR/gb_cna.dump"
fi

# ----- run -----
# Why --mpi-cmd mpirun (not srun) inside sbatch: the driver fires one LAMMPS
# per site, i.e. nested srun inside the outer SLURM allocation, which on
# recent SLURM can trigger "job step creation temporarily disabled".
# openmpi/4.1.6 mpirun picks up the SLURM allocation via PMI and avoids
# nested-step issues.
python "$DRIVER" \
    --annealed "$ANNEALED" \
    --gb-mask "$GB_MASK" \
    --potential "$POTENTIAL" \
    --n-gb "$N_GB" --n-bulk "$N_BULK" --seed "$SEED" \
    --etol "$ETOL" --ftol "$FTOL" \
    --lmp lmp --mpi-ranks "$SLURM_NTASKS" --mpi-cmd mpirun \
    --work-dir "$RUN_DIR/delta_e_run_n500" \
    --out-npz  "$RUN_DIR/delta_e_results_n500.npz" \
    --out-json "$RUN_DIR/delta_e_meta_n500.json"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
