#!/bin/bash
# Submit Phase 3 per-site ΔE_seg sampling on the Pt-Au prototype polycrystal.
#
# Assumes anneal + gb_identify have already produced:
#   <RUN_DIR>/poly_Pt_100A_8g_annealed.lmp
#   <RUN_DIR>/gb_mask_PtAu_100A.npy
#
# Usage:
#   sbatch /cluster/home/cainiu/Computational_modeling/project/PtAu/data/decks/submit_delta_e_PtAu_100A.sh
#
# Timing estimate: ~15 min on 16 cores for n=500 (tight CG) at ~64k atoms;
# 4 h budget for variance / first-run unknowns.

#SBATCH --job-name=delta_e_PtAu_100A_n500_tight
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/prototype_PtAu_100A/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/prototype_PtAu_100A/%x-%j.err

set -euo pipefail

# ----- environment -----
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4
source /cluster/home/cainiu/miniconda3/etc/profile.d/conda.sh
conda activate myenv

# ----- job parameters -----
RUN_DIR=/cluster/scratch/cainiu/prototype_PtAu_100A
ANNEALED=$RUN_DIR/poly_Pt_100A_8g_annealed.lmp
GB_MASK=$RUN_DIR/gb_mask_PtAu_100A.npy
POTENTIAL=/cluster/home/cainiu/Computational_modeling/project/PtAu/data/potentials/PtAu.eam.alloy
DRIVER=/cluster/home/cainiu/Computational_modeling/project/PtAu/scripts/sample_delta_e_PtAu.py

N_GB=500
N_BULK=10
SEED=42                # same as Al-Mg prototype for reproducibility across alloys
# CG tolerances inherit driver defaults (1e-25/1e-25/50000/5000000) — match
# Wagih's calculate_E_GB_solute.in.

# ----- bookkeeping -----
cd "$RUN_DIR"
echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Node           : $SLURM_JOB_NODELIST"
echo "MPI tasks      : $SLURM_NTASKS"
echo "Run dir        : $(pwd)"
echo "Annealed file  : $ANNEALED (~64k atoms, 10³ nm³, 8 grains, Pt host)"
echo "GB mask        : $GB_MASK"
echo "Potential      : $POTENTIAL"
echo "N_GB, N_BULK   : $N_GB, $N_BULK   (seed=$SEED)"
echo "Started        : $(date)"
echo "=========================================================="

# ----- sanity: did gb_identify run already? -----
if [[ ! -f "$GB_MASK" ]]; then
    echo "gb_mask missing — running gb_identify.py first" >&2
    python /cluster/home/cainiu/Computational_modeling/project/scripts/gb_identify.py \
        "$ANNEALED" \
        --parent fcc --lattice-a 3.9764 \
        --out-mask "$GB_MASK" \
        --out-report "$RUN_DIR/gb_info_PtAu_100A.json" \
        --out-dump "$RUN_DIR/gb_cna_PtAu_100A.dump"
fi

# ----- run -----
# mpirun (not srun) inside sbatch to avoid nested-step issues with openmpi/4.1.6.
python "$DRIVER" \
    --annealed "$ANNEALED" \
    --gb-mask "$GB_MASK" \
    --potential "$POTENTIAL" \
    --n-gb "$N_GB" --n-bulk "$N_BULK" --seed "$SEED" \
    --elements "Pt Au" --masses "195.0900 196.9665" \
    --lmp lmp --mpi-ranks "$SLURM_NTASKS" --mpi-cmd mpirun \
    --work-dir "$RUN_DIR/delta_e_run_n500_PtAu_100A_tight" \
    --out-npz  "$RUN_DIR/delta_e_results_n500_PtAu_100A_tight.npz" \
    --out-json "$RUN_DIR/delta_e_meta_n500_PtAu_100A_tight.json"

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
