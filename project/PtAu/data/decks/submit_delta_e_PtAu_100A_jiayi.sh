#!/bin/bash
# Jiayi-safe Pt(Au) per-site DeltaE submit script.
# Assumes submit_anneal_PtAu_100A_jiayi.sh already produced:
#   /cluster/scratch/jiayfu/prototype_PtAu_100A/poly_Pt_100A_8g_annealed.lmp

#SBATCH --job-name=delta_e_PtAu_100A_n500_tight
#SBATCH --time=04:00:00
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
ANNEALED=$RUN_DIR/poly_Pt_100A_8g_annealed.lmp
GB_MASK=$RUN_DIR/gb_mask_PtAu_100A.npy
POTENTIAL=$PROJECT/PtAu/data/potentials/PtAu.eam.alloy
DRIVER=$PROJECT/PtAu/scripts/sample_delta_e_PtAu.py

N_GB=500
N_BULK=10
SEED=42

if [[ ! -f "$ANNEALED" ]]; then
  echo "ERROR: annealed structure missing: $ANNEALED" >&2
  echo "Run submit_anneal_PtAu_100A_jiayi.sh first." >&2
  exit 1
fi

cd "$RUN_DIR"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "Run dir        : $RUN_DIR"
echo "Annealed file  : $ANNEALED"
echo "GB mask        : $GB_MASK"
echo "Potential      : $POTENTIAL"
echo "N_GB, N_BULK   : $N_GB, $N_BULK"
echo "Started        : $(date)"
echo "=========================================================="

if [[ ! -f "$GB_MASK" ]]; then
  echo "gb_mask missing; running gb_identify.py first"
  python3 "$PROJECT/scripts/gb_identify.py" \
    "$ANNEALED" \
    --parent fcc \
    --lattice-a 3.9764 \
    --out-mask "$GB_MASK" \
    --out-report "$RUN_DIR/gb_info_PtAu_100A.json" \
    --out-dump "$RUN_DIR/gb_cna_PtAu_100A.dump"
fi

python3 "$DRIVER" \
  --annealed "$ANNEALED" \
  --gb-mask "$GB_MASK" \
  --potential "$POTENTIAL" \
  --n-gb "$N_GB" \
  --n-bulk "$N_BULK" \
  --seed "$SEED" \
  --elements "Pt Au" \
  --masses "195.0900 196.9665" \
  --lmp lmp \
  --mpi-ranks "$SLURM_NTASKS" \
  --mpi-cmd mpirun \
  --work-dir "$RUN_DIR/delta_e_run_n500_PtAu_100A_tight" \
  --out-npz "$RUN_DIR/delta_e_results_n500_PtAu_100A_tight.npz" \
  --out-json "$RUN_DIR/delta_e_meta_n500_PtAu_100A_tight.json"

mkdir -p "$PROJECT/PtAu/output"
cp "$RUN_DIR"/delta_e_results_n500_PtAu_100A_tight.npz \
   "$RUN_DIR"/delta_e_meta_n500_PtAu_100A_tight.json \
   "$PROJECT/PtAu/output/"

echo "=========================================================="
echo "Finished       : $(date)"
echo "Copied compact outputs to $PROJECT/PtAu/output"
echo "=========================================================="
