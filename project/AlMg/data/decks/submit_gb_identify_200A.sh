#!/bin/bash
#SBATCH --job-name=gb_identify_AlMg_200A
#SBATCH --output=/cluster/scratch/cainiu/production_AlMg_200A/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/production_AlMg_200A/%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00

set -euo pipefail
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

RUN_DIR=/cluster/scratch/cainiu/production_AlMg_200A
cd "$RUN_DIR"

python3 /cluster/home/cainiu/Computational_modeling/project/scripts/gb_identify.py \
  "$RUN_DIR/poly_Al_200A_16g_annealed.lmp" \
  --parent fcc --lattice-a 4.05 \
  --lmp lmp \
  --out-mask   "$RUN_DIR/gb_mask_200A.npy" \
  --out-report "$RUN_DIR/gb_info_200A.json" \
  --out-dump   "$RUN_DIR/gb_cna_200A.dump"
