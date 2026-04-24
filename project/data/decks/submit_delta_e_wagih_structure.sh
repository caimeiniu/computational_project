#!/bin/bash
# Controlled experiment — run OUR Phase 3 pipeline on Wagih's OWN annealed
# structure, to separate "structure-generation variance" from "pipeline
# variance" in the 4.8 kJ/mol mean shift seen between our 500-site N=500
# production fit and Wagih's 82k Zenodo dataset.
#
# Uses Wagih's annealed dump0 (converted to a LAMMPS data file) and Wagih's
# own 82,646 GB site IDs (as our mask). Samples 500 of those sites with
# seed=42 and computes ΔE via our single-substitution + CG protocol.
#
# Expected outcome: if our ΔE distribution on Wagih's structure matches
# Wagih's 82k distribution on these same sites, our pipeline is fine and
# the residual is pure structure-realization variance.

#SBATCH --job-name=delta_e_wagih_struct_n500
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=/cluster/scratch/cainiu/wagih_pipeline_test/%x-%j.out
#SBATCH --error=/cluster/scratch/cainiu/wagih_pipeline_test/%x-%j.err

set -euo pipefail
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4
source /cluster/home/cainiu/miniconda3/etc/profile.d/conda.sh
conda activate myenv

RUN_DIR=/cluster/scratch/cainiu/wagih_pipeline_test
ANNEALED=$RUN_DIR/wagih_Al_200A.lmp
GB_MASK=$RUN_DIR/wagih_gb_mask.npy
POTENTIAL=/cluster/home/cainiu/Computational_modeling/project/data/potentials/Al-Mg.eam.fs
DRIVER=/cluster/home/cainiu/Computational_modeling/project/scripts/sample_delta_e.py

N_GB=500
N_BULK=10
SEED=42     # same as production run on our own structure
ETOL=1.0e-8
FTOL=1.0e-10

cd "$RUN_DIR"
echo "Job ID $SLURM_JOB_ID started $(date)"
echo "Structure : $ANNEALED (Wagih Zenodo, 483,425 atoms)"
echo "Mask      : $GB_MASK (82,646 GB atoms from seg_energies_Al_Mg.txt)"
echo "Sampling  : N_GB=$N_GB N_BULK=$N_BULK seed=$SEED"

python "$DRIVER" \
    --annealed "$ANNEALED" \
    --gb-mask "$GB_MASK" \
    --potential "$POTENTIAL" \
    --n-gb "$N_GB" --n-bulk "$N_BULK" --seed "$SEED" \
    --etol "$ETOL" --ftol "$FTOL" \
    --lmp lmp --mpi-ranks "$SLURM_NTASKS" --mpi-cmd mpirun \
    --work-dir "$RUN_DIR/delta_e_run" \
    --out-npz  "$RUN_DIR/delta_e_results_wagih_n500.npz" \
    --out-json "$RUN_DIR/delta_e_meta_wagih_n500.json"

echo "Finished $(date)"
