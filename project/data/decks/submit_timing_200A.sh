#!/bin/bash
#SBATCH --job-name=timing_AlMg_200A
#SBATCH --output=/cluster/scratch/cainiu/production_AlMg_200A/timing_200A-%j.out
#SBATCH --error=/cluster/scratch/cainiu/production_AlMg_200A/timing_200A-%j.err
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:30:00

set -e
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

cd /cluster/scratch/cainiu/production_AlMg_200A

mpirun lmp -in /cluster/home/cainiu/Computational_modeling/project/data/decks/timing_200A.lammps \
           -log timing_200A.log
