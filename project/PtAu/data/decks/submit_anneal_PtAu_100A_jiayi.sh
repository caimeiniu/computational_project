#!/bin/bash
# Jiayi-safe Pt(Au) anneal submit script.
# Generates the 100 A / 8-grain Pt polycrystal if missing, then runs the
# Wagih-style Pt(Au) anneal using Cainiu's PtAu deck.

#SBATCH --job-name=anneal_PtAu_100A
#SBATCH --time=06:00:00
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
DECK=$PROJECT/PtAu/data/decks/anneal_PtAu.lammps
POTENTIAL=$PROJECT/PtAu/data/potentials/PtAu.eam.alloy
DATAFILE=poly_Pt_100A_8g.lmp
OUTSTUB=poly_Pt_100A_8g
T_HOLD=816.0
T_HOLD_PS=250
COOL_RATE=3.0
VELOCITY_SEED=20240101

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

if [[ ! -f "$DATAFILE" ]]; then
  echo "Initial polycrystal missing; generating $RUN_DIR/$DATAFILE"
  python3 "$PROJECT/scripts/generate_polycrystal.py" \
    --structure fcc \
    --box 100 \
    --grains 8 \
    --lattice-a 3.9764 \
    --structure-seed 1 \
    --types 2 \
    --out "$DATAFILE"
fi

echo "=========================================================="
echo "Job ID       : $SLURM_JOB_ID"
echo "Run dir      : $RUN_DIR"
echo "Deck         : $DECK"
echo "Datafile     : $DATAFILE"
echo "Potential    : $POTENTIAL"
echo "T_hold       : $T_HOLD K"
echo "Started      : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var datafile      "$DATAFILE" \
  -var outstub       "$OUTSTUB" \
  -var T_HOLD        "$T_HOLD" \
  -var T_HOLD_PS     "$T_HOLD_PS" \
  -var COOL_RATE     "$COOL_RATE" \
  -var VELOCITY_SEED "$VELOCITY_SEED" \
  -var POTFILE       "$POTENTIAL" \
  -var EL1           Pt \
  -var EL2           Au \
  -var MASS1         195.0900 \
  -var MASS2         196.9665

echo "=========================================================="
echo "Finished     : $(date)"
echo "=========================================================="
