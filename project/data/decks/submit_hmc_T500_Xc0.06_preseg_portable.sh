#!/bin/bash
# Self-contained HMC job: T=500 K, X_c=0.06 (Al-Mg, preseg IC). Portable variant
# of submit_hmc_T500_Xc0.06_preseg_16rank.sh — no hard-coded paths to /cluster/
# home/cainiu/ or /cluster/scratch/cainiu/, so a teammate can run it from any
# folder on Euler (or any SLURM cluster with the same lammps build) without
# editing the script.
#
# Required files in the SAME folder as this script:
#   submit_hmc_T500_Xc0.06_preseg_portable.sh   (this file)
#   hmc_AlMg.lammps                              (LAMMPS deck, 5 KB)
#   Al-Mg.eam.fs                                 (Mendelev 2009 potential, 2.2 MB)
#   poly_AlMg_200A_preseg_Xc0.06.lmp             (pre-segregated initial config, 63 MB)
#
# How to run:
#   cd /path/to/folder/with/4/files
#   sbatch submit_hmc_T500_Xc0.06_preseg_portable.sh
#
# Wall ~17 h (16-rank); 20 h budget. Disk: ~700 MB peak (dump 600 MB + final 63 MB).
# Outputs land in $SLURM_SUBMIT_DIR (same folder).
#
# Send back to caimei (small files only, ~70 MB):
#   hmc_T500_Xc0.06_preseg.log
#   hmc_T500_Xc0.06_preseg_final.lmp
#   hmc_AlMg_T500_Xc0.06_preseg_16r-<jobid>.{out,err}

#SBATCH --job-name=hmc_AlMg_T500_Xc0.06_preseg_16r
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

# Euler module stack — adjust the `module load` line if running on a different
# cluster. Need: lammps built with MANYBODY (eam/fs) + MC (atom/swap) packages.
module purge
module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 lammps/20240829.4

cd "$SLURM_SUBMIT_DIR"

DECK=hmc_AlMg.lammps
PRESEG=poly_AlMg_200A_preseg_Xc0.06.lmp
POTFILE=Al-Mg.eam.fs

for f in "$DECK" "$PRESEG" "$POTFILE"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: missing $f in $(pwd)" >&2
    exit 1
  fi
done

T=500.0
XC=0.06
EQUIL_PS=10
PROD_PS=300
SWAPS_PER_CALL=100
OUTSTUB="hmc_T500_Xc0.06_preseg"

echo "=========================================================="
echo "Job ID         : $SLURM_JOB_ID"
echo "T              : $T K"
echo "X_c            : $XC"
echo "EQUIL_PS       : $EQUIL_PS ps"
echo "PROD_PS        : $PROD_PS ps"
echo "SWAPS_PER_CALL : $SWAPS_PER_CALL"
echo "MPI ranks      : $SLURM_NTASKS (16-rank variant)"
echo "Initial Mg     : PRE-SEGREGATED (X_GB(0) = 0.3206, all in GB)"
echo "canon-FD pred  : X_GB ≈ 0.2611 (ours, n=500)"
echo "Submit dir     : $SLURM_SUBMIT_DIR"
echo "Started        : $(date)"
echo "=========================================================="

srun -n "$SLURM_NTASKS" lmp \
  -in "$DECK" \
  -var annealed       "$PRESEG" \
  -var POTFILE        "$POTFILE" \
  -var outstub        "$OUTSTUB" \
  -var T              "$T" \
  -var XC             "$XC" \
  -var EQUIL_PS       "$EQUIL_PS" \
  -var PROD_PS        "$PROD_PS" \
  -var SWAPS_PER_CALL "$SWAPS_PER_CALL" \
  -var SKIP_PLACE     1

echo "=========================================================="
echo "Finished       : $(date)"
echo "=========================================================="
