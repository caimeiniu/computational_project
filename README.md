# Computational Project

## Setup

Load the environment (modules + venv):

```bash
source setup_env.sh
```

This loads:
- **Python 3.12.8** (gcc/12.2.0)
- **LAMMPS 20240829.4** (with OpenMPI 4.1.6)
- Virtual environment with pip

## Running LAMMPS

```bash
# Command line
lmp -in input.lammps

# MPI parallel
mpirun -np 4 lmp -in input.lammps

# Python interface
python3 -c "from lammps import lammps; lmp = lammps(); print(lmp.version())"
```
