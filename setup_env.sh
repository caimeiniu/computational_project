#!/bin/bash
# Load required modules and activate virtual environment
# Usage: source setup_env.sh

module load stack/2024-06 gcc/12.2.0 openmpi/4.1.6 python/3.12.8 lammps/20240829.4

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/venv/bin/activate"

echo "Environment ready: Python $(python3 --version 2>&1 | cut -d' ' -f2), LAMMPS $(lmp -h 2>&1 | head -1 | grep -oP '\d+ \w+ \d+')"
