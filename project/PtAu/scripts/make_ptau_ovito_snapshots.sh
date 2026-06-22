#!/bin/bash
# Generate OVITO-friendly Pt(Au) structure files for the 700 K IC checks.

set -euo pipefail

PROJECT=${PROJECT:-/cluster/home/jiayfu/computational_project/project}
RUN_DIR=${RUN_DIR:-/cluster/scratch/jiayfu/prototype_PtAu_100A}
GB_MASK=${GB_MASK:-$RUN_DIR/gb_mask_PtAu_100A.npy}
OUT_DIR=${OUT_DIR:-$PROJECT/PtAu/output/ovito_snapshots}
MARKER=$PROJECT/PtAu/scripts/mark_gb_solute_for_ovito_PtAu.py

mkdir -p "$OUT_DIR"

make_one() {
  local label=$1
  local input=$2
  local output=$OUT_DIR/${label}_ovito.lmp
  if [[ ! -f "$input" ]]; then
    echo "MISSING: $input" >&2
    return 1
  fi
  python3 "$MARKER" --input "$input" --gb-mask "$GB_MASK" --output "$output"
}

make_one X003_random_initial "$RUN_DIR/hmc_PtAu_T700_Xc0.03_random_initial.lmp"
make_one X003_random_final   "$RUN_DIR/hmc_PtAu_T700_Xc0.03_random_final.lmp"
make_one X003_gbseed_initial "$RUN_DIR/hmc_PtAu_T700_Xc0.03_gbseed_initial.lmp"
make_one X003_gbseed_final   "$RUN_DIR/hmc_PtAu_T700_Xc0.03_gbseed_final.lmp"

cat > "$OUT_DIR/README.txt" <<'EOF'
Pt(Au) OVITO type map:
  1 = Pt bulk
  2 = Au bulk
  3 = Pt GB
  4 = Au GB

Suggested OVITO colors:
  Pt bulk: light gray, low opacity
  Au bulk: orange/yellow
  Pt GB: blue or dark gray
  Au GB: red or bright yellow

Use a slab/slice through the 100 A box for readable before/after panels.
EOF

echo "OVITO snapshots written to $OUT_DIR"
