# Pt(Au) GB segregation pipeline

Mirror of the Al(Mg) pipeline at `project/data/decks/anneal_AlMg.lammps` +
`project/scripts/sample_delta_e.py`, adapted to the **Pt-host, Au-solute**
direction with the **O'Brien et al. 2017 PtAu EAM/alloy** potential.

Direction chosen because Wagih's polycrystal database
(`Pt/Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/Pt_Au_20nm_GB_segregation.dump`)
stores Pt(Au) segregation energies computed with this exact same EAM file —
so the resulting ΔE_seg spectrum can be KS-tested directly against per-site
Wagih data with no potential-mismatch confound.

## Files

```
project/PtAu/
├── data/
│   ├── potentials/PtAu.eam.alloy           # O'Brien 2017 (NIST), eam/alloy
│   └── decks/
│       ├── anneal_PtAu.lammps              # Wagih-style anneal, eam/alloy
│       ├── submit_anneal_PtAu_100A.sh      # SLURM submit for 10³ nm³ prototype
│       └── submit_delta_e_PtAu_100A.sh     # SLURM submit for n=500 sampling
├── scripts/
│   └── sample_delta_e_PtAu.py              # copy with pair_style eam/alloy
├── output/                                  # fits + figures (empty until run)
└── README.md                                # this file
```

Generic scripts used verbatim from `project/scripts/`:
`generate_polycrystal.py`, `gb_identify.py`, `fit_delta_e_spectrum.py`,
`fermi_dirac_predict.py`, `compare_vs_wagih.py`, `bootstrap_vs_wagih.py`.

## Parameters

| Quantity | Value | Source |
|---|---|---|
| Host (type 1) | Pt | O'Brien EAM file header |
| Solute (type 2) | Au | O'Brien EAM file header |
| Lattice a | 3.9764 Å | EAM equilibrium FCC line |
| Pt mass | 195.0900 g/mol | EAM header |
| Au mass | 196.9665 g/mol | EAM header |
| T_melt(Pt) | 2041 K | CRC handbook |
| T_hold | 816 K (= 0.4 T_melt) | matches Al(Mg) protocol choice |
| Box | 100 × 100 × 100 Å | 10³ nm³ prototype |
| Grains | 8 | seed=1 |
| Atoms | 62,096 | 2.4% deficit from GB close-pair removal |

## Pipeline

```bash
SCRATCH=/cluster/scratch/$USER/prototype_PtAu_100A
PROJECT=/cluster/home/$USER/Computational_modeling/project
mkdir -p "$SCRATCH" && cd "$SCRATCH"

# ----- Step 1: polycrystal (login node, ~1 min) — ALREADY DONE -----
# python "$PROJECT/scripts/generate_polycrystal.py" \
#   --structure fcc --box 100 --grains 8 --lattice-a 3.9764 \
#   --structure-seed 1 --types 2 \
#   --out poly_Pt_100A_8g.lmp

# ----- Step 2: anneal (~1 h on 16 cores, 6 h budget) -----
sbatch "$PROJECT/PtAu/data/decks/submit_anneal_PtAu_100A.sh"

# ----- Step 3: gb_identify (login node, ~1 min) — happens automatically inside step 4 -----
# python "$PROJECT/scripts/gb_identify.py" \
#   poly_Pt_100A_8g_annealed.lmp \
#   --parent fcc --lattice-a 3.9764 \
#   --out-mask gb_mask_PtAu_100A.npy \
#   --out-report gb_info_PtAu_100A.json

# ----- Step 4: per-site ΔE_seg (n=500, ~15 min on 16 cores, 4 h budget) -----
sbatch "$PROJECT/PtAu/data/decks/submit_delta_e_PtAu_100A.sh"

# ----- Step 5: fit + predict (login node, seconds) -----
python "$PROJECT/scripts/fit_delta_e_spectrum.py" \
    --npz "$SCRATCH/delta_e_results_n500_PtAu_100A_tight.npz" \
    --out-png "$PROJECT/PtAu/output/delta_e_spectrum_PtAu_100A.png" \
    --out-json "$PROJECT/PtAu/output/delta_e_fit_PtAu_100A.json"

python "$PROJECT/scripts/fermi_dirac_predict.py" \
    --ours-npz "$SCRATCH/delta_e_results_n500_PtAu_100A_tight.npz" \
    --T-list 500,700,900,1100 \
    --out-png "$PROJECT/PtAu/output/fd_curves_PtAu_100A.png" \
    --out-json "$PROJECT/PtAu/output/fd_curves_PtAu_100A.json"
```

## KS test against Wagih reference

Wagih reference is in the local tar at
`/cluster/scratch/cainiu/wagih_zenodo/learning_segregation_energies.tar.bz2`,
path inside:

```
learning_segregation_energies/segregation_spectra_database_accelerated_model/Pt/Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/Pt_Au_20nm_GB_segregation.dump
```

Extract once:

```bash
cd /cluster/scratch/$USER/wagih_zenodo
tar -xjf learning_segregation_energies.tar.bz2 \
    learning_segregation_energies/segregation_spectra_database_accelerated_model/Pt/Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/Pt_Au_20nm_GB_segregation.dump
```

Then run the Pt(Au) compare driver (parses the dump's `seg_kJ_per_mol`
column directly — no separate bulk-energy file needed because Wagih
already converted to ΔE_seg in kJ/mol):

```bash
WAGIH_DUMP=/cluster/scratch/$USER/wagih_zenodo/learning_segregation_energies/segregation_spectra_database_accelerated_model/Pt/Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/Pt_Au_20nm_GB_segregation.dump

python "$PROJECT/PtAu/scripts/compare_vs_wagih_PtAu.py" \
    --ours-npz "$SCRATCH/delta_e_results_n500_PtAu_100A_tight.npz" \
    --wagih-dump "$WAGIH_DUMP" \
    --out-png  "$PROJECT/PtAu/output/compare_vs_wagih_PtAu_100A.png" \
    --out-json "$PROJECT/PtAu/output/compare_vs_wagih_PtAu_100A.json"

python "$PROJECT/PtAu/scripts/bootstrap_vs_wagih_PtAu.py" \
    --ours-npz "$SCRATCH/delta_e_results_n500_PtAu_100A_tight.npz" \
    --wagih-dump "$WAGIH_DUMP" \
    --out-json "$PROJECT/PtAu/output/bootstrap_vs_wagih_PtAu_100A.json"
```

Project pass bar: KS p > 0.5 ("spectrum-level indistinguishable", per
`feedback_ks_test` memory). For Al(Mg) we got D=0.026, p=0.89 with
n=500 ours vs n=82,646 Wagih; expect a similar order for Pt(Au) if
the pipelines align.

## What's left for the teammate

The deck infrastructure is complete and validated by syntax/sanity checks.
The remaining steps are wall-time-bound execution:

1. **`sbatch submit_anneal_PtAu_100A.sh`** — produces
   `poly_Pt_100A_8g_annealed.lmp` (~1 h, 6 h budget).
2. **`sbatch submit_delta_e_PtAu_100A.sh`** — auto-runs `gb_identify.py` if
   the mask is missing, then n=500 sampling (~15 min, 4 h budget). Produces
   `delta_e_results_n500_PtAu_100A_tight.npz`.
3. **Fit + KS test** — run `fit_delta_e_spectrum.py` for the (μ, σ, α)
   skew-normal fit, then `compare_vs_wagih_PtAu.py` for the KS test
   against the Wagih dump (commands above).
4. **Sanity-check** the (μ, σ, α) against Wagih's SI Fig 4 panel for Pt(Au)
   (Wagih's superscript `Au^N` on the panel labels which SI reference number
   for the EAM; for O'Brien 2017 that is the same potential we used).

## Differences from Al(Mg) pipeline (audit trail)

| File | Change |
|------|--------|
| `anneal_PtAu.lammps` | `pair_style eam/fs` → `eam/alloy`; defaults `EL1=Al, MASS1=26.9815` → `EL1=Pt, MASS1=195.0900`; `EL2=Mg` → `EL2=Au`; `T_HOLD=373` → `T_HOLD=816`; `POTFILE` path |
| `sample_delta_e_PtAu.py` | `pair_style eam/fs` → `eam/alloy` (line 454 of the canonical script); `_DEFAULT_ELEMENTS = ("Al","Mg")` → `("Pt","Au")`; `_DEFAULT_MASSES = (26.9815, 24.3050)` → `(195.0900, 196.9665)`; docstring + CLI defaults updated |
| `submit_anneal_PtAu_100A.sh` | scratch dir `prototype_AlMg_100A` → `prototype_PtAu_100A`; deck path; T_HOLD; vars EL1/EL2/MASS1/MASS2 passed explicitly; 4 h → 6 h budget (higher T_hold → more cool ps) |
| `submit_delta_e_PtAu_100A.sh` | scratch dir; driver path; --elements/--masses passed explicitly; gb_identify auto-fallback uses --lattice-a 3.9764 (Pt) instead of 4.05 (Al) |
| `compare_vs_wagih_PtAu.py` | new file; `load_wagih(seg.txt + bulk.dat)` → `load_wagih_dump(dump)` reading `seg_kJ_per_mol` column directly. The Pt(Au) reference lives in the accelerated_model dump (already in kJ/mol per site, bulk atoms = 0), not in the smaller `seg_energies_*.txt + bulk_solute_*.dat` pair that exists only for Al(Mg) under `machine_learning_notebook/`. |
| `bootstrap_vs_wagih_PtAu.py` | new file; same `load_wagih_dump` adapter as `compare_vs_wagih_PtAu.py`; otherwise identical to the Al(Mg) bootstrap script. |

No canonical script (`project/scripts/sample_delta_e.py`,
`project/scripts/compare_vs_wagih.py`,
`project/scripts/bootstrap_vs_wagih.py`,
`project/data/decks/anneal_AlMg.lammps`) was modified, per the
no-in-place-script-edits rule.
