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
│       ├── submit_delta_e_PtAu_100A.sh     # SLURM submit for n=500 sampling
│       ├── hmc_PtAu.lammps                 # fresh HMC run deck
│       ├── hmc_PtAu_resume.lammps          # restart-based HMC continuation deck
│       ├── submit_hmc_PtAu_T700_bracket_jiayi.sh
│       │                                       # 700 K lower-Xc array bracket
│       ├── submit_hmc_PtAu_gbseed_check_jiayi.sh
│       │                                       # reverse IC check: Au starts on GB
│       ├── submit_hmc_PtAu_T700_Xc0.10_gbseed_resume_jiayi.sh
│       │                                       # continue slow Xc=0.10 GB-seeded run
│       └── submit_hmc_PtAu_T700_Xc0.10_random_resume66755862_jiayi.sh
│                                               # resumes the 700 K, Xc=0.10 job 66755862
├── scripts/
│   ├── sample_delta_e_PtAu.py              # copy with pair_style eam/alloy
│   ├── fit_delta_e_spectrum_PtAu.py        # copy with Pt(Au) title + Wagih Pt(Au) ref values
│   ├── fermi_dirac_predict_PtAu.py         # copy with Pt(Au) title + "100 Å" legend + --wagih-dump
│   ├── compare_vs_wagih_PtAu.py            # new; KS driver reading dump's seg_kJ_per_mol
│   ├── bootstrap_vs_wagih_PtAu.py          # new; bootstrap CI driver reading the same dump
│   ├── summarize_hmc_scan_PtAu.py          # HMC JSONs -> scan CSV with closed-box FD
│   ├── plot_hmc_scan_PtAu.py               # scan CSV -> final 700 K HMC/FD figure
│   ├── seed_gb_solute_PtAu.py              # make GB-seeded initial structures
│   ├── plot_hmc_initial_condition_check_PtAu.py
│                                               # random vs GB-seeded X_GB(t)
│   ├── mark_gb_solute_for_ovito_PtAu.py    # rewrite types for Pt/Au + GB coloring
│   └── make_ptau_ovito_snapshots.sh        # batch-export X=0.03 OVITO files
├── output/                                  # fits + figures (gitignored)
├── CHANGELOG.md                             # decisions + status log (reverse chronological)
└── README.md                                # this file
```

Generic scripts used verbatim from `project/scripts/` (alloy-agnostic via CLI flags):
`generate_polycrystal.py`, `gb_identify.py`.

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
WAGIH_DUMP=/cluster/scratch/$USER/wagih_zenodo/learning_segregation_energies/segregation_spectra_database_accelerated_model/Pt/Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/Pt_Au_20nm_GB_segregation.dump
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
# Pt(Au)-specific copies — canonical scripts in project/scripts/ are
# unchanged (they hardcode Al(Mg) title + Wagih Mg^15 reference).
python "$PROJECT/PtAu/scripts/fit_delta_e_spectrum_PtAu.py" \
    --npz "$SCRATCH/delta_e_results_n500_PtAu_100A_tight.npz" \
    --out-png "$PROJECT/PtAu/output/delta_e_spectrum_PtAu_100A.png" \
    --out-json "$PROJECT/PtAu/output/delta_e_fit_PtAu_100A.json"

python "$PROJECT/PtAu/scripts/fermi_dirac_predict_PtAu.py" \
    --ours-npz "$SCRATCH/delta_e_results_n500_PtAu_100A_tight.npz" \
    --wagih-dump "$WAGIH_DUMP" \
    --T-list 500,700,900,1100 \
    --n-total 62096 \
    --n-gb 23272 \
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

# ----- Step 6: HMC resume for T=700 K, Xc=0.10, job 66755862 -----
# Default continuation is 300 ps. Override with e.g. PROD_PS=600 sbatch ...
sbatch "$PROJECT/PtAu/data/decks/submit_hmc_PtAu_T700_Xc0.10_random_resume66755862_jiayi.sh"

# ----- Step 7: HMC lower-Xc bracket at 700 K -----
# Follows the Al(Mg) FD-first scan idea, now adjusted to Pt(Au):
# use the converged Xc=0.10 point as the high-concentration anchor and
# run lower total-Au fractions to locate the onset of HMC/FD divergence.
sbatch "$PROJECT/PtAu/data/decks/submit_hmc_PtAu_T700_bracket_jiayi.sh"

python "$PROJECT/PtAu/scripts/summarize_hmc_scan_PtAu.py" \
    "$PROJECT"/PtAu/output/hmc_PtAu_T700_Xc*_random_xgb_summary.json \
    "$PROJECT"/PtAu/output/hmc_PtAu_T700_Xc0.10_random_resume66755862_xgb_summary.json \
    --deltae-npz "$SCRATCH/delta_e_results_n500_PtAu_100A_tight.npz" \
    --n-total 62096 \
    --n-gb 23272 \
    --out-csv "$PROJECT/PtAu/output/hmc_PtAu_T700_scan_summary.csv"

python "$PROJECT/PtAu/scripts/plot_hmc_scan_PtAu.py" \
    --scan-csv "$PROJECT/PtAu/output/hmc_PtAu_T700_scan_summary.csv" \
    --out-png "$PROJECT/PtAu/output/hmc_PtAu_T700_scan.png" \
    --out-csv "$PROJECT/PtAu/output/hmc_PtAu_T700_plot_table.csv"

# ----- Step 8: reverse initial-condition check -----
# Starts with Au already preferentially on GB sites. If equilibrated, this
# over-segregated start and the homogeneous random start should converge to
# the same X_GB plateau.
sbatch "$PROJECT/PtAu/data/decks/submit_hmc_PtAu_gbseed_check_jiayi.sh"

python "$PROJECT/PtAu/scripts/plot_hmc_initial_condition_check_PtAu.py" \
    --random-csv "$PROJECT/PtAu/output/hmc_PtAu_T700_Xc0.02_random_xgb_timeseries.csv" \
    --gbseed-csv "$PROJECT/PtAu/output/hmc_PtAu_T700_Xc0.02_gbseed_xgb_timeseries.csv" \
    --out-png "$PROJECT/PtAu/output/hmc_PtAu_T700_Xc0.02_ic_check.png" \
    --title "Pt(Au) 700 K, X=0.02: random vs GB-seeded"

# If the high-concentration GB-seeded check remains above the random-start
# plateau, continue only that run:
sbatch "$PROJECT/PtAu/data/decks/submit_hmc_PtAu_T700_Xc0.10_gbseed_resume_jiayi.sh"

# OVITO-friendly structure snapshots for professor-facing structural checks.
bash "$PROJECT/PtAu/scripts/make_ptau_ovito_snapshots.sh"
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
   the mask is missing, then n=500 sampling (~80–90 min, 4 h budget;
   ~5× slower per site than Al(Mg) because Pt's stiffer EAM needs more
   CG iterations to hit the tight 1e-25 energy tolerance). Produces
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
| `fit_delta_e_spectrum_PtAu.py` | copy with three substantive label changes: plot title `"Al(Mg) per-site GB segregation spectrum"` → `"Pt(Au) per-site GB segregation spectrum"`; the `WAGIH_ALMG` constant (μ=9, σ=23, α=−2.3, source `Mg^15` SI Fig. 3) → `WAGIH_PTAU` (μ=3.65, σ=11.92, α=−1.42, source `Wagih 2020 Zenodo accelerated_model Pt_Au_20nm_GB_segregation.dump`); JSON key `wagih_alm_reference` → `wagih_ptau_reference`. Pt(Au) reference values come from `scipy.stats.skewnorm.fit` on the Wagih dump's 97,440 seg_kJ_per_mol entries (the same dump used by `compare_vs_wagih_PtAu.py`). |
| `fermi_dirac_predict_PtAu.py` | copy with two label changes: plot title `"Fermi-Dirac dilute-limit prediction — Al(Mg)"` → `"... Pt(Au)"`; legend `"ours 200 Å"` → `"ours 100 Å"` (we are on the 10³ nm³ prototype). Additionally swaps the dead `--wagih-seg`/`--wagih-bulk` flag pair (Al(Mg) seg.txt + bulk.dat format) for a single `--wagih-dump` flag that reads `Pt_Au_20nm_GB_segregation.dump` directly via the same loader as `bootstrap_vs_wagih_PtAu.py`. Default `--T-list` shifted `300,500,700,900` → `500,700,900,1100` to bracket Pt's higher T_melt. Pt(Au) also adds optional `--n-total/--n-gb` closed-box FD curves, because HMC fixes total Au while Wagih's reservoir formula fixes bulk Au. |

No canonical script (`project/scripts/sample_delta_e.py`,
`project/scripts/fit_delta_e_spectrum.py`,
`project/scripts/fermi_dirac_predict.py`,
`project/scripts/compare_vs_wagih.py`,
`project/scripts/bootstrap_vs_wagih.py`,
`project/data/decks/anneal_AlMg.lammps`) was modified, per the
no-in-place-script-edits rule.
