# Provenance

This file maps each report claim to source files. Use it to avoid mixing stale
presentation material with current Al(Mg) results.

## Geometry And Mask

| claim | source |
|---|---|
| N_total = 475715 | `assets/data/gb_mask_200A.npy` and HMC JSONs |
| N_GB = 89042 | `assets/data/gb_mask_200A.npy` and HMC JSONs |
| GB fraction = 0.187175 | HMC JSONs, `gb_mask_200A.npy` |
| 200 A Al polycrystal, 16 grains | `assets/notes/pipeline_walkthrough.md`, `assets/figures/01_polycrystal_geometry.png` |

## DeltaE Spectrum

| claim | source |
|---|---|
| n=500 sampled GB sites | `assets/data/delta_e_results_n500_200A_tight.npz` |
| ours mean DeltaE = -6.906 kJ/mol | `assets/data/compare_vs_wagih_200A_tight.json` |
| Wagih mean DeltaE = -6.814 kJ/mol | `assets/data/compare_vs_wagih_200A_tight.json` |
| KS p = 0.892 | `assets/data/compare_vs_wagih_200A_tight.json` |
| bootstrap stats inside 95% intervals | `assets/data/bootstrap_vs_wagih_200A_tight.json` |
| paired residual mean = -0.0346 meV | `assets/data/paired_pipeline_residual_n500_tight.json` |

## HMC Results

| claim | source |
|---|---|
| T500 Xc=0.10 latest upper bound X_GB = 0.18667 | `assets/data/hmc_T500_Xc010_fdseed_resume11.json` |
| T500 Xc=0.10 FD = 0.3519 and HMC-FD gap = -0.16523 | `assets/data/hmc_T500_Xc010_fdseed_resume11.json` |
| T700 Xc=0.10 upper branch X_GB = 0.22348 | `assets/data/hmc_T700_Xc010_fdseed_targeted_resume5.json` |
| T700 Xc=0.10 lower branch X_GB = 0.18722 | `assets/data/hmc_T700_Xc010_random_targeted_resume2.json` |
| T700 Xc=0.10 current bracket width = 0.03626 | `assets/data/almg_report_05_ergodicity_bracket_T700_T800_current.json` |
| T800 Xc=0.10 upper branch X_GB = 0.21575 | `assets/data/hmc_T800_Xc010_fdseed_targeted_resume5.json` |
| T800 Xc=0.10 lower branch X_GB = 0.19087 | `assets/data/hmc_T800_Xc010_random_targeted_resume3.json` |
| T800 Xc=0.10 current bracket width = 0.02489 | `assets/data/almg_report_05_ergodicity_bracket_T700_T800_current.json` |
| Current Fig. 2 HMC-vs-FD prediction panel | `assets/figures/02_prediction_hmc_vs_fd_current.png`, `assets/data/almg_report_02_prediction_hmc_vs_fd_current.json` |
| Current Fig. 5 ergodicity/bracket panel | `assets/figures/05_ergodicity_bracket_T700_T800_current.png`, `assets/data/almg_report_05_ergodicity_bracket_T700_T800_current.json` |

## Mechanism

| claim | source |
|---|---|
| Mg-Mg non-random spatial structure | `assets/data/solute_correlation_analysis.json`, `assets/figures/07_MgMg_clustering.png` |
| occupation below Wagih sigmoid at favorable DeltaE | `assets/figures/08_occupation_breakdown.png` |
| fixed-DeltaE favorable window [-30, -15] kJ/mol | `assets/data/03_mgmg_repulsion_fixed_dE.json` |
| slope vs local Mg count = -0.0804 per neighbor | `assets/data/03_mgmg_repulsion_fixed_dE.json` |

## Notes On Reliability

- The HMC concentration values above are not all equilibrium values. For the
  current report, distinguish practical anchors, descending upper bounds, and
  two-branch brackets.
- Jobs `1803476` and `1803474` completed. Job `1803475` timed out at the wall
  limit and was manually postprocessed from the scratch dump.
- The `/cluster/scratch/cainiu/...` path is currently broken with "Too many
  levels of symbolic links"; use repo-side snapshots and `/cluster/fs/scratchnv`
  paths for current work.
