# Result Snapshot: 2026-06-03

This supersedes `RESULT_SNAPSHOT_2026-06-02.md` after jobs `1803476`,
`1803474`, and `1803475` landed. Treat these values as the current report
source of truth.

## Job Landing

| job ID | run | scheduler outcome | report output |
|---:|---|---|---|
| 1803476 | `hmc_AlMg_T500_Xc0.10_fdseed_resume11` | `COMPLETED`, 21:17:45 | `results/hmc_T500_Xc0.10_fdseed_resume11.json` |
| 1803474 | `hmc_AlMg_T700_Xc0.10_fdseed_targeted_resume5` | `COMPLETED`, 22:54:34 | `results/hmc_T700_Xc0.10_fdseed_targeted_resume5.json` |
| 1803475 | `hmc_AlMg_T800_Xc0.10_fdseed_targeted_resume5` | `TIMEOUT`, manually postprocessed from scratch dump | `results/hmc_T800_Xc0.10_fdseed_targeted_resume5.json` |

The T800 timeout was a wall-time cutoff, not a simulation crash. The scratch
dump and restarts existed under `/cluster/fs/scratchnv/06/cainiu/hmc_AlMg`, and
`scripts/hmc_xgb_timeseries_targeted.py` produced the repo-side JSON/PNG.

## DeltaE Spectrum Validation

| quantity | ours | Wagih |
|---|---:|---:|
| n sites | 500 | 82646 |
| mean DeltaE (kJ/mol) | -6.906 | -6.814 |
| std DeltaE (kJ/mol) | 15.073 | 15.853 |
| skew | -0.213 | -0.224 |
| skew-normal mu (kJ/mol) | 6.343 | 6.718 |
| skew-normal sigma (kJ/mol) | 20.057 | 20.843 |
| skew-normal alpha | -1.465 | -1.395 |
| KS two-sample p | 0.892 |  |

Supporting files:

- `results/compare_vs_wagih_200A_tight.json`
- `results/bootstrap_vs_wagih_200A_tight.json`
- `results/paired_pipeline_residual_n500_tight.json`

Paired-pipeline residual on Wagih structure: mean PE residual
`-0.0346 meV`, Pearson r = `0.999999999996`.

## Current HMC Concentration Points

"UB" means descending upper-bound or bracket data, not equilibrium.

| run | T (K) | Xc | X_GB | CI95 HW | first -> last | FD | HMC - FD | fwd/rev | PE drift | status |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---|
| `hmc_T500_Xc0.01_fdseed_resume3` | 500 | 0.01 | 0.04053 | 0.000027 | 0.04021 -> 0.04052 | 0.0526 | -0.01207 | 1.019 | -82.4 eV | practical dilute anchor |
| `hmc_T500_Xc0.04_fdseed_resume` | 500 | 0.04 | 0.13109 | 0.000092 | 0.13192 -> 0.13051 | 0.1912 | -0.06011 | 0.536 | -56.0 eV | caveated; still has tail drift |
| `hmc_T500_Xc0.075_eq_cont` | 500 | 0.075 | 0.22893 | 0.002128 | 0.25398 -> 0.21772 | 0.3007 | -0.07177 | 0.292 | -325.5 eV | UB; first-breakdown region |
| `hmc_T500_Xc0.10_fdseed_resume11` | 500 | 0.10 | 0.18667 | 0.000108 | 0.18778 -> 0.18598 | 0.3519 | -0.16523 | 0.889 | -43.5 eV | latest headline UB |
| `hmc_T500_Xc0.10_multistart_xgb0.3` | 500 | 0.10 | 0.24591 | 0.003570 | 0.29961 -> 0.22778 | 0.3520 | -0.10609 | 0.270 | -2155.3 eV | older independent-IC UB |
| `hmc_T700_Xc0.10_random_targeted_resume2` | 700 | 0.10 | 0.18722 | 0.000456 | 0.18135 -> 0.18938 | 0.2956 | -0.10838 | 1.101 | -378.3 eV | lower branch, still rising |
| `hmc_T700_Xc0.10_fdseed_targeted_resume5` | 700 | 0.10 | 0.22348 | 0.000270 | 0.22425 -> 0.22167 | 0.2956 | -0.07212 | 0.935 | -75.9 eV | upper branch, still descending |
| `hmc_T800_Xc0.10_random_targeted_resume3` | 800 | 0.10 | 0.19087 | 0.000357 | 0.18556 -> 0.19286 | 0.2734 | -0.08253 | 1.069 | -68.7 eV | lower branch, still rising |
| `hmc_T800_Xc0.10_fdseed_targeted_resume5` | 800 | 0.10 | 0.21575 | 0.000368 | 0.21673 -> 0.21284 | 0.2734 | -0.05765 | 0.927 | -217.8 eV | upper branch, still descending |

Current bracket widths:

- T700 Xc=0.10 upper-lower mean gap: `0.03626`.
- T800 Xc=0.10 upper-lower mean gap: `0.02489`.

Relative to the 2026-06-02 snapshot, the new upper-branch means moved down by:

- T500 resume11 vs resume10: `-0.00200`.
- T700 resume5 vs resume4: `-0.00290`.
- T800 resume5 vs resume4: `-0.00357`.

## Regenerated Report Figures

- Fig. 2 candidate: `figures/02_prediction_hmc_vs_fd_current.png`
  with traceability JSON `results/almg_report_02_prediction_hmc_vs_fd_current.json`.
- Fig. 5 candidate: `figures/05_ergodicity_bracket_T700_T800_current.png`
  with traceability JSON `results/almg_report_05_ergodicity_bracket_T700_T800_current.json`.

## Mechanism Snapshot

Fixed-DeltaE local-density test at T=500 K, Xc=0.075:

- Favorable DeltaE window: `[-30, -15] kJ/mol`.
- Sites in window: `112`.
- Wagih FD predicted P range in this window: `0.750` to `0.991`.
- Empirical linear slope versus local Mg count: `-0.0804` per Mg neighbor.
- Occupancy drops from roughly `0.33-0.57` for local Mg count 2-5 to
  `0.13-0.15` for local Mg count 6-7.

Supporting files:

- `figures/01_MgMg_clustering.png`
- `figures/02_occupation_breakdown.png`
- `figures/03_mgmg_repulsion_fixed_dE.png`
- `results/solute_correlation_analysis.json`
- `results/03_mgmg_repulsion_fixed_dE.json`

## Writing Implications

- The T500 Xc=0.10 headline value is now `X_GB = 0.18667`, far below
  canonical FD `0.3519`; because it is still descending, call it an upper bound.
- The T700 and T800 two-branch checks remain open. They support the statement
  that both observed branches are below FD, but they do not certify equilibrium.
- Report wording should distinguish practical anchors, descending upper bounds,
  and two-branch bracket diagnostics.
