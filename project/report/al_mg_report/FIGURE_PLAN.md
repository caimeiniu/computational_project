# Figure Plan

This plan is for the Al(Mg) contribution only and respects the shared report
figure budget. The main text should use five Al(Mg)-related figures total:

1. Al(Mg) polycrystal / method structure figure.
2. Prediction result.
3. Mechanism explanation 1.
4. Mechanism explanation 2.
5. Ergodicity / bracket check.

Everything else in `assets/figures/` is provenance, backup, or supplement.

## Main-Text Figure Set

| slot | report role | current asset | status | final action |
|---:|---|---|---|---|
| Fig. 1 | Al(Mg) polycrystal structure / method basis | Already in `main.tex` through `methods_draft.tex`: `polycrystal_geometry.png` inside `Fig.~\\ref{fig:method}` | already allocated | Do not add a second standalone structure figure. Keep the Al(Mg) structure panel in the existing methods figure. |
| Fig. 2 | Prediction result: FD baseline vs HMC Al(Mg) | `assets/figures/02_prediction_hmc_vs_fd_current.png` | current candidate generated | Use this as the main result figure unless the caption needs a different point selection. |
| Fig. 3 | Mechanism 1: occupation probability vs DeltaE | `assets/figures/08_occupation_breakdown.png` | report-ready | Use single-panel version for clarity. Move 3-Xc version to supplement/internal backup. |
| Fig. 4 | Mechanism 2: fixed-DeltaE local Mg environment | `assets/figures/10_repulsion_fixed_deltaE.png` | report-ready | Use as the key mechanism proof: same DeltaE window, occupancy decreases with local Mg count. |
| Fig. 5 | Ergodicity / two-branch bracket check | `assets/figures/05_ergodicity_bracket_T700_T800_current.png` | current candidate generated | Use as the compact bracket figure; keep individual diagnostics as provenance. |

## What Not To Put In Main Text

| asset | use |
|---|---|
| `assets/figures/02_spectrum_match.png` | Important validation, but likely supplement or a small inset/table in methods, unless the report needs a separate validation figure. |
| `assets/figures/03_fd_curves_ours_vs_wagih_multiT.png` | Theory/method backup; not main text unless prediction figure needs an inset. |
| `assets/figures/07_MgMg_clustering.png` | Good exploratory mechanism figure, but less clean than fixed-DeltaE because it mixes geometry and interaction. Keep as supplement/Q&A. |
| `assets/figures/09_occupation_breakdown_3xc.png` | Backup to Fig. 3, likely supplement. |
| `assets/figures/11_repulsion_summary.png` | Talk-style summary; not needed if Fig. 4 is included. |
| `assets/figures/12_sampler_convergence_Xc0075.png` | Supplement / internal diagnostic. |
| `assets/figures/13_two_sided_verify_Xc005.png` | Supplement / internal diagnostic. |
| `assets/figures/14_panel_f_Xc004_gb_slab.png` | Optional visual; only use if Fig. 1 loses the polycrystal panel. |
| `assets/figures/15_ovito_segregation_legacy.png` | Old visual, do not use in main text. |
| workflow slides in `assets/figures/16_*` and `17_*` | Presentation backup, not report figures. |

## Regeneration Completed 2026-06-03

### Fig. 2: Prediction result

Current prediction figure generated from the latest JSONs:

`project/report/al_mg_report/assets/figures/02_prediction_hmc_vs_fd_current.png`

It includes:

- canon-FD prediction as the reference curve.
- T=500 K concentration series, with clear marker semantics:
  - practical anchors / caveated near-plateau points,
  - descending upper-bound points,
  - latest Xc=0.10 fdseed continuation.
- A concise legend. Avoid showing every diagnostic curve in the main figure.

### Fig. 5: Ergodicity / bracket check

Current bracket figure generated:

`project/report/al_mg_report/assets/figures/05_ergodicity_bracket_T700_T800_current.png`

It includes:

- x-axis: continuation index or cumulative production time.
- y-axis: `X_GB`.
- two panels or two colors for 700 K and 800 K.
- lower random-targeted branch and upper fdseed-targeted branch.
- FD reference line for each temperature.
- non-equilibrated bracket values; the current upper/lower gaps are still open.

## Stale Or Archive Only

- `project/report/figures/00_headline_hmc_vs_wagih_T500.png`: old defense
  state; predates later Xc=0.10 fdseed continuations.
- `project/report/figures/00_headline_hmc_vs_wagih_T500_9pt_2026-05-11.png`:
  useful archive only. It does not include `hmc_T500_Xc0.10_fdseed_resume10`
  or the current `resume11`.
- `project/report/README.md`: good explanatory background, but stale for final
  headline numbers.
- `project/report/69f62101465c6c925c0856e8/`: old embedded report/template
  directory. Ignore for current Al(Mg) drafting unless explicitly requested.
