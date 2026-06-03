# Figure Plan

This plan is for the Al(Mg) contribution only and respects the shared report
figure budget. The main text should use five Al(Mg)-related figures total:

1. Al(Mg) polycrystal / method structure figure.
2. Prediction result.
3. Mechanism explanation 1.
4. Mechanism explanation 2.
5. Ergodicity / bracket check.

Everything else in `figures/` is provenance, backup, or supplement.

## Main-Text Figure Set

| slot | report role | current asset | status | final action |
|---:|---|---|---|---|
| Fig. 1 | Al(Mg) polycrystal structure / method basis | Already in `main.tex` through `methods_draft.tex`: `polycrystal_geometry.png` inside `Fig.~\\ref{fig:method}` | already allocated | Do not add a second standalone structure figure. Keep the Al(Mg) structure panel in the existing methods figure. |
| Fig. 2 | Prediction result: FD baseline vs HMC Al(Mg) | `figures/02_prediction_hmc_vs_fd_current.png` | current candidate generated | Use this as the main result figure unless the caption needs a different point selection. |
| Fig. 3 | Mechanism 1: occupation probability vs DeltaE | `figures/02_occupation_breakdown.png` | report-ready | Use single-panel version for clarity. |
| Fig. 4 | Mechanism 2: fixed-DeltaE local Mg environment | `figures/03_mgmg_repulsion_fixed_dE.png` | report-ready | Use as the key mechanism proof: same DeltaE window, occupancy decreases with local Mg count. |
| Fig. 5 | Ergodicity / two-branch bracket check | `figures/05_ergodicity_bracket_T700_T800_current.png` | current candidate generated | Use as the compact bracket figure; keep individual diagnostics as provenance. |

## What Not To Put In Main Text

| asset | use |
|---|---|
| `figures/04_spectrum_match.png` | Important validation, but likely supplement or a small inset/table in methods, unless the report needs a separate validation figure. |
| `figures/01_MgMg_clustering.png` | Good exploratory mechanism figure, but less clean than fixed-DeltaE because it mixes geometry and interaction. Keep as supplement/Q&A. |

## Regeneration Completed 2026-06-03

### Fig. 2: Prediction result

Current prediction figure generated from the latest JSONs:

`project/AlMg/figures/02_prediction_hmc_vs_fd_current.png`

It includes:

- canon-FD prediction as the reference curve.
- T=500 K concentration series, with clear marker semantics:
  - practical anchors / caveated near-plateau points,
  - descending upper-bound points,
  - latest Xc=0.10 fdseed continuation.
- A concise legend. Avoid showing every diagnostic curve in the main figure.

### Fig. 5: Ergodicity / bracket check

Current bracket figure generated:

`project/AlMg/figures/05_ergodicity_bracket_T700_T800_current.png`

It includes:

- x-axis: continuation index or cumulative production time.
- y-axis: `X_GB`.
- two panels or two colors for 700 K and 800 K.
- lower random-targeted branch and upper fdseed-targeted branch.
- FD reference line for each temperature.
- non-equilibrated bracket values; the current upper/lower gaps are still open.

## Stale Or Archive Only

- Old defense headline plots and the embedded report/template directory are not
  part of this self-contained AlMg package.
