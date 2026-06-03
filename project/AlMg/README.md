# Al(Mg) Report Workspace

Clean Al(Mg)-only workspace for reproducing the report results and writing the
final Al(Mg) section. Everything needed for the report figures is inside this
directory; raw LAMMPS dumps, restarts, and final `.lmp` snapshots are excluded.

## Source of truth

- Latest numerical results: `results/*.json`.
- Canonical-FD input snapshots: `data/*.npz` and `data/*.npy`.
- Report-ready figure candidates: `figures/`.
- Methods and explanatory notes: `notes/`.
- Reproduction commands: `REPRODUCE.md`.

Use this workspace instead of browsing the full repo when drafting the report.
The full repo includes PtAu, older presentation figures, and cluster-specific
job decks that are not needed to reproduce the current Al(Mg) report figures.

## Figure budget

The Al(Mg) contribution should use five main-text figure slots:

1. Existing `main.tex` / methods figure with the Al(Mg) polycrystal structure.
2. Prediction result: current HMC vs FD.
3. Mechanism: occupation-vs-DeltaE breakdown.
4. Mechanism: fixed-DeltaE local Mg environment / repulsion.
5. Ergodicity check: 700 K / 800 K bracket figure.

All other linked figures are backup or supplement candidates.

## Current thesis

For Al(Mg) grain-boundary segregation, the independent-site Fermi-Dirac
baseline overpredicts GB Mg occupation once finite-Mg interactions matter.
The strongest current evidence is:

- Our n=500 Al(Mg) DeltaE spectrum matches the Wagih 2020 Al(Mg) spectrum:
  KS p = 0.892, with all six bootstrap statistics inside Wagih-like 95%
  intervals.
- At T=500 K, Xc=0.10, the latest fdseed continuation gives
  X_GB = 0.1867 versus canon-FD = 0.3519. This is still not equilibrium,
  but it is a descending upper-bound point already far below FD.
- At T=700 K and T=800 K, the same-protocol fdseed-targeted upper branches
  also remain below FD and still descend. Their lower random-targeted branches
  remain separated, so those temperatures are bracket data, not equilibrium.
- Mechanism figures show non-independent Mg placement through Mg-Mg spatial
  correlations, occupation-vs-DeltaE breakdown, and a fixed-DeltaE local-density
  dependence.

## Important caveat

Do not use the old `report/README.md` numbers verbatim for the headline
concentration plot. It was written for the defense/presentation state around
2026-05-05 to 2026-05-13 and predates the later Xc=0.10 continuation results.
It is useful for explanations and Q&A language, not as the final numerical
source of truth.

## Latest continuation status

These jobs were submitted on 2026-06-02 and checked on 2026-06-03:

| job ID | run | result |
|---:|---|---|
| 1803476 | T500 Xc=0.10 fdseed resume11 | completed; auto-postprocessed |
| 1803474 | T700 Xc=0.10 fdseed-targeted resume5 | completed; auto-postprocessed |
| 1803475 | T800 Xc=0.10 fdseed-targeted resume5 | timed out; manually postprocessed from scratch dump |

Use `RESULT_SNAPSHOT_2026-06-03.md` for current numbers. The June 2 snapshot is
kept as the pre-landing handoff record.

## Files here

- `FIGURE_PLAN.md`: proposed paper/report figures and their current status.
- `RESULT_SNAPSHOT_2026-06-03.md`: current numerical state after the
  resume11/resume5 batch landed.
- `REPRODUCE.md`: clean-clone commands for regenerating the current report
  figures and checking the headline numbers.
- `PROVENANCE.md`: claim-to-file mapping and stale/ignore list.
- `data/`: small canonical-FD input snapshots.
- `results/`: postprocessed JSON/PNG outputs used by the report.
- `figures/`: report-ready figure assets.
- `scripts/`: Python scripts needed to regenerate figures and inspect HMC JSONs.
- `notes/`: methods and explanatory notes.
