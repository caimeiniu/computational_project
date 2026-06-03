# Reproduce The Al(Mg) Report Results

This is the clean-clone path for reproducing the numerical report snapshot and
the current Fig. 2 / Fig. 5 candidates without rerunning the multi-day LAMMPS
jobs.

## Scope

Included in the minimal reproducibility package:

- Python scripts that regenerate the current report figures.
- Curated postprocessed HMC JSON files used by the report.
- DeltaE spectrum and GB-mask snapshots needed for the canonical FD curve.
- Report-ready figure PNGs and traceability JSONs.
- Provenance notes mapping each report claim to its source file.

Not included:

- Raw LAMMPS dumps, restart files, and final `.lmp` snapshots. These are tens of
  GB in aggregate and are not needed to reproduce the report plots from the
  postprocessed data.
- Local continuation `submit_*.sh` wrappers that depend on scratch restart
  paths. They are useful for the original Euler campaign, not for a clean clone.

## Environment

From the repository root:

```bash
cd project/AlMg
python -m pip install -r requirements-report.txt
```

The report figure path needs only `numpy`, `scipy`, and `matplotlib`. Full
simulation reruns additionally need LAMMPS with the MANYBODY, MC, and VORONOI
packages as described in `README.md`.

## Regenerate Current Report Figures

```bash
cd project/AlMg
python scripts/report_al_mg_current_figures.py
```

Expected outputs:

```text
results/almg_report_02_prediction_hmc_vs_fd_current.json
results/almg_report_02_prediction_hmc_vs_fd_current.png
results/almg_report_05_ergodicity_bracket_T700_T800_current.json
results/almg_report_05_ergodicity_bracket_T700_T800_current.png
```

The report-ready PNGs are also stored in:

```text
figures/02_prediction_hmc_vs_fd_current.png
figures/05_ergodicity_bracket_T700_T800_current.png
```

## Verify Headline Numbers

```bash
cd project/AlMg
python - <<'PY'
import json
from pathlib import Path

fig2 = json.loads(Path("results/almg_report_02_prediction_hmc_vs_fd_current.json").read_text())
fig5 = json.loads(Path("results/almg_report_05_ergodicity_bracket_T700_T800_current.json").read_text())

latest = fig2["descending_upper_bound_points"][-1]
print(f"T500 Xc=0.10 X_GB = {latest['X_GB_mean']:.5f}")
print(f"T500 Xc=0.10 HMC-FD gap = {latest['gap_hmc_minus_fd']:.5f}")
for panel in fig5["panels"]:
    print(f"T{panel['T']} latest bracket = {panel['latest_upper_minus_lower']:.5f}")
PY
```

Expected values:

```text
T500 Xc=0.10 X_GB = 0.18667
T500 Xc=0.10 HMC-FD gap = -0.16523
T700 latest bracket = 0.03626
T800 latest bracket = 0.02489
```

## Source Of Truth

- Current numerical snapshot:
  `RESULT_SNAPSHOT_2026-06-03.md`
- Claim-to-file mapping:
  `PROVENANCE.md`
- Main-text figure plan:
  `FIGURE_PLAN.md`
