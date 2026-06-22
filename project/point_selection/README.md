# Al(Mg) training-site selection

Reproducibility package for the training-site-selection part of the project.
It follows the accelerated spectrum-learning workflow of Wagih, Larsen, and
Schuh, *Nature Communications* **11**, 6376 (2020), using the associated
[Zenodo record](https://doi.org/10.5281/zenodo.4107058).

The controlled comparison keeps the descriptor, evaluation population, and
final predictor fixed. Only the indices chosen for atomistic labeling change.

- Data: 82,645 aligned Al(Mg) grain-boundary sites.
- Features: raw `pc1` ... `pc10`, matching the original notebook PCA output.
- Preprocessing: no `StandardScaler` and no new PCA in the main comparison.
- Target: `deltaE_kJmol`.
- Final predictor: `sklearn.linear_model.LinearRegression()` for every method.
- Evaluation: predictions over all 82,645 sites.

## Selection methods

| Method | Acquisition rule |
|---|---|
| Random | Uniform sampling without replacement |
| `kmeans_raw_pca` | Cluster raw PCA coordinates; choose the nearest real site to each centroid |
| Bootstrap active | Select the largest bootstrap prediction uncertainty, `sigma_i` |
| Tail-aware active | Select the largest `-mu_i + lambda * sigma_i` |
| XGB-aware active | Select the largest `sigma_i P_i (1-P_i)` |

The active methods use true segregation energies only for the currently
labeled set. Unlabeled targets are used only after selection for evaluation.

## Folder contents

```text
point_selection/
├── data/                         verified aligned feature-target table
├── scripts/                      reproducible analyses
├── results/reference/            compact result tables used in the report
└── outputs/                      generated locally; ignored by git
```

Important scripts:

- `scripts/compare_training_site_selection.py`: main 20-200-label comparison.
- `scripts/run_all_data_baseline.py`: all-data 10-PC LinearRegression reference.
- `scripts/run_extended_budget.py`: selection curves through 1,000 labels.
- `scripts/run_feature_dimension_sweep.py`: 10/20/50/100-PC all-data sweep.
- `scripts/resolve_al_mg_soap_energy_alignment.py`: optional, expensive proof
  of SOAP/energy row alignment from the original Zenodo files.

## Quick start

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r project/point_selection/requirements.txt

python project/point_selection/scripts/compare_training_site_selection.py
```

The main script reads `data/al_mg_pca_deltaE_verified.csv` and writes only to
`project/point_selection/outputs/main_comparison/`. A faster smoke run is:

```bash
python project/point_selection/scripts/compare_training_site_selection.py \
  --random-repeats 1 --active-repeats 1 --n-bootstrap 2
```

Other analyses:

```bash
python project/point_selection/scripts/run_all_data_baseline.py
python project/point_selection/scripts/run_extended_budget.py --budgets 100 500 1000
```

The extended active-learning experiment is computationally expensive and is
resumable. Existing method/budget/repeat rows are skipped unless
`--overwrite` is supplied.

## Feature-dimension sweep

The raw `GB_SOAP_Al_Mg.npy` matrix is about 671 MB and is intentionally not
stored in Git. Download the Zenodo archive, then run:

```bash
python project/point_selection/scripts/run_feature_dimension_sweep.py \
  --soap /path/to/GB_SOAP_Al_Mg.npy
```

The script fails if the SOAP row count is not exactly 82,645; it never silently
truncates or invents an alignment.

## Alignment provenance

The committed CSV contains 82,645 SOAP/PCA rows aligned with Mg-in-Al
segregation energies. The verification workflow found that current OVITO
identifies one additional candidate, site ID `50044`, which is absent from the
saved Zenodo SOAP matrix. Removing that candidate reproduces the saved matrix
row-for-row. Segregation energies are calculated as

```text
deltaE = (E_GB,solute - E_bulk,solute) * 96.485  [kJ/mol]
```

The optional alignment script requires the original Zenodo archive plus
OVITO, ASE, and `quippy-ase`; see `requirements-alignment.txt`.

## Reference results

At 100 labels, raw-PCA k-means gives global MAE 3.652 kJ/mol, tail MAE
5.669 kJ/mol, and absolute XGB error 0.006004. The all-data 10-PC linear fit
has MAE 3.558 kJ/mol. Increasing the all-data PCA representation from 10 to
100 components reduces MAE from 3.558 to 2.913 kJ/mol. These all-data values
are in-sample feature/model references, not independent generalization tests.

The compact source tables are in `results/reference/`; large predictions,
figures, slide decks, and intermediate SOAP arrays are deliberately excluded.
