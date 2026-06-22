# Data

`al_mg_pca_deltaE_verified.csv` contains 82,645 aligned Al(Mg)
grain-boundary sites with:

- `row_id`
- `site_id_if_verified`
- raw PCA coordinates `pc1` ... `pc10`
- `deltaE_kJmol`

The PCA coordinates reproduce the original Wagih notebook output. The target
is the Mg-in-Al segregation energy in kJ/mol. See the parent README and
`scripts/resolve_al_mg_soap_energy_alignment.py` for alignment provenance.

Source: Wagih, Larsen, and Schuh (2020), Zenodo
<https://doi.org/10.5281/zenodo.4107058>.

The 671 MB raw SOAP matrix is not committed to Git.
