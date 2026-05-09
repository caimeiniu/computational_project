# References used in `methods_draft.tex` — verification status

**Generated 2026-05-07; verification completed 2026-05-09.**

User's standing instruction: **"不可以编造文献!!!"** — every `\cite{}` key
in `methods_draft.tex` must be backed by a paper whose DOI we have personally
verified resolves to the correct article.

All 5 cite keys below were verified on 2026-05-09 via the Crossref REST API
(`https://api.crossref.org/works/<DOI>`). Crossref is the authoritative DOI
registry; resolution there is equivalent to the publisher landing page (which
returned HTTP 403 to automated fetches for Taylor & Francis, APS, and AIP
during this verification, but Crossref agreed with all suggested entries).

The bibitem block ready to paste into `main.tex`'s `thebibliography` is in
`report/methods_bibitems.tex`.

---

## Status table

| key            | resolved? | journal / vol / pages                   | DOI                                |
|----------------|-----------|------------------------------------------|------------------------------------|
| `Mendelev2009` | VERIFIED  | Phil. Mag. **89**, 3269–3285 (2009)      | `10.1080/14786430903260727`        |
| `Wagih2020`    | VERIFIED  | Nat. Commun. **11**, 6376 (2020)         | `10.1038/s41467-020-20083-6`       |
| `LAMMPS`       | VERIFIED  | Comput. Phys. Commun. **271**, 108171 (2022) | `10.1016/j.cpc.2021.108171`    |
| `Stukowski2012`| VERIFIED  | Modell. Simul. Mater. Sci. Eng. **20**, 045021 (2012) | `10.1088/0965-0393/20/4/045021` |
| `Sadigh2012`   | VERIFIED  | Phys. Rev. B **85**, 184203 (2012)       | `10.1103/PhysRevB.85.184203`       |

## Per-entry detail

### `Mendelev2009` — Al–Mg EAM/FS potential

Mendelev, M.I.; Asta, M.; Rahman, M.J.; Hoyt, J.J.,
"Development of interatomic potentials appropriate for simulation of
solid–liquid interface properties in Al–Mg alloys",
*Philosophical Magazine* **89**, 3269–3285 (2009).
DOI `10.1080/14786430903260727`. **VERIFIED.**

The project's `data/potentials/Al-Mg.eam.fs` corresponds to this paper.
Wagih et al. 2020 cites it. CHANGELOG entry 2026-04-22 (evening)
records the choice.

### `Wagih2020` — Wagih–Larsen–Schuh polycrystal segregation framework

Wagih, M.; Larsen, P.M.; Schuh, C.A., "Learning grain boundary
segregation energy spectra in polycrystals", *Nature Communications*
**11**, 6376 (2020). DOI `10.1038/s41467-020-20083-6`. **VERIFIED.**

This is the framework paper the project replicates and stress-tests.
The Zenodo dataset we KS-test against is the supplementary release of
this paper.

### `LAMMPS` — LAMMPS molecular dynamics package

Thompson, A.P. et al., "LAMMPS — a flexible simulation tool for
particle-based materials modeling at the atomic, meso, and continuum
scales", *Computer Physics Communications* **271**, 108171 (2022).
DOI `10.1016/j.cpc.2021.108171`. **VERIFIED.**

This is the current canonical LAMMPS citation requested by the
developers. The older Plimpton 1995 reference
(*J. Comput. Phys.* **117**, 1–19; DOI `10.1006/jcph.1995.1039`) is
acceptable in some groups; we use Thompson 2022 throughout.

### `Stukowski2012` — adaptive common-neighbour analysis (a-CNA)

Stukowski, A., "Structure identification methods for atomistic
simulations of crystalline materials", *Modelling and Simulation in
Materials Science and Engineering* **20**, 045021 (2012).
DOI `10.1088/0965-0393/20/4/045021`. **VERIFIED.**

Standard reference for the adaptive CNA algorithm OVITO uses by
default. The 2010 Stukowski OVITO-software paper
(*Modell. Simul. Mater. Sci. Eng.* **18**, 015012, DOI
`10.1088/0965-0393/18/1/015012`) is a separate citation, not used here
because we do not credit the OVITO software for a-CNA itself; OVITO
features only as the rendering tool for the geometry panel
(`report/figures/polycrystal_geometry.png`), which is generated via
`scripts/render_methods_polycrystal.py` and does not require its own
reference in Methods.

### `Sadigh2012` — semi-grand-canonical swap MC algorithm

Sadigh, B.; Erhart, P.; Stukowski, A.; Caro, A.; Martinez, E.;
Zepeda-Ruiz, L., "Scalable parallel Monte Carlo algorithm for
atomistic simulations of precipitation in alloys", *Physical Review
B* **85**, 184203 (2012). DOI `10.1103/PhysRevB.85.184203`. **VERIFIED.**

This is the standard reference for the semi-grand-canonical swap MC
algorithm that LAMMPS `fix atom/swap` implements. We use it in
canonical (closed-box) mode — `Sadigh2012` covers both the SGC and
canonical variants.

---

## Removed: `VoroPP` (no longer cited)

The 2026-05-07 candidate list included Rycroft 2009
(*Chaos* **19**, 041111; DOI `10.1063/1.3215722`) for "3D Voronoi
tessellation". On checking `scripts/generate_polycrystal.py`
(2026-05-09), the implementation does **not** use Voro++. Each lattice
site is assigned to its nearest grain seed under PBC via SciPy's
`cKDTree`, which produces the Voronoi cell partition implicitly
without computing the Voronoi diagram. This is a textbook construction
not requiring its own software reference, so the cite has been
dropped. The Methods narrative now describes the algorithm in one line
("k-d-tree nearest-neighbour lookup under periodic boundary
conditions") and does not cite a Voronoi library.

---

## Optional adds (only if you decide they belong)

- **Duane et al. 1987** is currently inline in the HMC footnote, not
  as a `\cite{}`. If you prefer a `\cite{Duane1987}` instead, the DOI
  `10.1016/0370-2693(87)91197-X` resolves and is verified.
- **OVITO software**: not currently cited (see Stukowski2012 entry
  above for rationale).
- **SciPy**: not cited; standard practice is to reserve bibliography
  for scientific results / methods papers, not Python toolkits.

---

## House-style decisions (apply uniformly across the bibliography)

These are still pending; pick once and `methods_bibitems.tex` already
follows the chosen style.

- **DOI link style**: chosen format in `methods_bibitems.tex` is
  `\href{https://doi.org/...}{doi:...}` (clickable). RevTeX 4-2 has
  `hyperref` loaded already.
- **Journal abbreviation**: chosen in `methods_bibitems.tex` is the
  short form (e.g. `Philos. Mag.`, not `Philosophical Magazine`).
  Match `main.tex`'s existing bibitems if those use a different
  convention.
