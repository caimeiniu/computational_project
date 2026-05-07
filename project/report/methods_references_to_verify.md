# References to verify before adding to bibliography

**Generated alongside `methods_draft.tex` on 2026-05-07.**

User's standing instruction: **"不可以编造文献!!!"** — every `\cite{}` key
in `methods_draft.tex` must be backed by a paper whose DOI you have personally
verified resolves to the correct article. The list below is my best
*candidate* for each cite key, with my confidence flagged honestly. **Do not
trust the DOI without resolving it yourself** (e.g. via `https://doi.org/<DOI>`
or a search on Google Scholar by author + year + journal).

For each entry I give: cite key → suggested paper → DOI candidate →
my confidence (high / medium / low) → what to check.

---

## `Mendelev2009` — Al–Mg EAM/FS potential

**Suggested**: Mendelev, M.I.; Asta, M.; Rahman, M.J.; Hoyt, J.J.,
"Development of interatomic potentials appropriate for simulation of
solid–liquid interface properties in Al–Mg alloys",
*Philosophical Magazine* **89**, 3269–3285 (2009).

**DOI candidate**: `10.1080/14786430903260727`

**Confidence**: HIGH — this is the canonical Al–Mg EAM/FS potential used
across the materials-modelling community; the project potential file
`Al-Mg.eam.fs` corresponds to this paper. Wagih et al. 2020 also cites it.

**Action**: confirm DOI resolves to the Phil. Mag. 89 paper with these
authors. The CHANGELOG entry 2026-04-22 (evening) of this project records
the choice and reference.

---

## `Wagih2020` — Wagih–Larsen–Schuh polycrystal segregation framework

**Suggested**: Wagih, M.; Larsen, P.M.; Schuh, C.A., "Learning grain
boundary segregation energy spectra in polycrystals", *Nature
Communications* **11**, 6376 (2020).

**DOI candidate**: `10.1038/s41467-020-20083-6`

**Confidence**: HIGH — this is the framework paper the entire project
replicates and stress-tests. The Zenodo dataset we KS-test against is
the supplementary release of this paper.

**Action**: confirm DOI resolves to Nat. Commun. 11:6376 (2020) by
Wagih, Larsen, Schuh. The Zenodo URL in the SI is also worth citing
separately if you reference the dataset directly.

---

## `LAMMPS` — LAMMPS molecular dynamics package

**Suggested (modern reference)**: Thompson, A.P. et al., "LAMMPS — a
flexible simulation tool for particle-based materials modeling at the
atomic, meso, and continuum scales", *Computer Physics Communications*
**271**, 108171 (2022).

**DOI candidate**: `10.1016/j.cpc.2021.108171`

**Confidence**: HIGH (Thompson 2022 is the current canonical LAMMPS
citation requested by the developers). Older code based on the original
1995 release sometimes cites Plimpton instead:
Plimpton, S., "Fast Parallel Algorithms for Short-Range Molecular
Dynamics", *J. Comput. Phys.* **117**, 1–19 (1995),
DOI candidate `10.1006/jcph.1995.1039`. Pick whichever matches your
group convention.

**Action**: pick one (Thompson 2022 vs Plimpton 1995); confirm DOI.

---

## `Stukowski2012` — adaptive common-neighbour analysis (a-CNA)

**Suggested**: Stukowski, A., "Structure identification methods for
atomistic simulations of crystalline materials", *Modelling and
Simulation in Materials Science and Engineering* **20**, 045021 (2012).

**DOI candidate**: `10.1088/0965-0393/20/4/045021`

**Confidence**: MEDIUM-HIGH — this is the standard reference for
adaptive CNA, the algorithm OVITO uses by default. There is also a
2010 Stukowski paper introducing OVITO itself
(*Modell. Simul. Mater. Sci. Eng.* **18**, 015012, DOI candidate
`10.1088/0965-0393/18/1/015012`); cite that separately if you decide
to credit the OVITO software in a-CNA / rendering.

**Action**: verify the 2012 paper DOI. If you also want to cite OVITO
the software, add `Stukowski2010` separately.

---

## `Sadigh2012` — semi-grand-canonical swap MC algorithm

**Suggested**: Sadigh, B.; Erhart, P.; Stukowski, A.; Caro, A.;
Martinez, E.; Zepeda-Ruiz, L., "Scalable parallel Monte Carlo algorithm
for atomistic simulations of precipitation in alloys",
*Physical Review B* **85**, 184203 (2012).

**DOI candidate**: `10.1103/PhysRevB.85.184203`

**Confidence**: MEDIUM — this is the standard cite for the SGC swap MC
algorithm that LAMMPS `fix atom/swap` implements. Verify author list,
volume, page, and DOI.

**Action**: verify DOI; check whether LAMMPS `fix atom/swap`
documentation cites this paper or a different one (sometimes
implementations cite the LAMMPS source code or a docs page in
addition).

---

## `VoroPP` — Voro++ Voronoi tessellation library

**Suggested**: Rycroft, C.H., "VORO++: A three-dimensional Voronoi cell
library in C++", *Chaos: An Interdisciplinary Journal of Nonlinear
Science* **19**, 041111 (2009).

**DOI candidate**: `10.1063/1.3215722`

**Confidence**: MEDIUM — this is the canonical Voro++ reference but
verify; depending on the Voronoi tool actually used in your
`scripts/generate_polycrystal.py` (could be Voro++, scipy's spatial
Voronoi, or a custom implementation), the cite might need to be
different. If you used a different library, replace this entry.

**Action**: open `scripts/generate_polycrystal.py`, check which Voronoi
implementation is imported, cite accordingly.

---

## Optional adds (only if relevant to your final draft)

- **Duane et al. 1987** (only in the HMC footnote, not as a `\cite{}`):
  Duane, S.; Kennedy, A.D.; Pendleton, B.J.; Roweth, D.,
  *Phys. Lett. B* **195**, 216 (1987). DOI candidate
  `10.1016/0370-2693(87)91197-X`. Used as inline disambiguation for
  the term "HMC"; if you want a `\cite{Duane1987}` instead of writing
  the reference inline in the footnote, verify and add.
- **Stukowski2010 OVITO software**: if you cite OVITO when describing
  Fig.~\ref{fig:method:geom} and the right-panel render, add
  Stukowski, A., "Visualization and analysis of atomistic simulation
  data with OVITO — the Open Visualization Tool",
  *Modell. Simul. Mater. Sci. Eng.* **18**, 015012 (2010), DOI candidate
  `10.1088/0965-0393/18/1/015012`.

---

## Once verified

Once each DOI resolves correctly, add one `\bibitem` per cite key to the
existing `\begin{thebibliography}` block in `main.tex`. Example skeleton:

```latex
\bibitem{Mendelev2009}
M.I.~Mendelev, M.~Asta, M.J.~Rahman, J.J.~Hoyt,
\emph{Development of interatomic potentials appropriate for simulation
of solid–liquid interface properties in Al–Mg alloys},
Philos. Mag. \textbf{89}, 3269–3285 (2009).
\href{https://doi.org/10.1080/14786430903260727}{doi:10.1080/14786430903260727}
```

Two house-style choices to make once and apply uniformly:
- whether to put DOI as a clickable `\href` (RevTeX has `hyperref` already
  loaded in main.tex)
- whether to include the journal abbreviation (`Philos. Mag.`) or full
  name; whichever, be consistent.
