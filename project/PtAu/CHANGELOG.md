# Pt(Au) pipeline changelog

Decisions and progress log for the Pt(Au) GB-segregation pipeline in
`project/PtAu/`. Reverse chronological (newest first). Companion to
`README.md`: the README is the operational "how to run" reference; this
file is the "what was done and why" record.

The end goal of this pipeline is to build a 3D Pt(Au) polycrystal whose
per-site segregation-energy spectrum (ΔE_seg) matches Wagih 2020's
reference data well enough to pass a two-sample Kolmogorov–Smirnov test.
The project's internal pass bar is `p > 0.5` (treated as "spectrum-level
indistinguishable" — strong enough that we accept the two distributions
as the same). Once that's validated, this folder is handed off for
downstream analysis.

This is a parallel scaffold of the existing Al(Mg) pipeline in
`project/`. Everything alloy-specific lives in `project/PtAu/`;
alloy-agnostic Python scripts in `project/scripts/` are reused verbatim.

---

## Current status (2026-05-13, 11:00 CEST)

| Stage | Status |
|---|---|
| Polycrystal generated | done — `poly_Pt_100A_8g.lmp` (62,096 atoms, 100³ Å, 8 grains) at `/cluster/scratch/cainiu/prototype_PtAu_100A/` |
| Anneal | done — SLURM job `66279426` finished 2026-05-13 00:14 (1 h wall, T_HOLD=816 K stable, force_max 3.34 eV/Å, 62,096 atoms preserved) |
| GB identification | done — 23,272 GB atoms (f_GB = 0.375; high vs Wagih's 19% because our 100³ Å box has higher GB-to-bulk surface ratio than Wagih's 200³ Å) |
| Per-site ΔE_seg sampling (n=500) | RUNNING — SLURM job `66391849` (started 2026-05-13 10:32 CEST, ~15 min ETA, 4 h budget) |
| Skew-normal (μ, σ, α) fit | pending — auto-runs after sampling lands |
| Wagih reference extraction from tar | done — `Pt_Au_20nm_GB_segregation.dump` (24 MB, 508,951 atoms, 97,440 GB sites) extracted at `/cluster/scratch/.../accelerated_model/Pt/Au_2017--OBrien.../` |
| KS test vs Wagih per-site array | pending — `compare_vs_wagih_PtAu.py` + `bootstrap_vs_wagih_PtAu.py` written and syntax-validated |

**Queue note.** The anneal cleared its queue at 2026-05-12 23:14 CEST after the Al(Mg) `X_c=0.30 fdseed_resume` job (66261260) released its 16-rank slot. Sampling job 66391849 launched immediately on 2026-05-13 10:32 CEST when CPU was free.

### Wagih database has BOTH Pt(Au) AND Au(Pt) for our EAM (correction)

A full `tar -tjf | grep -iE "Pt.*Au|Au.*Pt"` on 2026-05-13 returned
30+ paths covering BOTH directions for the O'Brien 2017 EAM (also for
the older Foiles/Adams EAMs). The 2026-05-12 conclusion that "the
database is one-directional, Pt(Au) only" was based on a partial grep
that missed the `Au/...` entries. **Practical impact: zero — Pt(Au)
remains a fine choice, and the KS comparison is on the matched
direction's reference. But the Au(Pt)-also-viable fact is recorded
here so future-me does not propagate the wrong premise.** See
[[reference_wagih_zenodo_layout]] memory (also corrected).

---

## 2026-05-12 — Pipeline scaffolded

### Trigger

Brief from the user (translated): build a parallel Pt(Au) pipeline on
top of the existing Al(Mg) work, produce a Wagih-comparable polycrystal,
KS-validate the per-site ΔE_seg distribution, then hand off to the
colleague for the next steps.

### Direction flipped: Au(Pt) → Pt(Au)

Initial plan was Au-host with Pt-solute, on the size-mismatch heuristic
(Au is the larger, softer atom — "big host, small solute" is the
conventional default for size-driven segregation). Two checks reversed
that decision:

1. **Potential.** The only EAM file available for this binary on NIST
   is O'Brien et al. 2017 `PtAu.eam.alloy` (downloaded into
   `data/potentials/`). Its header reads `2 Pt Au` with Pt equilibrium
   lattice 3.9764 Å and Au equilibrium lattice 4.1537 Å. The potential
   supports both directions.

2. **Wagih reference availability.** Wagih 2020's polycrystal database
   stores per-site ΔE_seg only in one host/solute direction per (alloy,
   potential) pair. A grep of the local archive
   `/cluster/scratch/cainiu/wagih_zenodo/learning_segregation_energies.tar.bz2`
   for any Pt or Au path returned 30 hits — **all** of them under
   `Pt/...`. The specific match for our EAM is:
   ```
   learning_segregation_energies/segregation_spectra_database_accelerated_model/Pt/Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/Pt_Au_20nm_GB_segregation.dump
   ```
   A confirming second sweep showed **no `Au/` host directory exists at
   all** in this archive — the database is genuinely one-directional
   for this potential.

If we had kept Au(Pt), we would have had no per-site Wagih reference
for the KS step. Flipping to Pt(Au) makes the KS comparison an
apples-to-apples per-site match using the same potential file on both
sides. After the flip:

- Folder name `AuPt/` → `PtAu/` (host-first, matching `AlMg/` style).
- Polycrystal lattice changed 4.1537 Å (Au) → 3.9764 Å (Pt).
- LAMMPS deck and Python sampler default elements/masses updated to
  `EL1=Pt, EL2=Au` (host first).

**Lesson for adding the next alloy:** before picking a host direction,
grep the Wagih tar for the specific (alloy, potential) match. The
database is one-directional per pair, and the host direction Wagih chose
is not always the size-mismatch default.

### Why `pair_style eam/alloy` here vs `eam/fs` in Al(Mg)

Both styles are EAM (Embedded Atom Method) in LAMMPS, with the same
total energy form:

    E_total = Σ_i F_i(ρ_i) + ½ Σ_{i,j} φ_ij(r_ij)

(embedding energy F_i of atom i in background electron density ρ_i,
plus a pair term φ_ij). The only structural difference is how ρ_i is
built up from neighbors:

| pair_style | ρ_j contribution depends on | density functions per system | file extension |
|---|---|---|---|
| `eam/alloy` (Daw–Baskes setfl) | element of neighbor j only | one ρ(r) per element | `.eam.alloy` |
| `eam/fs` (Finnis–Sinclair extended setfl) | elements of **both** i and j | one ρ(r) per (i, j) pair | `.eam.fs` |

`eam/fs` has an extra degree of freedom (asymmetric A→B vs B→A
densities); `eam/alloy` collapses that to a single per-element function.
For single-element systems they are mathematically equivalent — the
distinction matters only in alloys.

Etymology aside: Daw & Baskes proposed EAM in 1984 from effective-medium
theory; Finnis & Sinclair independently proposed the FS form the same
year from tight-binding 2nd-moment arguments, with the embedding
function constrained to `F(ρ) ∝ −√ρ`. LAMMPS exposes both as separate
`pair_style` strings because the file formats differ — `eam/fs` carries
an extra per-pair ρ column that `eam/alloy` does not.

**Consequence for this pipeline.** Mendelev 2009 Al(Mg) was fit in FS
form (`Al-Mg.eam.fs`); O'Brien 2017 Pt(Au) was fit in standard setfl
(`PtAu.eam.alloy`). LAMMPS will refuse to read one with the other's
pair_style — file column counts don't match. So the Pt(Au) LAMMPS deck
and the Pt(Au) Python sampler both need to declare `pair_style
eam/alloy`. The rest of the pipeline is alloy-agnostic.

The two lines where the pair_style is hardcoded:

- `data/decks/anneal_PtAu.lammps:70` — was line 61 in
  `project/data/decks/anneal_AlMg.lammps` with `eam/fs`.
- `scripts/sample_delta_e_PtAu.py:466` — was line 454 in
  `project/scripts/sample_delta_e.py` with `eam/fs`.

That is the entire substantive code delta between the two pipelines.
Everything else is parameter defaults and rename plumbing.

### Canonical scripts not modified in place

The project convention is that canonical generic scripts under
`project/scripts/` and `project/data/decks/` stay untouched when a
second alloy is added; behavior-extending changes go in renamed copies
inside the per-alloy folder. The reason is that the Al(Mg) pipeline is
still being actively run, and in-place edits to canonical files have
caused regressions before. So only the two files that hardcode the
pair_style were copied:

- `project/data/decks/anneal_AlMg.lammps` → `project/PtAu/data/decks/anneal_PtAu.lammps`
- `project/scripts/sample_delta_e.py`     → `project/PtAu/scripts/sample_delta_e_PtAu.py`

Every other Python script in `project/scripts/` already takes the alloy
parameters (elements, masses, lattice) through CLI flags, so it works
for Pt(Au) without modification.

### Files added (this folder)

```
project/PtAu/
├── README.md                                ← operational "how to run"
├── CHANGELOG.md                              ← this file
├── data/
│   ├── potentials/
│   │   └── PtAu.eam.alloy                    ← O'Brien 2017 (NIST), 875 KB
│   └── decks/
│       ├── anneal_PtAu.lammps                ← copy of anneal_AlMg.lammps
│       ├── submit_anneal_PtAu_100A.sh        ← 16-rank SLURM, 6 h budget
│       └── submit_delta_e_PtAu_100A.sh       ← n=500 sampling SLURM submit
├── scripts/
│   └── sample_delta_e_PtAu.py                ← copy of generic sampler, eam/alloy
└── output/                                   ← (empty; fits + figures land here)
```

Generic scripts in `project/scripts/` used verbatim (alloy-agnostic via
CLI flags — no copy needed):

- `generate_polycrystal.py` — Voronoi-based FCC polycrystal builder
- `gb_identify.py` — adaptive-CNA-based GB-atom mask
- `fit_delta_e_spectrum.py` — skew-normal (μ, σ, α) fit
- `fermi_dirac_predict.py` — temperature-dependent solute-occupation curves
- `compare_vs_wagih.py` — KS test driver (currently Al-Mg-specific; needs Pt-Au dump adapter)
- `bootstrap_vs_wagih.py` — bootstrap CI for the KS statistic

### File-by-file diff vs the Al(Mg) originals (audit trail)

| File | Change vs Al(Mg) original |
|---|---|
| `data/decks/anneal_PtAu.lammps` | `pair_style eam/fs` → `eam/alloy`; default elements `EL1=Al → Pt`, `EL2=Mg → Au`; masses `26.9815 → 195.0900` and `24.3050 → 196.9665`; `T_HOLD 373 → 816 K` (0.4 × T_melt(Pt) where T_melt(Pt) = 2041 K, CRC); potential file path |
| `scripts/sample_delta_e_PtAu.py` | `pair_style eam/fs` → `eam/alloy` (line 466); `_DEFAULT_ELEMENTS = ("Al","Mg") → ("Pt","Au")`; `_DEFAULT_MASSES = (26.9815, 24.3050) → (195.0900, 196.9665)`; docstring + CLI defaults updated. `--help` syntax-verified. |
| `data/decks/submit_anneal_PtAu_100A.sh` | scratch dir `prototype_AlMg_100A → prototype_PtAu_100A`; deck path; T_HOLD; EL1/EL2/MASS1/MASS2 passed explicitly; wall budget 4 h → 6 h (higher T_hold means more cool ps after the hold) |
| `data/decks/submit_delta_e_PtAu_100A.sh` | scratch dir; driver path; `--elements`/`--masses` passed explicitly; the auto-fallback `gb_identify.py` invocation uses `--lattice-a 3.9764` (Pt) instead of `4.05` (Al) |

### Polycrystal generated and validated

```
python project/scripts/generate_polycrystal.py \
    --structure fcc --box 100 --grains 8 --lattice-a 3.9764 \
    --structure-seed 1 --types 2 --out poly_Pt_100A_8g.lmp
```

Output at `/cluster/scratch/cainiu/prototype_PtAu_100A/poly_Pt_100A_8g.lmp`:

- 62,096 atoms in a 100 × 100 × 100 Å box, 8 grains, seed = 1.
- 2.4% deficit from GB close-pair removal (vs 1.6% for the Al(Mg)
  100 Å prototype — slightly higher because the smaller Pt lattice
  constant gives denser packing, so more close pairs get flagged).
- Type 1 only at this stage (host atoms; solute placement happens
  inside the ΔE_seg sampler, not the polycrystal generator).
- x/y/z all fully span the box; density 0.0621 Å⁻³ (ideal 0.0636 = 4/a³).

### Anneal submitted (job 66279426)

```
sbatch project/PtAu/data/decks/submit_anneal_PtAu_100A.sh
```

→ job `66279426`, name `anneal_PtAu_100A`, status PENDING. Waiting for
one of three RUNNING Al(Mg) HMC jobs to release a 16-rank slot (see
"Queue note" in the status table at the top).

### Pending follow-ups for the colleague

Steps 1–4 below are wall-time-bound execution. Operational commands for
each step are in `README.md`; only what to watch out for is repeated here.

1. **Anneal completion → ΔE_seg sampling.** Once job `66279426` writes
   `poly_Pt_100A_8g_annealed.lmp`, submit
   `submit_delta_e_PtAu_100A.sh`. The submit script auto-runs
   `gb_identify.py --lattice-a 3.9764` if the GB mask is missing, then
   does n=500 per-site sampling. Output:
   `delta_e_results_n500_PtAu_100A_tight.npz`.

2. **Fit + predict.** Run `project/scripts/fit_delta_e_spectrum.py` on
   the npz for the skew-normal `(μ, σ, α)` fit, then
   `project/scripts/fermi_dirac_predict.py` for the temperature curves.
   Outputs land in `project/PtAu/output/`.

3. **Wagih reference extraction.** Surgical extract from the local tar
   (~10–50 MB out of the 3.8 GB archive — fast):
   ```
   cd /cluster/scratch/$USER/wagih_zenodo
   tar -xjf learning_segregation_energies.tar.bz2 \
       learning_segregation_energies/segregation_spectra_database_accelerated_model/Pt/Au_2017--OBrien-C-J-Barr-C-M-Price-P-M-et-al--Pt-Au/Pt_Au_20nm_GB_segregation.dump
   ```

4. **KS test.** `project/scripts/compare_vs_wagih.py` currently parses
   the Al-Mg-specific dump format; the Pt-Au dump has the same Wagih
   schema (per-site `E_GB` column), so the adapter is small. Then
   `scipy.stats.ks_2samp` on the two per-site arrays. Pass bar:
   `p > 0.5`.

5. **Sanity-check `(μ, σ, α)` against Wagih's SI.** Wagih 2020 SI
   reports `(μ, σ, α)` per host-solute panel in Fig 4 and Fig 9; the
   superscript on each panel label indicates the SI reference number
   for the EAM used (`Au^N` = SI ref N). For O'Brien 2017 those labels
   should match the same potential we used here. Cross-check our fit
   against that panel — if `(μ, σ, α)` agree but the KS p-value fails,
   it points to a tail mismatch (likely a sampling/finite-size issue,
   not a physics issue).

### Git provenance

Two commits on branch `cainiu` (not pushed):

- `f42ab35` — initial scaffold: potential, deck, submit scripts,
  sampler copy, README, polycrystal generated.
- `bda3dcd` — expanded the Au(Pt) → Pt(Au) direction-flip narrative
  and added the eam/alloy vs eam/fs note. Originally written into
  the main project CHANGELOG; relocated into this file for handoff
  clarity.
