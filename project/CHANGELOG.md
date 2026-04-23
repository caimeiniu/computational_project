# Changelog

Entries in reverse chronological order (newest first).

## 2026-04-23 (evening) — A-pipeline end-to-end; Phase 3 sampling driver

Stood up the full "polycrystal → Wagih-anneal → GB mask → per-site ΔE"
chain on the 10³ nm³ / 8-grain Al prototype. First Phase 3 sample
completed end-to-end (job `64589990`, **3.5 min wall on 16 cores**).
Results at `/cluster/scratch/cainiu/prototype_AlMg_100A/delta_e_results.npz`
(not committed — `*.npz` gitignored).

**First ΔE_seg prototype spectrum (N_GB = 50, N_bulk_ref = 5, seed = 42)**:
- Bulk reference `E_bulk^Mg = −200223.604 eV`, std 0.009 eV across
  5 refs — baseline noise is ~1 meV, well below the segregation signal.
- CG convergence: all 55 sites hit `energy tolerance` (clean).
- **ΔE_seg range [−45.6, +19.6] kJ/mol**, mean −6.1, median −7.0.
- Sign and negative-skew match expectation (Mg segregates to Al GBs).
- Narrower than Wagih's Al(Mg) range [−60, +40] — expected for a
  50-site undersample vs Wagih's 10⁴; scaling to N_GB ≥ 500 should
  widen both tails. This is the first ΔE number on Al-Mg from our
  own pipeline and establishes the order of magnitude.

### Mendelev potential year correction (2014 → 2009)

Previous entries labelled the Wagih-default Al-Mg EAM as "Mendelev 2014".
The NIST entry actually cites Mendelev, Asta, Rahman & Hoyt, *Philos. Mag.*
89, 3269–3285 (**2009**). The downloaded file
`project/data/potentials/Al-Mg.eam.fs` (2.3 MB, `pair_style eam/fs`,
elements `Al Mg`) is the correct one — the year in the prior CHANGELOG
was a transcription error, not a different potential.

### Wagih-style anneal (`anneal_AlMg.lammps` + `submit_anneal.sh`)

Full protocol per `docs/paper_notes.md §1`:
1. CG #0 (absorb close-pair strain from fresh Voronoi construction)
2. NVT ramp 1 K → T_hold over 5 ps
3. **NPT hold at T_hold for 250 ps @ 0 bar**
4. NPT cool T_hold → 1 K at 3 K/ps
5. Final CG under `fix box/relax iso 0.0`

**T_hold = 373 K ≈ 0.4 × T_melt_Al (933 K)**. Wagih allows 0.3–0.5 ×
T_melt; 0.4 is the middle. 0.3 risks insufficient GB mobility; 0.5 risks
grain growth at the 100 Å prototype scale (finite-size). Deck documents
both extremes so 280 K / 467 K reruns are a one-`-var` change.

Tuning choices vs an earlier tighter draft:
- CG tolerances `1e-6/1e-8` (initial) and `1e-8/1e-10` (final) — loosened
  from `1e-10/1e-12` to stay inside the wall-time budget; Wagih is silent
  on tolerance and `1e-8` is industry standard for 0 K static relax.
- `neighbor 1.0 bin; neigh_modify every 10` — 5–10% faster at 300+ K MD
  than LAMMPS metal defaults, no correctness loss.
- `restart 50000 *.rst1 *.rst2` — double-buffered checkpoint every 50 ps
  so a killed job is resumable.

**Prototype result (job 64567957, 85 min wall on 16 cores)**:
- N_atoms = 59 224, box 100.77 Å at 373 K → 100.35 Å at 0 K (thermal
  contraction 0.4 %)
- PE −3.144 → −3.358 eV/atom (approaches Mendelev Al bulk −3.36; residual
  ~0.03 eV/atom is GB excess energy, expected)
- `poly_Al_100A_8g_annealed.lmp` is now the reusable artifact for every
  downstream phase.

### GB identification (`scripts/gb_identify.py`)

Pure numpy module, shells out to LAMMPS `compute cna/atom` once — no OVITO
/ ASE dependency so teammates can `import gb_identify` straight from
`myenv`. Returns `(mask, info)`; CLI writes `gb_mask.npy`, `gb_info.json`,
and a position-carrying `gb_cna.dump` (OVITO can colour-by-CNA without
reloading the data file).

**Not adaptive CNA**: LAMMPS' `compute cna/atom` is the conventional
fixed-cutoff variant (we use `0.854 a_fcc = 3.459 Å` for Al). For the
bulk/GB binary we need, this is equivalent to OVITO's a-CNA; for future
fine-grained GB-character analysis we'd switch to OVITO's
`CommonNeighborAnalysisModifier` or LAMMPS `compute ptm/atom`.

**Strict Wagih semantics**: bulk = parent structure only. In our FCC Al
system, any HCP-labelled atom (stacking fault, **SF**) gets flagged as
GB. On the annealed prototype this is 562 atoms (≈ 3% of n_gb), expected
given Al's low SF energy, and quantitatively negligible for downstream
ΔE_seg statistics. Docstring explicitly documents this semantic.

**Prototype result on annealed polycrystal**: `f_gb = 28.7 %`
(42 232 FCC / 16 430 Other / 562 HCP). Landed at the top of the expected
20–30 % window; 10 nm prototype grains give higher f_gb than Wagih's
20 nm / ~15 % due to finite-size effects, not a bug. Pre-anneal was
38.9 %; the 10-point drop confirms the anneal protocol relaxed the
Voronoi geometric distortion as intended.

### Phase 3 ΔE_seg sampling (`scripts/sample_delta_e.py` + `submit_delta_e.sh`)

Python driver + embedded LAMMPS deck template. For each sampled site:
one LAMMPS process does `read_data annealed → set atom <id> type 2 → CG
→ print pe`. One-process-per-site (startup overhead ~2 s × 55 sites
≈ 2 min, negligible next to CG) to keep state cleanly isolated between
substitutions — avoids snapshot/restore complexity of a single-session
loop. CG tolerances `1e-8/1e-10`, tighter than the anneal — ΔE_seg can
be < 1 kJ/mol so sub-meV numerics matter.

**Bulk reference divergence from Wagih**: Wagih uses "solute in a 6 nm
sphere of pure solvent". Our 100 Å prototype box cannot fit a 6 nm
sphere, so we take the **mean of 5 bulk-ref sites** chosen to sit ≥ 8 Å
(≈ 2 FCC NN shells) from any GB atom. σ/√5 of E_bulk gives ~0.4 kJ/mol
baseline uncertainty — small against Al(Mg) ΔE range [−60, +40] kJ/mol.
If/when we move to 20³ nm³ production, we can switch to the 6 nm-sphere
reference exactly.

**Checkpointing**: per-site results stream to `_results.csv` in work_dir
as soon as each LAMMPS finishes; `_run_meta.json` records
seed + sampled ids and verifies mismatches on resume. Prevents losing
40/55 sites if a SLURM job hits its time cap. Designed to scale to the
500+ sites we'll need for production skew-normal fit.

**`--mpi-cmd mpirun` (not srun) inside sbatch**: the driver fires one
LAMMPS per site, i.e. nested srun inside the outer SLURM allocation,
which on recent SLURM can trigger "job step creation temporarily
disabled". openmpi/4.1.6 `mpirun` integrates with SLURM via PMI and
avoids nested-step issues. The Python driver auto-detects (srun inside
SLURM, mpirun on login node) for interactive use.

### Housekeeping

- `project/data/potentials/` added with Mendelev 2009 `Al-Mg.eam.fs`
  (kept in-repo — stable, referenced exactly, teammates don't need to
  fetch from NIST).
- Stray SLURM `.out` / `.err` files (from the early anneal submit with
  relative `--output=%x-%j.out`) moved to scratch; both sbatch scripts
  now use absolute scratch paths; `.gitignore` gains
  `*-[0-9]*.{out,err}` + `slurm-*.{out,err}` as a safety net.
- No changes to the generator, UMA archive, or earlier decks.

## 2026-04-23 — 3D polycrystal generator implemented (FCC / BCC / HCP)

Added `scripts/generate_polycrystal.py` — parametric 3D Voronoi polycrystal
generator supporting all three lattice families covered by Wagih et al.
(FCC / BCC / HCP), reusable as a module (`build_polycrystal`, `write_lammps_data`)
or CLI. Replaces advisor's 2D columnar `create_nanocrystal.py`. Chose a pure
numpy/scipy implementation over Atomsk (next-steps item 1 from 2026-04-22
evening) so teammates can `import` it without an external tool dependency;
Atomsk remains a backup if we need exact Wagih reproducibility.

### Algorithm

1. Random grain centers in the periodic box + uniform SO(3) orientations per
   grain (scipy `Rotation.random`, which samples via quaternions with correct
   Haar measure — NOT naive Euler-angle sampling, which biases toward polar
   axes).
2. Per-structure crystal template (FCC conventional cell, BCC conventional
   cell, HCP primitive cell with hexagonal lattice vectors) sized to cover
   the box diagonal plus one lattice spacing after rotation.
3. Rotate + translate template to each grain center; keep atoms falling
   **inside** `[0, L)³` (no PBC wrap of the oversized template — wrapping
   would map multiple template atoms to the same position and double the
   density; an earlier draft hit this exact bug, 2× atom count on a 30 Å test).
4. Voronoi ownership via PBC-nearest grain center over 27 periodic images
   (KDTree built once outside the grain loop).
5. Close-pair removal at `NN_dist / 2` per structure (FCC `a/(2√2)`,
   BCC `a√3/4`, HCP `a/2`), matching the advisor's example logic.

### Safety checks

- **Minimum grain-center PBC separation** = `2 × lattice_a` (rejects literal
  overlaps; normal Poisson min separation is `~0.55 · L · n_grains^(-1/3)`,
  far looser).
- Perfect-lattice vs actual atom count printed so users can spot anomalies
  (> ~5% deficit = something wrong with Voronoi tessellation or close-pair
  cutoff).

### Validation

| System | Box | Grains | N_atoms | Ideal | Deficit | Mean NN | Ideal NN |
|--------|-----|--------|---------|-------|---------|---------|----------|
| Al/FCC | 60 Å | 8 | 12 691 | 13 006 | 2.4% | — | 2.864 |
| Al/FCC | 100 Å | 8 | 59 256 | 60 214 | 1.6% | 2.677 | 2.864 |
| Fe/BCC | 60 Å | 8 | 17 978 | 18 351 | 2.0% | 2.283 | 2.482 |
| Mg/HCP | 60 Å | 8 | 9 027 | 9 244 | 2.3% | 2.882 | 3.209 |

Prototype-scale FCC output matches the planned ~60 k atoms for 10³ nm³ Al.

### Known limitations / follow-ups

1. **No grain_id in LAMMPS output** — atoms carry only type 1. For 5D GB
   character analysis (misorientation × GB-plane normal) we'll need per-atom
   grain ownership. Options: sidecar `grain_ids.npy`, or switch to
   `atom_style molecular` and store grain_id in the molecule-ID field.
2. **No metadata sidecar** — should emit a JSON with `{box, structure,
   lattice_a/c, structure_seed, n_grains, centers, orientations_quat}` for
   reproducibility and downstream plotting.
3. **No misorientation-angle sanity check** — expected to follow the
   Mackenzie distribution (cubic symmetry, 0°–62.8°, peak at 45°) for FCC/BCC
   and Handscomb for HCP; small `n_grains` (= few `C(n, 2)` pairs) may
   visibly deviate and should be flagged.
4. **Close-pair removal is greedy by pair index**, not by inter-atom distance.
   Matches advisor's example; physically harmless but non-optimal. Could sort
   pairs by distance and drop the closer conflict first.
5. **`write_lammps_data` uses per-atom `f.write`** — ~5 s for 500 k atoms.
   Fine for one-shot use but could be vectorized with `np.savetxt`.
6. **`_crystal_template` sizing uses shortest lattice vector** — HCP over-
   generates ~2× along `c`; negligible at prototype scale.
7. **No unit tests yet.**

### Key concepts introduced (English terms for future reference)

- **Voronoi tessellation** — grain construction method; Wagih default.
- **SO(3) / Haar measure / Shoemake quaternion method** — correct uniform
  random rotation sampling.
- **Lattice parameter** `a` (and `c` for HCP); **conventional unit cell**
  (cubic for FCC/BCC) vs **primitive cell** (hexagonal for HCP).
- **Nearest-neighbor (NN) distance** and **coordination number** — FCC
  `a/√2` / 12, BCC `a√3/2` / 8, HCP `a` / 12.
- **Ideal c/a** = √(8/3) ≈ 1.633 for hard-sphere-packed HCP.
- **GB atom fraction** `f_gb ≈ 3t/d` (one-shell geometry estimate); higher
  orders come from triple-line/quadruple-point inclusion-exclusion
  corrections. Wagih measures ~15% via a-CNA on 20³ nm³ / 16-grain samples.
- **GB character (5D)** = 3 misorientation DOFs + 2 GB-plane-normal DOFs.
  2D columnar (advisor's original) only sampled a 2D slice of this space —
  one of the reasons we switched to 3D.
- **Mackenzie distribution** (cubic)/ **Handscomb** (hexagonal) —
  theoretical misorientation-angle distribution under random texture; used
  as sanity check on grain-orientation sampling.
- **Close pair** — construction artifact at GBs where two grains' rotated
  lattices place atoms < NN/2 apart; removed geometrically before any MD.
- **RNG seed separation**: `structure_seed` (this script) vs `solute_seed`
  (`set type/fraction`) vs `swap_seed` (`fix atom/swap`) vs `velocity_seed`
  (`velocity create`). Name them distinctly in all downstream LAMMPS decks.

### Impact on downstream work

- Equilibration deck (next-steps item 2) must handle a fresh Voronoi
  structure with large geometric distortion — the advisor's 10 ps NVT at
  300 K is insufficient; need Wagih's protocol (0.3–0.5 T_melt × 250 ps →
  3 K/ps cool → CG).
- `f_gb` target for prototype (10³ nm³, 8 grains, ~5 nm mean grain size) is
  ~20–30% (higher than Wagih's 15% at 20³ nm³ / 8 nm grains) — this is an
  expected finite-size effect, not a bug.
- Grain ID persistence becomes necessary once we move beyond scalar
  `X_GB(T, X_c)` curves to per-GB-character resolution (follow-up 1 above).

Covers next-steps items 1 and 3 from the 2026-04-22 evening entry.

## 2026-04-22 (evening) — Commit to 3D Voronoi and switch to Al-Mg

### Supersedes

- The 2026-04-22 afternoon decision to "start with 2D columnar Cu-Ni and scale to 3D later".
- HMC pipeline itself is unchanged (LAMMPS `fix atom/swap`, a-CNA GB identification,
  per-site ΔE, `(T, X_c)` scan). What's being replaced is the **structure generator**, the
  **equilibration protocol**, and the **alloy system**.

### 2D columnar → 3D Voronoi (committed)

The advisor's `create_nanocrystal.py` builds a 2D columnar polycrystal
(~5 × 200 × 173 Å, 4 grains rotated only around `[110]`, ~15 k atoms). Physical limits:

- GB character is restricted to a ~2D slice of the full 5D macroscopic GB character space
  — all boundaries are pure tilt around `[110]`; no twist, no mixed character.
- No triple points, no quadruple points (only triple *lines* along the ~5 Å thin x axis).
- ΔE histogram is artificially narrowed vs Wagih's 3D polycrystal; solute-solute g(r) at
  GB lacks triple-line / quadruple-point enrichment sites.
- Quantitative comparison to Wagih Fig 2 / Fig 5 is impossible because those are 3D
  (20³ nm³, 16 grains).

These all directly weaken the headline figure (HMC X_GB vs Fermi-Dirac X_GB): the
dilute-limit breakdown plausibly originates at GB heterogeneity / special sites absent from
columnar geometry.

**Decision:** skip the 2D pipeline-validation phase entirely and go straight to 3D. The
only things 2D would have validated that 3D would not (LAMMPS input-deck syntax,
`fix atom/swap` parameters) are trivial to check in isolation. The real risks — Voronoi
tessellation quality, Wagih-style annealing convergence, GB identification at scale — only
manifest in 3D and have to be debugged there anyway.

### Consequences of going 3D (what changes, what stays)

Reusable from advisor's example:
- `Cu_Ni_Fischer_2018.eam.alloy` potential file (only if we keep Cu-Ni — see below).
- HMC input-deck skeleton (`fix atom/swap 100 10 <seed> T ke yes types 1 2` is
  system-agnostic) and thermo/dump conventions.

Must be rebuilt:
- **Structure generator** — replace `create_nanocrystal.py` with Atomsk Voronoi
  tessellation (paper ref 71) or a custom `scipy.spatial.Voronoi` script. Atomsk is the
  default since it matches Wagih exactly.
- **Equilibration protocol** — replace the CG + 10 ps NVT at 300 K in
  `in_ncCuNi_equilibriate.lammps` with Wagih's protocol (paper_notes.md §1): anneal at
  0.3–0.5 T_melt for 250 ps, slow cool at 3 K/ps to 0 K, final CG. The 10 ps recipe is
  adequate for a columnar structure whose GBs have already been cleaned by `pdist`-based
  close-pair removal, but is insufficient for fresh Voronoi cells whose GBs start with
  large geometric distortion and high strain energy.

### Cu-Ni → Al-Mg

Cu-Ni was picked on 2026-04-22 afternoon specifically because the advisor shipped a
ready-to-run 2D structure generator + equilibration deck. With both being thrown out, the
only remaining Cu-Ni advantage is the Fischer 2018 EAM file — a single potential
downloadable from NIST in seconds. The tradeoff flips:

| Criterion | Cu-Ni | Al-Mg |
|-----------|-------|-------|
| Wagih paper coverage | Table 1 row + Fig 4/5 | Fig 2 headline + MAE benchmark |
| Spectrum params | μ ≈ −2, σ ≈ 8 kJ/mol | μ ≈ −2, σ ≈ 4, α ≈ −0.4 kJ/mol |
| ΔE range | ~60 kJ/mol | ~100 kJ/mol |
| X_GB at 5% X_c | ~15% | ~30% |
| Atoms in 20³ nm³ | ~680 k (a=3.61 Å) | ~481 k (a=4.05 Å, ~30% cheaper) |
| NIST EAM options | Fischer 2018 | Mendelev 2014, Mishin, Liu-Adams-Wolfer |
| Zenodo per-site ΔE | Ni(Cu) present | Al(Mg) present, directly matches Fig 2 |

**Decisive factor — signal strength:** the headline figure is X_GB^HMC vs X_GB^FD as a
function of X_c; the breakdown is a *deviation* of X_GB^HMC from the Fermi-Dirac prediction.
At 5% X_c, Al(Mg) reaches ~30% X_GB vs ~15% for Ni(Cu), giving 2× the baseline and 2× the
dynamic range for the deviation to be detected above HMC statistical noise. Al-Mg also
maps onto Wagih Fig 2 at the level of fitted skew-normal parameters `(μ, σ, α)`, not just
a scalar MAE — much more stringent comparison.

### Revised scale plan

| Stage | Box | Grains | Atoms (Al) | HMC / (T, X_c) (EAM, 16–32 core MPI) |
|-------|-----|--------|------------|-------------------------------------|
| Prototype | 10³ nm³ | 4–8 | ~60 k | ~10 min |
| Production | 20³ nm³ | 16 | ~481 k | ~1–2 h |

Prototype stage exists to debug Voronoi quality, the Wagih-style annealing protocol, and
GB identification on a system small enough that iteration is fast. Production mirrors
Wagih exactly for direct figure-by-figure comparison.

### Housekeeping

- `project/data/examples/cra_example/` removed (W CRA Frenkel-insertion template; its only
  pedagogical value was the hand-written MC loop pattern, which `fix atom/swap` replaces).

### Next steps (revised)

1. Install / locate Atomsk on Euler (check module availability first).
2. Download Al-Mg EAM potential from NIST (Mendelev 2014, the Wagih default).
3. Write `scripts/generate_polycrystal.py` — Atomsk Voronoi → LAMMPS data file,
   parameterized by box size, grain count, lattice parameter.
4. Extend equilibration deck to Wagih-style anneal
   (0.3–0.5 T_melt × 250 ps → 3 K/ps cool to 0 K → CG).
5. Run the prototype end-to-end (10³ nm³): structure → anneal → a-CNA GB ID →
   visualize in OVITO → confirm GB site count and GB fraction are sensible
   (expected f_gb ≈ 10–20% for 10 nm grain size).

## 2026-04-22 (afternoon) — Cu-Ni starting system + concrete HMC plan

### Example code review (from advisor)

Two reference archives received and extracted under `project/data/examples/`:

**`cra_example/`** — W bulk Frenkel-insertion (CRA) simulation, NOT GB segregation.
Irrelevant in physics, but useful as a LAMMPS scripting template:
- variable/loop/minimize structure for iterative MC-like operations
- `fix box/relax` for zero-stress relaxation
- random atom selection via LAMMPS `random()` variable

**`nc_swap_CuNi/`** — directly aligned with our project. Complete 3-step pipeline:
1. `create_nanocrystal.py` — 2D **columnar** FCC nanocrystal, 4 grains, `[110]×[112]×[111]`
   orientation, ~5 × 200 × 173 Å box, ~14,782 Cu atoms. Grains rotated by random
   angles around x. Removes too-close atom pairs at GBs via `scipy.spatial.pdist`.
2. `in_ncCuNi_equilibriate.lammps` — loads pure Cu NC, uses
   `set group all type/fraction 2 0.025` to convert 2.5% of atoms to Ni, CG-minimizes,
   then NVT 300 K for 10 ps. Outputs `initial_lattice_300K.lammps`.
3. `in_ncCuNi_hybrid_md-mc.lammps` — **HMC loop** implemented via LAMMPS built-in
   `fix atom/swap 100 10 <seed> 300.0 ke yes types 1 2` (every 100 MD steps, attempt
   10 Metropolis type 1↔2 swaps at T=300 K).

### Key simplifier: `fix atom/swap`

LAMMPS has a built-in semi-grand-canonical swap fix. We do **not** need to write our
own Python MC loop — `fix atom/swap` interleaves Metropolis swaps with NVT MD natively.
This collapses the previous "MC simulation script skeleton" step into configuring one
LAMMPS input line.

### Scale: columnar vs 3D

The example uses ~15 k atoms (2D columnar, 4 grains). Wagih's paper uses ~500 k atoms
(3D Voronoi, 16 grains, 20³ nm³). Columnar is ~30× cheaper per HMC step and sufficient
for validating the pipeline and mapping the X_GB vs X_c curve. 3D may be needed later
for publication-grade statistics on the ΔE spectrum, but is not required for the
dilute-limit breakdown result.

### System choice: Cu-Ni first, Al-Mg later (if time allows)

**Revised from earlier today**: start with **Cu-Ni** (not Al-Mg) because the advisor's
example includes:
- working structure generation script (columnar NC)
- validated `Cu_Ni_Fischer_2018.eam.alloy` potential (Fischer et al., Acta Mater. 2019)
- complete equilibration + HMC input decks

Al-Mg remains the preferred target for direct Wagih Fig 2/5 comparison, but Cu-Ni gets
us to a working HMC result fastest. Al-Mg becomes a Phase B goal after Cu-Ni validates.

**Trade-off accepted**: Wagih's Ni(Cu) data is less featured in the paper (μ ≈ -2,
σ ≈ 8, X_GB ≈ 15 % at 5 % total) than Al(Mg), but the Zenodo dataset should still have
per-site ΔE for Ni(Cu) for direct spectrum comparison.

### Execution plan (phases)

**Phase 0 — Euler environment** (day 0)
- Set up `/cluster/scratch/cainiu/nc_CuNi_HMC/` (home 50 GB limit)
- Load modules: `stack/2024-06`, `openmpi/4.1.6`, `lammps/<version>` (confirm on Euler)
- Verify `lmp -h` and `lmp -i in_ncCuNi_equilibriate.lammps -l log.test` runs

**Phase 1 — Reproduce example exactly** (day 1–2)
- Copy `nc_swap_CuNi/` to scratch, run all three steps unchanged
- Compare our `log_mcmd.lammps` against the supplied one (PE trajectory, swap stats)
- Visualize `dump_hybrid_md-mc.atom` in OVITO — confirm 4-grain columnar structure and
  visible Ni segregation to GBs

**Phase 2 — GB site identification** (day 2–3)
- Apply a-CNA via OVITO Python or `ase.spacegroup.symmetrize` to
  `initial_lattice_relaxed.lammps` (pure Cu, before adding Ni) — get boolean GB mask
- Record GB site fraction `f_gb` and number of GB sites N_GB

**Phase 3 — Per-site ΔE spectrum** (day 3–5)
- For each GB site (sample ≥200 sites if N_GB is large):
  - Substitute 1 Cu → Ni, CG relax, record ΔE_seg = E_GB^Ni − E_bulk^Ni
  - Reference E_bulk^Ni: 1 Ni substituted at an interior bulk site, CG relax
- Build ΔE histogram, fit skew-normal `F(ΔE) = (μ, σ, α)`
- Compare against Wagih's Ni(Cu) spectrum

**Phase 4 — HMC concentration scan** (week 2)
- X_c ∈ {0.005, 0.01, 0.02, 0.05, 0.10, 0.20} at T = 300 K
- Per run: equilibrate (Phase B style, 10 ps NVT) → HMC (100 ps, `atom/swap` every 100
  steps × 10 attempts) → measure X_GB from final frame averaged over last 10 ps
- Track swap acceptance rate and PE convergence to set run length adaptively

**Phase 5 — Compare to Fermi-Dirac prediction** (week 2–3)
- For each X_c, compute
  `X_GB^FD(X_c, T) = (1/N_GB) Σ_i 1/[1 + ((1-X_c)/X_c) exp(ΔE_i/kT)]`
  using Phase 3 ΔE spectrum
- Plot X_GB^HMC vs X_GB^FD — define breakdown concentration as first X_c where
  disagreement exceeds ML noise floor (~4 kJ/mol equivalent) or a fixed ratio threshold

**Phase 6 — Diagnostics** (week 3, if breakdown observed)
- Solute-solute g(r) at GB atoms
- Site occupation `P_i` vs local Ni density — look for cooperative / anti-cooperative
  signatures
- Temperature axis: repeat Phase 4 at T = 600 K to probe finite-T contribution

### Git housekeeping

- `*.atom` added to `.gitignore` (LAMMPS dump trajectories are large and regenerable)
- Original zip archives `cra_example.zip` and `nc_swap_simulation_CuNi.zip` removed —
  content is now unpacked under `project/data/examples/`, git history is authoritative
- `project/docs/project presentation.pptx` added — documents current scope and is a
  shared artifact for the team
- UMA-related files moved to `{docs,scripts}/archive/` (not deleted — kept as history)

## 2026-04-22 — Scope narrowed: drop UMA, single focus on dilute-limit breakdown

### Direction change

The project is narrowed to **one** central scientific question:

> At what solute concentration does Wagih's independent-site (dilute-limit) assumption
> break down, and what is the physical mechanism of the breakdown?

**Removed from scope:**
- UMA MLIP integration (previously Phase 2 / Extension ⑦ / fig3)
- "Systematic temperature effects" as a standalone extension

**Still in scope:**
- HMC scan on a `(T, X_c)` grid for the chosen binary system
- Comparison of `X_GB(T, X_c)` from HMC vs Fermi-Dirac prediction from 0 K ΔE spectrum
- Solute-solute spatial analysis (`g(r)`, local-density vs `P_i` correlation) to diagnose
  where and why the framework fails

### Why temperature is no longer a separate direction

Advisor's point: HMC is by construction a finite-T simulation — it samples the Boltzmann
distribution at temperature T with all degrees of freedom (vibrations, thermal expansion,
GB relaxation, solute-solute interactions). Comparing HMC against the paper's
Fermi-Dirac formula at multiple `(T, X_c)` automatically captures:
- vibrational / anharmonic corrections to ΔE (what the previous "temperature extension"
  was going to test)
- solute-solute interaction effects (the dilute-limit breakdown question)

The two effects can be separated cleanly at the experiment level: at `X_c → 0`, any
HMC/FD disagreement can only come from finite-T corrections, so the dilute end of the
concentration sweep isolates the temperature effect automatically. A separate
"temperature-only" sweep would reproduce work the grid scan already does.

**Conclusion:** Temperature is an axis of the core experiment, not a separate extension.

### Why UMA is dropped

Captured in the 2026-03-24 (evening) UMA feasibility analysis above, restated here:

- UMA + Wagih's accelerated ML (SOAP/PCA/k-means → 100 pts → 10⁵ sites) adds ~4 kJ/mol
  ML noise on top of any potential-accuracy difference, diluting UMA's benefit.
- Direct UMA relaxation on a full 20³ nm³ polycrystal (~480 K atoms) is ~7.5 days per
  100 sites on GPU — borderline infeasible at project scale.
- UMA's true value (alloys with no EAM, multi-component / HEA) is far from the core
  scientific question of this project.
- Keeping two loosely-coupled threads (HMC verification + UMA integration) dilutes the
  headline result. One sharp question is stronger than two weak ones.

### Choice of binary system (tentative): Al-Mg

**Starting system:** Al(Mg), pending advisor confirmation.

Reasoning:
1. Direct reproducibility of Wagih Fig 2 & Fig 5 (Al(Mg) MAE = 2.4 kJ/mol high-fidelity,
   4.2 kJ/mol accelerated) — gives an immediate sanity check for the pipeline.
2. Multiple well-validated EAM potentials on NIST (Mendelev, Mishin, Liu-Adams-Wolfer);
   LAMMPS `pair_style eam/alloy` reads them directly.
3. Moderate segregation strength (ΔE range [-60, +40] kJ/mol, X_GB ≈ 30% at X_tot = 5%)
   — strong enough to observe cleanly, not so strong that the GB saturates and hides the
   dilute→concentrated crossover.
4. Cheap: FCC, a = 4.05 Å, ~480 K atoms in 20³ nm³, fast CG relaxation on Euler.
5. Wagih's Zenodo data (doi:10.5281/zenodo.4107058) contains Al(Mg) → direct comparison
   of ΔE spectra.

**Backup systems** if Al-Mg runs into potential / convergence issues:
- Cu-Ag (FCC, strong segregation, classic test case)
- Ni-Cu (weaker signal but full Wagih data available)

### Revised project roadmap

1. **Baseline pipeline (Al-Mg)** — polycrystal generation → a-CNA GB ID → per-site ΔE
   for a subset of GB sites → compare to Wagih spectrum.
2. **HMC on (T, X_c) grid** — measure X_GB and site-resolved P_i.
3. **Dilute-limit failure analysis** — fit X_GB(X_c) from HMC vs Fermi-Dirac, define
   breakdown concentration, quantify divergence.
4. **Physical diagnosis** — solute-solute g(r) at GB, local density vs P_i correlation.
5. **(Optional)** replicate on a second alloy (Cu-Ag) to check universality.

### Figures status

- `fig_divergence_schematic.{py,png}` — **promoted to core figure** (no longer
  supplementary); illustrates HMC vs Wagih divergence with increasing X_c.
- `fig3_uma_integration.png` — retired. Keep in history but excluded from future writeups.
- `fig1_paper_pipeline.png`, `fig2_hmc_pipeline.png`, `method_overview.png`,
  `fig_procedure.png` — unchanged, still used.

### Next steps

1. Confirm Al-Mg choice with advisor (and which EAM potential to use)
2. Get the nanocrystalline structure generation code from advisor
3. Set up Euler environment (`module load stack/2024-06 openmpi/4.1.6 lammps/...`)
   and confirm available LAMMPS version
4. Move large simulation outputs to `/cluster/scratch/cainiu/` (home limit 50 GB)
5. Run 10-site sanity check: random GB sites → place Mg → CG relax → ΔE in
   [-60, +40] kJ/mol? If yes, pipeline is working.

## 2026-03-24 (evening) — Environment + UMA Feasibility + Project Directions

### LAMMPS installed

- `conda install -c conda-forge lammps` in gb-seg environment
- Version: 2024-08-29, Python interface working
- Key packages: MANYBODY (EAM), MEAM, MC, VORONOI, ML-IAP, USER-MLIP
- EAM pair styles: eam, eam/alloy, eam/fs
- Minimize styles: cg, fire, sd

### UMA CPU/GPU benchmark (Al FCC supercells, uma-s-1p1)

**CPU benchmark:**

| N_atoms | Box (Å) | Time (s) | ms/atom |
|---------|---------|----------|---------|
| 27      | 12.2    | 0.26     | 9.6     |
| 125     | 20.3    | 0.90     | 7.2     |
| 512     | 32.4    | 7.0      | 13.6    |
| 1,000   | 40.5    | 14.3     | 14.3    |

**GPU benchmark (NVIDIA RTX A6000):**

| N_atoms | Box (Å) | Time (s) | ms/atom |
|---------|---------|----------|---------|
| 125     | 20.3    | 0.07     | 0.58    |
| 512     | 32.4    | 0.24     | 0.47    |
| 1,000   | 40.5    | 0.45     | 0.45    |
| 3,375   | 60.8    | 1.50     | 0.45    |

GPU ~32× faster than CPU at 1000 atoms. Scaling ~N^1.05 (near-linear).

**Extrapolation to polycrystal sizes (GPU):**

| System           | N_atoms | 1 force eval | 1 relaxation (~30 CG) | 100 sites |
|------------------|---------|-------------|----------------------|-----------|
| Al 10³ nm³       | 60K     | ~27 s       | ~13 min              | ~22 hr    |
| Al 20³ nm³       | 480K    | ~216 s      | ~1.8 hr              | ~7.5 days |
| Ni 20³ nm³       | 730K    | ~336 s      | ~2.8 hr              | ~12 days  |

**Conclusion: CPU completely infeasible. GPU marginal for 20³ nm³ (7.5 days for 100 sites).**

### UMA + ML acceleration: marginal value argument

If using Wagih's accelerated model (SOAP + PCA + k-means → 100 training points):

**Error chain analysis:**
- Potential accuracy (EAM vs UMA vs DFT): ~5–20 kJ/mol difference
- ML prediction error (100 pts → 10⁵ pts): ~4.2 kJ/mol MAE

The ML prediction step introduces ~4 kJ/mol noise regardless of whether training data
comes from EAM or UMA. **UMA's accuracy improvement is diluted by ML prediction error.**

Why not use DFT for the 100 points? Because each relaxation is done on the FULL
supercell (~480K atoms) — DFT is impossible at this scale. So the choice is only
EAM vs UMA, and the difference may be smaller than the ML noise floor.

**Verdict: UMA + Wagih's accelerated ML pipeline = not worth the extra computation.**

UMA is only valuable when:
1. Computing ALL sites directly (no ML), on a small system (5³–10³ nm³)
2. Predicting alloys where NO EAM potential exists
But small systems can't capture polycrystalline GB statistics adequately.

### Revised project directions — beyond HMC verification

The project should focus on **stress-testing Wagih's framework**, not on replacing EAM.

**① Dilute-limit breakdown boundary (⭐⭐⭐ core contribution)**
- Wagih assumes independent sites (no solute-solute interaction)
- MC can test this directly by varying total concentration:
  X_tot = 1%, 5%, 10%, 20% → compare X_GB^MC vs X_GB^Wagih
- Define a "critical concentration" above which the framework fails
- This is the paper's potential key figure

**② Solute-solute spatial analysis (⭐⭐⭐ pairs with ①)**
- When MC disagrees with Wagih, WHY?
- Measure: solute clustering at GB, solute-solute g(r), P_i correlation with local density
- Provides physical basis for correcting Wagih's model

**③ Temperature effects (⭐⭐ easy addition)**
- Run MC at multiple temperatures (300K–1200K)
- Compare X_GB(T) curves with Wagih's Fermi-Dirac prediction
- 0K ΔE_seg may not capture vibrational entropy effects at high T

**④ Active learning for training point selection (⭐⭐ methodological)**
- Replace k-means (unsupervised) with uncertainty-guided active learning
- Could achieve same accuracy with 50 points instead of 100
- Or better accuracy (MAE < 4.2 kJ/mol) with same 100 points

**⑤ Multi-configuration statistics (⭐ supplementary)**
- Generate 5–10 different Voronoi polycrystal configurations
- Check if segregation spectra are robust across configurations
- Quantify configuration-to-configuration variance

**Recommended project story:**
1. MC verification of Wagih at standard conditions (baseline)
2. Concentration sweep → dilute-limit failure boundary (core result)
3. Spatial analysis of solute-solute interactions (physical insight)
4. Temperature sweep (supplementary)
5. (Optional) ML improvement / multi-config statistics

### Next steps

1. **Get advisor's nanocrystalline structure generation code**
2. **Pick alloy system** (likely Al-Mg, pending advisor input)
3. **Download EAM potential** from NIST repository
4. **Download Wagih's data** from Zenodo (doi:10.5281/zenodo.4107058) for comparison
5. **Prepare MC simulation script skeleton** (LAMMPS Python interface)

## 2026-03-24 — Project initialization

- Created project directory `/scratch/cainiu/UMA/Computational/`
- Cloned conda environment `gb-seg` from `tds-mlip` (Python 3.11, torch, ase, fairchem, phonopy)
- Initialized git repo with `.gitignore`
- Created project structure: `scripts/`, `data/`, `output/`, `docs/`

### Background

Based on Wagih, Larsen & Schuh, *Nature Communications* 11:6376 (2020):
"Learning grain boundary segregation energy spectra in polycrystals"

### Goals

1. **Hybrid Monte Carlo (HMC) verification**: Use MC + energy minimization to verify
   solute segregation behavior in a binary alloy polycrystal, reproducing key results
   from the paper (segregation enthalpy spectra, equilibrium GB solute concentration).
2. **UMA integration**: Replace classical interatomic potentials with UMA MLIP to
   evaluate whether universal ML potentials can improve segregation energy predictions
   without system-specific potential fitting.

### Method summary (from paper)

- Generate base-metal polycrystal: 20×20×20 nm³, 16 grains, Voronoi tessellation, thermally annealed
- Identify GB sites via common neighbor analysis (CNA)
- For each GB site: place solute atom, relax, compute ΔE_seg = E_GB^solute - E_bulk^solute
- Feature extraction: SOAP descriptors (r_cutoff=6 Å, F^SOAP=1015 features)
- Learning: linear regression on SOAP features → segregation energy
- Accelerated model: PCA (10 components) + k-means clustering (100 training points)

### Project positioning — MC as verification, not method innovation

Advisor's suggestion: *"you could test them by modelling atomistically segregation
in 2D or 3D nanocrystalline structure"* — the key word is **test** (verify).

**Project logic**:
1. Wagih predicts equilibrium GB solute concentration using per-site ΔE + Fermi-Dirac
2. We run MC simulation (more physically complete) to check if predictions hold
3. If they agree → validates Wagih's framework for that system
4. If they disagree → solute-solute interactions / concentration effects matter → new physics
5. UMA extension: can a universal MLIP replace system-specific EAM potentials?

### Why Wagih does NOT use MC (and why we do)

| Aspect                    | MC simulation (ours)            | Wagih per-site + ML           |
|---------------------------|---------------------------------|-------------------------------|
| Goal | Verify predictions | Build alloy-wide database       |
| Cost per (T, X_tot) point | Days of MC                      | Formula evaluation (instant)  |
| Change T or concentration | Re-run MC from scratch          | Plug into Fermi-Dirac (free)  |
| Change alloy system       | Re-run MC                       | Re-compute 100 training points|
| Physics completeness      | Full (solute-solute interaction)| Approximate (dilute limit)    |
| Throughput                | Low (1 alloy at a time)         | High (259 alloys scanned)     |

**Wagih needs throughput** → per-site + ML is the right choice for scanning 259 alloys.
**We need ground truth** → MC is the right choice for verifying specific predictions.

### Key concept: segregation energy spectrum

The spectrum is the distribution of ΔE_seg across all ~10⁵ GB sites in a polycrystal.
It is an **intrinsic material property** independent of T and concentration.

Once the spectrum is known, equilibrium GB occupation at site i follows Fermi-Dirac:

```
P_i(T, X_c) = 1 / [1 + ((1 - X_c) / X_c) × exp(ΔE_i / kT)]
```

- T = temperature, X_c = bulk solute concentration
- Negative ΔE_i → favorable segregation → P_i close to 1
- Average GB concentration: X_GB = (1/N_GB) × Σ P_i

This is analogous to electron filling energy levels — ΔE_i plays the role of energy
levels, and the Fermi-Dirac function determines occupation at given T.

**Compute the spectrum once → predict any (T, X_tot) instantly. This is the power of
Wagih's approach, but it relies on the dilute-limit (non-interacting) assumption.**

### Key concept: accelerated ML model (100 training points)

Full computation: 10⁵ GB sites × LAMMPS relaxation each = very expensive.

Accelerated pipeline:
1. Compute SOAP descriptors for all 10⁵ sites (cheap, no LAMMPS)
2. PCA: 1015-dim SOAP → 10-dim (captures >99% variance)
3. k-means clustering in 10-dim space → 100 clusters
4. Run LAMMPS only for the 100 cluster centroids (representative sites)
5. Train linear regression on these 100 (SOAP_10d, ΔE) pairs
6. Predict ΔE for all 10⁵ sites using the trained model

Result: 100 LAMMPS calculations instead of 100,000 = **1000× speedup**,
MAE increases only from ~2.5 to ~4.2 kJ/mol (Al-Mg example).

### Atom count estimates

For a 20×20×20 nm³ FCC polycrystal: N = (V_box / a³) × 4

| Metal | a (Å) | Atoms in 20³ nm³ |
|-------|-------|------------------|
| Al    | 4.05  | ~481,000         |
| Ni    | 3.52  | ~734,000         |
| Cu    | 3.61  | ~680,000         |

GB fraction depends on grain size: ~10–20% for 10–15 nm grains.

### Figures created

- `docs/fig1_paper_pipeline.png` — Wagih's per-site ΔE + SOAP + ML pipeline
- `docs/fig2_hmc_pipeline.png` — Our MC swap verification approach
- `docs/fig3_uma_integration.png` — UMA MLIP as drop-in replacement for EAM
- `docs/method_overview.png` — 4-panel overview (Voronoi, a-CNA, ΔE_seg, MC swap)
- Script: `scripts/pipeline_figures.py`, `scripts/demo_polycrystal_2d.py`

Note: Fig c revised — UMA now shown as replacement in per-site framework (not MC).
EAM box grayed out instead of strikethrough to avoid visual overlap.

### Project roadmap & priority assessment

**Main line (80% effort): MC verification (Fig b)**
- Goal: verify Wagih's per-site + Fermi-Dirac predictions with direct atomistic MC
- Clear scientific question, advisor-endorsed, standalone value
- If predictions agree → validates the framework
- If they disagree → solute-solute interaction effects → new physics

**Extension (20% effort): UMA as drop-in replacement (Fig c)**
- Goal: replace EAM with UMA in Wagih's per-site framework
- NOT using MC with UMA (too slow, no extra scientific value)
- Proof-of-concept on one binary alloy (e.g. Al-Mg), ~10–100 GB sites

### Critical assessment of UMA extension

**Concerns:**
- Speed: UMA (neural network, ~ms/call) much slower than EAM (~μs/call).
  10⁵ relaxations could go from hours to weeks.
- "Just swapping calculator": if UMA matches EAM → "UMA works" but limited novelty.
  If mismatch → who is right? Need DFT ground truth for a few sites.
- GB environments are highly disordered → may be out-of-distribution for UMA.
- Wagih already covered 259 alloys. Re-doing with UMA alone is not enough.

**Where UMA becomes truly valuable:**
- Alloys where NO EAM potential exists (beyond NIST's 259 alloys).
- Multi-component / high-entropy alloys where EAM fitting is impractical.
- If MC verification finds EAM predictions are inaccurate, and UMA is closer to DFT
  → proves UMA is a better energy calculator for segregation.

### Multi-component alloy extension — the long-term payoff

If UMA proves accurate on binary alloys, the natural extension is:

| Alloy type | Combinations | EAM coverage | UMA coverage |
|------------|-------------|-------------|-------------|
| Binary (A-B) | ~4,000+ | ~259 (6%) | All |
| Ternary (A-B-C) | ~50,000+ | Very few (<0.1%) | All |
| High-entropy (5+) | ~millions | ≈ 0 | All |

This is the strongest argument for UMA: **one model covers the entire alloy space**
without system-specific potential fitting.

Multi-component GB segregation is also more physically complex:
- Multiple elements compete for GB sites (co-segregation / site competition)
- Wagih's independent-site Fermi-Dirac model may need modification
- This itself is an open scientific question

**Three-step story for the project:**
1. MC validates Wagih framework on binary alloy (scientific rigor)
2. UMA matches EAM on binary alloy per-site ΔE_seg (proof-of-concept)
3. UMA predicts segregation in multi-component alloys (new predictions, the contribution)

Step 3 is the real contribution, but Steps 1–2 are necessary to establish credibility.

