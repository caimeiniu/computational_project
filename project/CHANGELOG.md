# Changelog

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
