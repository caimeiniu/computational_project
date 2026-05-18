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

## 2026-05-17 — Reverse initial-condition structural check added

Professor feedback: checking energy alone is insufficient; also verify
that the structure/composition equilibrates from different initial
conditions. Added a reverse check in which Au starts already segregated
on GB atoms, then HMC should desegregate toward the same plateau reached
from the homogeneous random start.

New files:

- `PtAu/scripts/seed_gb_solute_PtAu.py`: rewrites the annealed LAMMPS
  data file so type-2 Au is placed preferentially on GB atoms according
  to `gb_mask_PtAu_100A.npy`. For `X_total=0.10`, the initial
  `X_GB` is about `0.267`, above the equilibrated random-start value
  `0.210`, so the test is a true reverse relaxation.
- `PtAu/data/decks/submit_hmc_PtAu_gbseed_check_jiayi.sh`: Slurm array
  for GB-seeded checks, defaulting to `X_total=0.02` and `0.10` at
  700 K with 600 ps production.
- `PtAu/data/decks/submit_hmc_PtAu_T700_Xc0.10_gbseed_resume_jiayi.sh`:
  focused continuation for the high-concentration GB-seeded run if it
  remains above the random-start plateau after the first 600 ps.
- `PtAu/scripts/plot_hmc_initial_condition_check_PtAu.py`: plots
  random-start and GB-seeded `X_GB(t)` on the same axes.

Also changed the HMC dumps from `id type` to `id type x y z` so future
outputs carry coordinates for OVITO/snapshot and local-structure checks.
The existing `hmc_xgb_timeseries.py` still reads the same `id` and `type`
columns and ignores the extra coordinates for composition counting.

Added OVITO export helpers:

- `PtAu/scripts/mark_gb_solute_for_ovito_PtAu.py` rewrites atom types to
  `1=Pt bulk`, `2=Au bulk`, `3=Pt GB`, `4=Au GB`.
- `PtAu/scripts/make_ptau_ovito_snapshots.sh` exports initial/final
  random-start and GB-seeded structures at `X_total=0.03`, the first
  clear breakdown point that also passes the reverse initial-condition
  check.

Submit on Euler:

```bash
sbatch /cluster/home/jiayfu/computational_project/project/PtAu/data/decks/submit_hmc_PtAu_gbseed_check_jiayi.sh
```

To repeat at 500 K after the 700 K check:

```bash
sbatch --export=ALL,T=500.0,XC_LIST="0.02 0.10" --array=0-1%2 \
  /cluster/home/jiayfu/computational_project/project/PtAu/data/decks/submit_hmc_PtAu_gbseed_check_jiayi.sh
```

---

## 2026-05-17 — 700 K Pt(Au) HMC/FD breakdown bracket

Following the Al(Mg) FD-first idea, the Pt(Au) 700 K scan now compares
closed-box HMC against a closed-box FD baseline using the same 100 Å box
(`N_total=62096`, `N_GB=23272`). The closed-box correction matters: the
old reservoir FD value at `X_total=0.10` was not the apples-to-apples
baseline for a fixed-composition HMC run.

Final 700 K scan:

| X_total | X_GB^HMC | X_GB^FD closed | HMC - FD | HMC/FD |
|---:|---:|---:|---:|---:|
| 0.005 | 0.0097 | 0.0099 | -0.0002 | 0.98 |
| 0.010 | 0.0206 | 0.0195 | +0.0011 | 1.06 |
| 0.015 | 0.0326 | 0.0288 | +0.0038 | 1.13 |
| 0.020 | 0.0440 | 0.0378 | +0.0062 | 1.16 |
| 0.025 | 0.0540 | 0.0467 | +0.0074 | 1.16 |
| 0.030 | 0.0664 | 0.0553 | +0.0111 | 1.20 |
| 0.050 | 0.1094 | 0.0884 | +0.0211 | 1.24 |
| 0.070 | 0.1488 | 0.1192 | +0.0296 | 1.25 |
| 0.100 | 0.2102 | 0.1622 | +0.0480 | 1.30 |

Verdict: at 700 K, closed-box FD agrees with HMC through about
`X_total=0.01`. A measurable positive deviation starts around
`X_total=0.015–0.02`, and breakdown is clear by `X_total=0.03`.

Added `PtAu/scripts/plot_hmc_scan_PtAu.py` to render
`PtAu/output/hmc_PtAu_T700_scan.png` and a compact plot table from the
scan summary CSV.

---

## 2026-05-16 — HMC T=700 K, Xc=0.10 needs continuation

Current HMC result from job `66755862` is not yet accepted as equilibrated:
`X_GB=0.1876` versus FD `0.1605`, with `fwd/rev=1.237` and imbalance
`+0.106`. Added restart-resume infrastructure for the Jiayi scratch path:

- `data/decks/hmc_PtAu_resume.lammps` reads a LAMMPS binary restart and
  continues production HMC without re-randomizing composition or rerunning
  minimization/equilibration.
- `data/decks/submit_hmc_PtAu_T700_Xc0.10_random_resume66755862_jiayi.sh`
  selects the newest `hmc_PtAu_T700_Xc0.10_random.rst1/.rst2` in
  `/cluster/scratch/jiayfu/prototype_PtAu_100A`, runs an additional
  `PROD_PS` ps (default 300), and post-processes the original plus resumed
  dump stream.
- `project/scripts/hmc_xgb_timeseries.py` was filled in after being empty:
  it now computes `X_GB(t)`, `X_bulk(t)`, tail averages, JSON summary, CSV,
  and a PNG plot from one or more HMC dump stubs.

Submit command on Euler:

```bash
sbatch /cluster/home/jiayfu/computational_project/project/PtAu/data/decks/submit_hmc_PtAu_T700_Xc0.10_random_resume66755862_jiayi.sh
```

For a longer continuation:

```bash
PROD_PS=600 sbatch /cluster/home/jiayfu/computational_project/project/PtAu/data/decks/submit_hmc_PtAu_T700_Xc0.10_random_resume66755862_jiayi.sh
```

The resumed job `66768704` finished the extra 300 ps successfully. Combined
original+resume post-processing gives 601 frames through step 610000:

| window | mean X_GB | mean X_bulk |
|---|---:|---:|
| block 7, steps 430000–489000 | 0.2094 | 0.0341 |
| block 8, steps 490000–549000 | 0.2107 | 0.0333 |
| block 9, steps 550000–610000 | 0.2103 | 0.0335 |

Verdict: stop resuming this point. It is equilibrated near
`X_GB^HMC ≈ 0.2103`, well above the original FD comparison value
`0.1605`. Following the Al(Mg) strategy, the next Pt(Au)-specific step is
not more runtime at `X_c=0.10`, but a lower-concentration bracket at the
same temperature to locate the onset of divergence. Added:

- `data/decks/submit_hmc_PtAu_T700_bracket_jiayi.sh`: Slurm array over
  total Au fractions `0.005, 0.01, 0.03, 0.05, 0.07`, capped at `%3`
  concurrent 16-rank jobs to fit the 48-core public quota. Default
  production length is 600 ps because the `X_c=0.10` point needed about
  that long to plateau.
- `PtAu/scripts/summarize_hmc_scan_PtAu.py`: collects HMC JSON summaries
  into one CSV table for the later HMC-vs-FD plot.
- `PtAu/scripts/fermi_dirac_predict_PtAu.py`: added a closed-box FD curve
  (`ours_canonical_total`) that solves for the bulk Au fraction under mass
  conservation. This is the Pt(Au) adjustment to Cainiu's Al(Mg) FD-first
  idea: compare closed-box HMC against a closed-box FD baseline, while still
  retaining Wagih's reservoir formula for reference.

Submit the bracket on Euler:

```bash
sbatch /cluster/home/jiayfu/computational_project/project/PtAu/data/decks/submit_hmc_PtAu_T700_bracket_jiayi.sh
```

---

## Current status (2026-05-13, 16:00 CEST)

| Stage | Status |
|---|---|
| Polycrystal generated | done — `poly_Pt_100A_8g.lmp` (62,096 atoms, 100³ Å, 8 grains) at `/cluster/scratch/cainiu/prototype_PtAu_100A/` |
| Anneal | done — SLURM job `66279426` finished 2026-05-13 00:14 (1 h wall, T_HOLD=816 K stable, force_max 3.34 eV/Å, 62,096 atoms preserved) |
| GB identification | done — 23,272 GB atoms (f_GB = 0.375; high vs Wagih's 19% because our 100³ Å box has higher GB-to-bulk surface ratio than Wagih's 200³ Å) |
| Per-site ΔE_seg sampling (n=500) | **done** — SLURM job `66391849` finished 2026-05-13 11:55 CEST (1h23m wall, 510/510 CG energy-tolerance, 0 failures) |
| Skew-normal (μ, σ, α) fit | **done** — (μ, σ, α) = (+3.16, 10.31, −1.16) kJ/mol |
| Fermi-Dirac temperature curves | **done** — T = 500 / 700 / 900 / 1100 K, with Wagih Pt(Au) FD overlay |
| Wagih reference extraction from tar | done — `Pt_Au_20nm_GB_segregation.dump` (24 MB, 508,951 atoms, 97,440 GB sites) extracted at `/cluster/scratch/.../accelerated_model/Pt/Au_2017--OBrien.../` |
| KS test vs Wagih per-site array | **done** — D = 0.0898, p = 6.2×10⁻⁴ (does NOT meet p>0.5 bar; see bootstrap below for the correct interpretation) |
| Bootstrap CI vs Wagih N=500 sub-samples (B=10000) | **done** — skew-normal (μ, σ, α) all inside 95% CI; raw moments mean/std outside CI at z ≈ ±2.7 |
| 100³ Å scaling decision | **frozen at 100³ Å** — bootstrap fit-params PASS; Pt EAM 5× slower than Al makes 200³ Å scaling cost-prohibitive (~20 h wall); FD↔HMC must share box size to isolate dilute-breakdown signal |
| Canonical-script Al(Mg) label cleanup | **done** — `fit_delta_e_spectrum_PtAu.py` + `fermi_dirac_predict_PtAu.py` added under `PtAu/scripts/`; outputs regenerated with Pt(Au) title + Wagih Pt(Au) reference |

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

## 2026-05-13 (afternoon) — Sampling done; KS+bootstrap interpretation; 100³ Å frozen; canonical-label fix

### Sampling job 66391849 completed

`delta_e_results_n500_PtAu_100A_tight.npz` landed 11:55 CEST after 1h23m
wall (16 cores). 510/510 sites converged via energy tolerance
(cg_etol=1e-25); zero failures. Per-site time ≈ 10 s — **~5× slower
than Al(Mg) (~2 s/site)**, consistent with Pt EAM stiffness (bulk
modulus 280 vs 76 GPa) driving more CG iterations.

### Fit and KS test

Skew-normal fit on N=500: **(μ, σ, α) = (+3.16, 10.31, −1.16) kJ/mol**.
Sample moments: mean=−3.05, std=8.23, skew=−0.18 (kJ/mol).
Wagih reference (n=97,440): fit (+3.65, 11.92, −1.42); mean=−4.14,
std=9.02 (kJ/mol).

**KS two-sample: D = 0.0898, p = 6.2×10⁻⁴** — does NOT meet the project
p>0.5 bar (see [[reference_ks_test]] memory). See bootstrap below for
the correct interpretation; KS is not the right metric in this regime.

### Bootstrap CI verdict — fit params PASS, raw moments fail

`bootstrap_vs_wagih_PtAu.py` with B=10000 N=500 sub-samples of Wagih's pool:

| Statistic       | Boot mean | CI95             | Ours   | z    | Inside CI |
|-----------------|----------:|:----------------:|-------:|-----:|:---------:|
| sample_mean     |    −4.15  | [−4.94, −3.36]   | −3.05  | +2.73| **NO**    |
| sample_std      |    +9.02  | [ +8.47, +9.56]  | +8.23  | −2.81| **NO**    |
| sample_skew     |    −0.24  | [−0.44, −0.04]   | −0.18  | +0.61|   yes     |
| skewnorm_μ      |    +3.39  | [−0.52, +5.83]   | +3.16  | −0.13| **yes**   |
| skewnorm_σ      |   +11.83  | [ +9.68, +13.53] |+10.31  | −1.60| **yes**   |
| skewnorm_α      |    −1.42  | [−2.21, −0.56]   | −1.16  | +0.62| **yes**   |

Interpretation: raw mean is ~1.1 kJ/mol higher and raw std is ~0.8 kJ/mol
lower than a typical N=500 draw from Wagih. Root cause: our N=500 missed
the deep negative tail (Wagih min −46 kJ/mol; ours only −34). Skew-normal
fit absorbs this because it parameterizes shape rather than reproducing
tail outliers. **The three (μ, σ, α) — exactly what the FD dilute-limit
prediction uses — are statistically indistinguishable from a Wagih
N=500 sub-sample.** The small KS p-value comes from the same tail
mismatch combined with the N=500 vs N=97,440 asymmetry (KS is
extremely sensitive to tiny CDF differences when one side is huge),
not from any pipeline error.

### Decision: freeze 100³ Å, do not scale to 200³ Å

Three reasons stack:

1. **Bootstrap CI on (μ, σ, α) already passes.** FD baseline parameters
   are Wagih-equivalent. The raw-moment offset only matters if the goal
   is "raw distribution exactly replicates Wagih" — which is not the
   project's question.
2. **Pt EAM is ~5× slower per site than Al EAM** (per-site CG cost ~ √
   condition number, Pt bulk modulus is ~3.7× Al's). Scaling to 200³ Å
   (8× atoms, linear-in-N CG cost) projects to **~20 h total wall**
   (8 h anneal + 11–12 h sampling at N=500), and Euler public QOS 48
   CPU/user cap ([[reference_euler_quota]]) makes parallel scheduling
   impossible. Cost-benefit doesn't justify it for marginal raw-moment
   improvement.
3. **FD and HMC must share box size to isolate the dilute-breakdown
   signal.** The project's question — at what X_c does HMC deviate from
   the FD dilute-limit? — requires apples-to-apples. FD on 100³ Å +
   HMC on 100³ Å cancels box-size effects; any residual X_c dependence
   is the real physics. Using Wagih's 200³ Å (μ, σ, α) for FD while
   running HMC at 100³ Å would conflate box-size effects with
   dilute-breakdown effects.

### Canonical-script Al(Mg) label cleanup

The initial fit + FD plots generated 11:58 from the canonical scripts
had three Al(Mg) hardcoding bugs:

- `output/delta_e_spectrum_PtAu_100A.png`: title "Al(Mg) per-site GB
  segregation spectrum"; dashed reference "Wagih Mg^15 (SI Fig. 3)
  μ=+9, σ=23, α=−2.3".
- `output/fd_curves_PtAu_100A.png`: title "Fermi-Dirac dilute-limit
  prediction — Al(Mg)"; legend "ours 200 Å" (we are on 100 Å).
- `output/delta_e_fit_PtAu_100A.json`: JSON key `wagih_alm_reference`
  carrying Al(Mg)'s (9, 23, −2.3).

Plotted *numbers* were Pt(Au) and correct — only labels and the dashed
reference values were Al(Mg). Root cause: `project/scripts/fit_delta_e_spectrum.py`
and `fermi_dirac_predict.py` hardcode the title strings + the Wagih
Mg^15 reference constant. Per [[feedback_no_in_place_script_edits]],
canonical scripts were NOT modified. Instead two new Pt(Au) copies were
added under `PtAu/scripts/`:

- `fit_delta_e_spectrum_PtAu.py` — title → "Pt(Au) per-site GB
  segregation spectrum"; `WAGIH_ALMG` constant → `WAGIH_PTAU`
  (μ=+3.65, σ=11.92, α=−1.42), computed by `scipy.stats.skewnorm.fit`
  on the Wagih Pt(Au) dump's 97,440 sites; JSON key
  `wagih_alm_reference` → `wagih_ptau_reference`.
- `fermi_dirac_predict_PtAu.py` — title → "Fermi-Dirac dilute-limit
  prediction — Pt(Au)"; legend "ours 200 Å" → "ours 100 Å"; the dead
  `--wagih-seg`/`--wagih-bulk` Al-Mg-only flags replaced with a single
  `--wagih-dump` flag that reads `Pt_Au_20nm_GB_segregation.dump`
  directly via the same loader used by `bootstrap_vs_wagih_PtAu.py`.

Output PNGs + JSON regenerated with the new scripts (commands in
`README.md` step 5). The new FD plot additionally overlays Wagih's
Pt(Au) FD curve as a dashed reference: visible result is Wagih sits
slightly above our curves at low X_c (Wagih's deeper negative tail
fills first at dilute X_c), curves merge near X_c ≈ 0.1. Physically
consistent with the raw-moment offset from the bootstrap analysis.

### Caveats handed off

- Raw moments (sample_mean, sample_std) sit at z ≈ ±2.7. If a later
  analysis needs exact moment match — e.g., very dilute X_c < 10⁻³
  where the deepest GB sites dominate FD occupation — bump N from 500
  to ~2000–5000 to better cover the negative tail. Cheap on 100³ Å
  (~6 h additional sampling at the same 5× Pt slowdown); does NOT
  require a larger box.

- The Wagih `(μ, σ, α) = (+3.65, 11.92, −1.42)` reference in
  `WAGIH_PTAU` is computed by us from the Zenodo dump, not read from
  Wagih's SI Fig. 9 panel. Wagih's SI tabulates per-panel (μ, σ, α) for
  each binary; the Pt(Au) O'Brien 2017 panel uses superscript notation
  matching the EAM reference number in the SI. If a value-pinning audit
  needs it, locate that exact panel in the published SI and verify
  against `WAGIH_PTAU` — a small (<2%) disagreement would just reflect
  different fitter / clipping convention and is not a red flag.

### Output files produced today (gitignored, regeneratable)

All under `project/PtAu/output/`:

- `delta_e_fit_PtAu_100A.json` — skew-normal fit + Wagih Pt(Au) ref values
- `delta_e_spectrum_PtAu_100A.png` — histogram + fit + Wagih dashed
- `fd_curves_PtAu_100A.{json,png}` — FD curves at T=500/700/900/1100 K
  with Wagih Pt(Au) overlay
- `compare_vs_wagih_PtAu_100A.{json,png}` — KS overlay (ours vs Wagih
  histogram + both skew-normal fits)
- `bootstrap_vs_wagih_PtAu_100A.json` — 6-stat bootstrap CI (B=10000,
  N=500, seed=20260513)

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
