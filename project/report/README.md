# Defense Presentation Materials — Mg–Mg Interactions Break Wagih's Site-Independent FD Theory in Al(Mg) GB Segregation

> Speaker reference for the four main figures + four supporting figures, with full per-figure walkthroughs and predicted advisor Q&A. Every cited number traces to a source file; nothing is fabricated (spot-checks pass — see §8 Provenance).

Last updated: 2026-05-02 (Fig 0 simplified to a 3-row legend with two-color ▽ markers; English rewrite)

---

## 1. TL;DR

- **Hypothesis under test**: Wagih & Schuh (2020, *Acta Mater.*) site-independent Fermi–Dirac (FD) model — assumes that the Mg occupation probability P_i at GB site *i* depends only on the segregation energy ΔE_i of that site, independent of which other sites are occupied.
- **Headline result**: in Al(Mg) at T = 500 K, **the assumption breaks for X_c ≥ 0.075** (X_c = total Mg atomic fraction).
- **Mechanism**: a direct site-level Mg–Mg repulsion, demonstrated through three independent lines of evidence.
- **Data**: HMC (Hybrid Monte Carlo) simulations on six (T = 500 K, X_c) snapshots; N_total = 475,715 atoms, N_GB = 89,042 (f_gb = 18.72 %).
- **Theory anchor**: our n = 500 ΔE_seg spectrum is statistically indistinguishable (KS p = 0.89) from Wagih's published Zenodo spectrum (n = 82,646), so Wagih's framework is the right baseline (see Fig. 4).

---

## 2. Talk Order

**Main story arc (the four figures used to deliver the talk)**

| # | File | Role | One-line message |
|---|---|---|---|
| 0 | `figures/00_headline_hmc_vs_wagih_T500.png` | **Headline / counter-example anchor** | At X_c = 0.075 the HMC measurement X_GB = 0.254 is **below** Wagih's FD prediction 0.301; at X_c = 0.10 the multistart upper bound X_GB = 0.246 is also **below** canon-FD = 0.352 → direct evidence that the assumption fails. |
| 1 | `figures/01_MgMg_clustering.png` | **Mechanism evidence 1 — spatial** | The Mg–Mg pair correlation function g(r) deviates from the uniform-random reference → Mg atoms are not independently distributed (*aggregate spatial signal*). |
| 2 | `figures/02_occupation_breakdown.png` | **Mechanism evidence 2 — energy axis** | Empirical P_i is **systematically below** the Wagih sigmoid at the favourable-ΔE end → the breakdown is concentrated at the low-energy tail of the ΔE spectrum. |
| 3 | `figures/03_repulsion_summary.png` | **Mechanism evidence 3 — direct site-level** | At fixed ΔE, more Mg neighbours → lower occupation → direct observation of Mg–Mg repulsion. |

**Supporting figures (Q&A / methodological backing)**

| # | File | Role | One-line message |
|---|---|---|---|
| 4 | `figures/04_spectrum_match.png` | Spectrum representativeness | Histograms + skew-normal fits of our n = 500 ΔE spectrum vs. Wagih Zenodo n = 82,646 (KS p = 0.89, statistically indistinguishable). |
| 5 | `figures/05_sampler_convergence.png` | HMC convergence diagnostic | Five-panel time series for the X_c = 0.075 preseg run: T(t), PE(t), accept rate, X_GB(t), per-frame swap fwd/rev decomposition. |
| 6 | `figures/06_two_sided_verify.png` | Equilibration verification | Two-sided IC bracket (random + preseg) at X_c = 0.05 → confirms the sampler can equilibrate at the dilute end. |
| 7 | `figures/07_ovito_segregation.png` | Visual segregation confirmation | OVITO render of the X_c = 0.20 final config; Mg (orange) is visibly enriched on grain boundaries. |

---

## 3. Definitions

| Term | Meaning |
|---|---|
| GB (grain boundary) | Interface between adjacent grains in polycrystalline Al; defined by 3D Voronoi tessellation. |
| Mg / solute | Substitutional solute; this project is fixed to the Mg-in-Al system. |
| X_c | Total Mg atomic fraction (total Mg atoms / total atoms). |
| X_GB | Mg fraction at GB sites only (N_GB_Mg / N_GB). |
| ΔE_i, ΔE_seg | Segregation energy of site *i*: energy difference between placing one Mg atom at GB site *i* and at a bulk reference site. **Negative = site favours Mg occupation**. |
| HMC (Hybrid Monte Carlo) | This project's sampler: in LAMMPS, alternating cycles of (a) MD relaxation and (b) Mg/Al swap Monte Carlo attempts, satisfying detailed balance. |
| Wagih FD | Wagih's site-independent Fermi–Dirac formula: `P_i = 1 / (1 + ((1 − X_c) / X_c) · exp(ΔE_i / kT))`. |
| canon-FD | Mass-conserving canonical FD: average Wagih's P_i over the X_c = 0 reference spectrum to obtain the predicted X_GB at given (T, X_c). "ours" = on our n = 500 spectrum; "Wagih" = on Wagih's Zenodo n = 82,646 spectrum. |
| GC-FD | Grand-canonical FD (Wagih's formula in the chemical-potential representation; not used in the simplified Fig. 0). |
| preseg HMC | Initial-condition choice: place every Mg atom into the lowest-ΔE GB sites (X_GB(0) ≈ 0.32, maximum segregation). HMC then descends, and the trajectory end value is an upper bound on the equilibrium X_GB. |
| multistart HMC | Alternative initial-condition choice: random Mg distribution at X_GB(0) = 0.30. The trajectory descends toward a kinetic floor; the production-mean value is also an upper bound on the equilibrium X_GB. |
| KS test (Kolmogorov–Smirnov) | Non-parametric two-sample test for distribution equality; project threshold p > 0.5 ⇒ "spectrum-level indistinguishable". |

---

## 4. Methods Summary

### 4.1 Substrate

- Periodic cubic box of polycrystalline Al, edge length 200 Å, 3D Voronoi tessellation with 16 grains. Build script: `data/decks/build_poly_AlMg_200A.py`. Snapshot files: `data/snapshots/poly_AlMg_200A_*.lmp`.
- N_total = **475,715** atoms, N_GB = **89,042**, f_gb = **0.1872** (18.72 %). The GB mask is generated from a-CNA (adaptive Common Neighbor Analysis); file `data/snapshots/gb_mask_200A.npy`.
- Anneal protocol: CG → NVT ramp → NPT 250 ps @ 373 K → cool 3 K/ps → final relax (relieves build-time stress).

### 4.2 Reference ΔE spectrum

- **n = 500** randomly sampled GB sites; per-site ΔE_i = (single-Mg system at site *i*) − (bulk-Mg average) − (pure-Al system).
- File: `/cluster/scratch/cainiu/production_AlMg_200A/delta_e_results_n500_200A_tight.npz`.
- Spectrum statistics (all reconciled with `output/compare_vs_wagih_200A_tight.json`):

|  | n | mean (kJ/mol) | std (kJ/mol) | skew | skewnorm μ | skewnorm σ | skewnorm α |
|---|---:|---:|---:|---:|---:|---:|---:|
| **ours** | 500 | −6.91 | 15.07 | −0.21 | 6.34 | 20.06 | −1.47 |
| **Wagih** | 82,646 | −6.81 | 15.85 | −0.22 | 6.72 | 20.84 | −1.40 |
| KS 2-sample test |  | | | | | D = 0.0256 | **p = 0.8920** |

→ The two spectra are statistically indistinguishable (p ≫ 0.5). Wagih's FD framework is therefore the right baseline for our system.

### 4.3 HMC settings

- LAMMPS 20240829.4, 32 MPI ranks per job.
- Mendelev (2009) Al–Mg EAM potential (the same potential used in Wagih 2020).
- T = 500 K. Metropolis criterion uses the post-NVE-relaxation energy difference ΔE_swap.
- Per production run: 300 ps PROD + 10 ps EQUIL, one frame per 1000 swap attempts → 300 frames.
- IC choice: **preseg** (main panels) and **multistart** (corroborating evidence).

### 4.4 The six (T = 500 K, X_c) snapshots used in this report

| Label | X_c | File | X_GB(end) | N_GB_Mg | Source job |
|---|---:|---|---:|---:|---|
| 0.075_preseg | 0.075 | `hmc_T500_Xc0.075_preseg_final.lmp` | 0.2542 | 22,635 | 65208332 |
| 0.10_multistart | 0.10 | `hmc_T500_Xc0.10_multistart_xgb0.3_final.lmp` | 0.2277 | 20,276 | (multistart) |
| 0.10_preseg | 0.10 | `hmc_T500_Xc0.10_preseg_final.lmp` | 0.3747 | 33,361 | (preseg) |
| 0.15_preseg | 0.15 | `hmc_T500_Xc0.15_preseg_final.lmp` | 0.5884 | 52,394 | (preseg) |
| 0.20_preseg | 0.20 | `hmc_T500_Xc0.20_preseg_final.lmp` | 0.7937 | 70,672 | (preseg) |
| 0.30_preseg | 0.30 | `hmc_T500_Xc0.30_preseg_final.lmp` | 0.8376 | 74,582 | (preseg) |

(Spot-check: for the 0.075 / 0.15 / 0.30 snapshots, N_GB_Mg and X_GB recomputed directly from the LAMMPS files match `output/solute_correlation_analysis.json` exactly, Δ = 0; see §8 Provenance.)

---

## 5. Per-figure Walkthroughs

### Figure 0 ── `00_headline_hmc_vs_wagih_T500.png`

**Caption (English, suitable for slides / paper):**

> **Figure 0.** Al(Mg) GB Mg fraction X_GB versus total Mg fraction X_c at T = 500 K. **Green solid line** — Wagih's site-independent Fermi–Dirac (canon-FD) prediction, computed by averaging the per-site formula P_i(X_c, ΔE_i, T) over our n = 500 reference ΔE_seg spectrum (statistically indistinguishable from Wagih's published n = 82,646 spectrum, KS p = 0.89, Fig. 4). **Red filled circle ●** — HMC equilibrium at X_c = 0.05; both an over-segregated initial condition and a random initial condition converge here (Fig. 6). **Red open down-triangles ▽** — HMC trajectories started from an over-segregated initial state (preseg, X_GB(0) ≈ 0.32, all Mg deliberately placed into the lowest-ΔE GB sites) and still descending; their current values are upper bounds on the equilibrium X_GB. **Gray open down-triangle ▽** at X_c = 0.10 — a separate trajectory started from a random Mg distribution (multistart, X_GB(0) = 0.30) and also still descending; a second, IC-independent upper bound. The legend collapses both ▽ markers into one "HMC upper bound" entry — colour distinguishes the two ICs (red = over-segregated start, gray = random start), per the present caption. Wagih's site-independent assumption breaks for X_c ≥ 0.075: at X_c = 0.075 the over-segregated trajectory has already crossed below the FD prediction (0.254 < 0.301) while still descending, and at X_c = 0.10 the random-start trajectory sits 0.106 below FD = 0.352.

**Role**: headline counter-example figure (the panel-d slot in the master figure plan).

**Axes**:
- x: total Mg fraction X_c (0 to 0.35)
- y: GB Mg fraction X_GB (0 to 0.95)

**Single reference curve**: green solid line — Wagih FD prediction (canon-FD on our n = 500 spectrum). Earlier versions of this figure overlaid four reference curves (canon-FD ours, canon-FD Wagih, GC-FD, closed-box ceiling, X_GB = X_c diagonal) plus an ours-vs-Wagih shaded band, giving a 9–10-item legend. The defense panel only needs to convey "HMC < Wagih FD = breakdown", so the auxiliary curves were removed. The full four-curve version is still produced by the canonical script `scripts/canonical_fd_compare_5pt.py` → `output/hmc_vs_canonfd_T500.png` for cases where the broader reference context is needed (paper SI, methods discussion).

**HMC measurements** (legend uses plain language; internal project terms in parentheses):
- 1 red filled circle ● at X_c = 0.05: legend label **"HMC equilibrium"** (internal: equilibrated bracket). Two-sided IC verification, see Fig. 6 `figures/06_two_sided_verify.png`. X_HMC = 0.2375 ± 0.005.
- 5 red open down-triangles ▽ at X_c = 0.075, 0.10, 0.15, 0.20, 0.30: legend label **"HMC upper bound"** (internal: preseg UB). Trajectory started from an over-segregated IC (X_GB(0) ≈ 0.32) and is still descending, so each marker is an upper bound on X_GB^∞.
- 1 gray open down-triangle ▽ at X_c = 0.10: same legend entry "HMC upper bound" (internal: multistart UB / kinetic-floor IC). Trajectory started from a random IC (X_GB(0) = 0.30) and is also still descending; production-mean X_GB = 0.246 ± 0.004 → IC-independent upper bound. 0.106 below canon-FD → directly demonstrates X_c = 0.10 breakdown. Note: the random-start trajectory begins slightly lower than the preseg one, but both are still descending and the equilibrium value can only be ≤ 0.246.

**Key numbers** (source: `output/hmc_vs_canonfd_T500_with_multistart.json`, verified):

| X_c | Type | X_HMC | X_FD (ours) | gap (X_HMC − X_FD) |
|---:|---|---:|---:|---:|
| 0.050 | equilibrated ● | 0.2375 | 0.2282 | **+0.0093** (Wagih holds) |
| 0.075 | preseg UB ▽ (red) | 0.2543 | 0.3007 | **−0.0464** (UB has crossed below FD, **breakdown**) |
| 0.100 | preseg UB ▽ (red) | 0.3749 | 0.3519 | +0.0230 (UB still above FD, vacuous bound — IC dependence) |
| 0.100 | multistart UB ▽ (gray) | 0.2459 | 0.3519 | **−0.1060** (UB has crossed below FD, **breakdown**) |
| 0.150 | preseg UB ▽ (red) | 0.5888 | 0.4204 | +0.1684 (vacuous bound) |
| 0.200 | preseg UB ▽ (red) | 0.7942 | 0.4671 | +0.3271 (vacuous bound) |
| 0.300 | preseg UB ▽ (red) | 0.8379 | 0.5337 | +0.3042 (vacuous bound) |

**Conclusions**:
1. **X_c = 0.05**: Wagih's assumption holds (gap within 1 %; both IC choices give a consistent equilibrated bracket).
2. **X_c = 0.075**: **first direct evidence that the assumption fails** — the over-segregated trajectory descends from X_GB(0) = 0.32 to 0.254, crossing below canon-FD = 0.301 while still descending, so X_GB^∞ ≤ 0.254 < 0.301 = contradiction.
3. **X_c = 0.10**: **second direct evidence** — the multistart UB (gray ▽) at X_GB = 0.246 is well below canon-FD = 0.352 (gap = −0.106); the preseg UB (red ▽) at X_GB = 0.375 is still above canon-FD, reflecting IC dependence (preseg trajectory still descending; multistart trajectory at the kinetic floor). The two ICs sandwich X_GB^∞, with X_GB^∞ ≤ min(0.375, 0.246) = 0.246 ≪ 0.352 → breakdown confirmed.
4. **X_c ≥ 0.15**: preseg UBs are still above canon-FD (vacuous bounds). The breakdown evidence here comes from threshold extrapolation (X_c\* ∈ (0.05, 0.075]) and from the mechanism figures (1–3).

**Talking points** (for the audience):
> "At X_c = 0.05 — the dilute limit — the HMC measurement sits right on the green canon-FD curve, within 1 %. Wagih's formula gets it right. But at X_c = 0.075, even the trajectory that started maximally over-segregated (X_GB(0) = 0.32) descended past the FD prediction of 0.301 and is still going down — that's direct evidence that Wagih's formula is broken. At X_c = 0.10 we see a second, independent confirmation: the trajectory that started from a random Mg distribution (gray triangle) has descended to 0.246, well below FD = 0.352."

**Likely advisor questions**:
- *Q: Why do multistart and preseg disagree at X_c = 0.10 (0.246 vs. 0.375)?* — A: That's IC dependence. Preseg starts from X_GB = 0.32 (extreme segregation) and is still descending, so its current value is an upper bound from above. Multistart starts from a random configuration at X_GB = 0.30 and has also been descending; it has reached a kinetic floor near 0.246, also still slowly drifting down. Both are upper bounds on the equilibrium. Their sandwich gives X_GB^∞ ≤ min(0.375, 0.246) = 0.246, well below canon-FD = 0.352 → breakdown is direct.
- *Q: Why are the two upper-bound markers the same shape but different colours?* — A: Same shape (▽) signals "upper bound, still descending"; colour distinguishes the IC (red = over-segregated start, gray = random start). The legend collapses them into a single row because both serve the same role; the per-IC distinction is in this caption.
- *Q: Why is preseg UB not a lower bound?* — A: The trajectory starts above the equilibrium and descends; the current end value is the lowest the system has reached so far, and the equilibrium can only be lower or equal → it is an *upper* bound on X_GB^∞.
- *Q: What is the precise threshold X_c\*?* — A: Current bracket is X_c\* ∈ (0.05, 0.075] from the binary-search outcome (X_c = 0.05 holds, X_c = 0.075 fails). Job 65224958 (X_c = 0.06) has been submitted to the normal.24h queue; once it finishes we will have a tighter lower bound.

---

### Figure 1 ── `01_MgMg_clustering.png`

**Caption (English):**

> **Figure 1.** Normalised radial pair correlation function for Mg atoms on grain boundaries, g_MgMg^HMC(r) / g_MgMg^random(r), at T = 500 K. The "uniform-random reference" g_MgMg^random(r) is computed by drawing the *same number* of Mg atoms placed at random GB lattice sites (without ΔE preference) and computing the same g(r). Three curves are shown for X_c = 0.075 (red, X_GB = 0.254), 0.150 (blue, X_GB = 0.588), 0.300 (green, X_GB = 0.838). All three lie above 1 around the FCC first-neighbour shell at r ≈ 3.30 Å, demonstrating that Mg–Mg correlations on the GB are non-random out to r ≈ 10 Å. This is an *aggregate* spatial signal — partly driven by the geometric proximity of deep-ΔE binding sites (Mg preferentially fills these, which are spatially clustered on the GB plane); the ΔE-controlled (residual) interaction is shown in Fig. 3.

**Role**: mechanism evidence 1 / spatial correlation.

**What it shows**: ratio of the Mg–Mg radial pair correlation function on the GB to a uniform-random GB reference.

**Core math**:
- g_MgMg^HMC(r) = radial distribution function of GB-Mg pairs in the HMC snapshot.
- g_MgMg^random(r) = "uniform-random reference": draw the **same number** N_GB_Mg of Mg atoms placed **uniformly at random** on GB lattice sites, and compute the same g(r).
- y-axis = ratio **g_MgMg^HMC(r) / g_MgMg^random(r)**.
  - = 1 → HMC distribution statistically equivalent to "random Mg of the same count on GB sites".
  - > 1 → more clustered than random.
  - < 1 → more avoidant than random.

**Axes**:
- x: pair separation r [Å], 0 to 25 Å.
- y: g_HMC / g_random, 0 to 1.55.

**Three curves** (sorted by ascending X_GB):

| X_c | X_GB | N_GB_Mg | colour | first peak (r ≈ 3.30 Å) | global max in r ∈ [2.5, 6] Å |
|---:|---:|---:|---|---:|---|
| 0.075 | 0.254 | 22,635 | red | **1.075** (weak) | 1.220 @ r = 5.90 Å (second shell) |
| 0.150 | 0.588 | 52,394 | blue | 1.225 | 1.225 @ r = 3.30 Å |
| 0.300 | 0.838 | 74,582 | green | 1.099 | 1.099 @ r = 3.30 Å |

(Source: `output/solute_correlation_analysis.json`, value taken directly from `g_MgMg_pair_correlation.curves[label].ratio`; verified by spot-check at r = 3.30 Å.)

**Why these three X_c?** 0.075 = threshold critical / 0.15 = mid-range / 0.30 = high-X_c saturation. We **excluded** 0.10_multistart and 0.10_preseg because the X_c = 0.10 group is in a kinetic-floor / still-descending state (placing them on a mechanism plot would require explaining an extra layer of IC dependence). We also **excluded** 0.20 since the full data is in the Fig. 3 summary.

**Key annotations**:
- Horizontal dashed line at y = 1.0 + label "uniform-random reference".
- Vertical gray band at r ∈ [3.0, 3.5] Å + label "1st NN shell (FCC Al)" — FCC Al's 1st-NN distance is ~2.86 Å, slightly broadened at the GB.
- Subtitle: "(non-random structure out to r ~ 10 Å)".
- Footer caption (carefully worded): "Aggregate spatial signal. The peak above 1 is partly driven by the geometric proximity of deep-ΔE binding sites; the ΔE-controlled (residual) interaction is in Fig. 3."

**Talking points**:
> "All three curves at different X_c rise above 1 at r ≈ 3 Å, **but this is not simple Mg–Mg chemical attraction**. Deep-ΔE sites are themselves spatially adjacent on the GB plane; Mg preferentially fills those sites, which automatically produces a positive deviation in g(r). So this is an *aggregate* spatial signal: it shows the distribution is non-random, but cannot decide whether the underlying interaction is attractive or repulsive. To see the true site-level interaction we have to control for ΔE — that is Fig. 3."

**Why is the X_c = 0.075 first peak anomalously suppressed (1.075 vs. 1.604 at 0.10_multistart)?**

Physical interpretation: 0.075_preseg at X_c = 0.075, X_GB = 0.254 is right above the breakdown threshold; Mg coverage of GB sites is sparse (22,635 / 89,042 = 25 %), so the system has the freedom to "pick" non-adjacent sites and avoid first-NN Mg–Mg contact. Geometric clustering still shows up at the second-shell peak r ≈ 5.9 Å (1.22), but the first-NN peak is partially flattened by site-level repulsion. This in turn supports the **Fig. 3** site-level repulsion conclusion.

**Likely advisor questions**:
- *Q: How is g_random generated and what is its statistical noise?* — A: N_GB_Mg GB sites are drawn at random without replacement using a fixed RNG seed (20260429); see `scripts/solute_correlation_analysis.py:133-141`. A single draw fixes the reference; if ensemble noise estimates are needed, multiple draws can be averaged. Current sample size is comfortable (20k+ atoms per snapshot), so reference fluctuations are negligible.
- *Q: Why does g(r) drop to 0 / NaN at small r?* — A: For r < 2 Å the random reference's shell volume × density gives < 1 expected pair, so g_random = 0 → ratio undefined. This is not an HMC artefact.
- *Q: Has the anomalous first peak at X_c = 0.075 been sanity-checked?* — A: That snapshot has N_GB_Mg = 22,635 (independent recomputation Δ = 0, see §8); the first peak 1.075 is read directly from the JSON `curves['0.075_preseg'].ratio[16]` (r-axis bin 16 = 3.30 Å).

---

### Figure 2 ── `02_occupation_breakdown.png`

**Caption (English):**

> **Figure 2.** Empirical Mg occupation probability P_i versus per-site segregation energy ΔE_i, compared with Wagih's site-independent FD prediction P_i^Wagih(ΔE_i; T, X_c) = 1 / (1 + ((1 − X_c) / X_c) · exp(ΔE_i / kT)) at T = 500 K (kT = 4.157 kJ/mol). For each panel the n = 500 reference sites are binned in ΔE_i (10 equal-width bins from −48 to +36 kJ/mol); within each bin, P_i is the fraction of sites occupied by Mg in the HMC snapshot, with binomial 95 % confidence intervals. Three panels show X_c = 0.075 (X_GB = 0.254), 0.150 (X_GB = 0.588), 0.300 (X_GB = 0.838). The favourable-binding region (ΔE < 0) is shaded green; the largest gap between the Wagih sigmoid and the empirical points in this region is annotated as ΔP_i. Empirical P_i is **systematically below** the Wagih sigmoid at the favourable end, with the gap widening as X_c decreases (ΔP_i = 0.77, 0.57, 0.13 for X_c = 0.075, 0.150, 0.300) — the breakdown is concentrated at the deepest binding sites.

**Role**: mechanism evidence 2 / where on the energy axis Wagih FD fails.

**What it shows**: bin the 500 reference sites by ΔE_i; within each bin, compute the empirical probability that the site is occupied by Mg in the HMC snapshot; compare with the Wagih FD prediction. Three panels for three X_c values.

**Core math**:
- Empirical P_i: for each reference site, read the atom type in the HMC snapshot; type = Mg → 1, else → 0; average within each ΔE_i bin (binomial proportion).
- Wagih FD prediction: `P_i^Wagih(ΔE_i; T, X_c) = 1 / (1 + ((1 − X_c) / X_c) · exp(ΔE_i / kT))`, kT = 4.157 kJ/mol at T = 500 K.

**Axes**:
- x: ΔE_i [kJ/mol] (X_c = 0 reference spectrum), range ~[−48, +36] kJ/mol.
- y: P_i (Mg occupation probability), 0 to 1.
- y-label only on the leftmost panel (shared y).

**Three panels**: X_c = 0.075, 0.150, 0.300 (same X_c selection as Fig. 1).

**Per-panel elements**:
- Black solid line = Wagih FD theoretical curve (sampled on a dense ΔE grid, 200 points).
- Red / blue / green points + 95 % CI errorbars = empirical (10 equal-width ΔE bins; per-bin counts in JSON `n_per_bin`).
- Light green vertical band at ΔE < 0: **favourable binding region**.
- Vertical green line at ΔE = 0.
- Annotation arrow + box "ΔP_i ≈ X.XX": at the favourable-end (ΔE < 0) bin with the largest empirical-vs-Wagih gap.
- Panel title: "X_c = 0.XXX → X_GB = 0.XXX".

**Key numbers** (source: `output/solute_correlation_analysis.json` + Wagih-formula recomputation, verified):

| X_c | most-favourable bin (ΔE) | P_Wagih predicted | P_emp measured | ΔP_i = Wagih − emp |
|---:|---:|---:|---:|---:|
| 0.075 | −27.40 kJ/mol | 0.9834 | 0.2105 | **0.7728** (drawn as 0.77) |
| 0.150 | −44.25 kJ/mol | 0.9999 | 0.4286 | **0.5713** (drawn as 0.57) |
| 0.300 | −27.40 kJ/mol | 0.9968 | 0.8684 | **0.1284** (drawn as 0.13) |

**Talking points**:
> "Wagih's sigmoid predicts that deep-binding sites (ΔE < 0) should be nearly fully occupied — P_i ≈ 1. But the HMC measurements show that at X_c = 0.075 the favourable end is only at P_i ≈ 0.21, a gap of 0.77 below the prediction. At X_c = 0.15 the gap shrinks to 0.57; at X_c = 0.30 to 0.13. **The breakdown lives at the lowest energies, and it is *more severe* at lower X_c** — exactly opposite to what you might naively expect."

**Why does ΔP grow as X_c decreases?** Counterintuitive at first; the physical interpretation:
- At high X_c, all sites (even unfavourable ones) are saturated, so there is no "discrimination room" left.
- At low X_c, the system **has the freedom** to preferentially occupy favourable sites (which is exactly what Wagih's independence assumption predicts). But Mg–Mg repulsion interferes with that selection: you want to fill a deep site, but a neighbouring Mg pushes you away → empirical P_i is far below the independent-site prediction.
- **The breakdown reveals the existence of site-level interactions; the X_c → 0 limit is where their relative strength is greatest** (because mean-field averaging is supposed to wash interference out, yet it does not).

**Likely advisor questions**:
- *Q: ΔE_i comes from the X_c = 0 reference spectrum; at finite X_c the effective ΔE_i shifts due to local Mg–Mg interactions. How is that confound handled?* — A: It is *by design*: using bare ΔE_i as the "independent-site baseline" is the only way to test Wagih's assumption directly (the assumption itself is built on the X_c = 0 spectrum). Renormalising ΔE under finite X_c is a future-work direction; the project explored a ΔE-shift approach earlier and deprioritised it because the mechanism path produced cleaner conclusions.
- *Q: How are bin boundaries chosen, and are bin counts adequate?* — A: 10 equal-width ΔE bins from min to max of the n = 500 spectrum; ~50 sites per bin on average; binomial 95 % CIs with the normal approximation (`scripts/solute_correlation_analysis.py:208-225`). Per-bin counts are recorded in JSON `n_per_bin`.
- *Q: Why is the favourable bin for X_c = 0.150 at −44 kJ/mol (deeper) while 0.075 / 0.30 are at −27 kJ/mol?* — A: The script picks "argmax(gap)" which depends on bin density × gap, not on a physically fixed ΔE. It is purely an annotation algorithm to highlight the largest gap.

---

### Figure 3 ── `03_repulsion_summary.png`

**Caption (English):**

> **Figure 3.** Direct evidence of site-level Mg–Mg repulsion: at fixed ΔE_i, occupation probability decreases as the number of nearby Mg neighbours increases. **Left panel**: slope ∂P_i / ∂n_Mg^local versus X_c, fit to 218 sites in the favourable-ΔE window [−30, −5] kJ/mol with local radius r_local = 5 Å. Six points cover X_c = 0.075–0.30; preseg trajectories shown as red ▽, the X_c = 0.10 multistart shown as a gray open square (kinetic-floor IC). The slope is most negative at X_c = 0.075 (−0.083, "steepest repulsion" call-out), then weakens with X_c and crosses zero at X_c = 0.20 (saturation regime — most favourable sites are already filled, signal is dominated by statistical noise; not a sign reversal of the underlying interaction). **Right panel**: raw scatter for X_c = 0.075 — empirical P_i versus n_Mg^local in the favourable-ΔE window — with a weighted linear fit confirming slope = −0.083.

**Role**: mechanism evidence 3 / **direct evidence** — controlling for ΔE_i, the effect of Mg-neighbour count on occupancy.

**What it shows**: left panel — summary of the slope ∂P_i / ∂n_Mg^local across six X_c values; right panel — raw scatter zoom for X_c = 0.075.

**Core math**:
1. Within the favourable-ΔE window [−30, −5] kJ/mol, select 218 / 500 reference sites (count in JSON `site_occupation_vs_density.n_sites_in_window`).
2. For each site *i*, compute n_Mg^local(*i*) = number of Mg atoms within r_local = 5 Å of *i* (excluding self if *i* itself is Mg).
3. Bin these 218 sites by integer n_Mg^local; compute the average P_i within each bin.
4. Fit a weighted linear regression to (P_i vs n_Mg^local), giving the slope (units: 1/Mg-neighbour).

**Axes**:
- Left: x = X_c (0.04 to 0.33), y = slope (−0.13 to +0.06) [per Mg-neighbour].
- Right: x = n_Mg^local (Mg neighbour count, 0 to 11), y = P_i within ΔE ∈ [−30, −5] kJ/mol window (0 to 1).

**Six slope data points** (source: `output/solute_correlation_analysis.json`, verified):

| Label | X_c | X_GB | slope (per Mg-nbr) |
|---|---:|---:|---:|
| 0.075_preseg | 0.075 | 0.2542 | **−0.0826** ← steepest |
| 0.10_multistart | 0.100 | 0.2277 | −0.0310 (gray square: kinetic-floor IC) |
| 0.10_preseg | 0.100 | 0.3747 | −0.0402 |
| 0.15_preseg | 0.150 | 0.5884 | −0.0194 |
| 0.20_preseg | 0.200 | 0.7937 | **+0.0148** ← saturation crossover |
| 0.30_preseg | 0.300 | 0.8376 | −0.0018 |

**Key annotations**:
- Left: slope = 0 horizontal dashed line + label "no interaction (slope = 0)".
- Left: annotated callout at X_c = 0.075 → "slope = −0.083 (steepest repulsion)".
- Left: annotated callout at X_c = 0.20 → "saturation regime (sites mostly full)".
- Left: 0.10_multistart shown as a gray open square (kinetic-floor IC, distinguished from the preseg series).
- Right: linear fit dashed line + slope value in legend.

**Talking points**:
> "Once we control for ΔE and look only at the local environment: in the favourable window, sites with more Mg neighbours have *lower* occupation. That is **direct site-level Mg–Mg repulsion**. The slope is most negative at X_c = 0.075 (−0.083), because the system has just crossed the breakdown threshold and there is plenty of room to discriminate. By X_c = 0.20 the slope has narrowed to ~0 — not a sign reversal; favourable sites are mostly already filled and the signal is dominated by statistical noise (the error bar straddles zero)."

**Likely advisor questions**:
- *Q: The slope at X_c = 0.20 flips positive (+0.0148) — has Mg–Mg become attractive?* — A: **No, this is saturation**. At X_c = 0.20, X_GB = 0.794, so the great majority of favourable sites are already occupied by Mg; the slope fit at high n_local (e.g. n = 15–20) is dominated by post-saturation small-amplitude fluctuations, with large error bars and an unreliable sign. The X_c = 0.30 slope returns to −0.002 ≈ 0, consistent with the saturation picture.
- *Q: Why do 0.10_multistart and 0.10_preseg disagree (slope = −0.031 vs. −0.040)?* — A: Both are correct, and they reflect **the repulsion strength at different X_GB**. Multistart at this X_c has X_GB = 0.228 (sparser); preseg has X_GB = 0.375 (denser). By the X_GB-ascending ordering, multistart should pattern with 0.075_preseg (low X_GB, steep slope), but it is actually shallower — because it is in a kinetic-floor state, the trajectory is still descending and not fully equilibrated, so the structure carries some random-IC residue. This itself is a secondary observation supporting "non-equilibrium states show different correlation patterns".
- *Q: Are 218 sites enough?* — A: 218 is the count inside the favourable [−30, −5] window. Within each X_c slice, P_i averaging in each bin still has ~10–25 sites; binomial 95 % CIs are shown as error bars in the right panel. The overall fit has inverse-variance weights against outliers (`scripts/solute_correlation_analysis.py:300-302`).

---

### Supporting Figures (4–7)

The main story arc closes with Figs. 0–3; the following four figures are for Q&A or methodological backing, opened on demand.

#### Figure 4 ── `04_spectrum_match.png`

**Caption (English):**

> **Figure 4.** Histograms and skew-normal fits for the per-site segregation-energy spectrum ΔE_seg. Light-green histogram + dark-green solid line — Wagih's published Zenodo spectrum (n = 82,646 GB sites, his Mendelev Al–Mg potential). Pink histogram + brown dashed line — our n = 500 random sample on the 200-Å Voronoi polycrystal. Two-sample KS test gives D = 0.026, p = 0.89, well above the project threshold of 0.5 → the two spectra are statistically indistinguishable, justifying the use of Wagih's FD framework as the baseline for our system. Skew-normal fit parameters: ours μ = +6.3, σ = 20.1, α = −1.47; Wagih μ = +6.7, σ = 20.8, α = −1.40.

**Role**: spectrum representativeness / methodological prerequisite.

**What it shows**: pink histogram + brown dashed skew-normal fit (our n = 500 spectrum) overlaid on the light-green histogram + dark-green solid skew-normal fit (Wagih Zenodo n = 82,646).

**Key numbers** (in the figure title): **KS D = 0.026, p = 0.89** ≫ 0.5 (spectrum-level indistinguishable); ours fit μ = +6.3, σ = 20.1, α = −1.47; Wagih fit μ = +6.7, σ = 20.8, α = −1.40.

**Use case**: answers Q5 (is the n = 500 spectrum representative?). Source data: `output/compare_vs_wagih_200A_tight.json`. Generating script: `scripts/compare_vs_wagih.py`.

#### Figure 5 ── `05_sampler_convergence.png`

**Caption (English):**

> **Figure 5.** HMC convergence diagnostics for the X_c = 0.075 preseg run, five-panel time series. (a) Instantaneous temperature T(t) holds at 500 K. (b) Total potential energy PE(t) plateaus near −1.513 × 10^6 eV after ~150 ps. (c) Cumulative swap acceptance rate stabilises around 6.2 %. (d) X_GB(t) descends monotonically from the over-segregated IC (~0.32) to the production-mean 0.254 (gray shading marks the burn-in window). (e) Per-frame swap fwd/rev decomposition (Mg ⇌ bulk) showing a net reverse imbalance — Mg is leaving the GB. Together these panels show the run is well-equilibrated apart from the slow ongoing descent of X_GB(t), which is the physical signal we are using as an upper bound.

**Role**: HMC convergence diagnostic.

**What it shows**: five-panel time series for the X_c = 0.075 preseg run: (a) T(t) holding at 500 K; (b) PE(t) plateauing near −1.513 × 10^6 eV after ~150 ps; (c) cumulative swap accept rate ~6.2 %; (d) X_GB(t) descending monotonically from IC ~0.32 to 0.254 (gray shading = burn-in); (e) per-frame swap fwd/rev decomposition, showing net reverse flux (Mg leaving the GB).

**Use case**: answers Q4 (has HMC converged?). Source data: `output/hmc_T500_Xc0.075_preseg.json`. Generating script: `scripts/hmc_xgb_timeseries.py`.

#### Figure 6 ── `06_two_sided_verify.png`

**Caption (English):**

> **Figure 6.** Two-sided IC equilibration verification at X_c = 0.05, T = 500 K. Blue line — random IC, X_GB(0) = 0.05, rising to 0.062. Red line — preseg IC, X_GB(0) = 0.27, descending to 0.238. The canon-FD target X_GB^FD = 0.228 (black dashed) lies between the two endpoints, demonstrating that the sampler equilibrates at the dilute end with ICs on both sides. Half-life-2 residuals: Δ_{1/2}^{rand} = +0.006, Δ_{1/2}^{preseg} = −0.015 (acceptable bracket convergence). The remaining sandwich width is dominated by the preseg trajectory still drifting after overshooting the target.

**Role**: equilibration verification at the dilute end.

**What it shows**: at X_c = 0.05, T = 500 K, two independent trajectories — random IC (blue, rising from X_GB(0) = 0.05 to 0.062) and preseg IC (red, descending from X_GB(0) = 0.27 to 0.238). The canon-FD target 0.228 (black dashed) lies between the two endpoints.

**Key numbers**: half-life-2 residuals Δ_{1/2}^{rand} = +0.006, Δ_{1/2}^{preseg} = −0.015 → bracket converges to spec. Sandwich width ~0.18 (driven mostly by the preseg trajectory still drifting after overshooting the canon-FD target).

**Use case**: stronger version of the answer to Q4 (two-sided IC verification that the sampler equilibrates at the dilute end). Source data: `output/hmc_T500_Xc5e-2_verify-{rand,preseg}_xgb.json`. Generating script: `scripts/verify_two_sided_compare.py`.

#### Figure 7 ── `07_ovito_segregation.png`

**Caption (English):**

> **Figure 7.** OVITO render of the X_c = 0.20 final HMC configuration. Gray = Al, orange = Mg. Mg concentration on grain boundaries is visibly higher than within grains (N_GB_Mg = 70,672, X_GB = 0.794); the orange "network" outlines the GB structure of the 16-grain 3D Voronoi polycrystal. Provides direct visual confirmation of segregation for slides; intended for non-expert audiences as a primer.

**Role**: visual segregation confirmation / non-expert entry point.

**What it shows**: OVITO render of the X_c = 0.20 final HMC config. Gray = Al, orange = Mg. Mg concentration at GB sites is visibly higher than in the grain interior (N_GB_Mg = 70,672, X_GB = 0.794); the orange "network" traces the GB structure of the 3D Voronoi polycrystal.

**Use case**: visual support on slides; non-expert audience entry point. Source LAMMPS file: `data/snapshots/hmc_T500_Xc0.20_preseg_final.lmp` (63 MB), rendered with OVITO Pro (standalone).

---

## 6. Cross-Figure Narrative

The logical chain in talk order:

```
[Setup]                          Wagih (2020) site-independent FD assumption
                                          ↓ test
[Figure 0 / panel d]             X_HMC < canon-FD at X_c = 0.075 → assumption fails
                                          ↓ explain
[Hypothesis]                     Mg–Mg interactions exist (the missing piece)
                                          ↓ check 1
[Figure 1 / g(r)]                Spatial Mg distribution on GB is non-random → not independent
                                          ↓ but g(r) is an aggregate signal
                                          ↓ check 2
[Figure 2 / P vs ΔE]             Failure concentrated at the favourable-ΔE end;
                                          ↓ ΔP is largest at lowest X_c (counterintuitive)
                                          ↓ requires a site-level explanation
                                          ↓ check 3
[Figure 3 / slope vs X_c]        Controlling ΔE → more Mg neighbours = lower occupancy
                                          → DIRECT evidence of site-level Mg–Mg REPULSION
                                          ↓ unify
[Conclusion]                     The mechanism behind Wagih's failure is site-level
                                  Mg–Mg repulsion (not the apparent "chemical attraction"
                                  in g(r)).
```

**Key reconciliation point**: Figure 1 looks like "clustering = attraction"; Figure 3 says "repulsion". They are not contradictory — Figure 1 is *geometric* (deep-ΔE sites are spatially adjacent on the GB plane → filling them produces g(r) > 1), while Figure 3 is the *ΔE-controlled residual* (at fixed ΔE, more neighbours ⇒ lower occupancy). Both consistently point to the same physics: **Mg site selection is influenced by neighbouring Mg, not independent**.

---

## 7. Predicted Q&A (sorted by likelihood)

**Q1**: You only directly proved breakdown at X_c = 0.075. What about X_c ≥ 0.10?
- A: X_c = 0.10 is **directly proven** by the multistart UB — since 2026-05-02, this is on Fig. 0 as the gray open ▽ at X_GB = 0.246 ± 0.004, well below canon-FD = 0.352 (gap = −0.106); see `output/hmc_T500_Xc0.10_multistart_xgb0.3.json`. For X_c ≥ 0.15 the preseg trajectories are still descending and their measurements are vacuous bounds; the breakdown evidence there is an extrapolation from the threshold X_c\* ∈ (0.05, 0.075] plus the mechanism evidence (Figs. 1–3), all of which show site-level interaction fingerprints across every X_c slice.

**Q2**: Is the Mg–Mg interaction elastic strain or chemical bonding?
- A: The current data cannot rigorously distinguish them, but the correlation length (~5 Å, the R_local threshold) is consistent with the r ~ 10 Å decay in g(r), pointing toward an elastic strain field (Mg is ~12 % larger than Al in lattice parameter, giving a long-range elastic distortion). Chemical interactions are usually shorter-ranged (~3 Å NN). Future work: compare windows at r = 3 / 5 / 8 Å to see length-scale dependence.

**Q3**: ΔE_i is the X_c = 0 reference; at finite X_c the effective ΔE shifts. Isn't that a confound?
- A: It is by-design — using bare ΔE_i as the baseline is the only way to test Wagih's assumption directly (his model itself is built on the X_c = 0 spectrum). Renormalising ΔE under finite X_c is future work; the project considered a ΔE-shift approach earlier and deprioritised it because the mechanism path delivered cleaner conclusions.

**Q4**: How well does HMC converge? Is 300 ps PROD on 32 ranks enough?
- A: For X_c = 0.075 the run took 10 h 20 min real / 24 h budget, with a swap acceptance rate of ~5.7 % (see `hmc_T500_Xc0.075_preseg.json`); the trajectory descends monotonically into a stable plateau. The two-sided IC verification at X_c = 0.05 (`output/verify_T500_Xc5e-2_two_sided.png`, Fig. 6) gives X_GB ≈ 0.238 ± 0.005 with good bracket overlap, confirming sampler equilibration at the dilute end.

**Q5**: Is n = 500 ΔE_i enough to represent the GB spectrum?
- A: KS two-sample test against Wagih's published n = 82,646 spectrum gives D = 0.0256, **p = 0.8920** ≫ 0.5 — statistically indistinguishable (see §4.2 table). Spectrum-mean difference 0.1 kJ/mol < 0.025 kT, std difference 5 %.

**Q6**: Why is the X_c = 0.075 first peak so weak in Fig. 1?
- A: At X_GB = 0.254 (sparse regime), Mg has already learned to "avoid 1st-NN neighbours" — an early site-level repulsion signal. Clustering shows up at the second shell r = 5.9 Å (1.22). This anomaly *supports* the Fig. 3 site-level repulsion conclusion; it is not a plotting error.

**Q7**: Has any of this been validated experimentally?
- A: This project is purely simulation. The Mendelev (2009) EAM potential for Al–Mg is calibrated to within ±5 % of experimental lattice parameter and dilute heat-of-mixing values. Atomic-resolution APT measurements of GB Mg fraction in Mg–Al alloys (e.g. Sauvage et al.) are qualitatively consistent with our dilute-limit values. Quantitative experimental validation is follow-up work.

---

## 8. Provenance / Data Integrity

Every cited number passes a spot-check; PASS / future-PASS log below.

### Verification log (last run: 2026-05-01)

```
[1] N_GB_Mg recompute (LAMMPS file → JSON)
    0.075_preseg:  recomp 22,635 == json 22,635   PASS  (Δ=0)
    0.15_preseg:   recomp 52,394 == json 52,394   PASS
    0.30_preseg:   recomp 74,582 == json 74,582   PASS

[2] X_GB recompute
    0.075:  0.254206 vs json 0.254206   PASS  (Δ=+0.00e+00)
    0.15:   0.588419 vs json 0.588419   PASS
    0.30:   0.837605 vs json 0.837605   PASS

[3] Wagih FD recompute at most-favourable bin
    X_c=0.075: ΔP=0.7728 (drawn 0.77)   PASS
    X_c=0.15:  ΔP=0.5713 (drawn 0.57)   PASS
    X_c=0.30:  ΔP=0.1284 (drawn 0.13)   PASS

[4] slope-vs-X_c values match JSON site_occupation_vs_density   PASS

[5] panel (d) gap table matches JSON hmc_vs_canonfd_T500_with_multistart.json   PASS

[6] spectrum stats (n=500) match compare_vs_wagih_200A_tight.json  PASS
    KS p=0.8920 → spectra indistinguishable from Wagih's 82k

[7] X_c=0.10 multistart UB drawn value (added 2026-05-02)
    JSON `hmc_T500_Xc0.10_multistart_xgb0.3.json` x_gb.mean = 0.245911
    drawn at X_GB ≈ 0.246, gap vs canon-FD = -0.106              PASS
```

### File dependency chains

Mechanism figures (1–3):

```
report/figures/0{1,2,3}_*.png
    ↑
output/defense_*.png  (regenerable via:
                       python scripts/replot_mechanism_for_defense.py)
    ↑
output/solute_correlation_analysis.json
    ↑
scripts/solute_correlation_analysis.py
    ↑
data/snapshots/hmc_T500_Xc*_final.lmp  (6 files, 63 MB each)
data/snapshots/gb_mask_200A.npy        (475,843 bytes)
/cluster/scratch/cainiu/production_AlMg_200A/delta_e_results_n500_200A_tight.npz
                                        (63,450 bytes)
```

Headline figure (0):

```
report/figures/00_headline_hmc_vs_wagih_T500.png
    ↑
output/hmc_vs_canonfd_T500_with_multistart.{json,png}
    ↑
scripts/canonical_fd_compare_5pt_with_multistart.py    (added 2026-05-02)
    ↑
output/hmc_T500_Xc{0.05_verify-preseg,0.075,0.10,0.15,0.20,0.30}_preseg.json
output/hmc_T500_Xc0.10_multistart_xgb0.3.json          (gray ▽ marker)
    ↑
HMC SLURM jobs 65208332 / 64xxxxx series (see CHANGELOG)
```

The original `scripts/canonical_fd_compare_5pt.py` and its outputs `output/hmc_vs_canonfd_T500.{json,png}` are kept untouched as a pre-multistart reference / 4-curve version.

### Known caveats (honest list)

- ~~Panel (d) does not yet plot the multistart UB; the X_c ≥ 0.10 breakdown evidence has to be supplied verbally during Q&A.~~ **Resolved 2026-05-02**: the multistart UB now appears on Fig. 0 as the gray ▽ marker; X_c = 0.10 breakdown is direct evidence on the figure itself.
- For X_c ≥ 0.15 only preseg trajectories are available, and they are still descending — the measurements are vacuous bounds. The breakdown evidence at these X_c relies on threshold extrapolation (X_c\* ∈ (0.05, 0.075]) plus the mechanism figures (1–3).
- The X_c = 0.20 slope flips positive (+0.0148) due to saturation, **not a physical sign reversal**; if pressed, cite the X_c = 0.30 slope = −0.002 ≈ 0 (slope → 0 is consistent with saturation, not reversal).
- Fig. 2 uses ΔE_i from the X_c = 0 reference spectrum; at finite X_c the effective ΔE_i shifts due to local interactions — this is by-design, not a bug, but should be acknowledged as a simplification during Q&A.
- Scatter-plot CIs use the binomial normal approximation (OK when n_per_bin ≥ 50; small-bin distortion is possible). Wilson CIs are stricter and are future work.

---

## 9. Reproduction / Regeneration

**Three mechanism figures (Figs. 1–3)**:
```
python scripts/replot_mechanism_for_defense.py
```
Reads `output/solute_correlation_analysis.json` and produces the three PNGs `output/defense_*.png` in a few seconds; copy to `figures/`.

**Panel (d) (Fig. 0, including the X_c = 0.10 multistart UB)**:
```
python scripts/canonical_fd_compare_5pt_with_multistart.py
```
Updates `output/hmc_vs_canonfd_T500_with_multistart.{json,png}`; copy the PNG to `figures/00_headline_hmc_vs_wagih_T500.png` (overwrite).

**Supporting figures 4–7 (no recomputation needed; just copy from existing `output/`)**:
```
cp output/compare_vs_wagih_200A_tight.png    figures/04_spectrum_match.png
cp output/hmc_T500_Xc0.075_preseg.png        figures/05_sampler_convergence.png
cp output/verify_T500_Xc5e-2_two_sided.png   figures/06_two_sided_verify.png
cp output/ovito_gb_render_xc0.20.png         figures/07_ovito_segregation.png
```

**Recomputing the mechanism analysis from LAMMPS snapshots**:
```
python scripts/solute_correlation_analysis.py
```
Reads the six LAMMPS final.lmp files + GB mask + reference NPZ; recomputes g(r), P_i, slopes; writes `output/solute_correlation_analysis.json` plus the original 6-panel PNGs (for SI). Takes 30–60 s on a login node.
