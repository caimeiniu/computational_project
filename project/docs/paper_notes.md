# Technical Notes: Wagih et al. (2020)

**Paper**: "Learning grain boundary segregation energy spectra in polycrystals"
**Authors**: Malik Wagih, Peter M. Larsen, Christopher A. Schuh
**Journal**: Nature Communications 11:6376 (2020)
**DOI**: https://doi.org/10.1038/s41467-020-20083-6
**Code/Data**: https://doi.org/10.5281/zenodo.4107058

---

## 1. Polycrystal Generation

- **Software**: Atomsk (ref 71) for Voronoi tessellation
- **Box size**: 20 × 20 × 20 nm³
- **Grains**: 16 randomly oriented grains
- **Annealing procedure** (LAMMPS, ref 68-69):
  1. Thermally anneal at 0.3–0.5 × T_melt for 250 ps
  2. Nose-Hoover thermostat/barostat, timestep 1 fs
  3. Purpose: relax grain structure without exaggerated grain growth
  4. Slow cool to 0 K at 3 K/ps
  5. Final conjugate gradient energy minimization
- **Result**: ~10⁵ GB sites for a (50 nm)³ polycrystal with 10 nm grain size
- **Lattice correction**: scale cell to match equilibrium lattice parameter of each potential
  (e.g., Ni polycrystal: fitted to Ni lattice a = 3.518 Å,
   but Ni(Al) potential fitted to a = 3.520 Å → must rescale)

## 2. GB Site Identification

- **Method**: adaptive common neighbor analysis (a-CNA, ref 72)
- **Classification**: all atoms NOT identified as bulk crystal (FCC/BCC/HCP) are GB atoms
- **Software**: OVITO (ref 70) for visualization and structure identification

## 3. Segregation Energy Calculation

- **Definition**: ΔE_seg = E_GB^solute - E_bulk^solute
  - E_GB^solute: energy with solute at GB site i (relaxed)
  - E_bulk^solute: energy with solute at a bulk (intra-grain) site (relaxed)
  - Negative ΔE → favorable segregation (system lowers energy)
- **Procedure for each GB site**:
  1. Substitutionally place ONE solute atom at the GB site
  2. Conjugate gradient energy minimization
  3. Reference: solute in a 6 nm sphere of pure solvent (in the polycrystal)
     → avoids long-range interactions with other GB atoms
  4. All calculations at 0 K (static)
- **Note**: the reference bulk site E_bulk^solute is chosen far from any GB
  to isolate the enthalpic portion of segregation

## 4. Interatomic Potentials

- **Source**: NIST Interatomic Potentials Repository (refs 21-22)
- **Types**: EAM (Embedded Atom Method) for metallic alloys
- **Total**: 259 binary alloys surveyed
- **Important caveat**: for alloys with multiple available potentials,
  predictions can differ significantly (e.g., Al(Ni) system,
  Supplementary Fig. 3). The paper reports ALL potentials as a
  conservative choice.

## 5. Feature Extraction: SOAP

- **Method**: Smooth Overlap of Atomic Positions (SOAP, ref 30)
- **Software**: QUIP/GAP (ref 28)
- **Parameters**:
  - Cutoff radius: r_cutoff = 6 Å
  - n_max = l_max = 12 (max radial/angular quantum numbers)
  - σ_at = 1 Å (Gaussian width)
  - Feature vector dimension: F^SOAP = 1,015
- **Input**: pre-segregation LAE (local atomic environment) of PURE SOLVENT
  → the descriptor is computed on the undecorated GB site (before placing solute)
- **Key insight**: segregation energy depends on the LOCAL GEOMETRY of the site,
  not on the solute identity (solute info enters only through the potential)

## 6. High-Fidelity ML Model

- **Algorithm**: linear regression (min_w ||Xw - y||²)
- **Features**: full SOAP vector (1,015 features)
- **Training**: 50/50 holdout split, sample size > 10 × F^SOAP (~10⁴ sites)
- **Performance**: MAE typically < 6 kJ/mol, often < 1 kJ/mol
  - Al(Mg): MAE = 2.4 kJ/mol
  - Ag(Ni): 3.4, Cu(Zr): 5.5, Fe(Al): 1.4, Ni(Cu): 0.8, Pt(Au): 0.7, Zr(Ni): 12.8

## 7. Accelerated ML Model

- **Dimensionality reduction**: PCA on SOAP features → P = 10 principal components
  (captures > 99% variance)
- **Training point selection**: k-means clustering (k = P × 10 = 100 clusters)
  → train on cluster centroids only (100 data points)
- **Performance**: Al(Mg) MAE = 4.2 kJ/mol (vs 2.5 for high-fidelity)
  → 2 orders of magnitude fewer training points, minimal accuracy loss

## 8. Equilibrium Segregation (Thermodynamics)

- **Segregation isotherm** (Eq. 2-4):
  - Total solute: X^tot = (1 - f^gb) X^c + f^gb X̄^gb
  - X^c: bulk concentration, X̄^gb: average GB concentration
  - f^gb: GB site fraction (depends on grain size)
  - Fermi-Dirac at each site: occupation ∝ [1 + ((1-X^c)/X^c) exp(ΔE/kT)]^{-1}
- **Segregation spectrum** fitted to skew-normal function:
  - F(ΔE) = (1/√(2π)σ) exp(-(ΔE-μ)²/2σ²) × erfc(-α(ΔE-μ)/√2)
  - Three parameters: μ (characteristic energy), σ (width), α (shape/skewness)
- **Temperature**: T = 600 K used for equilibrium predictions (Fig. 5)
- **Grain size**: 15 nm (f^gb ≈ 10%)

## 9. Key Results to Reproduce

| Alloy | ΔE_seg range (kJ/mol) | μ | σ | α | GB solute at 5% total |
|-------|----------------------|---|---|---|----------------------|
| Al(Mg) | [-60, +40] | ~-2 | ~4 | ~-0.4 | ~30% |
| Ag(Ni) | [-75, +50] | - | - | - | ~55% |
| Fe(Al) | [-30, +30] | - | - | - | ~5% |
| Ni(Cu) | [-30, +30] | ~-2 | ~8 | ~0.0 | ~15% |
| Pt(Au) | [-25, +25] | - | - | - | ~10% |

(Values approximate from Figs. 2, 4, 5)

---

## 10. What WE Need to Do (Our Project)

### Phase 1: Atomistic MC segregation simulation (LAMMPS)
1. Generate nanocrystalline structure (advisor has codes)
2. Pick binary alloy + download EAM potential from NIST
3. Identify GB sites (a-CNA via OVITO or ASE)
4. Implement MC swap simulation:
   - Random swap of solute ↔ solvent atom
   - Energy minimization after each swap
   - Metropolis acceptance: accept if ΔE < 0, else accept with prob exp(-ΔE/kT)
   - Run until equilibrium → measure GB solute concentration
5. Compare equilibrium GB concentration to Wagih's predictions (Fig. 5)

### Phase 2: UMA MLIP as drop-in replacement
- Replace LAMMPS EAM calculator with UMA
- Same MC procedure
- Compare: does UMA give similar/better segregation energies?
- Advantage: UMA is universal → no need for system-specific EAM fitting

### Key questions to answer
- Does MC atomistic simulation agree with Wagih's spectral predictions?
- How many MC steps needed for convergence?
- Does grain size / polycrystal realization affect the result?
- Can UMA match EAM accuracy for segregation energies?

### Tools needed
- **Structure generation**: advisor's code (or Atomsk + LAMMPS annealing)
- **LAMMPS**: MC + energy minimization (or ASE with LAMMPS calculator)
- **OVITO / ASE**: GB site identification (a-CNA)
- **Python**: analysis, plotting, ML (scikit-learn for SOAP+regression if needed)
- **Zenodo data**: https://doi.org/10.5281/zenodo.4107058 (paper's code + database)
