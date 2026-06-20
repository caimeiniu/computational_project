# Grain-boundary segregation beyond the dilute limit

Spectrum-based models predict grain-boundary (GB) solute segregation from a
distribution of single-site segregation energies, converted to a
finite-temperature GB composition with an independent-site Fermi-Dirac (FD)
model. This project tests that pipeline at two stages: how well the energy
spectrum can be learned, and whether the independent-site assumption survives
at finite solute concentration, where solute-solute interactions matter.

## Components and folders

| Component | Folder |
|---|---|
| Accelerated spectrum learning / training-site selection | [`project/point_selection/`](project/point_selection/) |
| Pt(Au): same workflow on a second alloy | [`project/PtAu/`](project/PtAu/) |
| Al(Mg): independent-site FD vs finite-concentration hybrid MC/MD (HMC) | [`project/AlMg/`](project/AlMg/) |

## Main results

- **Al(Mg).** The n=500 segregation-energy spectrum matches the Wagih *et al.*
  reference (two-sample KS p = 0.89). Despite this, finite-concentration HMC at
  X_c = 0.10 gives GB Mg fractions **below** the canonical closed-box FD
  prediction at every temperature (X_GB = 0.186 / 0.206 / 0.201 vs FD 0.352 /
  0.296 / 0.273 at 500 / 700 / 800 K). Site-resolved analysis traces the gap to
  Mg-Mg interactions on the boundary that the independent-site model omits.
- **Pt(Au).** At 700 K the HMC matches closed-box FD in the dilute limit and
  rises **above** it at X_total = 0.03 (X_GB ~ 0.067 vs FD 0.055), the opposite
  sign to Al(Mg). The direction of the independent-site FD error is therefore
  system-dependent.

## Report

The written report (separate Overleaf document) covers all three components.
Each folder above holds the code and data for its part.
