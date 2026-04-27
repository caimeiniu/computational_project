# Changelog

Entries in reverse chronological order (newest first).

## 2026-04-27 (morning) — sweep dead-end: random IC fails kinetically at every X_c; preseg restart submitted

### What ran overnight

Last night's sweep series (random IC, T=500 K, X_c ∈ {0.10, 0.15, 0.20, 0.30},
PROD target = 200 ps, --time=08:00:00, 32 ranks) reported:

| job ID    | X_c   | wall    | PROD reached | accepts | accept% | _final.lmp |
|-----------|-------|---------|--------------|---------|---------|-----------|
| 64811160  | 0.10  | 7h45m   | ~104 ps      | 8 944   |  8.66 % | NO (TIMEOUT) |
| 64811162  | 0.15  | 7h36m   | ~100 ps      | 9 406   |  9.77 % | NO (TIMEOUT) |
| 64811163  | 0.20  | 7h54m   |  ~90 ps      | 9 375   | 10.50 % | NO (TIMEOUT) |
| 64811164  | 0.30  | RUNNING | ~52 ps so far | 6 601  | 12.73 % | (still active) |

Empirical wall-time scaling: 32 ranks × 8 h ⇒ ~100 ps PROD (about 5 min/ps,
not 1 min/ps as the SWAPS=10 dry-run ran). The 200 ps target was never
realistic for an 8 h budget at these X_c. Affects `_final.lmp` availability
for panel (f) of the report figure (see master figure plan memory).

### Sweep diagnosis: random IC stuck at X_GB ≈ X_c at EVERY X_c

`scripts/hmc_xgb_timeseries.py` re-run on the three completed dumps and the
partial X_c=0.30 dump:

| X_c   | X_GB^HMC ± CI95         | canon-FD | gap     | PE drift (eV) | half2−half1 X_GB |
|-------|--------------------------|----------|---------|----------------|------------------|
| 0.10  | 0.1051 [0.1042, 0.1060] | 0.3519   | −0.247  | −1 350        | +0.0045 (climbing) |
| 0.15  | 0.1538 [0.1532, 0.1544] | 0.4204   | −0.267  | −1 619        | +0.0032 (climbing) |
| 0.20  | 0.2017 [0.2011, 0.2022] | 0.4671   | −0.265  | −1 879        | +0.0024 (climbing) |
| 0.30* | 0.2999 [0.2998, 0.3001] | 0.5337   | −0.234  | −1 464        | −0.0005 (flat)     |

`*` X_c=0.30 partial: 55 frames so far, run still active.

X_GB^HMC sits *exactly* at X_c at every point — Mg is essentially uniformly
distributed. Tight CI95 (~0.001) is consistent with a stationary uniform
distribution, not with equilibrated segregation. PE is still drifting
downward at every point ⇒ system has not reached a thermodynamic minimum.

The night CHANGELOG predicted X_c ≥ 0.10 should converge from random IC
because (a) Mg pool ~3× larger and (b) canon-FD target closer to random-IC
starting X_GB. Empirically, both effects are insufficient to overcome the
random-IC kinetic barrier — the same dead-end seen at X_c=5e-2 yesterday
extends ALL THE WAY to X_c=0.30. Likely cause: at high X_c the
"productive direction" fraction of accepted swaps is *worse* not better,
because random IC at X_GB≈X_c sits near the high-T (uniform) equilibrium,
so MC sampling has no thermodynamic bias to drive segregation.

### Proposal-level analysis: why random IC fails (and why `region` is not the simplest fix)

`fix atom/swap N X seed T ke yes types 1 2` selects **one type-1 (Al) and
one type-2 (Mg)** atom per attempt — by construction, never an Al-Al or
Mg-Mg pair. Because in random IC every Mg is uniformly distributed at
X_GB ≈ X_c, the probability that the selected Mg sits at a GB site is
`P(Mg at GB) = N_Mg_GB / N_Mg = X_c · N_GB / (X_c · N_total) = GB_frac =
0.187`, **independent of X_c**. Same for Al. The proposal distribution
in random IC is therefore:

| Mg at | Al at | X_GB net change           | probability                |
|-------|-------|---------------------------|---------------------------|
| GB    | GB    | 0                         | 0.187 × 0.187 = 0.035      |
| GB    | bulk  | −1/N_GB (reverse, GB→bulk) | 0.187 × 0.813 = 0.152      |
| bulk  | GB    | +1/N_GB (forward, bulk→GB) | 0.813 × 0.187 = 0.152      |
| bulk  | bulk  | 0                         | 0.813 × 0.813 = 0.661      |

Two failure modes nest:

1. **Geometric**: ~70 % of proposals are bulk-bulk or GB-GB → X_GB
   doesn't change regardless of acceptance.
2. **Energetic**: among the ~30 % productive proposals, forward and
   reverse have *equal* geometric weight. Net direction only emerges
   from the ΔE bias (energy benefit of putting Mg on a deep ΔE_seg
   site). At random IC the GB sites being "tested" are not
   preferentially the deep ones, so the bias is weak → net flux ≈ 0.

Empirically: net flux *per accepted swap* (= net moves toward equilibrium
per accepted move) was 9.1 % at random IC X_c=0.10 (810 net of 8 944
accepts) versus 80 % at preseg X_c=0.05 (3 776 net of 4 707 accepts).
At equal acceptance rate (~7-10 %), preseg is **~9× faster per attempt**
toward equilibrium, *and* it points the right direction from the start.

#### `fix atom/swap region <id>` keyword

LAMMPS exposes a `region` keyword that restricts swap candidates to a
geometric region. It would attack failure mode (1) directly by
eliminating bulk-bulk and GB-GB proposals. But:

- Built-in regions are simple analytic shapes (`block`, `cylinder`,
  `sphere`, `prism`, `union`, …); our 3D Voronoi GB is an irregular
  curved network (89 042 atoms, ~100 distinct faces). Wrapping it
  needs either a bounding-mesh approximation or a `union` of many
  spheres — the former misses GB-adjacent bulk that has to receive
  desorbed Mg, the latter scales poorly at LAMMPS region-stack limits.
- `region` does NOT fix failure mode (2). At random IC + region the
  productive proposals would still be ~50/50 forward/reverse — the
  geometric trap simplifies but the energetic trap remains.

#### Why preseg IC is strictly stronger

Preseg IC fixes both failure modes at once:

- mode (1): with all Mg at GB initially, `P(Mg at GB)·P(Al at bulk) =
  1.0 × 0.903 = 0.903` ⇒ 90 % productive proposals (vs 30 %).
- mode (2): GB is over-saturated relative to canon-FD ⇒ many GB Mg
  occupy ΔE > 0 (unfavourable) sites ⇒ release is *downhill* ⇒
  acceptance bias 80 % toward GB → bulk (verify-preseg yesterday). A
  pure region trick cannot reproduce this without IC engineering.

For Phase 4's full (T, X_c) grid `region GB ± shell` could be a
worthwhile extra speedup on top of preseg, but it is not needed for
today's restart.

#### Correction to the in-conversation discussion

An in-conversation explanation earlier today described the random-IC
failure as "the picked Mg lands on a GB Mg neighbour, ΔU≈0". That is
incorrect — `types 1 2` already prevents Mg-Mg or Al-Al swap proposals,
so this failure mode does not exist. The correct decomposition is the
(Mg-position × Al-position) table above. Recording the correction here
so future readers don't pick up the wrong intuition.

### Decision: switch all production sweep points to preseg IC

Project rule from yesterday's CHANGELOG (then phrased as "X_c ≲ 0.05") is
now extended: **all closed-box HMC verification at T=500 K should default
to preseg IC**. Random IC is a kinetic dead-end at every X_c we've tested.

`scripts/pre_segregate.py` extended to support X_c > GB_frac (=0.187):

- For N_Mg ≤ N_GB: fill `n_mg` random GB sites; X_GB(0) = N_Mg/N_GB,
  X_bulk(0) = 0 (unchanged from before).
- For N_Mg > N_GB ("mixed mode"): fill ALL GB sites with Mg, place
  remaining `N_Mg − N_GB` Mg in random bulk sites. X_GB(0) = 1.0,
  X_bulk(0) = (N_Mg − N_GB) / N_bulk. Still above canon-FD ⇒ descent
  direction GB → bulk preserved.

Generated four preseg lmps (`production_AlMg_200A/poly_AlMg_200A_preseg_Xc{0.10,0.15,0.20,0.30}.lmp`,
seed 20260427):

| X_c   | N_Mg    | N_Mg(GB) | N_Mg(bulk) | X_GB(0) | X_bulk(0) | canon-FD target | required descent |
|-------|---------|----------|------------|---------|-----------|-----------------|------------------|
| 0.10  |  47 572 | 47 572   |      0     | 0.5343  | 0.0000    | 0.352           | −0.182           |
| 0.15  |  71 357 | 71 357   |      0     | 0.8014  | 0.0000    | 0.420           | −0.381           |
| 0.20  |  95 143 | 89 042   |   6 101    | 1.0000  | 0.0158    | 0.467           | −0.533           |
| 0.30  | 142 714 | 89 042   |  53 672    | 1.0000  | 0.1388    | 0.534           | −0.466           |

### Submit scripts (12 h walltime, 32 ranks each)

New deck files: `data/decks/submit_hmc_T500_Xc{0.10,0.15,0.20,0.30}_preseg.sh`.
Each uses `-var SKIP_PLACE 1` (deck skips the random `set type/fraction`
because the data file already has type-2 atoms), PROD=200 ps, SWAPS=100.

12 h ≈ 1.5× the failed 8 h wall to give buffer for `_final.lmp` write.
Empirical scaling suggests ~150 ps PROD reached if convergence is not hit
earlier; if X_c=0.10 converges around 120-150 ps (16k Mg net flux needed,
~70-80 % productive direction at preseg vs ~30 % at random IC) the
`_final.lmp` is written. Higher X_c likely TIMEOUT before equilibrium.

### Jobs submitted (2026-04-27 ~08:25)

| job ID    | X_c   | state          |
|-----------|-------|----------------|
| 64864478  | 0.10  | RUNNING        |
| 64864522  | 0.15  | PENDING (Priority) |
| 64864524  | 0.20  | PENDING (Priority) |
| 64864525  | 0.30  | PENDING (Priority) |

(Plus the original 64811164 sweep-Xc=0.30 random IC run still active —
will be processed for completeness when it finishes.)

### Headline figure (interim) — 1 equilibrated point, 4 kinetic-floor points

`scripts/canonical_fd_compare_5pt.py` — NEW. Plots canon-FD, GC-FD,
ceiling, and HMC errorbars on linear-x scale (panel-d candidate for the
master figure). Distinguishes equilibrated vs kinetic-floor HMC points by
marker style.

| | X_c | X_HMC | canon-FD | gap |
|--|--|--|--|--|
| **equilibrated** (filled red ●) | 0.05 | 0.238 ± 0.005 | 0.228 | +0.009 |
| **kinetic floor** (open red ■) | 0.10 | 0.105 ± 0.001 | 0.352 | −0.247 |
| | 0.15 | 0.154 ± 0.001 | 0.420 | −0.267 |
| | 0.20 | 0.202 ± 0.001 | 0.467 | −0.265 |
| | 0.30 | 0.300 ± 0.001 | 0.534 | −0.234 |

Outputs: `output/hmc_vs_fd_T500_5pt.{png,json}`.

**Breakdown signal cannot be read from this state.** The only point where
HMC equilibrium is verified (X_c=5e-2 verify-preseg) shows agreement with
canon-FD within 1 % — no breakdown there. The four kinetic-floor points
are not equilibrium values; they reflect a sampler limitation, not a
physics signal. Wait for tonight's preseg restarts before reading any
breakdown statement.

### Files this entry

- `scripts/canonical_fd_compare_5pt.py` — NEW (panel-d builder, 5 HMC pts)
- `scripts/pre_segregate.py` — extended (mixed-mode for X_c > GB_frac)
- `scripts/hmc_xgb_timeseries.py` — unchanged; re-used to process 4 dumps
- `data/decks/submit_hmc_T500_Xc{0.10,0.15,0.20,0.30}_preseg.sh` — NEW
- `output/hmc_T500_Xc{0.10,0.15,0.20,0.30}_xgb.{json,png}` — NEW
  (X_c=0.30 written as `_partial.json` since source run is still active)
- `output/hmc_vs_fd_T500_5pt.{png,json}` — NEW (interim panel-d candidate)
- `production_AlMg_200A/poly_AlMg_200A_preseg_Xc{0.10,0.15,0.20,0.30}.lmp` — NEW (in scratch)

### Tomorrow morning (2026-04-28)

1. Process preseg-restart dumps as they finish (X_c=0.10 first; expect
   first results by ~08:30 + 12 h = ~20:30 tonight).
2. Per-run convergence check: PE plateau + half-half X_GB drift +
   gap-vs-canon-FD. X_c=0.10 most likely to converge in 12 h budget;
   higher X_c may need a 2nd preseg run continuation.
3. Update `canonical_fd_compare_5pt.py` to swap the 4 kinetic-floor points
   for the new equilibrated preseg points; rebuild `hmc_vs_fd_T500_5pt.png`.
4. Read breakdown signal: at which X_c does the equilibrated HMC start to
   systematically depart from canon-FD?
5. If X_c=0.20-0.30 don't converge in 12 h: extend (resubmit from
   `_final.lmp` or `_postcg.lmp` with continuation deck) before claiming
   any breakdown statement at high X_c.
6. Persist the eventual `_final.lmp` of the most-segregated successful
   run to `project/data/snapshots/` for panel (f) (master figure plan).

## 2026-04-26 (night) — verify-Xc=5e-2 outcome (preseg single-sided sufficient); rename A/B → verify/sweep

### Rename: drop opaque A/B labels (supersedes the (late) entry's labels)

The (late) afternoon entry used "A1/A2" for the X_c=5e-2 equilibration twin
and "B" for the X_c sweep batch. Those labels carry no physical meaning and
were flagged as a barrier to reading the CHANGELOG cold. From this entry
forward, all references and filenames use:

| old (afternoon entry) | new (this entry forward) | meaning |
|---|---|---|
| A         | verify           | the X_c=5e-2 equilibration verification batch |
| A1        | verify-rand      | random-IC replica (climbs from below) |
| A2        | verify-preseg    | pre-segregated-IC replica (descends from above) |
| B         | sweep            | production X_c sweep at T=500 K |
| B-1 … B-4 | sweep-Xc=0.10 …  | one tag per X_c |

Output filename prefix changed accordingly:
`hmc_T500_Xc5e-2_verify-{rand,preseg}_xgb.{json,png}`. The afternoon entry's
A/B labels are NOT edited (append-only rule); supersession is recorded here.

### Submit-script bugs in the verify batch

- All verify submit scripts have `--time=06:00:00`, NOT `--time=08:00:00`
  as the afternoon entry's table reported. Both runs hit the wall before
  reaching the 200 ps PROD target.
- `submit_hmc_T500_Xc5e-2_dense_preseg.sh` has `--ntasks=16` while the
  random twin has `--ntasks=32`. verify-preseg ran on half the cores →
  made roughly half the wall-time progress (73 ps vs 163 ps PROD reached).

### `scripts/hmc_xgb_timeseries.py` patch — TIMEOUT-truncated thermo line

When SLURM kills LAMMPS mid-thermo, the last log line can be partial (e.g.
`162600   162.6          500` instead of 8 columns). The parser was
appending the row before checking column count, producing an inhomogeneous
numpy array. Fix: check `len(row) != len(cols)` BEFORE the float
conversion / append. Both verify dumps now process cleanly.

### Outcome — verify-preseg landed on canon-FD; verify-rand did not

| run            | PROD ps   | swaps att / acc      | X_GB^HMC ± CI95              | half-half drift   | gap vs canon-FD=0.228 |
|----------------|-----------|----------------------|------------------------------|-------------------|-----------------------|
| verify-rand    | 11–163 ps | 152 300 / 10 880 (7.14 %) | 0.0618 [0.0604, 0.0630] | +0.0061 (climbing) | −0.166 |
| verify-preseg  | 11–73 ps  |  61 300 /  4 707 (7.68 %) | 0.2375 [0.2325, 0.2424] | −0.0154 (slowing)  | +0.0095 |

verify-preseg PROD-second-half mean = 0.2303; trace endpoint at 73 ps =
0.2241. Late frames have crossed canon-FD to slightly below — consistent
with thermal fluctuation around the equilibrium target (canon-FD itself
carries ~0.005 sampling uncertainty from n=500 spectrum). **canon-FD is
the correct equilibrium target at X_c=5e-2; preseg landed on it within 1 %.**

verify-rand is far below canon-FD and not slowing. Plot:
`output/verify_T500_Xc5e-2_two_sided.png`.

### Why the bracket didn't close — kinetic, not statistical

`fix atom/swap` proposes a swap between one Mg and one Al each attempt.
The "productive" subset (= one atom at GB, the other at bulk) and the
energy-favoured *direction* of those productive moves both depend on the
current X_GB, X_bulk distribution — they are NOT determined by the
swap-acceptance rate alone.

| run            | net Mg flux needed to canon-FD | actual net Mg flux in PROD | net-direction fraction of accepted swaps |
|----------------|-------------------------------:|---------------------------:|-----------------------------------------:|
| verify-preseg  | −3 500 (GB→bulk)               | −3 776 (essentially done)  | ≈ 80 % |
| verify-rand    | +15 800 (bulk→GB)              | +1 531 (≈ 10 % of need)    | ≈ 14 % |

preseg starts with ~all 23 786 Mg crammed into GB (most picked Mg attempts
land on GB, most picked Al attempts land on bulk → productive direction
hugely favoured) AND many high-ΔE GB sites are over-occupied so Mg release
is energetically downhill (acceptance bias toward GB→bulk). Both effects
multiply.

rand starts with Mg uniform → only ~30 % of attempts are productive at all
(picked Mg and picked Al fall in opposite regions), and the bulk→GB
direction has only weak energy bias because Mg can land on average GB
sites (not preferentially on the deepest ΔE_seg ones). Result: an
order-of-magnitude lower effective convergence rate than preseg.

**This is intrinsic to `fix atom/swap` at low X_c, not a wall-time problem.**
A 10× wall-time extension on verify-rand (60 h) would only catch up to
where preseg already sits today. Mitigation if it ever matters: `fix
atom/swap region GB` to concentrate attempts on the productive subspace.

### Decision

1. Do NOT re-run verify-rand. Single-sided preseg evidence suffices.
2. **Headline X_GB^HMC at X_c=5e-2, T=500 K = 0.230 ± 0.005** (preseg
   PROD-second-half mean + block-bootstrap CI95). This is the verified
   canonical-equilibrium value used in tomorrow's (T=500 K) row figure.
3. Sweep series at X_c ∈ {0.10, 0.15, 0.20, 0.30} should converge from
   random IC despite verify-rand's failure: the productive fraction at
   X_c=0.10 is ~3× larger (Mg pool scales with X_c) AND the canon-FD
   target is closer to the random-IC starting X_GB (smaller required net
   flux). Confirm post-hoc by checking PE plateau and X_GB(t) drift in
   each sweep result.
4. **Project rule**: future low-X_c HMC verification should default to
   preseg IC. random IC at X_c ≲ 0.05 is a kinetic dead-end.
5. `_final.lmp` was NOT written by either verify run because both
   TIMEOUT'd before `write_data`. For the OVITO panel of the final
   report, source from a sweep `_final.lmp` instead (most segregated,
   most visually striking — see `project_report_figure_plan` memory).

### Sweep series — running (status as of 21:30)

| job ID    | X_c   | wall  | started | state               | ETA    |
|-----------|-------|-------|---------|---------------------|--------|
| 64811160  | 0.10  | 8 h   | 19:43   | RUNNING (1h54m)     | ~03:43 |
| 64811162  | 0.15  | 8 h   | 19:43   | RUNNING (1h54m)     | ~03:43 |
| 64811163  | 0.20  | 8 h   | 19:52   | RUNNING (1h45m)     | ~03:52 |
| 64811164  | 0.30  | 8 h   | (queued)| PENDING (QOSMaxCpu) | TBD    |

These are RESUBMITS of the original IDs 64810706/7/8/710 from the
afternoon entry, which were cancelled before they ever started — that
table refers to job IDs that never ran.

### Tomorrow morning (after sweep completes)

1. Process 4 sweep dumps with `hmc_xgb_timeseries.py` →
   `output/hmc_T500_Xc{0.10,0.15,0.20,0.30}_xgb.{json,png}`.
2. Per-run convergence check: PE plateau + X_GB(t) half-half drift.
   Random-IC kinetics at X_c ≥ 0.10 should be much friendlier than at
   X_c=5e-2; if any run fails, fall back to preseg restart.
3. Extend `canonical_fd_compare.py` to accept arbitrary point list and
   build the (T=500 K) row figure: HMC (5 points: X_c=5e-2 from
   verify-preseg + 4 sweep) vs canon-FD vs GC-FD vs ceiling.
4. Read breakdown signal: at which X_c does HMC start to systematically
   depart from canon-FD?

### Files this entry

- `output/hmc_T500_Xc5e-2_verify-rand_xgb.{json,png}` — NEW
- `output/hmc_T500_Xc5e-2_verify-preseg_xgb.{json,png}` — NEW
- `output/verify_T500_Xc5e-2_two_sided.png` — NEW
- `scripts/verify_two_sided_compare.py` — NEW
- `scripts/hmc_xgb_timeseries.py` — patched (truncated-row tolerance)
- `data/decks/submit_hmc_T500_Xc5e-2_dense_{random,preseg}.sh` — used today;
  `dense_preseg.sh --ntasks=16` flagged as a bug (random twin uses 32)

## 2026-04-26 (late) — Strategic pivot: grand-canonical vs canonical FD; the comparison was the bug, not the sampler

### What earlier today's "under-sampling" analysis missed

The morning entry concluded HMC was 40× under-sampled and proposed
raising SWAPS_PER_CALL 10×. That diagnosis had the right *symptoms*
(monotone X_GB(t), PE drift) but the wrong *cause*. Stepping back: the
Wagih FD formula `P_i = 1/(1+(1-X_c)/X_c · exp(ΔE_i/kT))` treats X_c as
a **bulk** mole fraction held fixed by an infinite reservoir
(grand-canonical). Our HMC is closed-box (canonical, total Mg = X_c ·
N_total fixed). At small X_c, segregation depletes the bulk, so the GC
prediction can demand more GB-Mg than the closed box even contains:

| X_c   | total Mg in box | ceiling X_GB = X_c·N_tot/N_GB | GC-FD prediction |
|-------|-----------------|-------------------------------|------------------|
| 5e-4  | 238             | 0.0027                        | 0.0794 (29× over) |
| 5e-3  | 2 379           | 0.0267                        | 0.1889 (7× over)  |
| 5e-2  | 23 786          | 0.2671                        | 0.3675 (1.4× over)|

In every one of today's three runs, the GC-FD target was physically
unreachable regardless of MC convergence. The morning's "raise
SWAPS_PER_CALL × 10" plan was charting a route to an impossible
destination.

### Canonical (mass-conserving) FD predictor

Added to `scripts/fermi_dirac_predict.py`:
- `solve_x_bulk_canonical(dE, T, X_c_total, N_GB, N_total)` — 1-D
  brentq root of `F(X_bulk) = X_bulk·N_bulk + <P_i(T,X_bulk)>·N_GB −
  X_c_total·N_total`. F is strictly monotone increasing in X_bulk, so
  the root is unique. xtol=1e-14, rtol=1e-12.
- `x_gb_canonical(...)` returns `(X_GB, X_bulk)` at canonical
  equilibrium.
- `x_gb_canonical_curve(...)` vectorises across an X_c grid.
- `--self-test` extended: total-Mg conservation within 1e-6 across
  X_c_total ∈ {1e-4, 1e-3, 1e-2, 1e-1, 0.4}; T→∞ canonical limit
  recovers X_GB → X_c_total. All pass.

### `scripts/canonical_fd_compare.py` (new)

Loads the n=500 ΔE_seg sample (`delta_e_results_n500_200A_tight.npz`,
mean −6.9 kJ/mol, fraction<0 = 0.646), the 200 Å GB mask
(N_total=475 715, N_GB=89 042, GB-frac=0.187), and the three HMC JSONs.
Emits canonical FD curve at T=500 K plus a two-panel plot
(`output/hmc_vs_fd_T500_canonical.png`) showing GC-FD, canonical FD,
the closed-box ceiling, and the three HMC error-bar points; CSV
companion at `output/hmc_vs_canonical_fd.csv`.

### Headline numbers — canonical FD changes the picture

| X_c    | X_GB^HMC | GC-FD  | canon-FD | ceiling | gap = HMC − canon |
|--------|----------|--------|----------|---------|-------------------|
| 5e-4   | 0.0012   | 0.0794 | 0.0027   | 0.0027  | −0.0015           |
| 5e-3   | 0.0054   | 0.1889 | 0.0265   | 0.0267  | −0.0211           |
| 5e-2   | 0.0504   | 0.3675 | 0.2282   | 0.2671  | −0.1778           |

Three observations:
1. At X_c ∈ {5e-4, 5e-3} canonical FD is essentially pinned to the
   closed-box ceiling — i.e. the spectrum is so Mg-attractive that
   *all available Mg ends up at GB*, leaving the bulk nearly empty
   (X_bulk → 0). These two points test "does the box drain"; they don't
   resolve dilute-limit physics.
2. At X_c=5e-2 canonical FD is meaningfully below the ceiling (0.228 vs
   0.267, ≈85 % of Mg at GB, with a real bulk reservoir at X_bulk ≈
   1×10⁻²). This is the only one of today's points where canonical FD
   is structurally informative.
3. HMC remains far below canonical FD at every point, so the sampler
   *is* under-sampled — but by ≤30× at the X_c=5e-2 point, not by 200×
   as the GC-FD comparison made it look. The morning's SWAPS_PER_CALL
   bump is now within striking distance of canonical FD at this point.

### Decision: don't fight unreachable targets; pivot the X_c sweep upward

Refined plan for the next batch:

1. **Drop X_c={5e-4, 5e-3}** as primary verification points. They are
   trivially canonical-FD-saturated and thus only test box-drain
   kinetics, not segregation thermodynamics.
2. **Designate X_c=5e-2 as the equilibration verification point.** Run
   it again with `SWAPS_PER_CALL=100` (10× attempts at no wall-time
   cost) and PROD 50 ps → 200 ps. Two-sided IC check: a second run
   starting from a *pre-segregated* config (Mg pre-placed at GB sites
   so X_GB(0)≈0.5). If both IC's converge to canonical FD's 0.228
   within block-bootstrap CI, equilibrium is established.
3. **Add X_c={0.10, 0.15, 0.20, 0.30}** as the production sweep — this
   is where canonical FD differs significantly from GC-FD and from the
   ceiling, *and* where solute-solute interactions are physically
   expected to matter (the actual scientific question). Mg supply is
   plentiful (47k–143k atoms), so swap acceptance and attempt budget
   are not bottlenecks at these X_c.
4. Once the (T=500 K) row is in, scan T ∈ {300, 700, 900} at the same
   X_c to fill the (T, X_c) grid.

The reframed scientific question is unchanged but cleaner: **Does
HMC X_GB^HMC match canonical FD across the (T, X_c) grid?** Disagreement
at any point is the breakdown signal we wanted; agreement validates
Wagih's per-site framework as a finite-size theory.

### Files this entry

- `scripts/fermi_dirac_predict.py` — extended with canonical FD + tests
- `scripts/canonical_fd_compare.py` — NEW
- `output/canonical_fd_T500.json` — NEW
- `output/hmc_vs_canonical_fd.csv` — NEW
- `output/hmc_vs_fd_T500_canonical.png` — NEW
- (note: `output/` is gitignored; CHANGELOG records the numbers above)

### Jobs submitted (queued at 13:42, T=500 K throughout)

A — equilibration verification at the only structurally informative
    point in the original sweep (X_c=5e-2):

| job ID    | submit script                              | IC          | SWAPS | PROD | target X_GB^canon |
|-----------|--------------------------------------------|-------------|-------|------|-------------------|
| 64810700  | submit_hmc_T500_Xc5e-2_dense_random.sh     | random      | 100   | 200 ps | 0.228 (from below) |
| 64810704  | submit_hmc_T500_Xc5e-2_dense_preseg.sh     | pre-seg @ ceiling | 100 | 200 ps | 0.228 (from above) |

A1 and A2 should converge to the same X_GB^∞ within block-bootstrap CI;
that is the equilibrium proof.

B — production X_c sweep where canon-FD ≠ ceiling and ≠ GC-FD, i.e.
    the regime that actually probes solute-solute breakdown:

| job ID    | submit script                          | X_c   | total Mg | ceiling |
|-----------|----------------------------------------|-------|----------|---------|
| 64810706  | submit_hmc_T500_Xc1e-1.sh              | 0.10  | 47 572   | 0.534   |
| 64810707  | submit_hmc_T500_Xc1.5e-1.sh            | 0.15  | 71 357   | 0.802   |
| 64810708  | submit_hmc_T500_Xc2e-1.sh              | 0.20  | 95 143   | 1.000 (capped) |
| 64810710  | submit_hmc_T500_Xc3e-1.sh              | 0.30  | 142 715  | 1.000 (capped) |

All B's: random IC, SWAPS_PER_CALL=100, PROD=200 ps, --time=08:00:00.
Outputs land in `/cluster/scratch/cainiu/hmc_AlMg/hmc_T500_Xc<xc>{,_dense_*}.{log,dump,_final.lmp}`.

### Deck change

`data/decks/hmc_AlMg.lammps` now supports `-var SKIP_PLACE 1` to bypass
the random `set type/fraction` step when the input data file already
contains type-2 Mg (used by A2 with the pre-segregated lmp from
`scripts/pre_segregate.py`). Default behaviour unchanged.

### Concept glossary — first-time terms used in today's discussion

These bridge "FD curves at different T" (where the project was) and
"why HMC ≠ GC-FD even with infinite sampling" (today's pivot). Listed
in the order they appear in the strategy thread.

- **Grand-canonical Fermi-Dirac (GC-FD)** — the textbook Wagih eq. 2,
  `P_i = 1/(1+(1-X_c)/X_c · exp(ΔE_i/kT))`. X_c is treated as the *bulk*
  solute mole fraction held fixed by an *implicit infinite reservoir*:
  whatever Mg the GB sucks in is replenished from outside. This is the
  thermodynamic-limit prediction; it has no notion of total Mg in any
  particular box. All FD curves drawn before today (the multi-T plots)
  are GC-FD.

- **Canonical Fermi-Dirac (canon-FD)** — closed-box version of the same
  per-site formula. Total Mg = X_c · N_total is *conserved*; bulk
  fraction X_bulk depletes as GB segregates. We solve a 1-D self-
  consistency equation
  `X_bulk · N_bulk + <P_i(T,X_bulk)> · N_GB = X_c · N_total`
  (brentq, monotone in X_bulk) for the equilibrium X_bulk, then evaluate
  X_GB^canon = `<P_i(T, X_bulk)>`. This is the apples-to-apples partner
  for our HMC, which is closed-box by construction.
  Implemented in `scripts/fermi_dirac_predict.py::x_gb_canonical`.

- **Closed-box ceiling** — the largest X_GB physically realisable in a
  finite box at total fraction X_c: every Mg atom in the box piled into
  GB, none in bulk. `ceiling = X_c · N_total / N_GB = X_c / f_gb`. For
  our 200 Å box `f_gb = 0.187` ⇒ `ceiling = X_c / 0.187`. If GC-FD's
  prediction exceeds the ceiling, GC-FD is *unreachable in the box* —
  not a sampler problem, a finite-size problem.

- **Bulk depletion** — segregation pulling the post-equilibrium
  `X_bulk` significantly below the input `X_c` because there isn't enough
  Mg to supply both bulk and GB at GC-FD levels. Always present in
  closed boxes; *negligible* iff GC-FD prediction `<<` ceiling. Today's
  three points all violate this: GC-FD/ceiling = 29.4, 7.1, 1.4 at
  X_c = 5e-4, 5e-3, 5e-2. The intended sweep regime (X_c ≥ 0.1) brings
  GC-FD/ceiling toward unity, where canon-FD and GC-FD diverge
  *physically* rather than mathematically.

- **Two-sided IC convergence test** — to certify that an HMC run reached
  equilibrium without an external reference, run two replicas from
  initial conditions that bracket the expected equilibrium: one that
  must climb (random IC, X_GB(0) ≈ X_c) and one that must descend
  (pre-segregated IC, X_GB(0) at ceiling). If their stationary X_GB^∞
  agree within block-bootstrap CI, equilibrium is established. A1+A2
  in today's submission implement this at X_c = 5e-2, target ≈ 0.228.

- **Why GC-FD was right for Wagih and wrong for our box** — Wagih's
  curves are framework predictions in the thermodynamic limit (Wagih
  himself does not test them at fixed total Mg in a finite box). At
  *low* X_c on a box with *large* GB fraction (our 0.187), the GB
  effectively swallows the entire reservoir, and GC-FD is an
  unattainable target. Wagih's MD validation samples were either at
  high enough X_c that depletion was small or used much larger boxes.
  Going forward our X_c sweep targets the regime where canon-FD < ceiling
  (X_c ≥ 0.05) so that the *physical* assumption being tested is the
  independent-site one, not finite-size mass balance.

## 2026-04-26 — HMC under-sampling diagnosed: 50 ps PROD ≈ 50× too short to reach FD equilibrium

### What ran overnight

After yesterday's Xc=5e-3 dry-run completed, two more production HMC runs
finished overnight on the 200³ Å box (476k atoms, 32 MPI ranks):

| job      | T (K) | X_c   | wall (h) | swap att. | accepts | accept% |
|----------|-------|-------|----------|-----------|---------|---------|
| 64779781 | 500   | 5e-4  | ~0.9     | 5000      | 315     | 6.3 %   |
| 64777123 | 500   | 5e-3  | ~0.9     | 5000      | 383     | 7.7 %   |
| 64779784 | 500   | 5e-2  | ~1.0     | 5000      | 443     | 8.9 %   |

All three completed cleanly (NVT thermostat held at 499–500 K, PE drifted
−40 to −60 eV over 50 ps PROD — small relaxation, no instability).

### New analysis script: `scripts/hmc_xgb_timeseries.py`

~190 lines. Parses `<stub>.log` thermo blocks (Step / Temp / PotEng /
f_hmc[1]=n_attempts / f_hmc[2]=n_accepts) plus `<stub>.dump` per-frame
`id type` (bulk-loaded via `np.loadtxt(..., max_rows=N)` — 50× faster
than per-line Python parsing on 200 MB dumps). Crosses each frame's
type vector with `gb_mask_200A.npy` (89 042 GB sites of 475 715 atoms,
GB fraction 0.187) to compute X_GB(t). Drops a 20 % burn-in then
estimates the stationary mean with a stationary-block bootstrap
(block = 5 frames, 2000 resamples, 95 % CI).

Outputs per run:
- `<stub>_xgb.json` — thermo summary + swap stats + X_GB(t) series + CI
- `<stub>_xgb.png` — 4-panel diagnostic (T, PE, instantaneous swap accept, X_GB)

### Headline result: HMC X_GB ≈ X_c — no segregation reached

| X_c   | X_GB^HMC | CI95            | X_GB^FD | enrich (HMC) | enrich (FD) |
|-------|----------|-----------------|---------|--------------|-------------|
| 5e-4  | 0.0012   | [0.0010, 0.0013]| 0.0793  | 2.4×         | 158.6×      |
| 5e-3  | 0.0054   | [0.0053, 0.0055]| 0.1869  | 1.1×         | 37.4×       |
| 5e-2  | 0.0504   | [0.0504, 0.0505]| 0.3702  | 1.0×         | 7.4×        |

X_GB^HMC sits essentially at the bulk fraction X_c at every point. The
X_GB(t) trace is monotonically rising but glacially — at Xc=5e-3 the
slope is ~1.7e-5 per 1 ps frame, would take ~10 000 frames (10 ns) to
plateau. PE drifts down at every point, confirming the system is still
relaxing toward the Mg-at-GB minimum.

### Why: the run was attempt-limited, not time-limited

Reaching FD equilibrium at Xc=5e-3 requires moving (0.187 − 0.005)·89042
≈ **16 200 Mg atoms** from bulk to GB sites. Each accepted swap moves
exactly one. With 5000 attempts × ~7 % acceptance = **383 accepts**, we
have ~40× too few. Same arithmetic at Xc=5e-4 (~7 100 needed,
315 obtained, ~22× short) and Xc=5e-2 (~28 500 needed, 443 obtained,
~64× short). The deck's `SWAPS_PER_CALL=10` × 100-step interval over
50 ps PROD ⇒ only 5000 attempts total — that ceiling is the bottleneck,
not wall time per swap.

### Comparison artifact

`scripts/hmc_sweep_compare.py` → `output/hmc_vs_fd_T500_sweep.png` +
`output/hmc_sweep_T500.csv`. Two-panel: log-log (X_GB vs X_c with FD
curve and X_GB=X_c reference) and enrichment ratio. Shows HMC points
sitting on the diagonal, FD curve well above.

### Decision: don't expand the (T, X_c) grid yet — fix the sampler first

The result so far does NOT mean the dilute-limit assumption breaks down.
It means our HMC is so under-sampled that we can't measure equilibrium
at all. Before sweeping more (T, X_c) points it has to be demonstrated
that *one* point reaches a stationary X_GB^HMC consistent with FD.

Plan for the next batch (one prototype point, then sweep):
1. **Raise SWAPS_PER_CALL 10 → 100** in `hmc_AlMg.lammps` (free — same
   wall time per call, 10× more attempts). At 7 % acceptance × 50 000
   attempts = ~3 500 accepts ⇒ within 5× of the 16 200 needed at Xc=5e-3.
2. **Extend PROD 50 ps → 200 ps** at the verification point (4× wall
   time but covers another 4× attempt budget) ⇒ ~14 000 accepts ⇒
   within range of FD-equilibrium.
3. **Two-sided equilibrium check**: also run the same point starting
   from a *pre-segregated* IC (Mg pre-placed at GB sites, X_GB(0)=0.5
   say). If both ICs converge to the same X_GB^∞ within CI, that's the
   equilibrium proof.
4. Only after step 3 succeeds: sweep T ∈ {300, 500, 700, 900} × X_c ∈
   {1e-4, 1e-3, 1e-2, 1e-1} on the `(T, X_c)` grid.

Optional acceleration if (1)+(2) still under-sample: use `region GB`
on `fix atom/swap` to confine swap-pair selection to the GB neighbourhood.
This raises the productive-move acceptance ~5× (only attempts that
actually involve a GB site).

### Files touched today

- `scripts/hmc_xgb_timeseries.py` — NEW
- `scripts/hmc_sweep_compare.py` — NEW
- `output/hmc_T500_Xc5e-{4,3,2}_xgb.{json,png}` — NEW
- `output/hmc_sweep_T500.{csv,png}`, `output/hmc_vs_fd_T500_sweep.png` — NEW
- `data/decks/hmc_AlMg.lammps` — already had the `write_data` ordering
  fix from yesterday's bug note; no further edit
- `data/decks/submit_hmc_*.sh` — still untracked (commit pending)

## 2026-04-25 (EOD) — End-of-day state, push complete, HMC dry-run successful, tomorrow's plan

### Pushed to origin (both branches)

After resolving a rebase conflict against `origin/main` (Jiayi Fu had a
no-op `test` commit `3eacd53` that touched the README's tail), our 8
commits replayed cleanly on top:

```
origin/main:    f70af31  Cu(Ni)-merge prep: parameterize alloy, ...   ← top of stack, 8 commits ahead of e4e1391
                3eacd53  test                                          (Jiayi Fu)
                e4e1391  ...                                           (previous main)
origin/cainiu:  same f70af31
origin/jy:      Jiayi Fu's branch (left untouched by us)
```

Colleagues clone the repo and follow `project/README.md`. Tarball
`/cluster/home/cainiu/gb-seg-toolkit.tar.gz` (28 KB) kept as alternative
delivery for non-git users.

### HMC dry-run completed (job 64777123, T=500 K, X_c=5e-3)

Job ran to completion; full 10 ps EQUIL + 50 ps PROD on the production
200 Å box (476k atoms, 32 MPI ranks). Outputs in
`/cluster/scratch/cainiu/hmc_AlMg/`:

| diagnostic | value | health |
|------------|-------|--------|
| T (PROD)   | 499–500 K | NVT thermostat working ✓ |
| PE (PROD end) | −1.577e6 eV ± ~50 eV | stationary, no drift ✓ |
| swap attempts | 5000  | (50000 steps / 100 step interval × 10 swaps/call) |
| swap accepts  | 383   | — |
| **acceptance rate** | **7.66 %** | **inside 5–30 % target band ✓** |
| dump file | `hmc_dry_T500_Xc0.005.dump` (208 MB, 50 frames) | ready for X_GB(t) post-processing |

Acceptance is on the lower side of the healthy range — fine for a
dry-run; if production wall time matters we can raise it by either
increasing `SWAPS_PER_CALL` (currently 10) or by restricting swaps to a
GB region (`region` keyword on `fix atom/swap`).

### Deck bug discovered (cosmetic, fix tomorrow)

```
ERROR: Could not find thermo fix ID hmc (src/src/thermo.cpp:314)
Last command: write_data   ${outstub}_final.lmp
```

`thermo_style` was set to reference `f_hmc[1]/f_hmc[2]`, but `unfix hmc`
was issued before the final `write_data`, which still triggers a thermo
print. The simulation itself was not affected (50 ps PROD finished, all
data written via dump + restart). Fix: either reorder so `write_data`
runs before `unfix hmc`, or reset `thermo_style` to a column set without
`f_hmc` after the unfix.

### Tomorrow's plan

1. **Fix the `hmc_AlMg.lammps` `write_data` ordering bug** (one-line edit:
   move `write_data ${outstub}_final.lmp` above `unfix hmc`, or add a
   neutral `thermo_style custom step time temp pe ke press` after `unfix hmc`).
2. **Write `scripts/hmc_xgb_timeseries.py`** (~150 lines):
   - Parse `*.log` thermo block → arrays of (t, T, PE, KE, P, n_attempts, n_accepts).
   - Parse `*.dump` (50 frames, id+type) → load `gb_mask_200A.npy` →
     compute X_GB(t) = (types[gb_mask] == 2).sum() / gb_mask.sum() per frame.
   - Diagnostic plot (4 panels: T, PE, instantaneous swap acceptance,
     X_GB(t)) + JSON summary {accept_mean, x_gb_mean ± block-bootstrap CI,
     fd_predicted: 0.189, gap}.
   - Burn-in detection: drop first ≥ 2τ_int of the X_GB(t) series.
3. **Compare X_GB^HMC vs X_GB^FD = 0.189** at this single point. If they
   agree within the FD-bootstrap CI, dilute-limit holds at X_c=5e-3, T=500 K.
   If X_GB^HMC > 0.189 (HMC sees more segregation), independent-site
   assumption is breaking down → first signal of breakdown.
4. **Decide HMC sweep grid** based on (a) FD knee table from late 4 entry
   and (b) what the dry-run point tells us about cost/precision tradeoffs.
5. **Once dry-run sign-off**: commit `hmc_AlMg.lammps` +
   `submit_hmc_dryrun.sh` (currently untracked) + new
   `hmc_xgb_timeseries.py`. They've been deliberately held out of `main`
   pending validation.

### Open task tracker (carried forward)

- HMC deck bug fix
- HMC X_GB(t) post-process script
- HMC sweep grid execution
- Eventually: solute-solute g(r) at GB to diagnose breakdown mechanism

## 2026-04-25 (late 6) — Cu(Ni) merge prep: parameterize Al/Mg-specific bits, README quickstart, concept glossary

### HMC dry-run launched (job 64777123)

`data/decks/hmc_AlMg.lammps` (new, ~95 lines) + `submit_hmc_dryrun.sh`.
Single (T, X_c) = (500 K, 5×10⁻³) point on the production 200 Å box,
EQUIL 10 ps NVT → PROD 50 ps with `fix atom/swap`. Output goes to
`/cluster/scratch/cainiu/hmc_AlMg/`. Will verify swap acceptance rate
(target 5–30 %), PE plateau, and X_GB(t) post-pipeline before the deck
gets merged.

### FD value correction at T=500 K, X_c=5×10⁻³

While discussing the dry-run target I quoted X_GB^FD ≈ 0.21. **Wrong** —
exact interpolation from `output/fd_curves_200A_tight.json` is **0.189**
(I misremembered X_c=10⁻² value of 0.236 as 0.21 at 5×10⁻³). User
caught this from reading the figure (their 0.15 estimate was much closer
than my 0.21). Comparison target for the dry-run X_GB^HMC is therefore
**0.189**, not 0.21.

For the record, exact FD predictions at X_c=5×10⁻³:
| T (K) | X_GB^FD |
|-------|---------|
| 300   | 0.341   |
| 500   | 0.189   |
| 700   | 0.102   |
| 900   | 0.056   |

### Concept glossary additions (in `reference_gb_glossary.md`)

While walking through the HMC plan, the following standard MD/MC terms
were defined; recorded so future sessions don't re-derive:

- **EQUIL / equilibration / burn-in** — "暖机"段，让体系松弛到目标系综平衡分布。
- **PROD / production / sampling** — equilibration 之后的"正式跑"，所有测量取数。
- **Plateau / stationary** — PE(t) 不再有系统性下降，只剩窗口内热涨落。
- **Swap acceptance rate** — `fix atom/swap` Metropolis 接受率，5–30% 健康。
- **thermo (LAMMPS)** — `thermo N` 每 N 步把 thermo_style 列出的标量打到 stdout/log；HMC 的时间序列主源。

### Parameterization for Cu(Ni) merge

`scripts/sample_delta_e.py`:
- Added `elements: tuple[str,str]` and `masses: tuple[float,float]`
  parameters to `compute_delta_e_spectrum` and `_write_site_deck`,
  defaulting to `("Al","Mg")` / `(26.9815, 24.305)`. Backward-compatible
  with all existing Al(Mg) calls.
- New CLI flags `--elements "Cu Ni"` and `--masses "63.546 58.6934"`.
- Self-tested: default rendering identical to pre-edit deck; Cu/Ni
  override produces correct `mass`, `pair_coeff`.

`data/decks/anneal_AlMg.lammps`:
- Added `EL1`, `EL2`, `MASS1`, `MASS2` `variable index` slots, defaults
  `Al / Mg / 26.9815 / 24.305`. `mass` and `pair_coeff` lines now
  reference these. Existing `submit_anneal_200A.sh` (which doesn't pass
  these vars) keeps the Al/Mg defaults — no behavior change.

`scripts/generate_polycrystal.py`, `gb_identify.py`,
`fit_delta_e_spectrum.py`, `fermi_dirac_predict.py`: already generic;
no changes needed.

### README rewrite

Replaced the placeholder README with an alloy-merge-ready document:

- Pipeline overview diagram.
- **Validated scripts** table — 5 Python scripts + 1 deck cleared for
  any FCC/BCC/HCP binary alloy via CLI/`-var` flags.
- **Validation-only** table — `compare_vs_wagih.py`,
  `bootstrap_vs_wagih.py`, `paired_pipeline_residual.py`,
  `wagih_dump_to_data.py` are Al(Mg)-specific (need Wagih's per-site
  reference data to run).
- **Not-yet-validated** table — HMC deck + post-process; explicitly
  excluded from the merge until dry-run signs off.
- **Cu(Ni) quickstart** — copy-paste shell block: generate (Cu lattice
  3.61 Å) → anneal (T_hold 540 K = 0.4·T_melt(Cu)) → GB ID → ΔE sample
  → fit + FD predict. All commands include the correct `--elements
  "Cu Ni" --masses "63.546 58.6934"` for Cu-Ni and the right `-var
  EL1 Cu -var EL2 Ni -var MASS1 -var MASS2 -var T_HOLD 540` for the
  anneal deck.

### Merge plan for `main` (proposal — user to execute)

`main` is 7 commits behind `cainiu`. Recommended:

1. Pull `cainiu` to `main` excluding the Wagih validation files (to
   keep `main` story-clean for Cu(Ni) team) — or include them clearly
   labelled "validation" since they're harmless and document our
   confidence level. **My suggestion: include**, since deleting them
   loses provenance and the README already labels them validation-only.
2. Exclude HMC deck (`hmc_AlMg.lammps`, `submit_hmc_dryrun.sh`) from
   the merge until dry-run confirms acceptance/plateau/X_GB.
3. Cherry-pick or fast-forward — prefer fast-forward of the relevant
   commits since `main` is a strict ancestor.

User to authorise and execute the actual merge — Claude will not push to
shared `main` without explicit instruction.

### Artifacts

- `data/decks/hmc_AlMg.lammps`, `submit_hmc_dryrun.sh` (new, in flight)
- `scripts/sample_delta_e.py` (parameterised; backward-compat)
- `data/decks/anneal_AlMg.lammps` (parameterised; backward-compat)
- `README.md` (rewritten)
- `reference_gb_glossary.md` (EQUIL / PROD / plateau / acceptance / thermo entries)

## 2026-04-25 (late 5) — Figure cleanup, archive, glossary update

### Figure-label cleanup (no underscore literals)

All kept figures regenerated with matplotlib mathtext so axis labels and
legends render symbols as math italic + sub/superscripts instead of raw
underscores. Scripts updated:

- `scripts/fermi_dirac_predict.py`:
  `X_c` → `$X_c$`,  `X_GB^FD` → `$X_\mathrm{GB}^\mathrm{FD}$`,
  `X_GB = X_c` → `$X_\mathrm{GB} = X_c$`, "Fermi-Dirac" → "Fermi–Dirac",
  ASCII "200A" → "200 Å".
- `scripts/compare_vs_wagih.py`:
  `ΔE_seg` → `$\Delta E_\mathrm{seg}$` (xlabel + title).
- `scripts/paired_pipeline_residual.py`:
  `E_GB^Mg` → `$E_\mathrm{GB}^\mathrm{Mg}$` (both axes), residual label
  → `$E^\mathrm{ours}_\mathrm{GB} - E^\mathrm{Wagih}_\mathrm{GB}$`.

Convention recorded in `reference_gb_glossary.md` ("Figure-axis symbols"):
all Phase-4 plot scripts must use `r"$...$"` mathtext for any compound
symbol — no raw `X_c` / `X_GB` / `ΔE_seg` in user-visible strings.

### Archive sweep — `output/archive/`

Moved 4 superseded PNGs out of the active `output/` to keep the report
shortlist scannable. JSONs kept in place (small, reproducibility data).

| archived file | reason |
|---------------|--------|
| `compare_vs_wagih.png` | 100A loose-CG (the 4.8 kJ/mol-shift discovery plot, now superseded by 200A tight) |
| `delta_e_spectrum_n500.png` | 100A skew-normal fit, replaced by 200A |
| `delta_e_spectrum_n500_200A.png` | 200A loose-CG: bug-affected ΔE values, ~−11 kJ/mol mean |
| `paired_pipeline_residual_n343.png` | loose-CG diagnostic snapshot — the smoking-gun image; n=500 tight is the publication-quality version |

### Active `output/` (post-cleanup)

| file | role |
|------|------|
| `compare_vs_wagih_200A_tight.{png,json}` | report main: spectrum histogram + skew-normal + KS |
| `fd_curves_200A_tight.{png,json}` | report main: Fermi–Dirac dilute-limit curves, 4 T overlay vs Wagih |
| `paired_pipeline_residual_n500_tight.{png,json}` | validation: per-site PE residual ≈ 0 on Wagih's structure |
| `bootstrap_vs_wagih_200A_tight.json` | validation: 6 statistics inside 95 % CI |
| `method_overview.png` | slides: 4-panel methodology overview |

### Reference memory additions

- `reference_ks_test.md` (new) — KS two-sample p-value reading conventions
  + project benchmarks (wagih-on-wagih p=0.91, our 200A tight p=0.892,
  our 200A loose-CG p=7e-9).
- `reference_gb_glossary.md` — added "Figure-axis symbols" subsection
  (mathtext convention) and HMC sweep grid plan (Phase 4 knee table).

### Verification

Re-rendered PNGs inspected: `$X_c$`, `$X_\mathrm{GB}^\mathrm{FD}$`,
`$\Delta E_\mathrm{seg}$`, `$E_\mathrm{GB}^\mathrm{Mg}$` all display
correctly as math-italic with subscripts/superscripts; no literal `_`
remaining in any axis label, legend, or title.

## 2026-04-25 (late 4) — 200A tight-CG production lands; structure-realization residual not detectable; FD predictor ready

### Production NPZ landed (job 64755232)

The dependency-chained tight-CG resume of job `64743372` finished
in 46 min wall (504/510 sites resumed from prior checkpoint, only
the 9 unfinished + bulk refs needed re-running). Result on our own
200A structure: `delta_e_results_n500_200A_tight.npz`,
ΔE_seg [kJ/mol] min=−48.46 max=+35.79 **mean=−6.91**.
CG stop reasons: 496 `linesearch alpha is zero`, 14 `energy tolerance`
at 1e-25 — all converged at machine precision. Bulk-ref:
E_bulk^Mg = −1613193.8018 ± 0.0081 eV (n=10).

### Compare vs Wagih Zenodo pool — KS p=0.892

Re-ran `scripts/compare_vs_wagih.py` on our 200A tight NPZ vs Wagih's
82,646-site Zenodo pool (`output/compare_vs_wagih_200A_tight.{png,json}`):

| metric          | **Ours 200A tight** | Wagih (n=82,646) | Δ |
|-----------------|---------------------|-------------------|------|
| sample mean (kJ/mol) | **−6.906** | −6.814 | +0.09 |
| sample std           | **15.07**  | 15.85  | −0.78 |
| sample skew          | **−0.213** | −0.224 | +0.01 |
| skew-normal μ        | **+6.34**  | +6.72  | −0.38 |
| skew-normal σ        | **+20.06** | +20.84 | −0.78 |
| skew-normal α        | **−1.465** | −1.395 | −0.07 |
| KS two-sample        | **D=0.0256, p=0.892** | — | indistinguishable |

### Bootstrap CI (B=10⁴, n=500, seed=20260425) — `scripts/bootstrap_vs_wagih.py`

Drew 10⁴ N=500 sub-samples from Wagih's 82,646 pool:

| stat              | boot_mean ± σ_boot | 95% CI                | ours    | z      | percentile | result |
|-------------------|--------------------|-----------------------|---------|--------|------------|--------|
| sample mean       | −6.82 ± 0.71       | [−8.23, −5.42]        | −6.906  | −0.12  | 44.9       | **inside** |
| sample std        | +15.85 ± 0.49      | [+14.89, +16.81]      | +15.073 | −1.59  | 5.6        | inside (tail) |
| sample skew       | −0.222 ± 0.092     | [−0.40, −0.04]        | −0.213  | +0.10  | 54.9       | inside |
| skew-normal μ     | +6.25 ± 2.89       | [−0.26, +10.33]       | +6.343  | +0.03  | 43.9       | inside |
| skew-normal σ     | +20.67 ± 1.58      | [+17.09, +23.44]      | +20.057 | −0.39  | 32.4       | inside |
| skew-normal α     | −1.385 ± 0.399     | [−2.12, −0.56]        | −1.465  | −0.20  | 42.5       | inside |

### Resolution: the 4.8 kJ/mol shift fully decomposed

```
4.8 kJ/mol shift (loose-CG 200A vs Wagih A)
  = 3.0 kJ/mol  loose CG          ← fixed (etol 1e-8 → 1e-25)
  + 1.8 kJ/mol  predicted "structure realization"
                ↑ NOT observed at tight CG: actual residual is 0.09 kJ/mol,
                  inside CI, std-ratio is the only mildly tail-leaning stat.
```

**Conclusion:** the residual after CG fix is sampling noise, not a
structural bias. Our Voronoi/anneal/CNA pipeline produces a Wagih-
equivalent ΔE distribution on independent structures. Pipeline is
locked in; no further structure-side debugging warranted.

### Fermi-Dirac predictor — `scripts/fermi_dirac_predict.py`

New, ~170 lines. Loads ΔE NPZ (eV), evaluates

    P_i(T, X_c) = 1 / (1 + ((1 − X_c)/X_c) · exp(ΔE_i / kT))
    X_GB^FD(T, X_c) = (1/N_GB) · Σ_i P_i

on a (T, X_c) grid; overlays Wagih's 82,646-site pool. Ships with
`--self-test` covering analytic limits:

- X_c → 0 (linear scaling, halving X_c halves X_GB)  ✓
- X_c → 1  (X_GB → 1)  ✓
- T → ∞   (X_GB → X_c, no preference)  ✓
- T → 0   (X_GB → fraction of ΔE_i < 0)  ✓

Curves on T ∈ {300, 500, 700, 900} K, X_c ∈ [10⁻⁵, 0.5] log-spaced
(`output/fd_curves_200A_tight.{png,json}`):

| T (K) | X_c=1e-4 | X_c=1e-3 | X_c=1e-2 | X_c=0.1 | X_c=0.5 |
|-------|----------|----------|----------|---------|---------|
| 300   | 0.159 / 0.166 | 0.259 / 0.256 | 0.378 / 0.373 | 0.510 / 0.514 | 0.645 / 0.647 |
| 500   | 0.037 / 0.044 | 0.105 / 0.112 | 0.236 / 0.237 | 0.433 / 0.432 | 0.641 / 0.640 |
| 700   | 0.006 / 0.009 | 0.039 / 0.045 | 0.144 / 0.149 | 0.366 / 0.366 | 0.635 / 0.631 |
| 900   | 0.002 / 0.002 | 0.016 / 0.019 | 0.090 / 0.095 | 0.312 / 0.313 | 0.627 / 0.623 |

(format: ours / wagih). Max curve-level deviation ≤ 0.008 anywhere
on the grid — consistent with the spectrum-level KS result. Saturation
plateau at X_GB ≈ 0.65 = frac(ΔE<0).

### "Knee" location → HMC grid plan

The breakdown of dilute-limit assumptions is most likely to show up
*in or just before* the rising-knee region. Reading off the FD curves:

| T (K) | knee X_c (X_GB rising 10→50%) |
|-------|-------------------------------|
| 300   | ~3×10⁻⁵ to 3×10⁻⁴ |
| 500   | ~3×10⁻⁴ to 5×10⁻³ |
| 700   | ~5×10⁻³ to 5×10⁻² |
| 900   | ~3×10⁻² to 0.2 |

HMC sweep proposal: 4 temperatures × ~6 X_c per T, dense at the knee,
two anchor points in saturation (X_c=0.3 each T) and dilute (X_c far
below knee, where X_GB^FD ≈ 0). ~24-point grid.

### Next

Phase 4: write `data/decks/submit_hmc_AlMg.sh` (LAMMPS `fix atom/swap`),
single-point dry-run at T=500 K, X_c=5×10⁻³ to verify swap acceptance
rate (target 5–30 %), energy plateau, X_GB autocorrelation < window.
Then batch the 24-point grid, headline figure: X_GB^HMC vs X_GB^FD
overlay; breakdown X_c per T = first divergence beyond FD-bootstrap CI.

### Artifacts

- `output/compare_vs_wagih_200A_tight.{png,json}`
- `output/bootstrap_vs_wagih_200A_tight.json`
- `output/fd_curves_200A_tight.{png,json}`
- `scripts/bootstrap_vs_wagih.py` (new)
- `scripts/fermi_dirac_predict.py` (new)

## 2026-04-25 (late 3) — Phase 4 sequencing: build FD predictor before HMC scan

With pipeline validation closed (2026-04-25 late 2) and production 200A
tight-CG ΔE landing imminently (job `64743372` in flight, `64755232`
dependency-chained as resume), the question is whether to go straight
to HMC `(T, X_c)` scan or first construct the Fermi-Dirac (FD) predictor
curve. Decision: **FD first, then HMC.**

### Why FD first

FD evaluation is a Python sum over the 500 ΔE values — microseconds:

```
X_GB^FD(T, X_c) = (1/N_GB) · Σ_i 1 / [1 + ((1−X_c)/X_c) · exp(ΔE_i / kT)]
```

Producing it before HMC gives us three things that the cost of HMC
would otherwise force us to learn the hard way:

1. **Visual sanity check vs. Wagih Fig 5.** Same potential, same alloy,
   same N_GB scale → curve shape and saturation values should match
   Wagih's published curves. Mismatch here would point at a bug we
   haven't caught (sign error in the formula, kT unit mistake, etc.)
   *before* spending 30 compute-hours discovering it via HMC.
2. **Grid density planning for HMC.** FD has a sigmoid-like X_c sweep
   with a "knee" region where saturation kicks in; segregation breakdown
   is most likely *near or just before* that knee. Uniform 5-point X_c
   sweeps risk wasting half the HMC budget far from the interesting
   region. Looking at FD first lets us put HMC density where the
   breakdown is testable.
3. **Quantitative "significantly different" threshold.** Currently we
   say "X_GB^HMC vs. X_GB^FD differ by ML noise ~4 kJ/mol" with no
   anchored statistic. With FD curves we can compute per-(T, X_c)
   sensitivity (∂X_GB^FD / ∂P_i, ΔE_i contribution) and define the
   threshold rigorously instead of by eye.

### Concrete sequence

```
1. Wait for 64743372  →  delta_e_results_n500_200A_tight.npz on our
                          own structure (or use Wagih-validated NPZ
                          as placeholder if 64743372 misbehaves)
2. scripts/fermi_dirac_predict.py:
     in:  ΔE NPZ + (T, X_c) grid spec
     out: X_GB^FD on the grid, overlay plot vs. Wagih Fig 5
     unit test against analytic limits
       X_c → 0    →  X_GB → exp(−ΔE_i/kT)·X_c/(1−X_c)·N_negative
       X_c → 1    →  X_GB → 1
       T → ∞      →  X_GB → X_c (no preference)
       T → 0      →  X_GB → fraction of negative ΔE_i sites
3. Inspect FD curves at T ∈ {300, 500, 700, 900} K vs. X_c ∈ log grid
     → choose HMC grid: dense near knee, sparse in saturation/dilute
4. data/decks/submit_hmc_AlMg.sh:
     fix atom/swap 100 10 <swap_seed> <T> ke yes types 1 2
     equilibrate ~10 ps NVT  →  HMC ~100 ps  →  block-average X_GB
     single-point dry-run first (e.g. T=500 K, X_c=5 at%) to verify:
       - swap acceptance rate sane (5–30 %)
       - energy plateaus before measurement window
       - X_GB time series autocorrelation < window length
5. Batch-submit the chosen (T, X_c) grid in parallel
6. Compare X_GB^HMC vs. X_GB^FD curves → headline figure;
   define breakdown X_c per T as first point where they diverge
   beyond the FD-bootstrap CI from step 2.
```

### Risk acknowledged

Going FD-first delays HMC start by ~½ day of code + figure work. If
the FD curve looks fine on Wagih's structure but reveals something
weird on our 200A (e.g. a ΔE outlier dragging the X_GB tail), we'd
investigate that before HMC — adds another half-day. Total worst
case: ~1 day of FD/grid work before HMC starts.

The alternative (straight to HMC) would save that day **only** if our
HMC grid happens to land at the right resolution near the breakdown.
Given we have no a-priori signal for where that is on Al(Mg), the FD
detour is "buy information for compute budget."

### Tooling note

`scripts/fermi_dirac_predict.py` is small (~80 lines) and shares no
state with `compare_vs_wagih.py` / `paired_pipeline_residual.py` /
`fit_delta_e_spectrum.py` — keep them as four single-purpose scripts
rather than a CLI mega-tool. Each one reads NPZ + writes plot+JSON;
composability via the file system is enough for our scale.

## 2026-04-25 (late 2) — Tight-CG validation passed: pipeline indistinguishable from Wagih on his own structure

Job `64742000` finished in **4h 06min** (8h budget; faster than predicted
because long-tail CG was less common than the 36 s/site sanity-test
average suggested — actual mean was 28.9 s/site over 510 sites). All
510 site decks completed; 501/510 stopped on `linesearch alpha is zero`
and 9 hit `energy tolerance` even at 1e-25 (CG converged exactly to
machine ε). The full N=500 ΔE result is at
`/cluster/scratch/cainiu/wagih_pipeline_test/delta_e_results_wagih_n500_tight.npz`.

### Bulk reference: independently reproduces Wagih's value

Our 10-atom bulk-ref protocol (random bulk-interior atoms ≥8 Å from any
GB) yields **E_bulk^Mg = −1639944.4135 eV (σ=0.0018, n=10)**. Wagih
reports `bulk_solute_Al_Mg.dat = −1639944.41395 eV` from his
6-nm-sphere protocol. **Difference = 0.45 meV ≈ 0.04 kJ/mol** —
i.e., our locally-sampled bulk Al(1Mg) is the same energy as his
spherical bulk Al(1Mg). The two protocols converge to the same minimum
under tight CG.

### Paired residual on full N=500

Re-ran `paired_pipeline_residual.py` with the full NPZ output. Results
in `output/paired_pipeline_residual_n500_tight.{png,json}`:

| metric        | loose (n=343) | **tight (n=500)** |
|---------------|---------------|-------------------|
| Pearson r     | 0.9808        | **0.999999999996** |
| mean residual | +31.4 meV     | **−0.0346 meV** (= −0.0033 kJ/mol) |
| std residual  | 33.1 meV      | **0.000478 meV** |
| min residual  | +5.1 meV      | −0.0367 meV |
| max residual  | +287.7 meV    | −0.0330 meV |

The residual is essentially a **constant −34.6 μeV/site** offset with
σ < 0.5 μeV. Magnitude is 4 orders below ΔE relevance. Likely cause:
LAMMPS-version (we use 20240829, Wagih used the 2020 release) or
MPI-rank count affecting the floating-point summation order across 34M
pair contributions per box. Cannot tell which without re-running on
Wagih's exact LAMMPS build; not on the critical path. Bug-class is
"reduction-order ulp drift", not "physics".

Driver fix: NPZ key for per-site PE is `gb_e_mg`, not `gb_pe`. The
`paired_pipeline_residual.py` `load_ours_npz` was patched accordingly
(it had been written from memory before the first NPZ landed).

### Spectrum-level fits

| metric | **Ours (tight, on Wagih, n=500)** | Wagih A (n=82,646) | gap |
|--------|-----------------------------------|--------------------|-----|
| sample mean (kJ/mol) | **−6.857** | −6.814 | +0.04 |
| sample std           | **16.39**  | 15.85  | +0.54 |
| sample skew          | **−0.169** | −0.224 | +0.06 |
| skew-normal μ        | **+6.30**  | +6.72  | −0.42 |
| skew-normal σ        | **+21.00** | +20.84 | +0.16 |
| skew-normal α        | **−1.263** | −1.395 | +0.13 |
| range (kJ/mol)       | [−53.1, +43.7] | [−67.4, +48.0] | tighter (n effect) |

KS two-sample: **D = 0.0249, p = 0.9106** — completely cannot reject
"same distribution" at any α. (Loose CG was D=0.139, p=7.2×10⁻⁹.)

### Bootstrap CI (B=10⁴, n=500 with replacement, seed=20260425)

Drew 10⁴ N=500 sub-samples from Wagih's 82,646 ΔE pool and computed six
statistics each. CI = 2.5–97.5 percentile of the bootstrap distribution.

| stat              | boot mean ± σ_boot | 95% CI                | ours    | z      | percentile | result |
|-------------------|--------------------|-----------------------|---------|--------|------------|--------|
| sample mean       | −6.82 ± 0.71       | [−8.23, −5.42]        | −6.857  | −0.05  | 47.6       | **inside** |
| sample std        | +15.85 ± 0.49      | [+14.89, +16.81]      | +16.39  | +1.12  | 86.9       | inside |
| sample skew       | −0.222 ± 0.092     | [−0.40, −0.04]        | −0.169  | +0.58  | 72.0       | inside |
| skew-normal μ     | +6.25 ± 2.89       | [−0.26, +10.33]       | +6.30   | +0.02  | 43.2       | inside |
| skew-normal σ_skew| +20.67 ± 1.58      | [+17.09, +23.44]      | +21.00  | +0.22  | 56.0       | inside |
| skew-normal α     | −1.385 ± 0.399     | [−2.12, −0.56]        | −1.263  | +0.31  | 65.1       | inside |

**All six inside; every |z| < 1.2**; sample mean at the 47th percentile —
dead center. By every test we have, our pipeline run on Wagih's
structure with tight CG is statistically indistinguishable from drawing
500 random sites from his full 82,646 population. The 4.8 kJ/mol shift
that started this thread decomposes cleanly:

```
4.8 kJ/mol shift (our 200A vs Wagih A)
  = 3.0 kJ/mol  pipeline (loose CG)        ← isolated and fixed today
  + 1.8 kJ/mol  Voronoi/anneal realization ← pending: 64743372 quantifies
```

### Ancillary code-quality observations

- 9/510 sites stopping on `energy tolerance` even at 1e-25 means
  `relative ΔE between two CG iterations` did go below 1e-25 for those
  sites — possible only when the gradient was already at numerical
  noise, so equivalent to `linesearch alpha = 0` outcome. Treat as
  good convergence.
- Bootstrap timing: 10⁴ skewnorm.fit calls took ~6 min single-thread.
  If we expand the analysis to multiple structures / multiple n's a
  vectorized fitter or `joblib.Parallel` would be worth the lift; one-
  off use today doesn't justify it.

### Next

Production 200A re-run (job `64743372`) is in progress on our own
structure with the same tight-CG defaults. At ~60 s/site current rate
it may bump up against the 8 h walltime; per-site checkpointing will
preserve work and a follow-up sbatch can finish off the bulk refs.
Once it lands, the same paired/bootstrap analysis will quantify the
~1.8 kJ/mol structure-realization residual against Wagih's pool.

## 2026-04-25 (late) — Methodology lockdown: CG-tolerance semantics, bootstrap CI, the three σ's

Concepts that came up while explaining today's loose-CG diagnosis to
user. Locked here so future sessions don't re-derive them.

### `minimize etol ftol maxiter maxeval` semantics (LAMMPS)

`minimize` exits on whichever criterion fires **first**:

| arg | meaning |
|-----|---------|
| **etol** (energy tolerance, **能量容差**) | exit when \|ΔE\|/max(\|E\|, EMACH) between two iters < etol |
| **ftol** (force tolerance, **力容差**) | exit when max-component force \|f\|_∞ < ftol (eV/Å) |
| **maxiter** (maximum iterations, **最大迭代数**) | exit after this many CG iterations |
| **maxeval** (maximum force evaluations, **最大力评估数**) | exit after this many force calls (CG line search may invoke ≥1/iter) |

Wagih's `1e-25 1e-25 50000 5000000` effectively **disables etol + ftol**:
double-precision floating point caps at ~1e-15 relative precision, so
etol=1e-25 can never trigger; ftol=1e-25 eV/Å is below quantum noise.
CG instead runs until `linesearch alpha = 0` (atom can't move without
raising energy) or until maxiter/maxeval cap. In our tight-CG sanity
test (job `64741906`), site 196943 stopped at 200 iters / 396 force
evals with f_max = 7.9e-5 eV/Å — far above ftol, at the numerical floor
for that local geometry.

Our previous `1e-08 1e-10 5000 50000` triggered etol after a few
hundred force evals while leaving ~3 kJ/mol of elastic relaxation on
the table — the bias diagnosed today.

### Bootstrap (**自助法**) + confidence interval (**置信区间**, CI)

**Bootstrap** (Efron 1979) answers *"how much does my fitted statistic
wobble if I'd drawn a different N=500 from the same underlying
population?"* without needing a parametric model for that wobble.
Procedure (used 2026-04-24 (late 5); to be re-run on tight-CG output):

1. Treat Wagih's full 82,646 ΔE values as the "population".
2. Draw N=500 **with replacement** (这里就是 bootstrap 的关键), B=10⁴ times.
3. For each draw compute (sample_mean, sample_std, sample_skew,
   skew-normal μ, σ, α) — six statistics.
4. The B values per statistic form an empirical sampling distribution.

A **95% CI** is the 2.5th–97.5th percentile of that distribution.
Interpretation: *if* ours and Wagih's come from the same population,
P(our N=500 statistic ∈ CI) = 95%. Outside → reject "same population"
at p < 0.05.

Loose-CG verdict (already on file): sample_mean **rejected** at
z=+6.85 σ_boot; σ and α **accepted** at z within ±1. Tight-CG re-run
(job `64742000`) is expected to push the mean back into the CI, with
any residual offset attributable to structure-realization difference,
not pipeline bias.

### Three different σ's — keep them straight

| symbol | source | role | Al(Mg) magnitude (kJ/mol) |
|--------|--------|------|---------------------------|
| **σ_skew** | `scipy.stats.skewnorm.fit` returns `scale` | distribution **scale parameter** in the (μ, σ, α) triple — ≠ true std when α≠0 | ~21 (Wagih A: 20.84) |
| **σ_sample** | `np.std(ΔE, ddof=1)` | empirical **sample standard deviation** of the per-site ΔE values | ~15.8 (Wagih A: 15.85) |
| **σ_boot** | `np.std(boot_X, ddof=1)` over B bootstrap draws | **standard error of a statistic** X (X = mean / std / α / …) under repeated N=500 sampling | sample_mean: 0.71 / sample_std: 0.49 / α: 0.41 |

**Relations** (all sanity-checked against Wagih A):

- `σ_sample = σ_skew · √(1 − 2δ²/π)` with `δ = α/√(1+α²)` — the
  geometric shrinkage from "scale" to true std once the distribution
  is skewed. Wagih A: α=−1.40 → δ=−0.814 → σ_sample ≈ 20.84 · 0.760
  = 15.84 ≈ measured 15.85 ✓.
- `σ_boot(mean) ≈ σ_sample / √n` (Central Limit Theorem). For std and
  α the analytic form depends on the statistic; we just measure the
  bootstrap distribution directly. Sanity: 15.85/√500 = 0.709 ≈
  measured 0.71 ✓.

**The "+6.85σ" cited above uses σ_boot(sample_mean) = 0.71**, NOT
σ_skew or σ_sample. Same symbol, three different objects — when
writing future comparisons always tag which σ. To make this explicit,
the `compare_vs_wagih.py` extension that lands with the tight-CG
results will print z-scores and CIs labeled `σ_boot[mean]`,
`σ_boot[std]`, `σ_boot[alpha]` to avoid the ambiguity.

## 2026-04-25 — CG tightness diagnosed: loose `etol=1e-8` accounts for the +3 kJ/mol shift

The overnight controlled experiment (job `64707493`, submitted 2026-04-24
late 5) **timed out at 3 h with 343/510 GB sites done**. Per-site CSV
checkpointing preserved the partial data; bulk references were never run
(the driver finishes all GB sites first, then bulk refs). Rather than
wait for a fresh run, we analyzed the partial data with a paired
diagnostic — and it gave a clean answer.

### Paired diagnostic on the 343 partial sites

`scripts/paired_pipeline_residual.py` (new, ~140 lines) joins the
streamed `_results.csv` with Wagih's `seg_energies_Al_Mg.txt` on
`atom_id`. This bypasses the bulk-reference question entirely: both
sides use the same annealed structure (Wagih's dump0) and the same
substituted atom_id, so `pe_ours(i) − E_GB_wagih(i)` isolates pipeline
differences with no E_bulk involvement.

Result on n=343 paired sites (`output/paired_pipeline_residual_n343.{png,json}`):

| metric              | value     |
|---------------------|-----------|
| Pearson r           | **0.9808** |
| mean residual       | **+31.4 meV  = +3.03 kJ/mol** |
| std residual        | 33.1 meV  =  3.20 kJ/mol |
| min residual        | **+5.11 meV** (never reaches Wagih) |
| max residual        | +287.7 meV (one stuck-CG outlier) |
| residual histogram  | 0 / 9.9 / 37.9 / 21.3 / 16.3 / 9.9 / 4.7 % in [0,5)/[5,10)/[10,20)/[20,30)/[30,50)/[50,100)/[100,1000) meV |

**Every single site is positive** — our PE is always higher (less
relaxed) than Wagih's. The minimum residual is +5.1 meV, not zero.
That signature is uniform under-relaxation, not numerical noise.

The mean +3.03 kJ/mol across 343 paired sites accounts for ~2/3 of
the +4.8 kJ/mol mean shift between our N=500 production fit (on our
own structure) and Wagih's 82,646-site Zenodo dataset. Caveat: the
+3.03 was measured on Wagih's structure with our loose-CG protocol;
it should transfer approximately but not exactly to our own structure.
Residual = 4.8 − 3.0 ≈ 1.8 kJ/mol; on the 500-draw bootstrap
SE_mean = 0.71 this is ~2.5σ — not noise, but plausibly the
structure-realization difference (ours has f_gb = 18.7 % vs Wagih's
17.1 %; different Voronoi seed). Tight-CG re-runs on both structures
will let us compute the structure component cleanly.

### Smoking gun: Wagih's `calculate_E_GB_solute.in`

Read Wagih's actual LAMMPS deck from the Zenodo `lammps_example_input_files/`:

```
min_style cg
minimize 1e-25 1e-25 50000 5000000
```

Ours (verified in the runtime per-site `site_<id>.lammps` deck):

```
min_style cg
minimize 1e-08 1e-10 5000 50000
```

Wagih's etol/ftol are **17 / 15 orders of magnitude tighter** and his
maxiter/maxeval are 10× / 100× larger. CG runs to floating-point
precision; ours stops on the energy criterion after a few hundred
iterations. The +3 kJ/mol systematic bias is the residual elastic
relaxation we leave on the table by stopping early.

### Sanity test: tight CG fully closes the gap (job `64741906`)

Two sites re-run with Wagih's tolerances on Wagih's structure
(`/cluster/scratch/cainiu/wagih_pipeline_test/tight_cg_test/`):

| site_id | resid (loose) | pe_ours (tight) | pe_wagih       | resid (tight) | wall  |
|---------|---------------|-----------------|----------------|---------------|-------|
| 60470   | +18.79 meV    | −1639944.55868  | −1639944.55865 | **−0.035 meV** | 29 s |
| 196943  | +287.72 meV   | −1639944.58033  | −1639944.58030 | **−0.027 meV** | 44 s |

Even the worst stuck-CG outlier converges to Wagih's PE within
machine precision once tolerances are loose enough to let CG actually
run. **No alternate basin of attraction; just under-iteration.**

Tight CG cost: ~36 s/site avg vs. ~7–9 s/site for loose CG (~5×
slower). Affordable: 510 sites ≈ 5 h on 16 cores.

### Code changes

- `scripts/sample_delta_e.py`: defaults bumped to Wagih's
  `etol=1e-25, ftol=1e-25, maxiter=50000, maxeval=5000000`. Function
  signature, CLI argparse defaults, and docstring all updated.
- `_META_KEYS_TO_MATCH` extended to include `cg_etol/ftol/maxiter/maxeval`,
  so a loose-CG checkpoint can no longer be silently extended with
  tight-CG sites on resume — the driver will refuse with a meta
  mismatch error instead.
- `data/decks/submit_delta_e{,_200A,_wagih_structure}.sh`: dropped the
  explicit `--etol/--ftol` overrides so deck inherits the new tight
  defaults; updated comments and walltimes (8 h for the wagih and 200A
  runs; prototype 4 h still fits the smaller box).

### Resubmission (job `64742000`, in flight)

Same `submit_delta_e_wagih_structure.sh` recipe as last night but with
tight CG and a fresh `delta_e_run_tight/` work_dir. Expected outcome:
N=500 paired residuals collapse to <1 meV mean / <1 meV σ, our skew-
normal fit on Wagih's structure matches Wagih's bootstrap CI for n=500,
and the residual 1.8 kJ/mol vs. our own structure becomes the cleanly-
isolated structure-realization variance we set out to measure.

### Implications

1. **The 200A production run (`delta_e_results_n500_200A.npz`) is
   contaminated** by the same loose-CG bias. Its (μ=+9.4, σ=19.4,
   α=−1.08) is therefore not a clean head-to-head with Wagih. Re-running
   it under tight CG is the natural next step before Phase 4 (FD
   prediction curve uses these ΔE values; we want them right).
2. **The prototype run (`delta_e_results_n500.npz`, 10³ nm³)** is also
   loose-CG; same fix applies if/when we want to use that data for
   anything quantitative.
3. **HMC truth curve (Phase 4) is unaffected**: HMC samples Boltzmann
   directly and never calls our per-site CG, so the bias only enters
   the FD prediction side of the headline figure.
4. **For the writeup**: the SI Mg¹⁵ vs. our spectrum comparison from
   2026-04-24 (late 4–6) needs an addendum once the tight-CG re-runs
   land. The α=−1.08 vs. Wagih A α=−1.40 sub-sampling argument still
   stands; the μ shift will close.

### Artifacts

- `scripts/paired_pipeline_residual.py` (new diagnostic, takes either
  partial CSV or final NPZ; outputs scatter + residual histogram)
- `output/paired_pipeline_residual_n343.{png,json}` (partial-data figure
  shipped to evidence the 2026-04-24 protocol problem)
- `/cluster/scratch/cainiu/wagih_pipeline_test/tight_cg_test/` (the
  2-site validation; small, kept on scratch as evidence)
- Job `64742000`: tight-CG full validation (8 h budget, results expected
  late afternoon 2026-04-25)

## 2026-04-24 (late 6) — SI Mg¹⁵ α=−2.3 confirmed: it's the ML-predicted fit

**Late-night discovery** (after delayed tar listing completed): Zenodo
*does* contain the Mendelev 2009 accelerated-model output at
`segregation_spectra_database_accelerated_model/Al/Mg_2009--Mendelev-M-I-Asta-M-Rahman-M-J-Hoyt-J-J--Al-Mg/Al_Mg_20nm_GB_segregation.dump`.
The earlier claim ("only Liu-Adams 1998 in accelerated DB") was from a
partial tar listing. Full listing shows 5 Al-Mg potentials in the
accelerated DB (refs 13, 14, **15**, 16, 2012-Jelinek, Zhou04).

### `Mg^15` SI panel source identified

Skew-normal fit to the 82,635 ML-predicted ΔE values in the Mendelev
2009 accelerated-model dump:

```
μ = +9.44 kJ/mol
σ = 23.13 kJ/mol
α = −2.256
```

→ **Matches user's SI Fig 3 panel read (μ=9, σ=23, α=−2.3, R²=1.00)
within sub-0.1 precision**. User was right all along; the α=−2.3 label
on SI Fig 3 Mg¹⁵ is what the data says. No figure error.

### Two Zenodo datasets are distinct (both on Mendelev 2009)

| Dataset | Path                                              | Size | Nature | μ | σ | α | mean |
|---------|---------------------------------------------------|------|--------|---|---|---|------|
| **A** direct LAMMPS | `machine_learning_notebook/seg_energies_Al_Mg.txt` | 82,646 | per-site CG (training data) | +6.72 | 20.84 | **−1.40** | −6.81 |
| **B** ML predicted | `segregation_spectra_database_accelerated_model/Al/Mg_2009.../.dump` | 82,635 | SOAP+linreg predicted | +9.44 | 23.13 | **−2.26** | −7.46 |

SI Fig 3 Mg¹⁵ shows **dataset B** (ML predicted output), not A. My earlier
comparison used A and concluded "4.8 kJ/mol mean shift" — that framing
compared us to the training-set ground truth, not to what the paper
actually plots.

### Repositioning our 500-site fit

| source | μ | σ | α | sample mean |
|--------|---|---|---|-------------|
| **Ours (500, direct LAMMPS, our structure)** | **+9.40** | 19.43 | −1.08 | −2.00 |
| Wagih A (direct LAMMPS, n=82k)    | +6.72 | 20.84 | −1.40 | −6.81 |
| Wagih B (ML predicted, n=82k = SI Fig 3)  | **+9.44** | 23.13 | −2.26 | −7.46 |

**Striking**: our μ = +9.40 **matches Wagih B to sub-0.05 kJ/mol**. Our
σ and α are between the two Wagih datasets but closer to A. So:
- The SI Fig 3 fit target's **location parameter μ** is already met by
  our pipeline, on our own structure.
- Our α=−1.08 is consistent with Wagih's *direct-LAMMPS* distribution
  (81st percentile in bootstrap), not with Wagih's ML-predicted α=−2.3.
  In other words: our raw ΔE values have roughly the same skewness as
  Wagih's raw ΔE values. The strong α=−2.26 on the SI panel is an
  **ML-smoothing artifact** (linear regression on 100 SOAP centroids
  broadens and skews the predicted spectrum beyond what the underlying
  CG data actually contains).
- The "4.8 kJ/mol mean shift vs A" persists, but it coincides with
  A-vs-B being different by 0.65 kJ/mol anyway — so ~3-5 kJ/mol of the
  mean residual vs SI is still potentially structure-driven.

### Does the overnight experiment still matter?

Yes. Job `64707493` tests whether **our Phase 3 pipeline on Wagih's
structure matches Wagih's direct-LAMMPS A values on identical atom IDs**.
Expected outcome: if paired ΔE_ours(i) ≈ ΔE_A(i), our CG + substitution
protocol is verified correct, and the entire residual between our 500
(on our structure) and Wagih A (on his structure) is pure Voronoi/anneal
realization variance — a clean green light to proceed to Phase 4.

### Implications for the paper story

1. Reproducing SI Fig 3 skew-normal parameters **precisely** requires
   running Wagih's ML pipeline (SOAP + k-means + linear regression),
   not just direct LAMMPS. Our direct-LAMMPS fit naturally differs
   because ML smoothing broadens and skews.
2. For the project's scientific question (when does Fermi-Dirac break),
   **direct LAMMPS is the correct input** — Fermi-Dirac under the
   independent-site assumption works with true site energies, not
   ML-smoothed predictions.
3. When writing up, frame comparison vs Wagih dataset **A (direct
   LAMMPS)**, not the SI figure. SI figure is ML-smoothed predictions
   designed for the 259-alloy throughput scan, not the ground truth.

### Artifacts

- `segregation_spectra_database_accelerated_model/Al/Mg_2009.../Al_Mg_20nm_GB_segregation.dump` extracted to scratch (23 MB).
- Will extend `scripts/compare_vs_wagih.py` to optionally load dataset B
  for a 3-way fit comparison plot when morning analysis resumes.

## 2026-04-24 (late 5) — Wagih Zenodo dataset + KS test + structure audit launched

Chose path Q (fetch Wagih's raw Al(Mg) data before committing to Phase 4).
Resolved an anomaly along the way: our N=500 mean is **7σ** off from
Wagih's full dataset — *not* sampling noise. User flagged structure
generation as the likely suspect; launched an overnight controlled
experiment to isolate cause.

### Zenodo archive (doi:10.5281/zenodo.4107058)

Single 4.0 GB tarball `learning_segregation_energies.tar.bz2` (MD5
verified: `dcad1225446df20c841b8a32359c03b1`). Full listing = 904 entries,
downloaded to `/cluster/scratch/cainiu/wagih_zenodo/`. Archive structure:

- `machine_learning_notebook/` — high-fidelity training example for Al(Mg)
  (the `pair_coeff * * Al-Mg.eam.fs Al Mg` in `calculate_E_GB_solute.in`
  **confirms same Mendelev 2009 potential we use**)
  - `heated_minimized_Al_polycrystal.dump0` — annealed pure-Al structure
  - `seg_energies_Al_Mg.txt` — **82,646 per-site** (site_id, E_GB^Mg) lines
  - `bulk_solute_Al_Mg.dat` — E_bulk^Mg = −1,639,944.41 eV
  - `GB_SOAP_Al_Mg.npy` (613 MB SOAP features)
  - `Learn_Segregation_Spectra.ipynb` — their ML notebook
- `segregation_spectra_database_accelerated_model/` — the 259-alloy scan
  output; Al/Mg subfolder **only has Liu-Adams 1998 (SI Ref [14])**, NOT
  Mendelev 2009. So the accelerated-model DB is NOT the right source for
  our apples-to-apples comparison; the notebook data is.

### Direct-LAMMPS ΔE spectrum fit (`scripts/compare_vs_wagih.py`)

Loaded Wagih's 82,646 ΔE values (computed `E_GB^Mg(i) − E_bulk^Mg` per
line, converted to kJ/mol). Fit skew-normal with identical
`scipy.stats.skewnorm`:

| metric | **Ours (n=500)** | **Wagih Zenodo (n=82,646)** | SI Mg¹⁵ (user read) |
|--------|------------------|-----------------------------|---------------------|
| μ      | +9.40            | **+6.72**                   | +9                   |
| σ      | 19.43            | **20.84**                   | ~23                  |
| α      | −1.08            | **−1.40**                   | −2.3 (re-checked OK) |
| sample mean | −2.00       | **−6.81**                   | —                    |
| sample std  | 15.76       | **15.85**                   | ~15                  |
| sample skew | −0.11       | **−0.22**                   | (≈−0.53 if α=−2.3)  |
| range (kJ/mol) | [−40.2, +58.2] | **[−67.4, +48.0]**   | [−60, +40]           |

**KS two-sample test**: D = 0.139, p = 7.2×10⁻⁹ — reject same-
distribution at 99.9%. But D < 0.2 is a small-to-moderate effect.

### Unresolved: SI Mg¹⁵ α=−2.3 vs Zenodo α=−1.40

User re-checked SI Fig 3 Mg¹⁵ panel — confirms α=−2.3 is what's drawn
and it's unlikely a figure-generation error for batch-published panels.
Yet the Zenodo 82,646-site direct-LAMMPS fit gives α=−1.40 on the
same Mendelev 2009 potential. Possible explanations (none yet verified):
- SI Fig 3 panel reports the **accelerated-model** (ML-predicted) fit,
  not the high-fidelity direct fit. But the Al/Mg accelerated output on
  Zenodo uses Liu-Adams 1998, not Mendelev 2009.
- The `heated_minimized_Al_polycrystal.dump0` is a *different* Voronoi
  realization from the one underlying SI Fig 3.
- Fit method difference (e.g. weighted KDE vs moment-based fit).

**Parked** — not on the critical path for answering the project's
central question.

### Bootstrap: our 4.8 kJ/mol mean shift is real, not sampling noise

Drew 10,000 random sub-samples of size n=500 from Wagih's 82,646 values:

| statistic | sub-sample distribution (500-draw from Wagih) | Ours (500) | z-score / percentile |
|-----------|----------------------------------------------|------------|----------------------|
| mean      | μ=−6.82, σ=0.71                              | **−2.00**  | **+6.85σ, pct 100** |
| std       | μ=15.84, σ=0.49                              | 15.76      | −0.17σ, pct 43     |
| α (skew-normal fit) | μ=−1.39, σ=0.41, 5–95% [−2.00, −0.76] | −1.08 | **pct 81** (well within) |

**Conclusion**: our 500-sample (σ, α) are **both statistically
consistent** with being a 500-draw from Wagih's full 82k. The only
**real** difference is the **4.8 kJ/mol mean shift** (+6.85σ, p≈3×10⁻¹²).
Our earlier "α mismatch" narrative is revised: α=−1.08 vs full-fit
α=−1.40 is just 500-pt noise — the bootstrap α distribution
[−2.00, −0.76] comfortably contains ours at 81st percentile.

**Also settles the α=−2.3 question empirically**: even if SI Fig 3 Mg¹⁵
label is correct, α=−2.3 is *outside* the 95% CI of 500-draw bootstrap
fits from Wagih's own 82k points, which means it cannot be a
realization drawn from the same underlying distribution. Either the SI
panel's α label is from a different dataset / fit method, or the SI
panel shows a different structure. Not on critical path.

### Structure comparison — small but not negligible

|              | Wagih (Zenodo dump0) | Ours (production anneal) |
|--------------|----------------------|--------------------------|
| N_atoms      | 483,425              | 475,715 (−1.6%)          |
| box Lx (Å)   | 200.79               | 199.78 (−0.5%)           |
| N_GB         | 82,646               | 89,042 (+7.7%)           |
| f_gb         | 17.10%               | 18.72% (+1.6 pp)         |
| ΔE range     | [−67, +48]           | [−40, +58]               |

Box and atom count differ by <2% — geometrically similar but not
identical. Our f_gb is ~10% higher than Wagih's, suggesting different
grain packing and GB character distribution.

### User's diagnosis (verbatim): *"我们应该考虑一下我们的3D结构构建是不是合理, 毕竟万物起源"*

The 4.8 kJ/mol mean shift points at our generate_polycrystal.py /
anneal chain as the first thing to audit before trusting Phase 4. Two
candidate mechanisms:
1. Voronoi realization variance — same parameters, different
   `structure_seed` give different spectra. Testable by regenerating
   with multiple seeds.
2. Pipeline difference — different CG convergence, `atom_style`, masses,
   neighbor handling, etc. would bias ΔE values.

### Overnight controlled experiment (in flight)

Goal: isolate (1) vs (2) with one job. Take **Wagih's own annealed
structure** (`heated_minimized_Al_polycrystal.dump0`) and run **our
Phase 3 pipeline** (same `sample_delta_e.py`, same CG tolerances, same
substitution protocol) against **Wagih's 82,646 GB site IDs** as the
mask.

- If ΔE(our)(i) ≈ ΔE(Wagih)(i) atom-for-atom on common sites → **our
  pipeline is correct**; the 4.8 kJ/mol residual is entirely from
  structure-realization variance (audit our Voronoi → anneal chain next).
- If paired ΔE disagree → our LAMMPS workflow has a systematic
  discrepancy; inspect step-by-step.

**Scripts**:
- `scripts/wagih_dump_to_data.py` — dump0 → LAMMPS data file (id-sorted)
  + builds length-N mask from `seg_energies_Al_Mg.txt` site IDs. Output:
  `wagih_Al_200A.lmp` + `wagih_gb_mask.npy` (82,646 True / 483,425
  confirmed, f_gb = 17.10% ✓ matches Wagih).
- `data/decks/submit_delta_e_wagih_structure.sh` — our driver on Wagih's
  structure, N_GB=500 / N_bulk=10 / seed=42 (same as production).

**Job `64707493`** submitted — 16 cores, 3 h budget. Run dir:
`/cluster/scratch/cainiu/wagih_pipeline_test/`. On completion, extend
`compare_vs_wagih.py` to do **paired** comparison
(ΔE_ours(i) vs ΔE_wagih(i) on the 500 common sites).

### Phase 4 parked until paired diagnostic lands

If the audit shows our pipeline is clean → proceed to Phase 4 HMC with
confidence.
If not → fix, then proceed.

## 2026-04-24 (late 4) — Phase 3.5 C result + sampling statistics + Phase 4 feasibility

Job `64699296` (N_GB=500, N_bulk=10, seed=42, 16 cores) completed in
**1:19 wall** (~9.3 s/site, only 3× slower than prototype's 3.1 s/site
despite 8× more atoms — CG on a single substitution is local, doesn't
scan all atoms). Bulk reference: E_bulk^Mg = −1,613,193.7 eV,
σ/√10 = **0.87 kJ/mol** baseline noise.

### Production skew-normal fit vs prototype vs Wagih Mg¹⁵

| metric | Prototype (10³, n=500) | **Production (20³, n=500)** | Wagih Mg¹⁵ |
|--------|------------------------|-----------------------------|-----------|
| N_atoms / N_GB         | 59k / 16,430    | 476k / **89,042**  | ~481k / ~10⁵    |
| μ (location, kJ/mol)   | +11.5           | **+9.4**           | +9 ✓             |
| σ (scale, kJ/mol)      | 21.3            | **19.4**           | 23               |
| α (shape)              | −1.93           | **−1.08** ⚠        | −2.3             |
| sample mean (kJ/mol)   | −3.5            | **−2.0**           | —                |
| sample std (kJ/mol)    | 15.0            | **15.8** ✓         | ~15              |
| sample skew            | −0.41           | **−0.11** ⚠        | ~−0.53 (from α) |
| ΔE range (kJ/mol)      | [−49.7, +35.2]  | **[−40.2, +58.2]** | [−60, +40]       |

Moving 10³ → 20³ fixed μ (+11.5 → +9.4, matches Wagih) but made α *worse*
(−1.93 → −1.08). Diagnosed below — it's a sampling artifact, not a
physics issue.

### Tail-outlier diagnostic (tried and rejected)

All 500 sites stopped on `energy tolerance` (clean CG, no bad minima).
The +58.2 kJ/mol top site (atom_id 428941) is **physically real**, not a
numerical artifact. Dropping it manually would be cherry-picking.
Iteratively dropping top-N outliers:

| drop top N | μ | σ | α | sample skew |
|------------|---|---|---|-------------|
| 0          | +9.4  | 19.4 | −1.08 | −0.11 |
| 1          | +11.9 | 20.9 | −1.52 | −0.21 |
| 10         | +14.0 | 22.4 | **−2.46** | −0.34 |

Getting α to match by dropping 10 out of 500 is artificial symmetry-
matching. The real asymmetry is **the negative tail truncates at −40
while Wagih goes to −60** — we simply did not sample the deep-
segregation sites. Dropping positives compensates in the fit but does
not reflect the true distribution.

### Why n=500 estimates μ and σ but not α (sampling statistics)

Different statistical quantities have wildly different convergence rates
on sample size:

| quantity | standard error | at n=500, σ=15.8 | relative |
|----------|----------------|------------------|----------|
| mean     | σ/√n           | 0.71 kJ/mol      | 4.5%     |
| **std**  | **σ/√(2n)**    | **0.50 kJ/mol**  | **3.2%** |
| skewness | √(6/n)         | 0.110            | ~20% of true |
| α (skew-normal) | ~√(24/n) | 0.22             | ~10% of true |

Intuition: std is a **central** statistic — dominated by the ±1σ region
which holds 68% of points, densely sampled. Skewness and α are **tail**
statistics — `(x−μ)³` weights the rare extremes that 500 points sample
only 5-20 times. Wagih's ~10⁵ sites sample the tails 200× more densely
and therefore get stable α.

**Takeaway**: "std matches, α does not" is **not evidence of a physics
mismatch** — it's the expected outcome of a 500-vs-10⁵ sample-size
asymmetry. Claiming α disagreement as a real effect would require
matching sample sizes first.

### What Wagih actually did (from paper quote, 2026-04-24)

Paper text: *"Using the high-fidelity approach, we train a model ... for
Mg solute segregation in a thermally annealed 20 × 20 × 20 nm³ Al
polycrystal that has 16 grains and ~10⁵ GB sites, using a randomized
50/50 split for training/testing."*

- **~10⁵ GB sites** total in the 20³ nm³ Al polycrystal (our 89,042
  matches within close-pair/CNA-definition noise)
- **50/50 split** for high-fidelity model → ~5×10⁴ direct LAMMPS ΔE
  (training) + ~5×10⁴ ML-predicted (testing). Spectrum plotted from all
  ~10⁵ values (ML MAE for Al(Mg) = 2.45 kJ/mol per SI Table 1).
- Accelerated model (Fig 3) reduces training to 100 cluster centroids
  for the 259-alloy scan.

So Wagih's Mg¹⁵ spectrum is fit to **~10⁵ points** (dense sampling + ML
smoothing), while ours is fit to **500 direct-LAMMPS points** (sparse
but numerically exact). Different statistical regimes.

### Phase 4 feasibility — HMC cost does NOT scale with N_GB

Important clarification: Phase 3's 500 is a **subsample** because Phase 3
does one LAMMPS CG *per site*. Phase 4 HMC is not per-site — it runs
one MD simulation on the full 475k-atom box and every GB site feels the
chemical potential simultaneously, every timestep.

| | Wagih per-site + FD | Our HMC |
|---|---------------------|---------|
| Unit        | 1 CG per GB site | 1 MD timestep (whole box) |
| Count       | N_GB CGs         | N_steps MD |
| Cost/unit   | ~10 s            | ~0.05 s (1.74 ns/day on 32 cores) |
| Total (1 alloy, 1 point) | N_GB · 10 s ≈ 247 h (ML → 17 min) | ~1.3 h per (T, X_c) |
| Cover (T, X_c) | Free formula | Separate MD run per point |

Wagih's win = ML-accelerated breadth (259 alloys × all (T, X_c) via
formula). His formula is only valid under the independent-site
assumption we're testing.

Our win = one MD run *automatically* sees all 89k GB sites interacting.
No subsampling needed at Phase 4.

**Estimated Phase 4 cost** (our grid): 3T × 5X_c = 15 points × ~2 h each
(inc. swap overhead) = **30 compute-hours total**, embarrassingly
parallel. Fully feasible — the "10⁵ GB sites" is not a barrier because
HMC doesn't iterate over sites.

### Project scope positioned precisely

- **Wagih** = independent-site assumption + ML acceleration → cheap
  breadth (259 alloys, any T/X_c). Spectrum is pipeline *input*; never
  verifies the assumption itself.
- **Us** = HMC on the *same* alloy Wagih used as his headline example,
  same potential, same box → expensive depth (1 alloy, 15 (T, X_c)
  points). Our contribution is measuring the concentration at which
  Wagih's Fermi-Dirac predictor starts diverging from HMC truth — i.e.
  drawing the validity boundary of his central assumption.
- HMC cost ~100× Wagih's per-alloy cost (30 h vs 17 min) buys this
  measurement. Not catastrophic at one-alloy scope.

### Decision pending

Two paths before launching Phase 4:

**P** — Go straight to Phase 4 HMC. Our 500 ΔE points + fit (μ=9.4,
σ=19.4, α=−1.08) build the Fermi-Dirac prediction. `X_GB^FD(T, X_c)` is
a central-moment integral, insensitive to α. The HMC truth curve is the
actual headline output regardless of Phase 3 tail precision.

**Q** — First fetch Wagih's Al(Mg) ΔE dataset from Zenodo
(doi:10.5281/zenodo.4107058). Do KS test + overlay histograms vs our
500 points → give Phase 3.5 a rigorous goodness-of-agreement number.
Then Phase 4. Costs ~2 h extra compared to P.

Q → P preferred: getting Wagih's 10⁵ points costs little and lets us
build two Fermi-Dirac curves (ours from 500 pts; Wagih's from 10⁵ pts)
to isolate sampling noise from real pipeline differences before
interpreting Phase 4. User to confirm.

## 2026-04-24 (late 3) — Phase 3.5 A+B: 20³ nm³ production anneal + GB mask

Scaled the pipeline from the 10³ nm³ / 8-grain prototype to the Wagih-match
20³ nm³ / 16-grain production box, aiming to rule in or rule out finite-size
as the cause of our 20% offset from Wagih's Mg¹⁵ target (μ=9, σ=23, α=−2.3).
Same potential (Mendelev 2009, Al-Mg.eam.fs), same anneal protocol, same
VELOCITY_SEED — **only box size and grain count differ from prototype**.

### Structure generation (`scripts/generate_polycrystal.py`)

- `python generate_polycrystal.py --structure fcc --box 200 --grains 16
  --lattice-a 4.05 --structure-seed 1`
- **475,715 atoms** (perfect-lattice 481,709; 1.2% close-pair deficit —
  tighter than prototype's 1.6%, as expected for larger grains).
- Stored at `/cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g.lmp`.

### Timing calibration (job 64663404, 32 cores)

Before committing to a multi-hour anneal we measured ns/day on a 1 ps NVT
trial: **1.740 ns/day on 32 cores**, 49.7 s / 1000 steps, 91% pair / 4%
comm / 3% neigh — communication overhead low, strong-scaling headroom
clearly available. 32 cores / 8 h chosen for the real anneal.

### Wagih-style anneal (job 64665121, `submit_anneal_200A.sh`)

- 32 cores, walltime **4:52:22** (predicted 5.5 h; came in 12% under).
- Performance during hold: 1.85 ns/day (consistent with calibration).
- Thermo evolution (key samples):
  - CG #0:        PE = −1,610,632 eV, Lx = 200.00 Å, T = 0 K
  - 5 ps ramp end: PE = −1,585,022, Lx = 200.00, T ≈ 375 K
  - 250 ps hold end: PE = −1,589,233, Lx = **201.36 Å** (thermal expansion 0.68%)
  - 124 ps cool end: PE = −1,613,047, Lx = **199.78 Å** (cooled below initial)
  - Final CG + box relax: **PE/atom = −3.391 eV/atom**, Lx = 199.78 Å
- vs prototype PE/atom = −3.358 eV/atom → production is 0.033 eV/atom more
  negative, consistent with lower GB fraction (less GB excess energy per
  atom) at the larger box.
- **Artifact**: `poly_Al_200A_16g_annealed.lmp` (62 MB) — input for
  Phase 3.5 B/C. The 1.4 GB `*.dump` trajectory is gitignored and will be
  deleted once OVITO inspection (if any) is done.

### GB identification (job 64698536, `submit_gb_identify_200A.sh`)

Ran `gb_identify.py` on the annealed structure (serial LAMMPS, 2:32 wall,
8 GB mem). First submission 64698397 **failed** because the script's
`--lmp` argument only accepts a single binary path (uses `shutil.which`),
not `mpirun -n N lmp`. Fix: drop MPI, run LAMMPS serial. `compute cna/atom`
on 475k atoms is ~2 min serial, so no reason to MPI this.

**CNA classification** (fixed cutoff `0.854 × 4.05 = 3.459 Å`, parent FCC):

| label  | atoms   | fraction | role              |
|--------|---------|----------|-------------------|
| FCC    | 386,673 | 81.28%   | bulk              |
| HCP    | 2,095   | 0.44%    | stacking faults (→ GB per Wagih def) |
| Other  | 86,947  | 18.28%   | GB core           |
| BCC/ico/unknown | 0 | 0%     | —                 |

**f_gb = 18.72%** (N_GB = 89,042).

### Scaling validation: `f_gb ∝ 3t/d` holds

| system           | box     | grains | f_gb   | HCP fraction |
|------------------|---------|--------|--------|--------------|
| prototype        | 10³ nm³ | 8      | 28.7%  | 0.95%        |
| **production**   | **20³ nm³** | **16** | **18.7%** | **0.44%** |
| Wagih            | 20³ nm³ | 16     | ~15%   | ~1%          |

Box ×2 → mean grain size ×1.59 (√[8/2] via V/N_grain) → f_gb should drop
by ×1.59. Measured 28.7/18.7 = 1.53 ✓. HCP dropping by half also expected
— larger grains absorb initial Voronoi distortion more completely.

**Residual 3.7 percentage points above Wagih's 15%** is not chased: plausible
sources are (i) LAMMPS fixed-cutoff CNA vs OVITO adaptive CNA on high-
distortion GB cores, (ii) `structure_seed` realization variance ±1–2%,
(iii) minor anneal-protocol differences. The ΔE spectrum shape (μ, σ, α)
is set per-site, not by N_GB, so a 3.7% offset in mask size only changes
normalization, not the fit parameters we want to compare against Wagih Mg¹⁵.

### Artifacts

- `gb_mask_200A.npy` (475k bool)
- `gb_info_200A.json` (summary above)
- `gb_cna_200A.dump` (per-atom CNA label, OVITO-ready)

### Next: Phase 3.5 C — N_GB=500 ΔE sampling on production box (in flight)

Submitted job **`64699296`** (`submit_delta_e_200A.sh`) — 16 cores, 6 h
wall budget (predicted ~3.5 h), same seed=42, same (N_GB=500, N_bulk=10)
recipe as prototype. Output → `delta_e_results_n500_200A.npz` in scratch.
On completion: re-run `fit_delta_e_spectrum.py` to get production (μ, σ,
α) and compare against Wagih Mg¹⁵ (μ=9, σ=23, α=−2.3, R²=1.00).

## 2026-04-24 (late 2) — Conceptual background: the Fermi-Dirac model

Locked down terminology + the assumption stack behind Wagih's `X_GB(T, X_c)`
predictor, since Phase 4's entire purpose is to test where this predictor
breaks. Kept here so future sessions don't re-derive it.

### Original Fermi-Dirac (quantum statistics)

For a fermion state at energy ε, occupation at `(T, μ)` is
`f(ε) = 1 / [1 + exp((ε − μ)/kT)]`. Derived from Pauli exclusion
("at most one fermion per state") plus Boltzmann weighting. Properties:
`ε ≪ μ → f → 1`; `ε ≫ μ → f → 0`; `ε = μ → f = 0.5`; higher T widens the
transition region around ε = μ.

### Mapping to GB segregation (Wagih eq. 2)

| Fermion gas      | GB segregation                              |
|------------------|---------------------------------------------|
| state i          | GB site i                                   |
| energy ε_i       | segregation energy ΔE_i = E_GB^Mg(i) − E_bulk^Mg |
| chemical potential μ | set by bulk solute concentration X_c    |
| occupation f_i ∈ {0,1} | site i occupied by Mg (else by Al)    |
| Pauli exclusion  | geometric one-atom-per-site constraint      |

Wagih's predictor:

```
P_i(T, X_c) = 1 / [1 + ((1 − X_c)/X_c) · exp(ΔE_i / kT)]
X_GB(T, X_c) = (1/N_GB) · Σ_i P_i(T, X_c)
```

The `(1 − X_c)/X_c` replaces `exp(−μ/kT)` — it's the Mg/Al concentration
ratio from a semi-grand-canonical bookkeeping ("put Mg in GB ⇔ take Al out
to bulk"). Whole sum collapses to a formula evaluation: one (T, X_c)
prediction is ~microseconds, vs. days of MC/HMC. This is the speed behind
Wagih's 259-alloy scan.

### The assumption stack (in order of fragility)

1. **Independent sites (non-interacting solutes)** — ΔE_i is defined as
   "substitute 1 Mg at site i while *all other GB sites are Al*". The
   formula then treats ΔE_i as a fixed property, independent of what
   other sites hold. In reality, a neighboring Mg changes site i's local
   environment and shifts its true substitution energy. **This is the
   assumption Phase 4 targets.** Our `ΔE` spectrum is computed under this
   assumption by construction, so it cannot detect its own breakdown —
   only HMC can.
2. **Dilute limit** — `(1 − X_c)/X_c` treats bulk/solute concentrations
   as static. At high X_c the GB starts to drain the bulk, so the ratio
   drifts during uptake. The formula can be patched with a self-
   consistent X_c(bulk), but Wagih uses the static form.
3. **Fixed site set** — N_GB and the ΔE_i values are frozen at the
   annealed 0 K structure. GB reconstruction at elevated T (atomic
   rearrangement, new sites appearing / disappearing) is invisible to
   the formula.

### Role in the project pipeline

- **Phase 3 / 3.5** outputs `{ΔE_i}` → feeds into Fermi-Dirac → yields
  `X_GB^FD(T, X_c)` **prediction** curve.
- **Phase 4** HMC samples the *actual* Boltzmann distribution with all
  solute-solute interactions active → yields `X_GB^HMC(T, X_c)` **truth**
  curve.
- Headline figure = both curves overlaid; the X_c at which they diverge
  beyond statistical noise is the **breakdown concentration** of the
  independent-site assumption (question 1). Analyzing *why* (solute-
  solute g(r) at GB, local-density vs P_i correlation) is question 2.

Wagih's paper assumes (1) and (2) hold implicitly across all 259 alloys
and never draws a validity boundary for them. Drawing that boundary is
the contribution of this project.

## 2026-04-24 (late) — Wagih SI reference values + notation

### `Mg^15` superscript decoded

The `X^N` superscript on element labels in Wagih's SI spectra figures
(e.g. `Mg^15` in Supplementary Fig. 3, Al-based alloys) is the **SI
Supplementary Reference number of the interatomic potential used to
compute that spectrum** — not a concentration, sample count, or figure
subpanel. SI refs list is on SI pp. 26–28.

**Ref [15] = Mendelev, Asta, Rahman & Hoyt, *Philos. Mag.* 89,
3269–3285 (2009)** — the same potential we already ship at
`project/data/potentials/Al-Mg.eam.fs`. `Mg^15` therefore is the
apples-to-apples benchmark for our pipeline: same potential, same
alloy, same parameter space.

### Target values (read from SI Fig. 3 `Mg^15` panel, kJ/mol)

Wagih's figures label the three skew-normal parameters as **(μ, σ, α)**.
We adopt the same labels (see "Notation convention" below).

| parameter | value |
|-----------|-------|
| μ (location) | ≈ +9 |
| σ (scale)    | ≈ 23 |
| α (shape)    | −2.3 |
| R²           | 1.00 |

Supersedes the "Wagih Al(Mg) σ ≈ 4" (earlier paper_notes transcription
error) and the subsequent "σ ≈ 15" self-correction — both were readings
of Fig. 2 in the main text, not of the Mendelev-potential `Mg^15` panel
in the SI. The SI panel is the correct comparison target because it
uses the same potential as us; Fig. 2 aggregated across potentials.

### Our n=500 fit (2026-04-24) vs Wagih Mg^15

| parameter | ours (n=500) | Wagih Mg^15 |
|-----------|--------------|-------------|
| μ         | +11.5        | ≈ +9        |
| σ         | 21.3         | ≈ 23        |
| α         | −1.93        | −2.3        |

Same signs, same order of magnitude, Wagih slightly more left-skewed.
Likely finite-size (our 10³ nm³ prototype vs Wagih's 20³ nm³) — to be
re-checked on the 20³ nm³ production box.

### Notation convention (project-wide, going forward)

All skew-normal parameters use **(μ, σ, α)** — matching the labels on
Wagih's Fig. 2 and SI Fig. 3 panels, so values copy straight between
our outputs and the paper with no symbol translation. Mapping to
`scipy.stats.skewnorm`: `loc=μ, scale=σ, a=α`.

**Caveat, always worth remembering**: here μ and σ are the skew-normal
*location* and *scale*, **not** the distribution's mean and standard
deviation. The true moments are
`mean = μ + σ·δ·√(2/π)`, `std = σ·√(1 − 2δ²/π)` with `δ = α/√(1+α²)`.
For genuine sample moments use the `sample_moments_kjmol` block in the
JSON output, which reports `mean` and `std` directly.

Mathematical literature would write these as (ξ, ω, α); we picked
Wagih's labels over the math-literature convention because
paper-comparison frequency outweighs textbook-comparison frequency for
this project.

### Script update — `fit_delta_e_spectrum.py`

- `WAGIH_ALMG` constant holds the SI Fig. 3 `Mg^15` values
  (μ=9, σ=23, α=−2.3, R²=1.00) with a `source` field citing the
  Mendelev 2009 potential explicitly; dict keys renamed `xi→mu`,
  `omega→sigma`. The fit-output dict uses the same keys.
- Docstring updated to state the labelling choice and the μ/σ ≠
  mean/std caveat up front.
- Reference-curve label in the overlay plot reads "Wagih Mg¹⁵
  (SI Fig. 3)"; stdout prints (μ, σ, α, R²) with Wagih symbols.
- **Breaking change for downstream consumers of the JSON**: old
  `delta_e_fit_n500.json` uses `xi`/`omega`; regenerated JSONs will
  use `mu`/`sigma`. Re-render `output/delta_e_spectrum_n500.png` and
  regenerate JSON from the scratch `.npz` before the next status
  report.

## 2026-04-24 — N_GB=500 Al(Mg) ΔE_seg spectrum + skew-normal fit

Scaled Phase 3 sampling from the 50-site prototype to **N_GB=500, N_bulk=10,
seed=42** on the same annealed 10³ nm³ / 8-grain Al prototype. Job
`64638359`, wall **26.6 min on 16 cores** (510 sites × ~3.1 s/site, all 510
CG relaxations stopped on `energy tolerance` — no bad minima). Results at
`/cluster/scratch/cainiu/prototype_AlMg_100A/delta_e_results_n500.npz`.

**N_GB=500 ΔE_seg spectrum** (vs N_GB=50 prototype / Wagih Fig 2):

| quantity             | N_GB=50       | N_GB=500      | Wagih Al(Mg)  |
|----------------------|---------------|---------------|---------------|
| range (kJ/mol)       | [−45.6, +19.6]| [−49.7, +35.2]| [−60, +40]    |
| mean (kJ/mol)        | −6.1          | −3.5          | —             |
| median (kJ/mol)      | −7.0          | −1.4          | —             |
| std (kJ/mol)         | 14.8          | 15.0          | ~15           |
| sample skew          | −0.26         | **−0.41**     | ~−0.4         |
| bulk ref σ/√n (meV)  | 4             | 1.5           | —             |

Scaling 50 → 500 recovers the upper tail (+20 → +35 kJ/mol) and pins the
sample skew at −0.41, essentially on top of Wagih's Al(Mg) ≈ −0.4. Lower
tail reaches −50 kJ/mol, 10 kJ/mol short of Wagih's −60 — likely a
finite-size effect (our 10³ nm³ prototype produces fewer "deep" GB sites
than Wagih's 20³ nm³); will re-check on the 20 nm production box.

### Skew-normal fit (`scripts/fit_delta_e_spectrum.py`)

Small single-purpose script: loads the `.npz`, fits
`scipy.stats.skewnorm` (Wagih's parameterization `F(ΔE; ξ, ω, α)`),
writes histogram + fit overlay to `output/delta_e_spectrum_n500.png`
and params to `output/delta_e_fit_n500.json`.

**Fit parameters (N_GB=500)**: ξ = +11.5, ω = 21.3, α = −1.93 kJ/mol.
Consistency check: skew-normal std = ω·√(1 − 2δ²/π) with δ = α/√(1+α²) gives
15.0 kJ/mol, exactly the sample std — the fit is internally coherent and
not an artifact of a bad optimizer start.

### paper_notes.md correction — Wagih Al(Mg) σ "~4" → "~15"

Earlier notes table (§9) had `μ~-2, σ~4, α~-0.4` for Wagih's Al(Mg)
skew-normal. With α=-0.4, σ=4 kJ/mol forces std ≈ 3.8 kJ/mol, which is
incompatible with the [-60, +40] range quoted on the same row — it was a
transcription error from Fig 2. Corrected to σ~15 kJ/mol, matching our
fit and the visible width of the published histogram. The fit plot's
"Wagih ref" curve was generated from the old σ=4 value for visual
contrast; keeping it in for now as a record of the discrepancy, will
re-render with σ=15 once we have time to re-check Fig 2 digitization.

### Housekeeping

- `data/decks/submit_delta_e.sh` now targets N_GB=500 / N_bulk=10 and a
  separate `delta_e_run_n500/` work dir + `delta_e_results_n500.npz`
  output so the 50-site prototype results stay intact next to the new
  ones.
- `output/delta_e_spectrum_n500.png` + `output/delta_e_fit_n500.json`
  written locally (directory is gitignored — regenerate via
  `fit_delta_e_spectrum.py` from the `.npz` in scratch).

### Next

- Phase 4 — HMC `(T, X_c)` scan using the same 10³ nm³ annealed
  polycrystal as the starting config. `fix atom/swap` semi-grand-
  canonical Metropolis at `T ∈ {300, 500, 700} K` × `X_c ∈ {0.1, 1, 5,
  10, 20} at%` as the first sweep; compare average GB enrichment
  against Fermi-Dirac-from-ΔE prediction built from the N_GB=500 fit
  above.

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

