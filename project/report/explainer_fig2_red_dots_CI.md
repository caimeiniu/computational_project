# Explainer: Fig 2 — what are the red dots and the error bars?

**Status**: PENDING re-explanation. User said "还是没有搞的很清楚" on
2026-05-02 and asked me to write the answer down here so we can pick it
up next session. **Do not delete after re-explaining; convert this file
into the §5 Figure 2 walkthrough subsection if the explanation lands.**

---

## User's questions (verbatim)

1. "红色远点是我们实际算出来的均值吗"
   = Are the red dots the actual mean we calculated?
2. "置信区间是怎么算的呢"
   = How is the confidence interval calculated?

## Why the previous explanation did not land (hypothesis)

I delivered the answer in formula form (Wald CI = p̂ ± 1.96·SE) with one
numerical example. User responded "还是没有搞的很清楚". Likely reasons:

- The "binomial proportion" framing was abstract — the user may not
  have recognised that "fraction of Mg in a bin" *is* the same thing as
  a "sample mean of {0,1} variables".
- The Wald formula appeared without geometric / coin-flip intuition; it
  felt like a recipe rather than a derived expression.
- The CI was introduced as a definition rather than as something that
  *follows* from the question "if I redid this with a fresh snapshot,
  how much would the dot move?".

**Re-explanation strategy to try next session:**

1. Start with one bin (the −27.40 kJ/mol bin). Forget the rest of the
   figure for a moment.
2. Show the 38 sites in that bin as 38 boxes. Mark which 8 are Mg
   (●) and which 30 are Al (○). The red dot height is just
   `count(●) / 38 = 0.211`.
3. Coin-flip analogy: imagine a biased coin with unknown probability `p`
   of heads. We flip it 38 times, get 8 heads. p̂ = 8/38 is our best
   guess at `p`, but we are not certain. The CI is: "values of `p` that
   are consistent with what we observed". For Wald, the math reduces to
   p̂ ± 1.96 √[p̂(1−p̂)/n].
4. THEN pull back and say: each red dot in the figure is one such
   experiment, run on a different ΔE bin. The bins differ in `n`, which
   is why the error bars differ in width (small n at the tails → wide
   bars).
5. Only after the intuition lands, show the formula and the table.

---

## Answer 1: Are the red dots the empirical mean?

**Yes.** Concretely, for the most-favourable bin at ΔE = −27.4 kJ/mol
in the X_c = 0.075 snapshot:

```
1. Take the 500 reference sites, group by ΔE_i value, 10 equal-width
   bins from −48 to +36 kJ/mol.
2. The −27.4 bin contains n = 38 sites.
3. In the HMC snapshot, look at each of those 38 sites: is the atom
   there Mg or Al? Result: 8 are Mg, 30 are Al.
4. The red dot's height is p̂ = 8 / 38 = 0.211.
```

This is mathematically identical to averaging 38 binary values
(Mg=1, Al=0). So "empirical fraction" / "binomial proportion" /
"sample mean of binary outcomes" all describe the same number.

## Answer 2: How is the CI calculated?

The script `scripts/solute_correlation_analysis.py:207-225` computes a
**Wald 95 % normal-approximation CI**, *not* Wilson (the code comment
mislabels it; the actual formula is Wald):

```
SE   = √[p̂ (1 − p̂) / n]
CI   = [p̂ − 1.96·SE,  p̂ + 1.96·SE], clipped to [0, 1]
```

For the −27.4 bin (n = 38, p̂ = 0.211):

```
SE       = √[0.211 × 0.789 / 38]  = √0.00438       ≈ 0.0662
1.96·SE  ≈ 0.130
CI       = [0.211 − 0.130, 0.211 + 0.130] = [0.081, 0.341]
                                              ↑ matches JSON [0.081, 0.340]
```

**Intuition**: "If I redid the experiment with a fresh independent
snapshot, p̂ would land in this range about 95 of every 100 reruns."

## Per-bin table for X_c = 0.075 (verified from JSON)

```
  bin (ΔE)    n       p̂    low95   high95    width
    -44.25    7    0.429    0.062    0.795    0.733
    -35.83   25    0.440    0.245    0.635    0.389
    -27.40   38    0.211    0.081    0.340    0.259  ← Fig 2 main-panel ΔP_i
    -18.97   85    0.388    0.285    0.492    0.207
    -10.55   88    0.295    0.200    0.391    0.191
     -2.12  105    0.229    0.148    0.309    0.161
      6.30   92    0.174    0.096    0.251    0.155
     14.73   43    0.233    0.106    0.359    0.253
     23.16   13    0.308    0.057    0.559    0.502
     31.58    4    0.250    0.000    0.674    0.674
```

The error bar width tracks `n`: small n at the tails → wide bars; large
n in the middle → narrow bars. This is the visual reason the leftmost
and rightmost dots have "exaggerated" error bars in the figure — it is
not a bug, it is small-sample uncertainty.

## Caveat (advisor-question level, currently not in caption)

Wald CI assumes the 38 site outcomes in a bin are *independent*. They
are not strictly independent because:

1. Mg–Mg interactions create local correlation — the very signal the
   figure is meant to expose.
2. Global mass conservation (∑ Mg = const) couples all sites.

So Wald CI **slightly under-estimates** the true uncertainty (treats
correlated samples as independent). Strict treatment would use a block
bootstrap on the trajectory, or thin samples by autocorrelation length.
We chose Wald for simplicity; for the favourable end the gap ΔP = 0.77
is large enough that no realistic CI inflation flips the conclusion.

## File pointers

- Source script: `scripts/solute_correlation_analysis.py:207-225`
  (computes p_low95 / p_high95, mislabelled as "Wilson approx" but
  actually Wald).
- Cached values: `output/solute_correlation_analysis.json` →
  `site_occupation_vs_energy.snapshots["0.075_preseg"]` →
  `p_empirical`, `p_low95`, `p_high95`, `n_per_bin`.
- Figure caption: `report/README.md` §5 Figure 2 (already includes the
  Wald formula and the worked example for the −27.4 bin).
