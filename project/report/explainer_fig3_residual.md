# Fig 3 — Residual vs Wagih (X_c=0.075) — talk version

## Status

Test version of Fig 3 redesign. Source: `scripts/residual_vs_wagih_test.py`;
output `output/residual_vs_wagih_test_Xc0.075.{png,json}`. NOT yet copied
to `report/figures/` and `report/README.md` §5 Fig 3 still describes the
old window-slope version. Decision deferred to next session after
re-explanation.

## Why we replaced the old Fig 3

Old: P̄ vs n_local inside the [−30, −5] kJ/mol window, slope = −0.083.

Methodological problem: P_Wagih varies 0.99 → 0.21 across that window
(at X_c=0.075, T=500 K, kT=4.16 kJ/mol). So binning only by n_local
inside the window does NOT fix ΔE — the slope conflates the Mg-Mg
interaction signal with omitted-variable bias from ΔE varying within
the window. Calling it "controlling for ΔE" was wrong.

New (residual approach): for each of the 500 reference sites, compute
r_i = I_i^HMC − P_Wagih(ΔE_i, X_c, T), then bin by integer n_local.
Tower property: E[r_i | n_local_i = n] = 0 in Wagih's null hypothesis,
regardless of how (ΔE_i, n_local_i) are jointly distributed
geometrically. So slope ∂r̄/∂n_local ≠ 0 is unbiased evidence of
site-level Mg-Mg interaction, and the sign of the slope tells you
attraction vs repulsion.

## 30-second talk version (deliver slowly tomorrow — three steps, 10s each)

### Step 1 — what's on the axes

> "Y 轴是 **残差** —— 每个 GB site 在 HMC 里的实际占据（Mg=1, Al=0），
> 减去 Wagih 独立位点公式给这个 site 的预测概率 P_Wagih(ΔE_i)，再按
> n_local 取平均。所以 **y=0 (灰虚线) 就是 Wagih 完全预测正确的情况** ——
> 如果 Mg 之间没有相互作用，所有点都应该贴在这条线上。"

### Step 2 — what we see

> "n=2 个邻居的 site 残差 +0.28（Wagih 反而低估），n=9 个的残差 −0.37
> （Wagih 严重高估），slope = −0.069 / Mg-邻居 —— **每多一个 Mg 邻居，
> 这个 site 比 Wagih 预测多 under-fill 7 个百分点**。"

### Step 3 — physics

> "ΔE 这个变量已经被 Wagih 公式吃进去了，所以残差对 n_local 的依赖
> **只能来自 site 之间的相互作用**。负 slope 说明：一个 site 周围已经
> 有 Mg 时，再放 Mg 进去要付能量代价 —— 就是 site-level Mg-Mg 排斥。"

### Optional Q&A handle

> "和 Fig 1 g(r) > 1 的互补关系：g(r) 看到了 Mg 在空间上扎堆，但分不清
> 是化学吸引还是深 ΔE 位点几何相邻；这张残差图把 ΔE 控制掉，方向
> 拍板成排斥。"

## Key numbers (from `output/residual_vs_wagih_test_Xc0.075.json`)

| Quantity | Value | Meaning |
|---|---|---|
| N_sites total | 500 | All reference sites used (no window) |
| N_sites in valid bins (n_sel ≥ 10) | 488 | n_local ∈ [2, 9] |
| Global residual mean ⟨r⟩ | −0.133 | "Uniform breakdown" — Wagih over-predicts globally by 13 pp |
| Slope ∂r̄/∂n_local | −0.069 / Mg-nbr | Per Mg-neighbour, this site is 6.9 pp more under-filled vs Wagih |
| Intercept | +0.208 | Linear fit's y-value at n_local=0 |

Sign and order of magnitude consistent with the old window slope
−0.083; quantitative shift attributable to (a) full-spectrum vs
[−30,−5] window, (b) ΔE confound removed.

## Pending decisions (resolve next session)

1. **Replace `report/figures/03_repulsion_summary.png`** with the new
   single-panel residual PNG?
2. **Filename**: keep `03_repulsion_summary.png` or rename to
   `03_residual_vs_wagih_Xc0.075.png` (more descriptive)?
3. **Rewrite `report/README.md` §5 Fig 3** caption + likely-Qs to
   switch from window-slope narrative to residual narrative?
4. **X_c-sweep saturation crossover** (currently in README §5 as a
   key observation): drop entirely, or keep as a text remark only?
