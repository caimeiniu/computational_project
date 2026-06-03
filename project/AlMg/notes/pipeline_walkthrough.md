# Al(Mg) HMC GB-segregation pipeline — end-to-end walkthrough

*Written 2026-05-13 as a reference for future-me or new teammates joining the
project. Uses the 200³ Å Al(Mg) production run as the worked example; Pt(Au)
follows the same 8 stages with system-specific parameters substituted per
Table 1 in `report/methods_formula_patches.tex`.*

## DAG of the 8 pipeline stages

```
[Stage 1]  generate_polycrystal.py     login node    1 min
              ↓ poly_*.lmp
[Stage 2]  anneal_AlMg.lammps          SLURM         ~1–2 h     (NVT/NPT 退火)
              ↓ poly_*_annealed.lmp
[Stage 3]  gb_identify.py              login node    1 min      (a-CNA GB 识别)
              ↓ gb_mask_*.npy + gb_info_*.json
[Stage 4]  sample_delta_e.py           SLURM         ~15-30 min (per-site CG, 500 sites)
              ↓ delta_e_results_*.npz + meta JSON
[Stage 5]  fit_delta_e_spectrum.py     login node    sec        (skew-normal fit)
           fermi_dirac_predict.py      login node    sec        ((μ,σ,α) → FD curves)
              ↓ output/{delta_e_fit,fd_curves}_*.json + .png
[Stage 6]  hmc_AlMg_v2.lammps          SLURM         12-24 h    (每个 (T,X_c) 一个作业)
              ↓ hmc_*.{log,dump,rst1,rst2,_final.lmp}
[Stage 7]  hmc_xgb_timeseries.py       login node    sec        (X_GB(t) + block bootstrap)
              ↓ output/hmc_T*_Xc*.{json,png}
[Stage 8]  report_headline_T500_9pt.py login node    sec        (aggregate 9 点 → panel d)
              ↓ output/panel_d_*.{json,png}
```

## Two top-level storage locations

- **`/cluster/scratch/cainiu/production_AlMg_200A/`** (Pt(Au) 平行目录:
  `/cluster/scratch/cainiu/prototype_PtAu_100A/`) — 大文件(数百 MB ~ 数 GB):
  polycrystal `.lmp`, annealed `.lmp`, dump 轨迹, restart 文件, ΔE npz 结果。
  Euler scratch 15-day TTL,见 `~/scripts/backup_scratch.sh` 维护增量备份。
- **`/cluster/home/cainiu/Computational_modeling/project/output/`** — 小文件
  (KB 量级):summary JSON, figure PNG, panel-d aggregated 数据。git-tracked
  via 例外条目;`output/` 默认 gitignore 但 individual JSON 大多保留。

---

## Stage 1: Polycrystal generation

**脚本**:`scripts/generate_polycrystal.py`(login node,~1 min)

**输入**(命令行参数,无外部文件):
```bash
python scripts/generate_polycrystal.py \
  --structure fcc --box 200 --grains 16 --lattice-a 4.05 \
  --structure-seed 1 --types 2 \
  --out /cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g.lmp
```

关键参数:
- `--box 200`:200 Å 立方边长(Pt(Au) 用 100)
- `--grains 16`:Voronoi seed 数量(Pt(Au) 用 8)
- `--lattice-a 4.05`:FCC Al 晶格常数(Pt 用 3.98)
- `--types 2`:LAMMPS atom types(host=type 1, solute=type 2)
- `--structure-seed 1`:伪随机种子(锁定 GB 拓扑)
- 内部 `r_cut = a/(2√2) ≈ 1.43 Å` 自动剔除 close-pair

**输出**:1 个 LAMMPS data 文件 `poly_Al_200A_16g.lmp`

内容(LAMMPS native format,475{,}715 atoms,2 atom types,box 200³ Å,
全部 type=1)。这阶段还没退火,close-pair 已剔除。

---

## Stage 2: Anneal

**Submit script**:`AlMg/data/decks/submit_anneal_200A.sh` → sbatch
**LAMMPS deck**:`AlMg/data/decks/anneal_AlMg.lammps`(被 submit 调用)

```bash
sbatch project/AlMg/data/decks/submit_anneal_200A.sh
```

SLURM 参数:32 ranks × 8 h budget,实际 ~1–2 h。

**LAMMPS deck 输入**(`-var` 传):
- `annealed`:输入 `poly_Al_200A_16g.lmp`
- `outstub`:输出文件前缀
- `T_HOLD=373`:hold 温度(K),= 0.4 × T_melt(Al)
- `EL1=Al MASS1=26.9815 EL2=Mg MASS2=24.305`:元素 + 质量
- `POTFILE`:`AlMg/data/potentials/Al-Mg.eam.fs`

**LAMMPS 5 阶段**(`anneal_AlMg.lammps`):

1. CG minimize(吸收 close-pair 残应力,`etol/ftol = 1e-6/1e-8`)
2. NVT 1→373 K 5 ps(τ_T=100 fs)
3. NPT hold @ 373 K, 250 ps(τ_T=100 fs, τ_P=1000 fs, P=0)
4. NPT cool 373→1 K, 3 K/ps(同 barostat)
5. `fix box/relax` + 终 CG(`etol/ftol = 1e-8/1e-10`)

**输出**(全部到 `/cluster/scratch/cainiu/production_AlMg_200A/`):

| 文件 | 内容 | 用途 |
|---|---|---|
| `poly_Al_200A_16g_annealed.lmp` | 退火后最终 atomic structure | ↓ Stage 3/4/6 的核心输入 |
| `poly_Al_200A_16g_anneal.dump` | 退火轨迹(per-step positions + per-atom energy) | OVITO 可视化 / 检查热应变 |
| `poly_Al_200A_16g_anneal.log` | LAMMPS thermo log:T(t), PE(t), box volume(t)... | 检查温度/压力/能量是否收敛 |
| `poly_Al_200A_16g.rst1, rst2` | LAMMPS restart 二进制(轮换) | 接续跑用 |
| `anneal_AlMg_200A-<jobid>.{out,err}` | SLURM stdout/stderr | debug |

**验证 checklist**(从 .out / log 找):`T_HOLD = 373 K ± 5 K`,
`force_max < 5 eV/Å`,atom count 保留。

---

## Stage 3: GB identification(adaptive CNA)

**脚本**:`scripts/gb_identify.py`(login node,~1 min)

```bash
python scripts/gb_identify.py \
  /cluster/scratch/cainiu/production_AlMg_200A/poly_Al_200A_16g_annealed.lmp \
  --parent fcc --lattice-a 4.05 \
  --out-mask  /cluster/scratch/cainiu/production_AlMg_200A/gb_mask_200A.npy \
  --out-report /cluster/scratch/cainiu/production_AlMg_200A/gb_info_200A.json \
  --out-dump   /cluster/scratch/cainiu/production_AlMg_200A/gb_cna_200A.dump
```

**输入**:`poly_Al_200A_16g_annealed.lmp`(LAMMPS data)

**输出**:

| 文件 | 数据结构 |
|---|---|
| `gb_mask_200A.npy` | numpy bool array,长度 N_tot=475{,}715。True = 该 atom_id 是 GB site,False = bulk FCC |
| `gb_info_200A.json` | summary 统计(下面) |
| `gb_cna_200A.dump` | CNA-annotated LAMMPS dump(每个原子带 CNA classification),OVITO 可视化 |

`gb_info_200A.json` 实际结构:

```json
{
  "n_atoms": 475715,
  "n_gb": 89042,
  "f_gb": 0.18717509...,
  "cna_counts": {"unknown": ..., "fcc": ..., "hcp": ..., "bcc": ..., "ico": ...},
  "parent_structure": "fcc",
  "bulk_cna_int": 1,
  "cna_cutoff_angstrom": 3.4587,
  "data_file": ".../poly_Al_200A_16g_annealed.lmp"
}
```

a-CNA 把每个 atom 分类为 `{fcc, hcp, bcc, ico, unknown}`。**GB site = 不是 fcc 的**。

---

## Stage 4: 每点 ΔE_seg 采样(production)

**Submit script**:`AlMg/data/decks/submit_delta_e_200A.sh` → sbatch
**驱动**:`scripts/sample_delta_e.py`(被 submit 调用)

```bash
sbatch project/AlMg/data/decks/submit_delta_e_200A.sh
```

16–32 ranks × 4 h budget,实际 ~15-30 min for n=500 (Al-Mg);Pt(Au) ~1h23m
因 Pt EAM ~5× 慢(CG 公差严+刚度大)。

**输入**:
- `--annealed`:`poly_Al_200A_16g_annealed.lmp`
- `--gb-mask`:`gb_mask_200A.npy`
- `--potential`:`Al-Mg.eam.fs`
- `--n-gb 500`:采 500 个 GB sites
- `--n-bulk 10`:10 个 bulk 参考位
- `--seed 42`:随机种子
- `--elements "Al Mg" --masses "26.9815 24.305"`

**内部循环**(每个 GB site):

1. Copy annealed config
2. `type[i]`:1 (Al) → 2 (Mg)
3. CG minimize(`etol = ftol = 1e-25`)
4. Record total PE
5. Restore site → 下一个

Bulk reference:同样 10 次 single Al→Mg substitution,在 bulk 位点。

**输出**:

| 文件 | 数据结构 |
|---|---|
| `delta_e_results_n500_200A_tight.npz` | per-site arrays(下面详) |
| `delta_e_meta_n500_200A_tight.json` | meta info(potential, tolerances, wall time, n_atoms...) |
| `delta_e_AlMg_200A_n500_tight-<jobid>.{out,err}` | SLURM logs |

`.npz` 结构:

```python
{
  "gb_site_ids":      int64[500],   # 选中的 500 个 GB site 的 atom_id (1-based, LAMMPS conv)
  "gb_e_mg":          float64[500], # 每个 site 装 Mg 后的总 PE (eV)
  "gb_delta_e":       float64[500], # 关键! ΔE_i in eV = gb_e_mg - bulk_e_mean
  "gb_stop_reasons":  <U24[500],    # CG 收敛原因字符串 (e.g. "energy tolerance")
  "bulk_ref_ids":     int64[10],    # 10 个 bulk reference atom_ids
  "bulk_e_mg":        float64[10],  # 10 个 bulk Mg 后 PE
  "bulk_stop_reasons": <U24[10],
  "bulk_e_mean":      scalar,       # <E^bulk-Mg>,Stage 5 用
  "bulk_e_std":       scalar        # σ_bulk(我们 Methods 报告的 8.1 meV)
}
```

**单位换算**:`gb_delta_e * 96.485` → kJ/mol(per Eq. `unit_kjmol` in Methods)。

---

## Stage 5: skew-normal fit + FD prediction

**两个脚本,login node,秒级**:

### 5a. `scripts/fit_delta_e_spectrum.py`

```bash
python scripts/fit_delta_e_spectrum.py \
  --npz /cluster/scratch/cainiu/production_AlMg_200A/delta_e_results_n500_200A_tight.npz \
  --out-png project/output/delta_e_spectrum_n500_200A.png \
  --out-json project/output/delta_e_fit_n500_200A.json
```

**输入**:Stage 4 的 npz
**输出**(都到 `project/output/`):

`delta_e_fit_n500_200A.json`:

```json
{
  "sample_moments_kjmol": {"n": 500, "min", "max", "mean", "median", "std", "skew"},
  "skewnorm_fit_kjmol":   {"alpha", "mu", "sigma"},      // ← 喂给 FD
  "wagih_alm_reference":  {"mu": 9.0, "sigma": 23.0, "alpha": -2.3, ...},
  "source_npz":           "..."
}
```

`delta_e_spectrum_n500_200A.png`:直方图 + skew-normal 拟合 + Wagih dashed 参考。

### 5b. `scripts/fermi_dirac_predict.py`

```bash
python scripts/fermi_dirac_predict.py \
  --ours-npz .../delta_e_results_n500_200A_tight.npz \
  --T-list 300,500,700,900 \
  --out-png project/output/fd_curves_200A_tight.png \
  --out-json project/output/fd_curves_200A_tight.json
```

**输入**:同 npz + T 列表 + X_c grid 参数(默认 1e-5 到 0.5,120 点)

**输出**:`fd_curves_*.json` 含 `X_c` grid + 每个 T 的 `X_GB^FD(X_c)` curve。
这是后面 panel(d) 画 hypothesis 线的数据源。

---

## Stage 6: HMC,一个 (T, X_c) 一个 job

**Submit script** 例:`AlMg/data/decks/submit_hmc_T500_Xc0.10_fdseed.sh`(每个 X_c/T 一个 .sh)
**LAMMPS deck**:`AlMg/data/decks/hmc_AlMg_v2.lammps`(所有 X_c/T 共用)

每个 SLURM 作业 16 ranks × 24 h,实际跑满或 plateau。

**LAMMPS deck 内部**(`hmc_AlMg_v2.lammps`,关键行):

```
# Inputs (via -var):  annealed, outstub, T, XC, EXTRA_PS, RESTART_PS, ...
read_data       ${annealed}            # fdseed 改造后的 .lmp:X_GB(0) = X_FD
velocity        all create T ...
fix nvt         all nvt temp T T 0.1
fix hmc         all atom/swap 100 100 SEED T ke yes types 1 2
                                ^      ^
                                每 100 timestep 试 100 次 Mg↔Al swap
thermo_style    custom step time temp pe ke press f_hmc[1] f_hmc[2]
dump d1         all custom DUMP_EVERY ${outstub}.dump id type
restart         N_RESTART ${outstub}.rst1 ${outstub}.rst2
run             N_PROD                  # PROD_PS=300 → N=300,000 steps @ Δt=1fs
write_data      ${outstub}_final.lmp
```

**关键预处理**(submit script 内,fdseed IC):
`scripts/build_fdseed_inits.py` 把 annealed config 改成
`poly_AlMg_200A_fdseed_T*K_Xc*.lmp` —— X_GB(0) ≈ X_GB^FD(T, X_c)。

**输入**:
- `poly_AlMg_200A_fdseed_T500K_Xc0.1.lmp`(fdseed IC,Stage 6 起点)
- `Al-Mg.eam.fs` 势函数

**输出**(全部到 `/cluster/scratch/cainiu/hmc_AlMg/`):

| 文件 | 内容 | 大小量级 |
|---|---|---|
| `hmc_T500_Xc0.10_fdseed.log` | thermo log:每 ~100 步一行 step/time/T/PE/KE/Pressure/f_hmc[1]/f_hmc[2] | KB |
| `hmc_T500_Xc0.10_fdseed.dump` | per-frame `id type` 转储,每 1 ps 一帧,300 帧 | ~1-2 GB |
| `hmc_T500_Xc0.10_fdseed.rst1, rst2` | LAMMPS restart(轮换) | ~40 MB 每个 |
| `hmc_T500_Xc0.10_fdseed_final.lmp` | 最后 step 的 atomic structure(post-CG 写出) | ~60 MB |
| `hmc_AlMg_T500_Xc0.10_fdseed-<jobid>.{out,err}` | SLURM logs | KB |

`.dump` 是 X_GB(t) 时间序列的**唯一**信息源(每帧只存 id+type,无 positions)。

`f_hmc[1]` = swap attempts 累计;`f_hmc[2]` = accepts 累计。
**swap acceptance = f_hmc[2]/f_hmc[1]**,健康范围 5-30%。

---

## Stage 7: HMC post-process

**脚本**:`scripts/hmc_xgb_timeseries.py`(login node,秒级,**已在 submit script 末尾自动调用**)

```bash
python scripts/hmc_xgb_timeseries.py \
  --stub /cluster/scratch/cainiu/hmc_AlMg/hmc_T500_Xc0.10_fdseed \
  --gb-mask project/AlMg/data/gb_mask_200A.npy \
  --xc 0.10 --temp 500 --fd-pred 0.3519 \
  --out-prefix project/output/hmc_T500_Xc0.10_fdseed
```

**做什么**:
1. 读 `.dump`(每帧 `id type`),计算 X_GB(t) = (1/N_GB)·#{i∈GB | type[i]=2}
2. 读 `.log`,提取 PE(t)、swap accept rate
3. 丢前 20% 帧 burn-in
4. moving block bootstrap 算 ⟨X_GB⟩ + 95% CI(L=5, B=2000)
5. 生成 4-panel 诊断图

**输入**:`.dump` + `.log` + `gb_mask.npy`

**输出**(`project/output/`):

`hmc_T500_Xc0.10_fdseed.json` 关键 keys:

```json
{
  "stub": "...", "T": 500.0, "X_c": 0.10,
  "n_atoms": 475715, "n_gb_atoms": 89042, "gb_fraction": 0.187,
  "swap": {"total_attempts", "total_accepts", "accept_rate_overall"},
  "thermo_prod": {"T_mean", "T_std", "PE_initial", "PE_final", "PE_drift_eV"},
  "x_gb": {                                  // ← X_GB 主结果
    "n_frames_total": 300,
    "n_frames_burnin_dropped": 60,
    "block_bootstrap_block": 5,
    "mean": 0.2785,                          // ⟨X_GB⟩
    "ci95_lo": 0.27,
    "ci95_hi": 0.28,
    "x_total_mean": 0.10,                    // X_c 校验
    "enrichment_vs_bulk": 2.78               // X_GB / X_c
  },
  "fd_predicted": 0.3519,                    // X_GB^FD(T, X_c)
  "gap_hmc_minus_fd": -0.073,                // ← 主信号
  "swap_decomposition": {...},
  "series": {                                // X_GB(t) 时间序列(画图用)
    "timestep", "x_gb", "x_total",
    "fwd_gb_per_frame", "rev_gb_per_frame"
  }
}
```

`hmc_T500_Xc0.10_fdseed.png` 4 panel:X_GB(t), PE(t), swap accept rate,
ΔE histogram。

---

## Stage 8: panel (d) headline aggregation

**脚本**:`scripts/report_headline_T500_9pt.py`(login node,秒级)

```bash
python scripts/report_headline_T500_9pt.py
```

**输入**:**9 个 Stage 7 的 JSON**(每个 X_c 一个)+ FD curves JSON + gb_mask
+ spectrum npz。脚本里硬编码了文件列表(`hmc_T500_Xc0.04_fdseed_resume.json`,
`hmc_T500_Xc0.05_fdseed.json`, ..., `hmc_T500_Xc0.40_fdseed.json`)。

**做什么**:
- 收每个点的 ⟨X_GB⟩, CI95, swap imbalance, drift
- 按 imbalance/drift 阈值标 ● plateau 或 ▽ UB
- 画 X_c → X_GB,叠 FD canon curve + ceiling line + diagonal

**输出**:
- `project/output/00_headline_hmc_vs_wagih_T500_9pt_2026-05-11.json`(每点
  source / X_HMC / X_FD / gap…)
- `project/report/figures/00_headline_hmc_vs_wagih_T500_9pt_2026-05-11.png`
  (panel-d 图)

---

## 数据流总结表

| 数据类型 | 文件扩展 | 位置 | 大小 | git 状态 |
|---|---|---|---|---|
| LAMMPS structure | `.lmp` | scratch | 60-500 MB | gitignored |
| LAMMPS dump 轨迹 | `.dump` | scratch | 100 MB - 2 GB | gitignored |
| LAMMPS restart | `.rst1`/`.rst2` | scratch | 40 MB | gitignored |
| LAMMPS log | `.log` | scratch | KB | gitignored |
| GB mask | `.npy` | scratch + `project/AlMg/data/` | KB | tracked(snapshot copy) |
| Per-site 能量 | `.npz` | scratch | KB | gitignored |
| Summary stats | `.json` | scratch (meta) + `project/output/` (report) | KB | tracked |
| Figure | `.png` | `project/output/` 或 `project/report/figures/` | KB | tracked |
| Script | `.py` | `project/scripts/` | KB | tracked |
| LAMMPS deck | `.lammps` | `project/AlMg/data/decks/` | KB | tracked |
| Submit script | `.sh` | `project/AlMg/data/decks/` | KB | tracked |

---

## 最常用的 read-out 模式

**"我想看 X_c=0.10 T=500 K 的 X_GB^HMC 是多少?"**

```python
import json
d = json.load(open('project/output/hmc_T500_Xc0.10_fdseed.json'))
print(d['x_gb']['mean'], d['x_gb']['ci95_lo'], d['x_gb']['ci95_hi'])
```

**"我想看 ΔE_seg 分布的 (μ, σ, α)?"**

```python
import json
d = json.load(open('project/output/delta_e_fit_n500_200A.json'))
print(d['skewnorm_fit_kjmol'])
```

**"我想看 X_GB(t) 时间序列?"**

```python
d = json.load(open('project/output/hmc_T500_Xc0.10_fdseed.json'))
ts = d['series']['x_gb']  # length 300, X_GB(t) at 1 ps spacing
```

**"重画 panel (d) 的数据从哪儿?"**

`project/output/00_headline_hmc_vs_wagih_T500_9pt_2026-05-11.json`,里面
`equilibrated[]` 和 `upper_bound[]` 两个 list,每项有 `Xc`, `X_HMC`,
`X_FD_canon_ours`, `gap_HMC_minus_canon_ours`。

---

## Pt(Au) 平行管线的差异(2026-05-13 之前)

Pt(Au) 走的是同样 8 个 stage,但参数从 `report/methods_formula_patches.tex`
的 Table 1 substituted。当前 Pt(Au) 只跑到 Stage 5(spectrum + Wagih 比对
+ bootstrap CI 通过 fit-参数判据),Stage 6 HMC 在同事手上,还没开始。

具体差异:

- 输出根目录 `/cluster/scratch/cainiu/prototype_PtAu_100A/`
- 100³ Å 而非 200³ Å,N_tot = 62{,}096
- EAM `eam/alloy`(O'Brien 2017)而非 `eam/fs`(Mendelev 2009)
- T_hold = 816 K(0.4 T_melt(Pt)=2041 K)
- Stage 4 single-site CG ~ 5× 慢(Pt 刚度大)
- 详细 see `project/PtAu/CHANGELOG.md`
