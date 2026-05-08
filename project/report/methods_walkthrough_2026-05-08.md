# Methods 草稿 walk-through — team meeting 2026-05-08

会议用 1 页提纲。源文件:`report/methods_draft.tex` (375 行)、
`report/methods_references_to_verify.md` (169 行)。Methods 部分和附录均
是独立 standalone 文件,**还没**粘进 Overleaf `main.tex` —— 等会上对齐
完再操作,避免和 Cu(Ni) 合作者改动冲突。

---

## 1. 章节地图 (methods_draft.tex 里有什么)

| 段落 | 标签 | 内容 |
|---|---|---|
| §Methodology 引子 | `sec:meth` | 一段话讲清"hypothesis vs ground-truth"框架:FD 给预测、HMC 给 ground truth、gap 为可测量 |
| §Atomistic base + ΔE 谱 | `sec:meth:base` | 200³ Å³ 16-grain Voronoi 多晶 / Mendelev EAM-FS / a-CNA 鉴别 GB / N_tot=475715, N_GB=89042, f_GB=0.187 / n=500 ΔE_i 子样本 |
| §FD 预测 | `sec:meth:FD` | GC-FD vs canon-FD 两个变体;canon-FD 是闭盒 N_Mg 守恒下的可比量 |
| §HMC 采样 | `sec:meth:HMC` | LAMMPS `fix atom/swap` + NVT;**fdseed IC 是关键方法学选择**;production 10 ps eq + 300 ps prod, 16 MPI ranks |
| §$(T, X_c)$ 网格 + 比较 | `sec:meth:grid` | T ∈ {300, 500, 700} K,X_c ∈ [0.04, 0.30];稀薄端 0.04 由 ceiling compression 限定 |
| §统计分析 | `sec:meth:stats` | 两套 CI:trajectory block-bootstrap (5-frame block, 2000 resamples) + per-site Wald binomial CI |
| §方法图 | `fig:method` | 双副面板:左 TikZ 流程图(yellow callout 突出 fdseed)、右 OVITO 多晶渲染 (`polycrystal_geometry.png`) |
| App. A 原子系统 | `app:atomistic` | 完整退火协议:CG → NVT 100→373 K (50 ps) → NPT 373 K (250 ps) → 3 K/ps quench → CG + box relax |
| App. B 关键公式 | `app:eq` | Eq. (1) ΔE_i 定义;Eq. (2) Wagih 单点 FD;Eq. (3) canonical FD self-consistent;Eq. (4) Wald binomial CI |
| App. C HMC 算法 | `app:alg` | 完整伪代码 box(`algorithm` 包) |
| App. D ΔE 谱收敛 + Wagih 比较 | `app:dE_compare` | n=500 选取理由 + 与 Wagih Zenodo n=82646 KS test (p > 0.5) |

---

## 2. 会议决策钩子(方法草稿头部已标注)

**A. Cu(Ni) 合作者的原子学协议是否与我们一致?**
- 如果一致 → 把 `sec:meth:base` 合并成共享一段,Cu(Ni) 与 Al(Mg) 共用退火协议、GB 鉴别、ΔE 谱采样描述
- 如果不一致 → 保持 base 段独立 (各自一小段),`sec:meth:FD` 和 `sec:meth:HMC` 仍是共享的(框架与采样代码完全相同)
- 需要 Cu(Ni) 同事现场对齐:退火温度斜率 / 250 ps NPT hold / 3 K/ps quench / 是否一样的 a-CNA 参数

**B. 选点优化合作者(refining Wagih's training point selection)的成果**
- 我的建议:放在新 appendix(不进 main Methods),因为他们的工作"细化谱输入"而非"测试 FD 假设" — 与项目主线交叉但不在主测试链路上
- 现场要确认:成果成熟度、是否已经有可放图、是否需要 Methods 引子里加一句话指向附录

**C. HMC 命名方案(已草稿落地了 (ii))**
- 目前选项:首次出现展开成 "hybrid Monte Carlo / molecular dynamics",脚注里说明与 Hamiltonian Monte Carlo (Duane et al. 1987) 区分
- 如果合作者更倾向其它方案(比如 "MC-MD" 简写、或干脆叫 "swap MC"),可以现场敲定再统一改
- 影响范围:methods_draft.tex 里 ~6 处 "HMC" 字眼

**D. 方法图布局(已草稿落地双副面板)**
- 左:TikZ 流程图(`base → FD predictor / HMC → grid comparison`,fdseed 黄色 callout)
- 右:OVITO 多晶渲染(已生成,**今天换为 steel-blue GB**;之前是暖红被反馈太刺眼)
- 现场可确认:双副面板宽度 (0.36\textwidth) 是否合适、左右对调还是当前左流程右几何更好

**E. 报告主线"broken X_c band"措辞放在哪里**
- 现状:Methods 不放 critical-X_c 措辞 — 那是 Results / Discussion 内容
- Methods 里只描述网格范围 [0.04, 0.30] 与稀薄端边界由 ceiling compression 限定
- 现场可对齐:Results 章节里谁负责写 broken band 这个 take-away

---

## 3. 文献待验证清单(细节见 methods_references_to_verify.md)

**6 个 \\cite{} key**(每个都标了候选 DOI 与置信度;**用户规则:不可以编造文献**):

| key | 文献 | 候选 DOI | 置信度 |
|---|---|---|---|
| `Mendelev2009` | Mendelev–Asta–Rahman–Hoyt, Phil. Mag. 89, 3269 (2009) | `10.1080/14786430903260727` | HIGH |
| `Wagih2020` | Wagih–Larsen–Schuh, Nat. Commun. 11:6376 (2020) | `10.1038/s41467-020-20083-6` | HIGH |
| `LAMMPS` | Thompson et al., Comput. Phys. Commun. 271, 108171 (2022) **OR** Plimpton 1995 | `10.1016/j.cpc.2021.108171` 或 `10.1006/jcph.1995.1039` | HIGH (二选一) |
| `Stukowski2012` | Stukowski, MSMSE 20, 045021 (2012) | `10.1088/0965-0393/20/4/045021` | MEDIUM-HIGH |
| `Sadigh2012` | Sadigh et al., PRB 85, 184203 (2012) | `10.1103/PhysRevB.85.184203` | MEDIUM |
| `VoroPP` | Rycroft, Chaos 19, 041111 (2009) | `10.1063/1.3215722` | MEDIUM (依赖于实际所用 Voronoi 库) |

**可能需要现场分配**:让 Cu(Ni) 同事/选点优化同事各自承担 1-2 篇 verify(打开 DOI 链接确认条目正确)。

**两条房屋样式选择**(讨论一次,统一应用):
- DOI 是否做成 `\href` 可点击链接(主 main.tex 已 `\usepackage{hyperref}`)
- 期刊缩写还是全名(选 Phil. Mag. 还是 Philosophical Magazine,然后全文一致)

---

## 4. Overleaf 操作的前置条件(粘贴前要做的事)

1. **Preamble 加包**(main.tex 当前没有,要先加上):
   ```
   \usepackage{tikz}
   \usetikzlibrary{arrows.meta, positioning, shapes}
   \usepackage{algorithm}
   \usepackage{algpseudocode}
   ```
2. **\\bibitem 块逐条加**到现有 `\\begin{thebibliography}` 段(skeleton 在 references-to-verify.md 末尾)
3. **OVITO PNG 上传**:今天生成的 `report/figures/polycrystal_geometry.png`(steel-blue 版本)上传到 Overleaf,`\\includegraphics{polycrystal_geometry.png}` 自动 resolve
4. **粘贴顺序**:Methodology section → Method figure(放在 Methodology 之后)→ Appendices 全部放在已有 `\\appendix` 行下面;不要动 `app1` `app2` placeholders

---

## 5. 风险 / 仍要小心的点

- **`report/69f6210146...` 是 Overleaf checkout** — `.gitignore` 里(CHANGELOG 2026-05-07 entry),其 `.git` 不要被 main repo 当 submodule 处理;会议后所有改动**走 Overleaf 端**而非这个本地 checkout
- **Methods 草稿提到了 Tachyon 渲染需要 OVITO ≥ 3.0** — 同事如要本地复现 OVITO 渲染,确认 OVITO 版本(我们用的是 3.15.4)
- **Algorithm box 排版**:`algorithm` + `algpseudocode` 在 RevTeX 4-2 的 twocolumn 模式下偶尔会跨栏排版,排版结果若不好可降为单栏 figure*

---

## 6. 我能在会上即时做的事

- 现场打开 `methods_draft.tex` 任何章节做 ad-hoc 编辑,提交到 cainiu 分支
- 现场打开 `polycrystal_geometry.png`、TikZ 流程图初稿,做布局示意
- 关于"Cu(Ni) 是否合并 base 段"的两种结果都准备好了 diff 模板
