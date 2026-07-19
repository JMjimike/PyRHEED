# DyBa₂Cu₃O₇ Growth Intermediate Phase CIF Library

**Language:** [English](#english) | [中文](#中文)

<a id="english"></a>

## English

CIFs for **secondary / intermediate phases** commonly seen during **DyBCO (Dy-123)** growth by PLD, MOCVD, or melt-textured growth (MTG), for PyRHEED diffraction simulation.

### Simplified growth phase diagram

```
Precursor powders
  BaO  +  Dy₂O₃  +  CuO  →  mix / pre-react
         ↓
High T / oxygen-poor              Low T / oxygen-rich
────────────────────────────────────────────────────
  BaCuO₂ (liq/solid)  +  CuO
         ↓  peritectic / reaction
  Dy₂BaCuO₅ (211 green phase)  ← solid inclusion, common in MTG
         ↓
  DyBa₂Cu₃O₆ (tetragonal O₆)  ← before OHT; no Dy entry in MP/COD*
         ↓  oxygenation + OHT (~650–700 °C)
  DyBa₂Cu₃O₇ (orthorhombic O₇)  ← final SC phase → ../orthorhombic_O7/DyBa2Cu3O7_mp-622105.cif
         ‖ coexist
  DyBa₂Cu₄O₈ (124)      ← high-P or specific thermal history
         ‖ coexist
  REBa₄Cu₃O₈.₅ (143)    ← non-SC cubic; Dy143 / Y143, …
  RE₃Ba₅Cu₈O₁₈ (385/358)← high-Tc family (5×CuO₂ planes + 3×CuO chains)
```

\* Materials Project and COD do not list **Ba₂Dy(CuO₂)₃** (tetragonal O₆). This folder provides an **isostructural Y-123 tetragonal** CIF as an O₆ reference; tune lattice constants for the Dy ionic radius if needed.

### Phase naming

| Common name | Formula | Notes |
|-------------|---------|-------|
| **DBCO143** / RE143 | REBa₄Cu₃O₈.₅₊δ | 1:4:3 cation ratio; **non-superconducting** (Zhu et al. 1999) |
| **DBCO385** / RE358 | RE₃Ba₅Cu₈O₁₈₋δ | 3:5:8 ratio; also written **358**; orthorhombic Pmm2, Tc ~80–102 K |
| DBCO | DyBa₂Cu₃O₇₋δ | Dy-123 (final O₇ phase) |

RE may be Y, Gd, Sm, Ho, etc. This folder provides **Dy / Y / Gd** representative structures.

### 143 phase (REBa₄Cu₃O₈.₅)

| File | Phase | Space group | Source |
|------|-------|-------------|--------|
| `DyBa4Cu3O8.5_143_COD1522334.cif` | Dy143 (DBCO143) | P23 | [COD #1522334](https://www.crystallography.net/1522334.html) |
| `YBa4Cu3O8.5_143_COD1539855.cif` | Y143 | Pm-3 | [COD #1539855](https://www.crystallography.net/1539855.html) |

Dy143 and Y143 have **different space groups** (P23 vs Pm-3); Gd143 is Pm-3 but not in COD.

### 385 phase (RE₃Ba₅Cu₈O₁₈ / RE358)

COD/MP have **no** complete experimental CIF. Files below are generated from Landínez-Téllez et al. (*Physica C* **590**, 2021) **Table 2** Rietveld model; lattice constants from **Table 3** (and KIT compilations):

| File | Phase | a, b, c (Å) | Space group |
|------|-------|-------------|-------------|
| `Dy3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` | Dy385 (DBCO385) | 3.833, 3.867, 31.025 | Pmm2 |
| `Y3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` | Y358 | 3.888, 3.823, 31.013 | Pmm2 |
| `Gd3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` | Gd358 | 3.851, 3.877, 31.120 | Pmm2 |

385 models use **simplified atomic coordinates** (34 atoms/cell), suitable for qualitative RHEED; see the paper for refined structures.

### Other intermediate-phase files

| File | Phase | Space group | Role in growth | Source |
|------|-------|-------------|----------------|--------|
| `Dy2BaCuO5_211_COD1526528.cif` | Dy₂BaCuO₅ (211) | Pnma | **Green phase**; peritectic solid, MTG inclusion | [COD #1526528](https://www.crystallography.net/1526528.html) |
| `Dy2BaCuO5_211_mp22550.cif` | BaDy₂CuO₅ (211) | Pnma | Same (DFT, isostructural with COD) | [MP mp-22550](https://legacy.materialsproject.org/materials/mp-22550/) |
| `DyBa2Cu4O8_124_mp6691.cif` | Ba₂Dy(CuO₂)₄ (124) | Cmmm | Double CuO-chain phase; can coexist with 123 | [MP mp-6691](https://legacy.materialsproject.org/materials/mp-6691/) |
| `BaCuO2_tetragonal_mp752398.cif` | BaCuO₂ | P4/mmm | Ba–Cu–O-rich liquid/solid by-phase | [MP mp-752398](https://legacy.materialsproject.org/materials/mp-752398/) |
| `BaCuO2_orthorhombic_mp997034.cif` | BaCuO₂ | Cmcm | Another stable polymorph | [MP mp-997034](https://legacy.materialsproject.org/materials/mp-997034/) |
| `CuO_mp14549.cif` | CuO | C2/c | Binary oxide by-phase / feedstock | [MP mp-14549](https://legacy.materialsproject.org/materials/mp-14549/) |
| `BaO_mp1342.cif` | BaO | Fm-3m | **Barium oxide**; cubic rock-salt (MP stable) | [MP mp-1342](https://legacy.materialsproject.org/materials/mp-1342/) |
| `BaO_COD1527735.cif` | BaO | P4/nmm | Dense tetragonal polymorph (experimental) | [COD #1527735](https://www.crystallography.net/1527735.html) |
| `Dy2O3_mp2345.cif` | Dy₂O₃ | Ia-3 | **Dysprosium oxide**; C-type cubic sesquioxide | [MP mp-2345](https://legacy.materialsproject.org/materials/mp-2345/) |
| `Dy2O3_COD1539592.cif` | Dy₂O₃ | Ia-3 | Same (neutron/powder experimental) | [COD #1539592](https://www.crystallography.net/1539592.html) |
| `YBa2Cu3O7_tetragonal_O6_reference_mp22215.cif` | Ba₂Y(CuO₂)₃ | P4/mmm | **O₆ tetragonal 123 structure reference** (not Dy) | [MP mp-22215](https://legacy.materialsproject.org/materials/mp-22215/) |

**Final product (O₇):** `../orthorhombic_O7/DyBa2Cu3O7_mp-622105.cif`

### Use in PyRHEED

```
Simulation → Diffraction Pattern → Add CIF
```

Suggested combinations:

| Scenario | Suggested CIFs |
|----------|----------------|
| Unreacted precursors | `BaO_mp1342.cif`, `Dy2O3_mp2345.cif`, `CuO_mp14549.cif` |
| MTG green-phase inclusions | `Dy2BaCuO5_211_COD1526528.cif` |
| 124 coexist with 123 | `DyBa2Cu4O8_124_mp6691.cif` + O₇ |
| BaCuO₂ / CuO by-phases | `BaCuO2_*.cif`, `CuO_mp14549.cif` |
| 143 / 385 coexist | `DyBa4Cu3O8.5_143_COD1522334.cif`, `Dy3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` |
| Y / Gd analogues | `YBa4Cu3O8.5_143_*`, `Y3Ba5Cu8O18_385_*`, `Gd3Ba5Cu8O18_385_*` |
| High-T tetragonal O₆ (structure ref.) | `YBa2Cu3O7_tetragonal_O6_reference_mp22215.cif` |

Prefer **COD experimental** 211 when available; use `mp22550` when aligning with the MP DFT series.

### Citations

- Materials Project: A. Jain et al., *APL Materials* **1**, 011002 (2013)
- COD #1526528: Kan et al., *J. Eur. Ceram. Soc.* **21**, 2593 (2001)
- 211 review: Michel et al., *J. Solid State Chem.* **68**, 271 (1987)
- 143: Zhu et al., *J. Mater. Res.* **14**, 334 (1999); COD #1522334, #1539855
- 385/358: Landínez-Téllez et al., *Physica C* **590**, 2021; Aliabadi et al., *Physica C* **469**, 2012

---

<a id="中文"></a>

## 中文

本目录收录 **DyBCO（Dy-123）** 在 PLD、MOCVD 或熔融织构（MTG）生长过程中常见的 **伴生相 / 中间相** 结构文件，供 PyRHEED 衍射模拟使用。

### 生长相图（简化）

```
原料粉末
  BaO  +  Dy₂O₃  +  CuO  →  混合 / 预反应
         ↓
高温 / 贫氧                    低温 / 富氧
─────────────────────────────────────────────
  BaCuO₂ (液/固)  +  CuO
         ↓  peritectic / 反应
  Dy₂BaCuO₅ (211 绿相)  ← 固相夹杂物，MTG 中常见
         ↓
  DyBa₂Cu₃O₆ (四方 O₆)  ← OHT 转变前；MP/COD 无 Dy 条目*
         ↓  吸氧 + OHT (~650–700 °C)
  DyBa₂Cu₃O₇ (正交 O₇)  ← 最终超导相 → ../orthorhombic_O7/DyBa2Cu3O7_mp-622105.cif
         ‖ 共存
  DyBa₂Cu₄O₈ (124)      ← 高压或特定热历史下出现
         ‖ 共存
  REBa₄Cu₃O₈.₅ (143)    ← 非超导立方相；Dy143 / Y143 等
  RE₃Ba₅Cu₈O₁₈ (385/358)← 高 Tc 家族相（5×CuO₂ 面 + 3×CuO 链）
```

\* Materials Project 与 COD 均未收录 **Ba₂Dy(CuO₂)₃**（四方 O₆）。目录中提供 **同型 Y-123 四方相** 作为 O₆ 结构参考；晶格常数需按 Dy 离子半径微调。

### 相命名说明

| 俗称 | 化学式 | 备注 |
|------|--------|------|
| **DBCO143** / RE143 | REBa₄Cu₃O₈.₅₊δ | 1:4:3 阳离子比；**非超导**（Zhu et al. 1999） |
| **DBCO385** / RE358 | RE₃Ba₅Cu₈O₁₈₋δ | 3:5:8 比；文献亦写 **358**；正交 Pmm2，Tc 约 80–102 K |
| DBCO | DyBa₂Cu₃O₇₋δ | 即 Dy-123（O₇ 终相） |

RE 可替换为 Y、Gd、Sm、Ho 等；本目录提供 **Dy / Y / Gd** 代表结构。

### 143 相 / 385 相 / 其他中间相

表格与文件列表见上方英文部分（文件名与来源相同）。

### 在 PyRHEED 中使用

```
Simulation → Diffraction Pattern → Add CIF
```

推荐组合模拟场景见英文表。211 相优先使用 **COD 实验结构**。

### 引用

- Materials Project: A. Jain et al., *APL Materials* **1**, 011002 (2013)
- COD #1526528: Kan et al., *J. Eur. Ceram. Soc.* **21**, 2593 (2001)
- 211 相综述: Michel et al., *J. Solid State Chem.* **68**, 271 (1987)
- 143 相: Zhu et al., *J. Mater. Res.* **14**, 334 (1999); COD #1522334, #1539855
- 385/358 相: Landínez-Téllez et al., *Physica C* **590**, 2021; Aliabadi et al., *Physica C* **469**, 2012
