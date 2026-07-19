# DyBa₂Cu₃O₇ 生长过程中间相 CIF 库

本目录收录 **DyBCO（Dy-123）** 在 PLD、MOCVD 或熔融织构（MTG）生长过程中常见的 **伴生相 / 中间相** 结构文件，供 PyRHEED 衍射模拟使用。

## 生长相图（简化）

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

## 相命名说明

| 俗称 | 化学式 | 备注 |
|------|--------|------|
| **DBCO143** / RE143 | REBa₄Cu₃O₈.₅₊δ | 1:4:3 阳离子比；**非超导**（Zhu et al. 1999） |
| **DBCO385** / RE358 | RE₃Ba₅Cu₈O₁₈₋δ | 3:5:8 比；文献亦写 **358**；正交 Pmm2，Tc 约 80–102 K |
| DBCO | DyBa₂Cu₃O₇₋δ | 即 Dy-123（O₇ 终相） |

RE 可替换为 Y、Gd、Sm、Ho 等；本目录提供 **Dy / Y / Gd** 代表结构。

## 143 相（REBa₄Cu₃O₈.₅）

| 文件 | 相 | 空间群 | 来源 |
|------|-----|--------|------|
| `DyBa4Cu3O8.5_143_COD1522334.cif` | Dy143 (DBCO143) | P23 | [COD #1522334](https://www.crystallography.net/1522334.html) |
| `YBa4Cu3O8.5_143_COD1539855.cif` | Y143 | Pm-3 | [COD #1539855](https://www.crystallography.net/1539855.html) |

Dy143 与 Y143 **空间群不同**（P23 vs Pm-3）；Gd143 为 Pm-3 但 COD 未收录。

## 385 相（RE₃Ba₅Cu₈O₁₈ / RE358）

COD/MP **无**完整实验 CIF。下列文件由 Landínez-Téllez et al. (*Physica C* **590**, 2021) **Table 2** Rietveld 模型生成，晶格常数取自该文 **Table 3**（及 KIT 汇编）：

| 文件 | 相 | a, b, c (Å) | 空间群 |
|------|-----|-------------|--------|
| `Dy3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` | Dy385 (DBCO385) | 3.833, 3.867, 31.025 | Pmm2 |
| `Y3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` | Y358 | 3.888, 3.823, 31.013 | Pmm2 |
| `Gd3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` | Gd358 | 3.851, 3.877, 31.120 | Pmm2 |

385 模型为 **简化原子坐标**（34 原子/胞），适用于 RHEED 定性模拟；精修结构见原文。

## 其他中间相文件列表

| 文件 | 相 | 空间群 | 生长中的角色 | 来源 |
|------|-----|--------|--------------|------|
| `Dy2BaCuO5_211_COD1526528.cif` | Dy₂BaCuO₅ (211) | Pnma | **绿相**；peritectic 固相、MTG 夹杂物 | [COD #1526528](https://www.crystallography.net/1526528.html) |
| `Dy2BaCuO5_211_mp22550.cif` | BaDy₂CuO₅ (211) | Pnma | 同上（DFT 优化结构，与 COD 同型） | [MP mp-22550](https://legacy.materialsproject.org/materials/mp-22550/) |
| `DyBa2Cu4O8_124_mp6691.cif` | Ba₂Dy(CuO₂)₄ (124) | Cmmm | 双 CuO 链相，可与 123 共存 | [MP mp-6691](https://legacy.materialsproject.org/materials/mp-6691/) |
| `BaCuO2_tetragonal_mp752398.cif` | BaCuO₂ | P4/mmm | Ba–Cu–O 富液相 / 副相 | [MP mp-752398](https://legacy.materialsproject.org/materials/mp-752398/) |
| `BaCuO2_orthorhombic_mp997034.cif` | BaCuO₂ | Cmcm | 另一稳定多型 | [MP mp-997034](https://legacy.materialsproject.org/materials/mp-997034/) |
| `CuO_mp14549.cif` | CuO | C2/c | 二元氧化物副相 / 原料 | [MP mp-14549](https://legacy.materialsproject.org/materials/mp-14549/) |
| `BaO_mp1342.cif` | BaO | Fm-3m | **氧化钡**；立方岩盐相（MP 稳定相） | [MP mp-1342](https://legacy.materialsproject.org/materials/mp-1342/) |
| `BaO_COD1527735.cif` | BaO | P4/nmm | 氧化钡致密四方多型（实验结构） | [COD #1527735](https://www.crystallography.net/1527735.html) |
| `Dy2O3_mp2345.cif` | Dy₂O₃ | Ia-3 | **氧化镝**；C 型立方 sesquioxide | [MP mp-2345](https://legacy.materialsproject.org/materials/mp-2345/) |
| `Dy2O3_COD1539592.cif` | Dy₂O₃ | Ia-3 | 同上（中子/粉末衍射实验结构） | [COD #1539592](https://www.crystallography.net/1539592.html) |
| `YBa2Cu3O7_tetragonal_O6_reference_mp22215.cif` | Ba₂Y(CuO₂)₃ | P4/mmm | **O₆ 四方 123 结构参考**（非 Dy，同型） | [MP mp-22215](https://legacy.materialsproject.org/materials/mp-22215/) |

**最终产物（O₇）** 见：`../orthorhombic_O7/DyBa2Cu3O7_mp-622105.cif`

## 在 PyRHEED 中使用

```
Simulation → Diffraction Pattern → Add CIF
```

推荐组合模拟场景：

| 场景 | 建议 CIF |
|------|----------|
| 原料粉末（未反应） | `BaO_mp1342.cif`、`Dy2O3_mp2345.cif`、`CuO_mp14549.cif` |
| MTG 绿相夹杂 | `Dy2BaCuO5_211_COD1526528.cif` |
| 124 相共存 | `DyBa2Cu4O8_124_mp6691.cif` + O₇ |
| 副相 BaCuO₂ / CuO | 对应 `BaCuO2_*.cif`、`CuO_mp14549.cif` |
| 143 / 385 相共存 | `DyBa4Cu3O8.5_143_COD1522334.cif`、`Dy3Ba5Cu8O18_385_Pmm2_Landinez2021.cif` |
| 换 Y / Gd 模拟 | `YBa4Cu3O8.5_143_*`、`Y3Ba5Cu8O18_385_*`、`Gd3Ba5Cu8O18_385_*` |
| 高温四方 O₆（结构参考） | `YBa2Cu3O7_tetragonal_O6_reference_mp22215.cif` |

211 相优先使用 **COD 实验结构**；需要与 MP 系列 DFT 结构统一时可用 `mp22550` 版本。

## 引用

- Materials Project: A. Jain et al., *APL Materials* **1**, 011002 (2013)
- COD #1526528: Kan et al., *J. Eur. Ceram. Soc.* **21**, 2593 (2001)（Dy₂BaCuO₅ 结构参数）
- 211 相综述: Michel et al., *J. Solid State Chem.* **68**, 271 (1987)
- 143 相: Zhu et al., *J. Mater. Res.* **14**, 334 (1999); COD #1522334, #1539855
- 385/358 相: Landínez-Téllez et al., *Physica C* **590**, 2021; Aliabadi et al., *Physica C* **469**, 2012
