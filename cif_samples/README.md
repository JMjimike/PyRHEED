# REBCO (LnBa₂Cu₃O₇) CIF Library

**Language:** [English](#english) | [中文](#中文)

<a id="english"></a>

## English

Structure files for **123-phase** high-Tc cuprates (REBCO / Ln-123) for PyRHEED diffraction simulation.

### Directory layout

| Subdirectory | Content | Count |
|--------------|---------|-------|
| `orthorhombic_O7/` | Orthorhombic Ba₂LnCu₃O₇ (Materials Project DFT) | 12 |
| `tetragonal_O6/` | Tetragonal Ba₂Ln(CuO₂)₃ / O₆ (oxygen-deficient 123) | 5 |
| `experimental_COD/` | Experimental YBCO structures (Crystallography Open Database) | 2 |
| `dybco_growth_intermediates/` | DyBCO growth intermediates (143, 385/358, 211, 124, BaO, Dy₂O₃, …) | 16 |
| `substrates/` | Common substrates (SrTiO₃, etc.) | 3 |

### Orthorhombic O₇ (`orthorhombic_O7`)

Space group **Pmmm**, isostructural with classic YBCO. Source: [Materials Project](https://materialsproject.org).

| File | Formula | mp-id |
|------|---------|-------|
| YBa2Cu3O7_mp-20674.cif | YBCO | mp-20674 |
| LaBa2Cu3O7_mp-622210.cif | La-123 | mp-622210 |
| PrBa2Cu3O7_mp-20936.cif | Pr-123 | mp-20936 |
| NdBa2Cu3O7_mp-22719.cif | Nd-123 | mp-22719 |
| SmBa2Cu3O7_mp-21451.cif | Sm-123 | mp-21451 |
| EuBa2Cu3O7_mp-622211.cif | Eu-123 | mp-622211 |
| GdBa2Cu3O7_mp-19813.cif | Gd-123 | mp-19813 |
| DyBa2Cu3O7_mp-622105.cif | DyBCO | mp-622105 |
| HoBa2Cu3O7_mp-6616.cif | Ho-123 | mp-6616 |
| ErBa2Cu3O7_mp-622110.cif | Er-123 | mp-622110 |
| TmBa2Cu3O7_mp-622108.cif | Tm-123 | mp-622108 |
| LuBa2Cu3O7_mp-20324.cif | Lu-123 | mp-20324 |

**Ln-123 not in Materials Project:** Ce, Tb, Yb (Ce/Tb 123 unstable or non-superconducting; Yb ion too small).

### Tetragonal O₆ (`tetragonal_O6`)

Space group **P4/mmm**, low-oxygen 123 variants.

| File | Formula | mp-id |
|------|---------|-------|
| YBa2Cu3O7_mp-22215.cif | Ba₂Y(CuO₂)₃ | mp-22215 |
| NdBa2Cu3O7_mp-614981.cif | Ba₂Nd(CuO₂)₃ | mp-614981 |
| SmBa2Cu3O7_mp-622576.cif | Ba₂Sm(CuO₂)₃ | mp-622576 |
| HoBa2Cu3O7_mp-616166.cif | Ba₂Ho(CuO₂)₃ | mp-616166 |
| LuBa2Cu3O7_mp-21868.cif | Ba₂Lu(CuO₂)₃ | mp-21868 |

### Experimental structures (`experimental_COD`)

| File | Notes | Source |
|------|-------|--------|
| YBa2Cu3O7_COD1001453.cif | Superconducting phase at 100 K, neutron powder diffraction | Capponi et al., Europhys. Lett. (1987) |
| YBa2Cu3O6.9_COD1000030.cif | Joint X-ray/neutron refinement | Williams et al., Phys. Rev. B (1988) |

### DyBCO growth intermediates

Common secondary/intermediate phases during DyBa₂Cu₃O₇ growth (PLD / MTG, etc.). See [`dybco_growth_intermediates/README.md`](dybco_growth_intermediates/README.md).

### Use in PyRHEED

```
Simulation → Diffraction Pattern → Add CIF
```

Prefer `orthorhombic_O7/DyBa2Cu3O7_mp-622105.cif` or `experimental_COD/YBa2Cu3O7_COD1001453.cif`. For growth by-phases, use files under `dybco_growth_intermediates/`.

### Citations

- Materials Project: A. Jain et al., *APL Materials* **1**, 011002 (2013)
- COD entries: see `_journal_*` fields inside each CIF

---

<a id="中文"></a>

## 中文

本目录包含 **123 相** 高温超导铜氧化物（REBCO / Ln-123）的结构文件，供 PyRHEED 衍射模拟使用。

### 目录结构

| 子目录 | 内容 | 数量 |
|--------|------|------|
| `orthorhombic_O7/` | 正交相 Ba₂LnCu₃O₇（Materials Project DFT 结构） | 12 |
| `tetragonal_O6/` | 四方相 Ba₂Ln(CuO₂)₃ / O₆（氧欠缺 123 相） | 5 |
| `experimental_COD/` | 实验测定 YBCO 结构（Crystallography Open Database） | 2 |
| `dybco_growth_intermediates/` | DyBCO 生长中间相（143、385/358、211、124、BaO、Dy₂O₃ 等） | 16 |
| `substrates/` | 常用衬底（SrTiO₃ 等） | 3 |

### 正交相 O₇（orthorhombic_O7）

空间群 **Pmmm**，与经典 YBCO 同型。来源：[Materials Project](https://materialsproject.org)。

| 文件 | 化学式 | mp-id |
|------|--------|-------|
| YBa2Cu3O7_mp-20674.cif | YBCO | mp-20674 |
| LaBa2Cu3O7_mp-622210.cif | La-123 | mp-622210 |
| PrBa2Cu3O7_mp-20936.cif | Pr-123 | mp-20936 |
| NdBa2Cu3O7_mp-22719.cif | Nd-123 | mp-22719 |
| SmBa2Cu3O7_mp-21451.cif | Sm-123 | mp-21451 |
| EuBa2Cu3O7_mp-622211.cif | Eu-123 | mp-622211 |
| GdBa2Cu3O7_mp-19813.cif | Gd-123 | mp-19813 |
| DyBa2Cu3O7_mp-622105.cif | DyBCO | mp-622105 |
| HoBa2Cu3O7_mp-6616.cif | Ho-123 | mp-6616 |
| ErBa2Cu3O7_mp-622110.cif | Er-123 | mp-622110 |
| TmBa2Cu3O7_mp-622108.cif | Tm-123 | mp-622108 |
| LuBa2Cu3O7_mp-20324.cif | Lu-123 | mp-20324 |

**Materials Project 中未收录的 Ln-123：** Ce、Tb、Yb（Ce/Tb 的 123 相不稳定或非超导；Yb 离子过小）。

### 四方相 O₆（tetragonal_O6）

空间群 **P4/mmm**，氧含量较低时的 123 相变体。

| 文件 | 化学式 | mp-id |
|------|--------|-------|
| YBa2Cu3O7_mp-22215.cif | Ba₂Y(CuO₂)₃ | mp-22215 |
| NdBa2Cu3O7_mp-614981.cif | Ba₂Nd(CuO₂)₃ | mp-614981 |
| SmBa2Cu3O7_mp-622576.cif | Ba₂Sm(CuO₂)₃ | mp-622576 |
| HoBa2Cu3O7_mp-616166.cif | Ba₂Ho(CuO₂)₃ | mp-616166 |
| LuBa2Cu3O7_mp-21868.cif | Ba₂Lu(CuO₂)₃ | mp-21868 |

### 实验结构（experimental_COD）

| 文件 | 说明 | 来源 |
|------|------|------|
| YBa2Cu3O7_COD1001453.cif | 100 K 超导相，中子粉末衍射 | Capponi et al., Europhys. Lett. (1987) |
| YBa2Cu3O6.9_COD1000030.cif | 联合 X 射线/中子精修 | Williams et al., Phys. Rev. B (1988) |

### DyBCO 生长中间相（dybco_growth_intermediates）

详见 [`dybco_growth_intermediates/README.md`](dybco_growth_intermediates/README.md)。

### 在 PyRHEED 中使用

```
Simulation → Diffraction Pattern → Add CIF
```

推荐优先使用 `orthorhombic_O7/DyBa2Cu3O7_mp-622105.cif` 或 `experimental_COD/YBa2Cu3O7_COD1001453.cif`。

### 引用

- Materials Project: A. Jain et al., *APL Materials* **1**, 011002 (2013)
- COD entries: 见各 CIF 文件内 `_journal_*` 字段
