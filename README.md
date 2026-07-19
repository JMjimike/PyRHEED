# <img src="https://github.com/yux1991/PyRHEED/blob/master/src/icons/icon.png" width="48"/> PyRHEED
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/yux1991/PyRHEED/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/yux1991/PyRHEED.svg)](https://github.com/yux1991/PyRHEED/blob/master/LICENSE) [![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](mailto:yux1991@gmail.com)

**Language:** [English](#english) | [中文](#中文)

---

<a id="english"></a>

## English

### Table of Contents
1. [Upstream & Changes in This Fork](#upstream--changes-in-this-fork)
2. [Description](#description)
3. [Requirements](#requirements)
4. [Installation & Usage](#installation--usage)
5. [Modules](#modules)
6. [Demonstrations](#demonstrations)
7. [Tutorial](#tutorial)
8. [Contact](#contact)

### Upstream & Changes in This Fork

This repository is derived from [yux1991/PyRHEED](https://github.com/yux1991/PyRHEED) (MIT License). It keeps Yu Xiang’s core RHEED analysis and simulation features, with the following fixes and extensions:

| Category | Changes |
|----------|---------|
| **Bug fix** | NumPy 2.x compatibility: replace removed `ndarray.ptp()` with `np.ptp()`, fixing TAPD “Add Structure” crashes |
| **Docs** | Chinese README and detailed workflow tutorial [`TUTORIAL_zh.md`](TUTORIAL_zh.md) |
| **Tooling** | Dependency management with [uv](https://docs.astral.sh/uv/) (`pyproject.toml` / `uv.lock`) |
| **CIF library** | [`cif_samples/`](cif_samples/): STO substrates, REBCO O₇/O₆, DyBCO growth intermediates (211/124/143/385, etc.) |
| **XRD extensions** | Four-circle geometry, Ewald sphere, continuous scans (`xrd_*.py`, `ewald.py`, `four_circle_geometry.py`, …) and tests under `scripts/` |

Upstream install: `pip install pyrheed` / `python -m pyrheed`. **This fork recommends:**

```bash
uv sync
uv run src/main.py
```

### Description

PyRHEED is for **reflection high-energy electron diffraction (RHEED)** data analysis and theoretical simulation.

RHEED uses a relatively high-energy (5–30 keV) electron beam at grazing incidence. It is highly surface-sensitive (penetration of only a few nanometers). Because the electron scattering factor is about four orders of magnitude larger than that of X-rays, RHEED is especially suited to 2D materials (e.g. graphene) that are hard to detect with XRD. The large beam footprint (~1 cm) also enables wafer-scale averages of lattice constants, grain orientation, and defect density.

Developed and tested with **Python ≥ 3.12 (64-bit)** and a **PyQt6** GUI. The *simulate_RHEED* module uses [pymatgen](http://pymatgen.org/) to read CIF files and build structures. Main features:

1. Process RHEED images with [rawpy](https://pypi.org/project/rawpy/); extract intensity profiles with [numpy](https://www.numpy.org/) vectorization; build 2D reciprocal-space maps and pole figures; export 3D data as `*.vtp` for [ParaView](https://www.paraview.org).
2. Batch-fit RHEED line profiles with predefined peak shapes (Gaussian, Voigt, …) and save formatted results.
3. Visualize fit results and generate reports.
4. Simulate a statistical factor from a Markov process for hexagonal surfaces with given step-density distributions. See [Spadacini et al.](https://www.sciencedirect.com/science/article/pii/0039602883904922).
5. Read crystal structures from CIF, stack materials into custom heterostructures, and compute diffraction patterns with kinematic theory.
6. Simulate 3D reciprocal space; build atomic models in-app, especially for 2D translational antiphase domain (TAPD) diffraction. See [Lu et al.](https://www.sciencedirect.com/science/article/pii/0039602881905410).

### Requirements

- Python ≥ 3.12
- [uv](https://docs.astral.sh/uv/getting-started/installation/) package manager
- Linux and Windows supported; macOS should work but is less tested

### Installation & Usage

```bash
git clone https://github.com/JMjimike/PyRHEED.git
cd PyRHEED
uv sync
uv run src/main.py
```

1. **Load data** — Open RAW or common image formats via the file dialog. Images are converted to grayscale. Prefer TIFF/PNG/RAW over JPEG for higher dynamic range.
2. **Analysis** — Choose modules by research goal. Pattern / structure-factor simulation does **not** require experimental images.
3. **Run Scenario** — Load `default_scenario.ini` or a custom scenario and run to generate a full result set. See the [Chinese tutorial](TUTORIAL_zh.md) for workflows.

### Modules

| Module | Role |
|--------|------|
| `broadening` | Batch auto-fit of RHEED line profiles |
| `browser` | Browse files in the working directory |
| `canvas` | Display images, draw shapes, take user input |
| `configuration` | Default configuration |
| `cursor` | Cursor-related information |
| `generate_report` | Visual report of fit results |
| `graph_3D_surface` | 3D reciprocal-space surface view |
| `kikuchi` | Kikuchi pattern for non-reconstructed surfaces |
| `main` | Application entry point |
| `manual_fit` | Manual line-profile fit (initial parameters) |
| `my_widgets` | Custom widgets |
| `plot_chart` | QChart widget for fit visualization |
| `preference` | User preferences |
| `process` | Backend processing |
| `process_monitor` | Memory monitoring (incomplete) |
| `profile_chart` | QChart widget for line profiles |
| `properties` | Dynamic program parameters |
| `reciprocal_space_mapping` | 2D/3D reciprocal maps and pole figures |
| `scenario` | Select/customize and run scenarios |
| `simulate_RHEED` | Diffraction from atomic structures; TAPD/APD models |
| `statistical_factor` | Markov-process statistical factor |
| `translational_antiphase_domain` | 1D/2D TAPD analytic profiles |
| `window` | Main window |
| `write_scenario` | Write default scenario config |
| `xrd_*` / `ewald` / `four_circle_geometry` | XRD / geometry extensions (this fork) |

### Demonstrations

1. Extract line profiles:

![line](https://user-images.githubusercontent.com/38077812/111377405-9a688e00-866e-11eb-8ef2-b25386f10d27.gif)

2. Extract features on a sphere:

![tilt](https://user-images.githubusercontent.com/38077812/111377452-aa806d80-866e-11eb-91eb-8a7f103c2077.gif)

3. Construct azimuthal RHEED:

![Azimuth](https://user-images.githubusercontent.com/38077812/111377562-cbe15980-866e-11eb-8c64-5fa6137a0d96.gif)

4. Vertical scan:

![Vertical](https://user-images.githubusercontent.com/38077812/111377572-ce43b380-866e-11eb-8b8f-e6ccd2e74a68.gif)

5. 3D surface view:

![surface](https://user-images.githubusercontent.com/38077812/111377787-1236b880-866f-11eb-8e52-60f3235085df.gif)

6. Regression analysis:

![fit](https://user-images.githubusercontent.com/38077812/111377799-16fb6c80-866f-11eb-96cd-f01dff3425ab.gif)

7. Interactive data visualization:

![report](https://user-images.githubusercontent.com/38077812/111377803-18c53000-866f-11eb-94ff-4ea16daaef3e.gif)

8. Kikuchi line simulation:

![kikuchi](https://user-images.githubusercontent.com/38077812/111377813-1bc02080-866f-11eb-8043-28bc199f8cd5.gif)

9. Domain boundary statistics:

![simulation](https://user-images.githubusercontent.com/38077812/111377823-1f53a780-866f-11eb-8b26-4638de0200c0.gif)

### Tutorial

A detailed Chinese guide (UI, workflows, scenarios, FAQ) is available here:

**[→ TUTORIAL_zh.md](TUTORIAL_zh.md)**

### Contact

Questions or suggestions about the **original** project: [yux1991@gmail.com](mailto:yux1991@gmail.com)

This fork: [https://github.com/JMjimike/PyRHEED](https://github.com/JMjimike/PyRHEED)

---

<a id="中文"></a>

## 中文

### 目录
1. [上游与本仓库改动](#上游与本仓库改动)
2. [项目简介](#项目简介)
3. [环境要求](#环境要求)
4. [安装与启动](#安装与启动)
5. [模块说明](#模块说明)
6. [功能演示](#功能演示)
7. [详细教程](#详细教程)
8. [联系方式](#联系方式)

### 上游与本仓库改动

本仓库基于 [yux1991/PyRHEED](https://github.com/yux1991/PyRHEED)（MIT License）衍生开发，保留原作者 Yu Xiang 的核心 RHEED 分析与模拟能力，并在此基础上做了以下扩展与修复：

| 类别 | 改动说明 |
|------|----------|
| **Bug 修复** | NumPy 2.x 兼容：将已移除的 `ndarray.ptp()` 改为 `np.ptp()`，修复 TAPD「Add Structure」崩溃 |
| **文档** | 中文 README、[`TUTORIAL_zh.md`](TUTORIAL_zh.md) 详细教程与工作流 |
| **工程化** | 使用 [uv](https://docs.astral.sh/uv/) 管理依赖（`pyproject.toml` / `uv.lock`） |
| **CIF 样品库** | [`cif_samples/`](cif_samples/)：STO 衬底、REBCO O₇/O₆、DyBCO 生长中间相（211/124/143/385 等） |
| **XRD 扩展** | 新增四圆几何、Ewald、连续扫描等模块（`xrd_*.py`、`ewald.py`、`four_circle_geometry.py` 等）及 `scripts/` 测试脚本 |

上游官方安装方式为 `pip install pyrheed` / `python -m pyrheed`；本仓库推荐：

```bash
uv sync
uv run src/main.py
```

### 项目简介

本项目用于**反射高能电子衍射（RHEED）**数据的分析与理论模拟。

RHEED 是一种电子衍射技术，使用相对高能量（5~30 keV）的电子束，以掠入射角照射样品表面。它对表面极为敏感，穿透深度通常只有数纳米。由于电子的散射因子比 X 射线高约四个数量级，RHEED 特别适合表征 XRD 难以检测的二维材料（如石墨烯）。RHEED 的另一优势是束斑尺寸很大（约 1 cm），可在晶圆尺度上统计材料的晶格常数、晶粒取向分布乃至缺陷密度等信息。

本项目基于 **Python 3.12（64 位）** 编写与测试，图形界面使用 **PyQt6** 构建。*simulate_RHEED* 模块借助 [pymatgen](http://pymatgen.org/) 读取 CIF 文件并构建晶体结构。主要功能包括：

1. 使用 [rawpy](https://pypi.org/project/rawpy/) 处理 RHEED 原始图像，并通过 [numpy](https://www.numpy.org/) 向量化加速强度剖面提取；自动构建二维倒易空间图与极图；三维数据可导出为 *.vtp* 格式，供 [ParaView](https://www.paraview.org) 进一步处理。
2. 批量拟合 RHEED 线剖面，支持预定义峰形函数（高斯、Voigt 等），并保存格式化结果。
3. 可视化拟合结果并生成分析报告。
4. 基于马尔可夫过程模拟统计因子，假设六方表面具有特定台阶密度分布。详见 [Spadacini 等](https://www.sciencedirect.com/science/article/pii/0039602883904922) 的论文。
5. 从 CIF 文件读取晶体结构，通过堆叠不同晶态材料构建自定义结构，并基于运动学衍射理论计算衍射图样。
6. 模拟给定结构的三维倒易空间；可在程序内构建原子模型，专门用于模拟二维平移反相畴（TAPD）模型的衍射，详见 [Lu 等](https://www.sciencedirect.com/science/article/pii/0039602881905410) 的论文。

### 环境要求

- Python ≥ 3.12
- [uv](https://docs.astral.sh/uv/getting-started/installation/) 包管理器
- 支持 Linux 与 Windows；macOS 尚未充分测试，但理论上可直接运行

### 安装与启动

```bash
git clone https://github.com/JMjimike/PyRHEED.git
cd PyRHEED
uv sync
uv run src/main.py
```

1. **加载数据** — 通过「打开文件」对话框加载 RAW 或常见图像格式；自动转灰度。高动态范围请用 TIFF / PNG / RAW，避免 JPEG。
2. **数据分析** — 按研究目的选用模块；图样/结构因子模拟**不依赖**实验图像。
3. **场景批处理** — 加载 `default_scenario.ini` 或自定义场景后运行。详见 [详细使用教程与工作流](TUTORIAL_zh.md)。

### 模块说明

| 模块 | 功能 |
|------|------|
| `broadening` | RHEED 线剖面批量自动拟合 |
| `browser` | 浏览工作目录中的文件 |
| `canvas` | 显示图像、绘制图形并接收用户输入 |
| `configuration` | 修改默认配置 |
| `cursor` | 显示光标相关信息 |
| `generate_report` | 生成拟合结果可视化报告 |
| `graph_3D_surface` | 倒易空间三维曲面可视化 |
| `kikuchi` | 模拟非重构晶体表面的 Kikuchi 图样 |
| `main` | 主程序入口 |
| `manual_fit` | 手动拟合线剖面，用于初始化拟合参数 |
| `my_widgets` | 自定义控件 |
| `plot_chart` | 基于 QChart 的拟合结果可视化控件 |
| `preference` | 修改默认设置 |
| `process` | 后端处理逻辑 |
| `process_monitor` | 内存占用监控（尚未完善） |
| `profile_chart` | 基于 QChart 的线扫描剖面可视化控件 |
| `properties` | 控制程序动态参数 |
| `reciprocal_space_mapping` | 构建二维/三维倒易空间图与极图 |
| `scenario` | 选择或自定义场景并批量运行 |
| `simulate_RHEED` | 从给定原子结构模拟衍射图样；支持平移反相畴（APD）模型 |
| `statistical_factor` | 基于马尔可夫过程计算统计因子 |
| `translational_antiphase_domain` | 计算平移 APD 模型的一维/二维剖面 |
| `window` | 主窗口 |
| `write_scenario` | 写入默认场景配置 |
| `xrd_*` / `ewald` / `four_circle_geometry` | XRD / 几何扩展（本仓库新增） |

### 功能演示

（演示动图与上方英文部分相同，见 [Demonstrations](#demonstrations)。）

### 详细教程

完整的界面说明、典型工作流、场景配置与常见问题，请参阅：

**[→ 详细使用教程与工作流（TUTORIAL_zh.md）](TUTORIAL_zh.md)**

### 联系方式

原项目问题或建议：yux1991@gmail.com

本仓库：https://github.com/JMjimike/PyRHEED
