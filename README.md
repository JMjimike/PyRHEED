# <img src="https://github.com/yux1991/PyRHEED/blob/master/src/icons/icon.png" width="48"/> PyRHEED
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/yux1991/PyRHEED/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/yux1991/PyRHEED.svg)](https://github.com/yux1991/PyRHEED/blob/master/LICENSE) [![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](mailto:yux1991@gmail.com)

## 目录
1. [项目简介](README.md#项目简介)
2. [环境要求](README.md#环境要求)
3. [安装与启动](README.md#安装与启动)
4. [模块说明](README.md#模块说明)
5. [功能演示](README.md#功能演示)
6. [详细教程](TUTORIAL_zh.md)
7. [联系方式](README.md#联系方式)

## 上游与本仓库改动

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

## 项目简介

本项目用于**反射高能电子衍射（RHEED）**数据的分析与理论模拟。

RHEED 是一种电子衍射技术，使用相对高能量（5~30 keV）的电子束，以掠入射角照射样品表面。它对表面极为敏感，穿透深度通常只有数纳米。由于电子的散射因子比 X 射线高约四个数量级，RHEED 特别适合表征 XRD 难以检测的二维材料（如石墨烯）。RHEED 的另一优势是束斑尺寸很大（约 1 cm），可在晶圆尺度上统计材料的晶格常数、晶粒取向分布乃至缺陷密度等信息。

本项目基于 **Python 3.12.2（64 位）** 编写与测试，图形界面使用 **PyQt6** 构建。*simulate_RHEED* 模块借助 [pymatgen](http://pymatgen.org/) 读取 CIF 文件并构建晶体结构。主要功能包括：

1. 使用 [rawpy](https://pypi.org/project/rawpy/) 处理 RHEED 原始图像，并通过 [numpy](https://www.numpy.org/) 向量化加速强度剖面提取；自动构建二维倒易空间图与极图；三维数据可导出为 *.vtp* 格式，供 [ParaView](https://www.paraview.org) 进一步处理。
2. 批量拟合 RHEED 线剖面，支持预定义峰形函数（高斯、Voigt 等），并保存格式化结果。
3. 可视化拟合结果并生成分析报告。
4. 基于马尔可夫过程模拟统计因子，假设六方表面具有特定台阶密度分布。详见 [Spadacini 等](https://www.sciencedirect.com/science/article/pii/0039602883904922) 的论文。
5. 从 CIF 文件读取晶体结构，通过堆叠不同晶态材料构建自定义结构，并基于运动学衍射理论计算衍射图样。
6. 模拟给定结构的三维倒易空间；可在程序内构建原子模型，专门用于模拟二维平移反相畴（TAPD）模型的衍射，详见 [Lu 等](https://www.sciencedirect.com/science/article/pii/0039602881905410) 的论文。

## 环境要求

- Python ≥ 3.12
- [uv](https://docs.astral.sh/uv/getting-started/installation/) 包管理器
- 支持 Linux 与 Windows；macOS 尚未充分测试，但理论上可直接运行

## 安装与启动

### 1. 安装依赖

确保已安装 [uv](https://docs.astral.sh/uv/getting-started/installation/)。

### 2. 克隆并安装

```bash
git clone https://zivgitlab.uni-muenster.de/ag-salinga/pyrheed.git
cd pyrheed
uv sync
```

### 3. 启动程序

```bash
uv run src/main.py
```

### 4. 加载数据

RHEED 数据通常为图像文件。可通过「打开文件」对话框直接加载 RAW 原始图或常见压缩格式。程序会自动转换为灰度图。**注意：** JPEG 仅支持最高 8 位 RGB，若需要更高动态范围，请使用 TIFF、PNG 或 RAW 等格式。

### 5. 数据分析

根据研究目的，可使用不同分析模块。RHEED 图样模拟、结构因子模拟等模块**不依赖**实验数据，可独立运行。

### 6. 场景批处理

加载预定义场景（`default_scenario.ini`）或自定义场景，点击运行即可自动生成整套结果。详见 [详细使用教程与工作流](TUTORIAL_zh.md)。

## 模块说明

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

## 功能演示

1. 提取线剖面：

![line](https://user-images.githubusercontent.com/38077812/111377405-9a688e00-866e-11eb-8ef2-b25386f10d27.gif)

2. 球面特征提取：

![tilt](https://user-images.githubusercontent.com/38077812/111377452-aa806d80-866e-11eb-91eb-8a7f103c2077.gif)

3. 方位角 RHEED 构建：

![Azimuth](https://user-images.githubusercontent.com/38077812/111377562-cbe15980-866e-11eb-8c64-5fa6137a0d96.gif)

4. 垂直扫描：

![Vertical](https://user-images.githubusercontent.com/38077812/111377572-ce43b380-866e-11eb-8b8f-e6ccd2e74a68.gif)

5. 三维曲面视图：

![surface](https://user-images.githubusercontent.com/38077812/111377787-1236b880-866f-11eb-8e52-60f3235085df.gif)

6. 回归分析：

![fit](https://user-images.githubusercontent.com/38077812/111377799-16fb6c80-866f-11eb-96cd-f01dff3425ab.gif)

7. 交互式数据可视化：

![report](https://user-images.githubusercontent.com/38077812/111377803-18c53000-866f-11eb-94ff-4ea16daaef3e.gif)

8. Kikuchi 线模拟：

![kikuchi](https://user-images.githubusercontent.com/38077812/111377813-1bc02080-866f-11eb-8043-28bc199f8cd5.gif)

9. 畴界统计模拟：

![simulation](https://user-images.githubusercontent.com/38077812/111377823-1f53a780-866f-11eb-8b26-4638de0200c0.gif)

## 详细教程

完整的界面说明、典型工作流、场景配置与常见问题，请参阅：

**[→ 详细使用教程与工作流（TUTORIAL_zh.md）](TUTORIAL_zh.md)**

## 联系方式

如有问题或建议，请联系：yux1991@gmail.com
