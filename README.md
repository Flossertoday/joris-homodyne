# Joris homodyne

本项目研究非共振泵浦、optical-Stark 触发的 polariton Josephson dimer，并为后续单次扫描轨迹、相位和强度关联分析建立可核查的数据流程。

当前优先级不是直接宣称 squeezing 或把 Cauchy–Schwarz inequality（CSI）违反设为成败标准，而是：

1. 确认积累测量中看到的强度振荡能否在 single-sweep 数据中恢复；
2. 检查左右阱是否呈反相交换，并排除 reservoir/common-mode 瞬态；
3. 在此基础上再测相对相位和完整的 \(g^{(2)}_{ij}(t_1,t_2)\) 矩阵。

## 目录

- `DATA/`：原始实验数据及与实验一一对应的分析副本；不进入 Git。
- `analysis/`：可复现的数据分析代码和可审阅的结果图。
- `docs/`：用户可审阅的项目判断、实验路线、文献依据和数据约定。
- `scripts/`：可重复运行的数据转换或分析程序。
- `tests/`：脚本测试。
- `.codex/`：项目级 Codex 配置。
- `MEMORY.md`：简短的项目状态和细节索引，不承载长篇论证。

## 当前数据

- 原始数据：`DATA/26-07-31-single-sweep/`
- 无损分析副本：`DATA/26-07-31-derived-single-sweep/26-07-31_single_sweep.h5`
- 转换器：`scripts/convert_matlab_v73_to_analysis_hdf5.py`

原始 MATLAB 7.3 文件保持不变；日常 Python 分析读取实验专属的 HDF5 副本。当前副本只包含 streak-camera 图像，两个条件由 `ZeroOrderPowerRatio` 标识为 0.05 和 0.18。

当前求和示例选择 0.18 条件；1000 帧的逐像素和呈现四条亮带。见 [分析说明与结果](analysis/README.md)。

## 文档入口

- [实验路线](docs/experiment-plan.md)
- [测量与物理主张](docs/measurement-claims.md)
- [文献判断和来源](docs/literature.md)
- [数据与环境约定](docs/data.md)

