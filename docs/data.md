# 数据与环境约定

## 原则

- 原始 MATLAB 7.3 文件是不可变的 source of truth，不删除、不覆盖、不原地转换。
- 每次实验使用一对同级目录：原始目录采用 `YY-MM-DD-<experiment>`，分析目录采用 `YY-MM-DD-derived-<experiment>`。
- 分析副本必须保留条件、采集顺序、源文件、帧号、轴状态和哈希 provenance。
- 机器专属解释器路径、包版本和一次性检查输出不写入长期项目文档。

## 2026-07-31 single sweep

- 原始目录：`DATA/26-07-31-single-sweep/`
- 分析目录：`DATA/26-07-31-derived-single-sweep/`
- 分析文件：`26-07-31_single_sweep.h5`

本批实验已确认为 streak-only。像素 payload 只取 `structVar/StreakCamImage`；`MainCamImage`、`SecondCamImage` 和 Cam2 状态不进入分析条件或验证。`ZeroOrderPowerRatio` 是 central-dot strength 的权威字段，对应 0.05 和 0.18 两个各 1000 帧的条件。

HDF5 中每帧都能追溯到源 MAT 文件和源文件内帧号，并记录源文件 SHA-256。根属性 `source_root` 是转换时的绝对路径快照，目录重命名或移动后可能过时；定位现有源文件时，应将当前原始目录与 HDF5 内的 `source_files/relative_path` 组合。精确文件大小、输出哈希和逐像素验证证据只维护在 `DATA/26-07-31-derived-single-sweep/README.md`，避免在多份文档中复制后过时。

## Python

本地科学代码使用从机器 `sci` 环境克隆出的、位于非 OneDrive 路径的项目环境。核心库由 `requirements-analysis.txt` 声明；不要在仓库中提交环境目录。受限沙盒统一从 `scripts/init_sandbox.ps1` 进入：脚本优先复用工作树 `.venv`，否则从 `sci` 建立带系统包的 `.venv`，然后只安装缺失依赖。

初始化并运行完整测试：

```powershell
.\scripts\init_sandbox.ps1
& $env:PROJECT_PYTHON -m unittest discover -s tests -p "test_*.py" -v
```

初始化后，`PROJECT_PYTHON` 是权威解释器，`JORIS_DATA_ROOT` 是实际数据根。`-CheckOnly` 不创建环境或安装包；`-AsJson` 输出适合自动化保存的预检摘要。详见 [沙盒初始化方法与反馈](sandbox-initialization.md)。

转换器拒绝覆盖既有输出，并先写 `.partial` 文件后原子改名。需要重新生成时，应指定新的输出文件名；不要覆盖已经验证的分析副本。

## 轴标定

当前 HDF5 的两个 detector axes 仍标记为 `uncalibrated`。在获得实验记录或标定图前，不应擅自把它们写成 time、space、energy、left/right 或确定方向的 x/y。

