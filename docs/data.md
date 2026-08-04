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

本地科学代码使用以机器 `sci` 解释器创建的隔离 `.venv`（不共享其 site-packages，也不向 `sci` 安装任何包）。核心库由钉版本的 `requirements.txt` 声明；不要在仓库中提交环境目录。Codex Desktop 激活 `.codex/environments/environment.toml` 时会自动调用 `scripts/init_sandbox.ps1`；若未激活或 setup 失败，代理按 `AGENTS.md` 自行调用。脚本优先复用隔离的工作树 `.venv`，否则从 `sci` 重建不带系统包的隔离 `.venv`，并在每次运行时按钉版清单 pip 同步（全部满足时 pip 本地 no-op）。

自动初始化后运行完整测试：

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -p "test_*.py" -v
```

初始化结果保存在忽略 Git 的 `tmp/sandbox-initialization.json`，其中 `python` 是权威解释器、`data_root` 是实际数据根。用户无需维护这份状态；若状态缺失或失败，`AGENTS.md` 要求 Codex 代理自行恢复。详见 [沙盒初始化方法与反馈](sandbox-initialization.md)。

转换器拒绝覆盖既有输出，并先写 `.partial` 文件后原子改名。需要重新生成时，应指定新的输出文件名；不要覆盖已经验证的分析副本。

## 轴标定

纵轴（detector axis 1）已确认为 streak 时间轴：整帧高度对应采集时间窗，窗口由源 metadata 的 `StreakRange` 字段（HPD-TA "Time Range" 设备档位，枚举值）配合采集 README 的记录确定；26-07-31 批次为 1 ns（1016 px，约 0.984 ps/px）。约定映射集中在根目录 `metadata-conventions.toml`。横轴（detector axis 2）仍是未标定的探测器像素，不应写成 space、energy、left/right 或确定方向。派生 HDF5 中的 `axis_semantics_status="uncalibrated"` 是转换时的保守标签，轴语义以本节和 `metadata-conventions.toml` 为准。

