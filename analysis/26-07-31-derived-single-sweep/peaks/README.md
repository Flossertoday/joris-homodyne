# central_dot_0p18 峰数据集

本目录保存 0.18 条件全量 1000 帧的稀疏峰检测结果，作为后续分析（x 剖面拟合、关联等）的原始数据。

## `26-07-31_central_dot_0p18_peaks.h5`（全量峰表）

由 `../extract_condition_peaks.py` 生成：对 `DATA/26-07-31-derived-single-sweep/26-07-31_single_sweep.h5` 的 `central_dot_0p18` 条件逐帧运行 `find_frame_peaks.py` 的生产峰检测（trackpy Crocker–Grier，`diameter=7, minmass=300, percentile=0`）。每个峰一行：

- `frame`（int32）：0–999 帧号；
- `y`, `x`（float64）：亚像素坐标，对应帧数组的 `[detector_axis_1, detector_axis_2]`；
- `mass`（float64）：带通积分亮度，仅作诊断。

根属性记录完整 provenance：条件名、`central_dot_strength`、检测参数、`frame_count`、`peak_count`、`frame_shape`、源 HDF5 文件名与 SHA-256、生成时间与脚本名。

校验（2026-08-03）：共 30542 峰；逐帧计数中位 30、均值 30.54、范围 16–53（p1–p99：18–45），与 `find_frame_peaks.py` 文档记录的试运行统计一致；帧 230 的 22 行与已审阅 CSV `results/26-07-31_central_dot_0p18_frame0230_peaks.csv` 逐行一致（最大差 5e-5，即 CSV 的四位小数舍入）。

复现（在项目根目录，拒绝覆盖既有输出，重跑需先改名）：

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/extract_condition_peaks.py `
  (Join-Path $sandboxState.data_root "26-07-31-derived-single-sweep\26-07-31_single_sweep.h5") `
  .\analysis/26-07-31-derived-single-sweep/peaks/26-07-31_central_dot_0p18_peaks.h5 `
  --central-dot-strength 0.18
```

## `26-07-31_central_dot_0p18_peaks_xrange3p0sigma.h5`（x 范围筛选峰表）

由 `../fit_spot_x_ranges.py` 生成：把全量峰表按 x 做 1 px 直方图（跨帧求和、压缩纵轴），用"常数背景 + 四高斯"拟合；四个高斯从左到右对应 condensate A arm A、condensate B arm A、condensate A arm B、condensate B arm B。保守判据取每高斯 μ ± 3σ，接受域为四区间并集；只保留横坐标落在接受域内的峰（schema 与全量峰表相同，另附筛选描述与区间属性）。

拟合结果（30542 峰，背景 10.94 峰/bin；完整不确定度见 `../results/26-07-31_central_dot_0p18_four_gaussians_fit.csv`，诊断图见 `../results/26-07-31_central_dot_0p18_xprofile_four_gaussians.png`）：

| 光斑 | 物理对应 | μ [px] | σ [px] | 接受区间 [px] |
|---|---|---|---|---|
| 1 | condensate A arm A | 425.93 | 12.69 | [387.86, 464.00] |
| 2 | condensate B arm A | 489.83 | 14.44 | [446.50, 533.16] |
| 3 | condensate A arm B | 657.03 | 17.34 | [605.02, 709.04] |
| 4 | condensate B arm B | 734.24 | 16.28 | [685.39, 783.09] |

成对区间相互重叠，并集实际为两段连续区域：约 [387.9, 533.2]（arm A）与 [605.0, 783.1]（arm B）。保留 17481 / 30542 峰（57.2%）；被排除的主要是全探测器均匀分布的背景检出（约 11 峰/bin）。

复现（同样拒绝覆盖）：

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/fit_spot_x_ranges.py `
  .\analysis/26-07-31-derived-single-sweep/peaks/26-07-31_central_dot_0p18_peaks.h5 `
  --figure .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_xprofile_four_gaussians.png `
  --fit-csv .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_four_gaussians_fit.csv `
  --filtered-h5 .\analysis/26-07-31-derived-single-sweep/peaks/26-07-31_central_dot_0p18_peaks_xrange3p0sigma.h5 `
  --sigma-multiple 3.0
```

`--sigma-multiple` 可调（如更严的 2.5 或更宽的 3.5），改判据倍数不需要改代码；筛选文件名应随之改名以保留判据记录。
