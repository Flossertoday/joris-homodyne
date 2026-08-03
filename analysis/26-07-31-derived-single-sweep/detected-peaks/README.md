# detected-peaks：central_dot_0p18 原始帧的逐帧峰检测数据集

## 这里面对 DATA 做了什么

输入是 `DATA/26-07-31-derived-single-sweep/26-07-31_single_sweep.h5` 中 `central_dot_0p18` 条件的 1000 张原始探测器帧。本目录的数据是把这 1000 张帧**逐帧做稀疏峰检测**得到的峰坐标表：每张帧不再按像素存储，而是压缩成其中亮峰的亚像素坐标（每帧中位 30 个峰），作为后续分析（x 剖面拟合、关联等）的输入数据。处理分两步，各产生一个 HDF5 文件：

1. `../extract_condition_peaks.py` 对全部 1000 帧逐帧运行 `find_frame_peaks.py` 的生产峰检测（trackpy Crocker–Grier），把每帧的 `(frame, y, x, mass)` 汇总成 `26-07-31_central_dot_0p18_peaks.h5`（全量峰表，30542 峰）；
2. `../fit_spot_x_ranges.py` 把全量峰表按 x 做直方图并拟合"常数背景 + 四高斯"，只保留横坐标落在四个光斑接受区间内的峰，得到 `26-07-31_central_dot_0p18_peaks_xrange2p0sigma.h5`（x 范围筛选峰表，15958 峰）。

两个文件的逐字段 schema、provenance、校验与复现命令如下。

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
  .\analysis/26-07-31-derived-single-sweep/detected-peaks/26-07-31_central_dot_0p18_peaks.h5 `
  --central-dot-strength 0.18
```

## `26-07-31_central_dot_0p18_peaks_xrange2p0sigma.h5`（x 范围筛选峰表）

由 `../fit_spot_x_ranges.py` 生成：把全量峰表按 x 做 1 px 直方图（跨帧求和、压缩纵轴），用"常数背景 + 四高斯"拟合；四个高斯从左到右对应 condensate A arm A、condensate B arm A、condensate A arm B、condensate B arm B。只有横坐标落在四个高斯接受区间内的峰才可能是真实数据；本文件只保留这些峰（schema 与全量峰表相同，另附筛选描述与区间属性）。

判据与重叠处理：保守起点取每高斯 μ ± 3σ。但 3σ 下同一 arm 内两个光斑的区间重叠（spot 1–2 共享 [446.50, 464.00]，spot 3–4 共享 [685.39, 709.04]），而重叠区内一个峰属于哪个光斑目前**没有解决办法**（脚本 `disjoint_intervals` 的注释记录了这一点）。因此脚本自动把倍数按 0.5 步长收缩到区间互不重叠为止：本数据集落在 **μ ± 2σ**。`--allow-overlap` 可关闭收缩（保留重叠区间）。

拟合结果（30542 峰，背景 10.94 峰/bin；完整不确定度见 `../results/26-07-31_central_dot_0p18_four_gaussians_fit.csv`，诊断图见 `../results/26-07-31_central_dot_0p18_xprofile_four_gaussians.png`）：

| 光斑 | 物理对应 | μ [px] | σ [px] | 接受区间 (2σ) [px] |
|---|---|---|---|---|
| 1 | condensate A arm A | 425.93 | 12.69 | [400.55, 451.31] |
| 2 | condensate B arm A | 489.83 | 14.44 | [460.95, 518.72] |
| 3 | condensate A arm B | 657.03 | 17.34 | [622.35, 691.70] |
| 4 | condensate B arm B | 734.24 | 16.28 | [701.67, 766.81] |

四个区间互不重叠，保留 15958 / 30542 峰（52.2%）；被排除的主要是全探测器均匀分布的背景检出（约 11 峰/bin）和各高斯 2σ 以外的尾部。

复现（同样拒绝覆盖；收缩是自动的，起点倍数可用 `--sigma-multiple` 调）：

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/fit_spot_x_ranges.py `
  .\analysis/26-07-31-derived-single-sweep/detected-peaks/26-07-31_central_dot_0p18_peaks.h5 `
  --figure .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_xprofile_four_gaussians.png `
  --fit-csv .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_four_gaussians_fit.csv `
  --filtered-h5 .\analysis/26-07-31-derived-single-sweep/detected-peaks/26-07-31_central_dot_0p18_peaks_xrange2p0sigma.h5
```

筛选文件名记录了最终使用的判据倍数；若未来解决重叠归属问题并采用更宽区间，应以新文件名重新生成。
