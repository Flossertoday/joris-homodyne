# 数据分析

这个目录保存可复现的数据分析代码和可审阅的小型结果图；原始帧与大型中间数组仍留在忽略 Git 的 `DATA/` 中。

## 2026-07-31 central-dot 条件求和

`DATA/26-07-31-single-sweep/` 中的第 2、3 组不是两个需要分别转换的同类文件，而是同一次实验批次里的两个 central-dot 条件。规范化 HDF5 因而用一个文件承载两组数据：

- `conditions/central_dot_0p05/frames`：原始第 2 组，central-dot strength 0.05，1000 帧；
- `conditions/central_dot_0p18/frames`：原始第 3 组，central-dot strength 0.18，1000 帧。

本次选择 0.18 条件，因为它的逐像素帧和清楚呈现四条亮带。结果见 [`results/26-07-31_central_dot_0p18_sum.png`](results/26-07-31_central_dot_0p18_sum.png)。

计算保持原始 detector 网格，对 1000 帧执行逐像素 `uint64` 求和：

\[
S(i,j)=\sum_{n=1}^{1000} I_n(i,j).
\]

PNG 使用线性色阶，并把显示范围裁剪到求和图的第 1–99.9 百分位；这只影响显示，不修改求和数组。两个 detector axes 尚未标定，因此图中只标为像素轴。

在项目根目录运行：

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/sum_streak_condition.py `
  (Join-Path $sandboxState.data_root "26-07-31-derived-single-sweep\26-07-31_single_sweep.h5") `
  .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_sum.png `
  --central-dot-strength 0.18
```

脚本通过 HDF5 condition 的 `central_dot_strength` 属性选择数据，不依赖 group 名称，并按批次读取，避免把全部三维帧栈一次载入内存。

## 2026-08-03 central-dot 0.18 单帧渲染

`render_condition_frames.py` 把 0.18 条件下均匀间隔的 100 帧（每第 10 帧，索引 0–990）写成 16-bit 灰度 PNG，结果在 `drafts/26-07-31_central_dot_0p18_frames/`（约 143 MB 的草稿输出，已被 Git 忽略，不随仓库分发）。PNG 的 `(y, x)` 即帧数组的 `[detector_axis_1, detector_axis_2]`，不做翻转。

数据在 16-bit 容器里只用到 0–6499 counts（全条件 1000 帧的全局范围），精确值 PNG 在普通查看器里接近全黑，因此渲染统一经过一条对该全局范围定义的平方根 LUT：

\[
\mathrm{PNG}(v)=\Big\lfloor 65535\sqrt{\tfrac{v-0}{6499-0}}+0.5\Big\rceil,
\quad v\in[0,6499].
\]

该 LUT 严格单调递增（已逐值验证无碰撞），因此像素强度与探测器计数保持一一对应；背景（约 105 counts）映射到约 13% 灰度，噪声纹理可见，稀疏亮事件映射到高亮。同一 LUT 应用于全部 100 张图，帧间亮度可互相比较。已逐张回读校验，100 张全部等于 LUT(源帧)。

单帧接近单光子计数 regime：背景约 105 counts，只有约 0.03% 像素显著发亮（峰值数千 counts）；四条亮带是 1000 帧求和后才显现的统计结构。

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/render_condition_frames.py `
  (Join-Path $sandboxState.data_root "26-07-31-derived-single-sweep\26-07-31_single_sweep.h5") `
  .\drafts\26-07-31_central_dot_0p18_frames `
  --central-dot-strength 0.18
```

## 2026-08-03 central-dot 0.18 单帧峰检测

`find_frame_peaks.py` 定位单帧中的稀疏类高斯峰。核心是 `trackpy.locate`（Crocker–Grier）：带通滤波后在滤波域取局部极大，再在带通图上以质心精修到亚像素。`find_peaks(frame)` 返回每峰一行 `(y, x, mass)`；`mass` 是带通后的积分亮度，只作诊断，峰的亮度与直径不进入结果。

默认参数按 0.18 条件的实测特性标定（背景全局平坦 ~105 counts、读出噪声 σ≈3、峰 FWHM≈6–7 px、每帧 ~500 个单像素脉冲噪声、脉冲最大仅 ~110 counts 超出背景）：

- `diameter=7`：匹配实测峰宽；
- `percentile=0`：关闭 trackpy 默认的 64 百分位亮度截断。该截断在这类稀疏帧上会误删峰顶仅超背景 ~110–150 counts 的真实暗峰；
- `minmass=300`：带通后积分亮度下限。实测脉冲噪声的 mass ≲85，余量充足。

灵敏度下限约为峰顶超出背景 ~110–150 counts；比这更暗的斑不计入。抽验帧（0/230/500/700/999）按"7×7 窗口内超阈像素数 ≥5"的独立判据核对，检出零误检。

验证：合成数据单元测试在 `tests/test_find_frame_peaks.py`（含暗峰、贴边峰、空帧与实测强度的脉冲噪声）。对 0.18 条件全部 1000 帧运行：每帧检出数中位 30、均值 30.5、范围 16–53（p1–p99：18–45），与单帧 ~20 个量级一致。帧 230 的叠加图（22 峰）与 `drafts/26-07-31_central_dot_0p18_frames/frame_0230.png` 目检一致：圈全部落在可见峰上、无空圈；仅一个峰顶 ~108 counts 的极暗斑低于灵敏度下限未检出。

在项目根目录运行：

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/find_frame_peaks.py `
  (Join-Path $sandboxState.data_root "26-07-31-derived-single-sweep\26-07-31_single_sweep.h5") `
  .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_frame0230_peaks.csv `
  --central-dot-strength 0.18 --frame-index 230 `
  --overlay-png .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_frame0230_peaks.png
```

CSV 每行一个峰 `y,x,mass`；`--diameter` 与 `--minmass` 可调，其余行为与条件选择逻辑同 `sum_streak_condition.py`。

## 2026-08-03 峰定位方法学审计：流程可视化与高斯拟合对照

两个审计脚本为方法学报告服务（结论与引用见 [`docs/peak-localization.md`](../../docs/peak-localization.md)），默认参数直接 import 自 `find_frame_peaks.py`，不会与被审计的生产参数发生漂移：

- `visualize_peak_pipeline.py` 逐步复现 `tp.locate` 在单帧上的操作并渲染审查图：全帧六联版（原始帧 → 高斯低通 → boxcar 背景 → 带通 → 候选检出与 minmass 筛选 → 精修质心叠加）加两个峰的逐像素质心迭代放大。frame 230：29 个候选、7 个被 minmass 拒绝、22 个保留。注意精修在带通图上进行；脚本为可视化用浮点带通图，而 `locate` 内部用 uint8 量化副本（对质心影响 <0.02 px）。
- `compare_gauss_refinement.py` 对可复现抽样（`--n-frames 50 --seed 20260803`，脚本会打印帧号清单）逐峰对照 Crocker–Grier 质心与 `tp.refine_leastsq` 二维高斯拟合位置，输出逐峰 CSV、全部审计统计（终端打印）、四联总结图与最大离群峰放大图。1518 峰：92.2% 差值 ≤0.2 px（中位 0.052 px，无系统偏差）；离群全部是高斯拟合撞到 ±3 px 位置边界（非高斯形态斑），质心始终留在亮核上。

在项目根目录运行：

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/visualize_peak_pipeline.py `
  (Join-Path $sandboxState.data_root "26-07-31-derived-single-sweep\26-07-31_single_sweep.h5") `
  --central-dot-strength 0.18 --frame-index 230 `
  --steps-png .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_frame0230_pipeline_steps.png `
  --zoom-png .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_frame0230_refinement_zoom.png

& $sandboxState.python `
  .\analysis/26-07-31-derived-single-sweep/compare_gauss_refinement.py `
  (Join-Path $sandboxState.data_root "26-07-31-derived-single-sweep\26-07-31_single_sweep.h5") `
  .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_gauss_vs_centroid_peaks.csv `
  --central-dot-strength 0.18 --n-frames 50 --seed 20260803 `
  --figure .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_gauss_vs_centroid.png `
  --outliers-figure .\analysis/26-07-31-derived-single-sweep/results/26-07-31_central_dot_0p18_gauss_outliers_zoom.png
```
