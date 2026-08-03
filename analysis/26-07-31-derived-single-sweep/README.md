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
