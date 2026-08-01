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
.\.venv\python.exe `
  .\analysis\sum_streak_condition.py `
  .\DATA\26-07-31-derived-single-sweep\26-07-31_single_sweep.h5 `
  .\analysis\results\26-07-31_central_dot_0p18_sum.png `
  --central-dot-strength 0.18
```

脚本通过 HDF5 condition 的 `central_dot_strength` 属性选择数据，不依赖 group 名称，并按批次读取，避免把全部三维帧栈一次载入内存。
