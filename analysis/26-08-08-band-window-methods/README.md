# 带窗定位方法对比（26-08-08）

本目录承载 26-08-04 HBT 数据带窗（ROI）定位方法对比调查的全部可复现交付物。
调查结论、方法定义、定量结果与文献依据见
[`docs/26-08-08-band-window-methods.md`](../../docs/26-08-08-band-window-methods.md)；
背景见 `docs/26-08-04-fig2a.md` 的带窗定位方法节。

- `compare_band_window_methods.py` —— 8 种候选方法的实现、公共局部地板与评分；
- `compare_band_window_figures.py` —— 方法栅格/参数平台/增长曲线/合成基准四图；
- `compute_full_column_profiles.py` —— 两批次全列 analog 廓线（复核检测列界用）；
- `render_roi_morphology.py` —— 发射支撑区行走 + ROI 正确性指标 + 形态图；
- `band_windows_method_grid.png` —— 方法栅格对比（sqrt 压缩纵轴，两批次）；
- `band_windows_robustness.png` —— 阈值参数平台；
- `band_windows_snr_growth.png` —— S/√(S+N) 增长曲线与各方法窗宽落点；
- `band_windows_synthetic.png` —— 合成恢复基准；
- `roi_fullrange_overview.png` —— 全列 analog 廓线：检测列界 350–1150 排除/切到什么；
- `roi_band_morphology.png` —— 逐带发射形态与各方法窗缘位置；
- `results/` —— `band_windows_summary.json`、`roi_morphology_summary.json`。

测试：`tests/test_compare_band_window_methods.py`（17 项，自含合成数据）。
