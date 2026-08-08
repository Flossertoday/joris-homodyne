# 带窗定位方法对比（26-08-08）

本目录承载 26-08-04 HBT 数据带窗（ROI）定位方法对比调查的全部可复现交付物。
调查结论、方法定义、定量结果与文献依据见
[`docs/26-08-08-band-window-methods.md`](../../docs/26-08-08-band-window-methods.md)；
背景见 `docs/26-08-04-fig2a.md` 的带窗定位方法节。

- `compare_band_window_methods.py` —— 8 种候选方法的实现、公共局部地板与评分；
- `compare_band_window_figures.py` —— 四张定稿图与合成基准；
- `band_windows_method_grid.png` —— 方法栅格对比（sqrt 压缩纵轴，两批次）；
- `band_windows_robustness.png` —— 阈值参数平台；
- `band_windows_snr_growth.png` —— S/√(S+N) 增长曲线与各方法窗宽落点；
- `band_windows_synthetic.png` —— 合成恢复基准；
- `results/band_windows_summary.json` —— 机读评分汇总。

测试：`tests/test_compare_band_window_methods.py`（14 项，自含合成数据）。
