# streak 纵轴时间窗：metadata 调查与约定

本文记录 2026-08-04 对 26-07-31 single-sweep 批次源 metadata 的调查，结论是根 `AGENTS.md` 中 streak 渲染规则（纵轴为时间并正确缩放）和根目录 `metadata-conventions.toml` 约定值的依据。

## 结论速览

- 纵轴（detector axis 1）= 时间；横轴（detector axis 2）= 探测器像素，无坐标变换。
- 26-07-31 批次时间窗 = **1 ns**，覆盖整帧 1016 px，即 ≈ 0.984 ps/px。
- 时间范围字段：`structVar/StreakRange`（HPD-TA "Time Range" 设备档位，本批次枚举值 5）；窗口实际值以采集 `README.txt` 的 "Scan of 1 ns window" 为准。
- 约定值集中在根目录 `metadata-conventions.toml`，渲染规则在根 `AGENTS.md`。

## 调查对象与方法

抽查文件（MATLAB 7.3 = HDF5，每个文件一个 `structVar` group，约 300–400 个字段；字符串字段是 uint16 字符数组）：

- `2 1000 images single sweep/100_image_savedStructure.mat`（条件 0.05）；
- `3 1000 images single sweep/112_image_savedStructure.mat` 与 `211_image_savedStructure.mat`（条件 0.18）；
- `1 Experiment structure parameters/26-07-31_17-54-10_FourDots_ConjugateGradient_savedStructure.mat`（批次结构参数快照）。

复核方法（任意 image MAT 文件）：

```python
import h5py
with h5py.File("DATA/26-07-31-single-sweep/2 1000 images single sweep/100_image_savedStructure.mat", "r") as h:
    s = h["structVar"]
    print(float(s["StreakRange"][0, 0]))          # 5.0
    print(float(s["ZeroOrderPowerRatio"][0, 0]))  # 0.05
    print(s["StreakCamImage"].shape)              # (10, 1016, 1344)
```

## 发现

1. **条件字段**：`ZeroOrderPowerRatio` 在 block 2 各文件为 0.05、block 3 为 0.18，与既有 central-dot 条件约定一致。
2. **streak 控制链路**：`StreakCamMode = 'Acquire'`，`StreakCmd = 'Appinfo(HPDTA)'`，经 TCP 端口 1001/1002 远程控制——采集软件是 Hamamatsu HPD-TA。
3. **时间范围字段**：`StreakRange` 在全部抽查文件（两个条件）中一致为 5.0。HPD-TA 的 "Time Range" 是设备枚举参数，允许值由扫描单元型号决定（[HPD-TA RemoteEx Programmers Handbook](http://www-bl20.spring8.or.jp/detectors/manual/RemoteExProgrammersHandbook91.pdf) 中 `DevParamInfoEx` 的应答示例即列出 "5 ns, 10 ns, 20 ns, …" 的允许值表），因此数字 5 本身不直接给出窗口。窗口实际值来自采集笔记：`DATA/26-07-31-single-sweep/README.txt` 记录 "Scan of 1 ns window"。
4. **帧布局**：`StreakCamImage` 形状 `(10, 1016, 1344)`——纵轴 1016 px 为时间，横轴 1344 px 为探测器像素。
5. **相关但未采用**：`StreakCamDelay = 24270`（触发延迟），单位未在 metadata 中记录（疑为 ps），不影响窗口缩放，未纳入约定。

## 边界与未决

- `StreakRange` 档位 → 时间窗的映射依赖具体扫描单元，未从数据本身获得；本批次按采集 README 定为 1 ns。将来改档位或换批次时，须由实验者确认新窗口并更新 `metadata-conventions.toml`。
- 横轴保持未标定像素，不应写成 space、energy、left/right 或确定方向。
- 派生 HDF5 中的 `axis_semantics_status="uncalibrated"` 是转换时的保守标签；轴语义以本文与 `metadata-conventions.toml` 为准。
