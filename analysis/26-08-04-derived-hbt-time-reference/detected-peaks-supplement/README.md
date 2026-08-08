# detected-peaks-supplement：26-08-04 HBT 补充批次(72000 帧)逐光子事件表

## 这里面对 DATA 做了什么

输入是 `DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference_supplement.h5`
(7200 个 MAT 文件 × 10 帧 = 72000 张原始 streak 帧;转换与采集条件确认见
`docs/26-08-06-hbt-supplement.md` §1–§2)。本目录的数据是把这 72000 帧用与
`../detected-peaks/`(原 20000 帧批次)**完全相同的单遍检测**逐帧做光子事件检测得到的
逐事件坐标表:实测背景(自原批次复制,适用性已验证)扣除、参考脉冲 CFD50 定时间零点后,
把 8-连通亮斑(≥2 px、积分幅度 > 25 counts)压缩为亚像素坐标。与 fig2a 迹线由
`../extract_hbt_traces.py` 同一遍检测产出,末尾自动执行同遍自洽验证。

注意:补充批次与原批次**不作为合并数据集**(带内信号率 ~65%、L/R 比不同,见
`docs/26-08-06-hbt-supplement.md` §1),本表独立于原批次的 `../detected-peaks/` 事件表。

文件:

- `26-08-04_hbt_time_reference_supplement_events.h5`:事件表(661090 事件,
  71971/72000 帧保留);
- `26-08-04_hbt_time_reference_supplement_events_validation.json`:同遍自洽验证结果;
- `26-08-04_hbt_time_reference_supplement_events_sample_overlay.png`:3 帧抽样叠加图
  (事件点叠在扣背景帧上,白带边线、绿 t0 线),供目视检查检测质量。

## `26-08-04_hbt_time_reference_supplement_events.h5`(事件表)

schema 与原批次逐字一致:实测背景扣除 → 参考脉冲前沿恒比定时 CFD50 定 t0
(`t0_method=cfd50, cfd_fraction=0.5, t0_roi=(110,210,170,215)`)→ 8-连通组分检测
(阈值 25 counts 上背景、min 2 px、仅列 350–1150)→ 计数加权亚像素质心。阈值由灵敏度
扫描拐点自动选定(本遍:`knee: relative drop 0.060 < 0.25`,见 `threshold_method` 与
`threshold_study_json` 属性)。

每个事件一行,存于 `events/` 组:

- `frame`(int32):0–71999 帧号;
- `x`, `y`(float64):亚像素质心(探测器列、绝对行,px);
- `tau`(float64):`y - frame_t0[frame]`,抖动校正后的对齐坐标(px);
- `amplitude`(float64):组分积分幅度(counts 上背景)。

表级数据集:

- `frame_t0`(float64,72000 行):逐帧 CFD50 t0 行,t0 失败为 NaN;
- `failed_frames`(int32,29 行):被排除的帧号,这些帧无任何事件。

根属性记录完整 provenance(同原批次):源 HDF5 路径与 SHA-256(非默认文件现算)、
检测参数、t0 方法/ROI、`ps_per_px`、帧数/保留帧数/失败数/事件数、生成时间与脚本名。

## 横轴 ROI 与带窗

- **检测时**:只对列 350–1150 建掩膜;质心 `x` 实测落在
  [349.99999999999994, 1149.0000000000002],越界仅 ±2e-13(float64 加权均值边界舍入),
  无实质越界。
- **分析时**:`band_windows` 属性是本遍**自动带检测结果,仅 3 条部分带**
  ((550,660)/(790,850)/(980,1080))——自动检测在补充数据上失效(带外结构化辉光,
  72000 帧统计量下 baseline+5σ 只抓到峰顶子窗,见 `docs/26-08-06-hbt-supplement.md`
  §4)。**g2 定案带窗不以此属性为准**,下游统一用
  `FIG2C_BAND_WINDOWS="350,440,540,680,770,860,970,1100"`(B1 350–440、B2 540–680、
  B3 770–860、B4 970–1100)按 x 指派。验证 JSON 中的 `band_event_totals`
  (B1 140957 / B2 55997 / B3 115680)对应的是自动检出的 3 条部分带。

## 脉冲光校准(时间抖动校正)

每个事件的 `tau` 都以本帧参考脉冲的 CFD50 穿越行为原点,全表逐事件验证
`tau = y - frame_t0[frame]` 最大偏差 0.0;t0 失败的 29 帧被排除并记录在
`failed_frames`。t0 统计:中位行 159.43,RMS 抖动 9.26 ps
(`ps_per_px` = 0.984252,1 ns 单扫窗)。

## 验证(2026-08-07,单遍全量提取,累积 6220 s)

1. **同遍自洽**(`..._validation.json`):由表重建的 3 条自动带 photon 迹线与同遍 npz
   在 0–820 ps 窗口内最大差 [4,1,3] counts/bin(容差 6,按事件数为原批次 2.4 倍放宽),
   PASS。
2. **与旧双遍产物逐位一致**:本表与废止脚本第二遍的
   `drafts/26-08-06_fig2c_supplement/hbt_event_table.h5`(gitignored,保留作对照)逐列
   (frame/x/y/tau/amplitude)、逐帧 t0(含 NaN)、失败帧列表全部 `array_equal`;
   事件数 661090、保留帧 71971、失败帧 29 全部一致。因此
   `docs/26-08-06-hbt-supplement.md` §6 基于旧表的全部 g2 结论对本表继续有效。
3. **叠加图目视**:3 帧抽样(100/5000/12000)事件点均落在亮斑上,参考脉冲列无事件,
   t0 线穿过参考脉冲高度(见 `..._sample_overlay.png`)。

## 复现(拒绝覆盖既有输出,重跑需先改名或 `--events-h5` 另指路径)

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-08-04-derived-hbt-time-reference/extract_hbt_traces.py `
  --h5 "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference_supplement.h5" `
  --out-dir .\drafts\26-08-06_fig2a_supplement_rerun `
  --events-h5 .\analysis\26-08-04-derived-hbt-time-reference\detected-peaks-supplement\26-08-04_hbt_time_reference_supplement_events.h5 `
  --validation-tolerance 6
# 迹线 npz 副产物写到 --out-dir;只要迹线、不写事件表时加 --no-events
```

下游消费者:`compute_g2_fig2c.py --table <本文件>`(g2 带窗以 `FIG2C_BAND_WINDOWS`
覆盖为准,见上);与旧表逐位一致,故 `drafts/26-08-06_fig2c_supplement/` 的既有 g2
产物即为本表结果。
