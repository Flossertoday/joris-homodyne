# detected-peaks：26-08-04 HBT 条件逐光子事件表

## 这里面对 DATA 做了什么

输入是 `DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5` 中
`hbt_condensate_josephson_regime` 条件的 20000 张原始 streak 帧。本目录的数据是把这
20000 帧**逐帧做光子事件检测**得到的逐事件坐标表：每帧经实测背景扣除、参考脉冲
CFD50 定时间零点后，把 8-连通亮斑（≥2 px、积分幅度 > 25 counts）压缩为亚像素坐标，
作为下游分析（Fig. 2(c) g2、未来的二维 g2/cross-mesa）的输入数据。

**2026-08-07 起，事件表与 Fig. 2(a) 迹线由 `../extract_hbt_traces.py` 的同一遍检测
产出**（原第二遍脚本 `extract_hbt_event_table.py` 废止删除）：fig2a 累积迹线与事件表
不可能再相互漂移，同遍自洽验证（表重建带 photon 迹线 vs 同遍 npz）在每次全量提取
末尾自动执行。

文件：

- `26-08-04_hbt_time_reference_events.h5`：事件表（275218 事件，19992/20000 帧保留）；
- `26-08-04_hbt_time_reference_events_validation.json`：同遍自洽验证结果；
- `26-08-04_hbt_time_reference_events_sample_overlay_frameNNNNN.png`：5 张抽样
  叠加图，**每图一个种子随机帧**（种子 20260807，从 t0 保留帧中抽取），**全帧
  不裁剪**：事件点叠在扣背景帧上，白带边线（四个数据 ROI 的边界）、点线检测列
  范围 350–1150、虚线 t0 参考区（reference-pulse 估计区，非 ROI）、绿 t0 线（穿过参考脉冲前沿——只有全帧图才能看到参考脉冲列）；图内
  右下角标注 bin size（1 行 = 0.9843 ps，1 列 = 未标定探测器像素）。色标上限
  取检测列范围的 p99.9，参考脉冲饱和显示以保全带光子可见。供目视检查检测
  质量；`--overlays-only` 无需重跑提取即可重生成。

## `26-08-04_hbt_time_reference_events.h5`（事件表）

由 `../extract_hbt_traces.py` 生成（与 fig2a `hbt_traces.npz` 同一遍检测；检测参数与
fig2a 定案逐字一致）：实测 `StreakCamBackG` 背景扣除 → 参考脉冲前沿恒比定时 CFD50
定 t0（`t0_method=cfd50, cfd_fraction=0.5, t0_roi=(110,210,170,215)`）→ 8-连通组分
检测（阈值 25 counts 上背景、min 2 px、仅列 350–1150）→ 计数加权亚像素质心。
阈值由灵敏度扫描拐点自动选定（`threshold_method` 与 `threshold_study_json` 属性）。

每个事件一行，存于 `events/` 组：

- `frame`（int32）：0–19999 帧号；
- `x`, `y`（float64）：亚像素质心（探测器列、绝对行，px）；
- `tau`（float64）：`y - frame_t0[frame]`，抖动校正后的对齐坐标（px）；
- `amplitude`（float64）：组分积分幅度（counts 上背景）。

表级数据集：

- `frame_t0`（float64，20000 行）：逐帧 CFD50 t0 行，t0 失败为 NaN；
- `failed_frames`（int32，8 行）：被排除的帧号
  （3543, 4458, 4589, 6614, 7935, 8190, 9765, 11225），这些帧无任何事件。

根属性记录完整 provenance：源 HDF5 路径与 SHA-256、frames/background 数据集路径、
检测阈值与拐点方法、min_size、检测列范围、t0 方法/比例/ROI、ps_per_px、tau 栅原点、
子窗边、带窗、帧数/保留帧数/失败数/事件数、生成时间与脚本名。

## 横轴 ROI（列范围与带窗）

- **检测时**：`detect_events` 只对列 350–1150（`col_range` 属性）建掩膜，参考脉冲列
  （170–215）与探测器边缘在检测外；全表事件中参考脉冲列计数为 0。所有检出斑的
  **像素**按构造都在该范围内；质心 `x` 实测落在
  [349.99999999999994, 1149.0000000000002]，越界仅 ±2e-13（float64 加权均值的
  边界舍入，12 个事件），无实质越界（判据容差 1e-9）。
- **分析时**：四个发射带带窗存于 `band_windows` 属性（B1 380–450、B2 570–680、
  B3 810–880、B4 1010–1110，与本遍自动带检测结果一致）；带内事件
  B1 31893 / B2 56056 / B3 29063 / B4 48551，带外 109655（g2 的残余背景 β 即由
  带外事件估计）。带指派在下游按 x 进行（`fig2c_g2_estimator.assign_bands`；
  补充批次用 `FIG2C_BAND_WINDOWS` 环境变量覆盖）。

## 脉冲光校准（时间抖动校正）

每个事件的 `tau` 都以**本帧**参考脉冲的 CFD50 穿越行为原点：`tau = y - t0(frame)`，
全表逐事件验证一致（atol 1e-9）；t0 失败的 8 帧被排除且记录在 `failed_frames`。
t0 统计：中位行 154.78，RMS 抖动 13.17 ps（`ps_per_px` 属性 = 0.984252，即
1 ns / 1016 行）。
样品脉冲与参考点延迟实测 0.32 ps（见 `docs/08-04_sync_oscillation.md`）。

## 验证（2026-08-07，合并单遍全量提取，4609 s）

1. **同遍自洽**（`..._validation.json`）：由表重建的带 photon 迹线与同遍 npz 在
   0–820 ps 窗口内最大差 [1,1,1,2] counts/bin（整数计数判据 ≤2，残差为 float64
   求和次序的 bin 边 ties），PASS。
2. **与旧双遍产物逐位一致**：本表与废止脚本第二遍的
   `drafts/26-08-06_fig2c/hbt_event_table.h5`（gitignored，保留作对照）逐列
   （frame/x/y/tau/amplitude）、逐帧 t0、失败帧列表 `array_equal`；本遍 npz 与
   已审计 `drafts/26-08-06_fig2a_cfd50/hbt_traces.npz` 全部数组逐位一致。
3. **g2 复现**：`compute_g2_fig2c.py` 对本表重跑，`g2_summary.json` 全部 66 个
   数值叶与旧结果一致（4 px pooled g2(0)：L 0.3104±0.3795、R 1.1138±0.1884，
   jackknife）。
4. **叠加图目视**：旧 3 帧合版图（100/5000/12000，仅显示检测列范围裁剪区）确认
   事件点落在亮斑上、参考脉冲列无事件；2026-08-07 起改为 5 张全帧独立图
   （种子随机帧 2414/4546/5339/10880/16085），逐张复核以上两项，且"t0 线穿过
   参考脉冲前沿"一项在全帧图上才真实可见（裁剪图看不到参考脉冲列，该检查
   此前无法在此图上执行）。
5. 旧 `drafts/26-08-06_fig2c/event_table_validation.json` 的 `pass:false` 为陈旧
   记录：它写于对账容差放宽（1e-12 → ≤2 计数）之前，放宽后该文件未被重写
   （`--validate-only` 只打印不写盘）；以其自身记录的差值按最终判据折算实为 PASS
   （详见 `docs/26-08-04-fig2c.md` 方法说明第 2 条）。

合成数据单元测试：`tests/test_hbt_event_table.py`（13 项：schema/dtype、拒绝覆盖与
`.partial` 清理、检测列范围、τ 与 frame_t0 一致、失败帧无事件、同遍重建零差、
叠加图种子抽样排除失败帧、逐帧独立 PNG 渲染）。

## 复现（拒绝覆盖既有输出，重跑需先改名或 `--events-h5` 另指路径）

```powershell
$sandboxState = Get-Content .\tmp\sandbox-initialization.json -Raw | ConvertFrom-Json
& $sandboxState.python `
  .\analysis/26-08-04-derived-hbt-time-reference/extract_hbt_traces.py `
  --out-dir .\drafts\<新的 fig2a 输出目录>
# 事件表默认写到本目录 <h5 主名>_events.h5；只要迹线、不写事件表时加 --no-events
# 更密数据集（如 72000 帧补充批次）按事件数比例放宽同遍判据：--validation-tolerance 6
# 只重生成抽样叠加图（不重跑提取）：--overlays-only（另可选 --overlay-count/--overlay-seed）
```

下游消费者：`compute_g2_fig2c.py`（`DEFAULT_TABLE` 即本文件）；`fig2c_g2_estimator.py`
的带窗常量/环境变量覆盖与本表 `band_windows` 属性独立，按 x 指派。
