# 26-08-04 HBT 补充批次（72000 帧）g2 与 intensity 分析

目标：对补充数据集 `DATA/26-08-04-hbt-time-reference/New folder (10)/New folder (10)/`
（7200 个 MAT 文件、编号 1–7200 连续、每文件 10 帧，共 72000 帧）用与
[26-08-04-fig2a.md](26-08-04-fig2a.md) / [26-08-04-fig2c.md](26-08-04-fig2c.md) 相同的管线做
intensity（Fig. 2(a) 式）与 g⁽²⁾(0)(t)（Fig. 2(c) 式）分析。本文件为执行文档；
论文侧事实与判定框架见 [adiyatullin2017-fig2a.md](adiyatullin2017-fig2a.md)、
[adiyatullin2017-fig2c.md](adiyatullin2017-fig2c.md)。

## 1. 采集条件确认（2026-08-06，先于一切分析）

逐字段对比补充文件（`1_image_savedStructure.mat`）与原批次
`01_HBT_condensate_Josephson_regime/2_image_savedStructure.mat` 的 `structVar`：

- 字段集合完全相同；标量字段仅 2 处不同：
  - `SweepNbStepsVar1`：7200 vs 2000 —— 即扫描长度本身（7200 文件 vs 2000 文件）；
  - `lastBufferIndexStreakCam`：1 vs 6 —— 缓冲记账字段，非采集参数。
- 关键参数全部一致：`ZeroOrderPowerRatio` 0.26、`StreakRange` 5（1 ns 窗）、
  `StreakCamImage` (10, 1016, 1344) uint16、streak 曝光/延迟、SLM 图样参数。
- 补充文件不含 `StreakCamBackG`（原批次也只有条件 09 文件含）。原批次实测背景
  对补充数据的适用性验证：补充 50 帧（跨全程 5 个采样点）逐像素中位数与实测背景
  在非信号区差 +0.02 counts（背景 ~105 counts/px），与原批次自身对照（+0.07）同水平
  → **背景沿用原批次实测值**。

结论：采集条件与原 HBT 条件一致；补充批次是同一条件的 3.6 倍加长采集。
**但图像内容有差别**（metadata 不记录的物理漂移）：补充批次带内信号约为原批次的
65%（photon 4.6 vs 7.0 events/帧·带），L/R 比 0.55 vs 0.578，因此本批次按
**独立数据集**分析，不与原 20000 帧合并。

## 2. 派生 HDF5

`analysis/26-08-04-derived-hbt-time-reference/convert_hbt_supplement.py`（本 worktree
新增）生成 `DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference_supplement.h5`：

- schema-2 布局；`conditions/hbt_condensate_josephson_regime/frames`
  uint16[72000,1016,1344]（gzip-4、shuffle、fletcher32、逐文件 provenance + SHA-256）；
- `conditions/background_and_streak_pulses_delay_calibration/background`
  自原批次 HDF5 复制（属性中注明来源与 §1 的验证）；
- 分层逐帧校验：每 29 个文件抽首末帧共 498 帧与源 MAT 逐像素比对，0 不符。

## 3. 管线适配（仅两处，均不改逻辑）

1. `extract_hbt_traces.py` / `extract_hbt_event_table.py`：`H5_SHA256` 硬编码常量改为
   `resolve_h5_sha256()`（默认文件用常量，其他文件现算），溯源元数据修正。
2. `fig2c_g2_estimator.py`：`BAND_WINDOWS` 支持环境变量 `FIG2C_BAND_WINDOWS` 覆盖
   （默认原值）。补充批次带窗与原批次不同（见 §4），g2 各脚本统一以
   `FIG2C_BAND_WINDOWS="350,440,540,680,770,860,970,1100"` 运行。

相关单元测试：适配改动后 39 项（test_extract_hbt_traces / test_fig2c_pipeline /
test_fig2c_g2_estimator）通过；收尾时全量 118 项通过（工作树自带 .venv）。

## 4. 带定案（补充批次）

fig2a 管线的自动带检测在补充数据上失效：带外存在结构化辉光，子窗计数的
baseline+5σ 阈值在 72000 帧统计量下只抓到峰顶子窗（全量运行仅检出
(550,660)/(790,850)/(980,1080) 三条部分带）。改用全量对齐子窗廓线（71971 帧，
analog 与 photon 双源一致）手动定案：

- **B1 = 350–440、B2 = 540–680、B3 = 770–860、B4 = 970–1100**（10 px 网格；
  中心 395/610/815/1035；mesa 间距 215/220、臂间距 420/425，与原批次
  415/625/845/1060（mesa 间距 210/215、臂间距 430/435）同几何、整体平移 ≤20 px）。
- 配对不变：L = {B1, B3}、R = {B2, B4}。臂内 L/R 比：analog 0.535/0.556、
  photon 0.543/0.571，两臂互洽（交叉配对给不出一致比），配对成立；
  窗口 L/R ≈ 0.55（原批次 0.578）。
- 静态遮挡斑（两批次共有的探测器尘埃，cols ~760–790 × rows ~300–470）部分落入
  补充 B3 的 770–800 段：该段 analog 计数在 300–470 ps 被压低于基线；photon 廓线
  平滑，影响为 B3 单臂局部少量信号损失，不剔除（原管线同样未剔除）。
- 带重聚合：`analysis/26-08-04-derived-hbt-time-reference/reaggregate_supplement_bands.py`（复用
  `extract_hbt_traces.aggregate_band_traces` 的 mean/SEM 语义）从全量 npz 的
  80 子窗迹线重聚合 → `drafts/26-08-06_fig2a_supplement/hbt_traces_bands.npz`，
  下游渲染/检验脚本原样消费（`--bands-l 0,2 --bands-r 1,3`）。

## 5. Fig. 2(a) 式 intensity 分析（71971 帧，29 帧 t0 失败 = 0.04%）

管线：CFD50 t0、实测背景扣除、事件检测（8-连通、25 counts、min 2 px）、
子窗对齐累计（全量提取 94 分钟）→ 带重聚合 → `render_fig2a.py`
（photon 源、rebin 8 ≈ 7.9 ps/bin、0–820 ps、峰值归一）。

产物：`drafts/26-08-06_fig2a_supplement/`：`hbt_traces.npz`（子窗级）、
`hbt_traces_bands.npz`（带级）、`fig2a_replica_supplement.png`（+`_matched.png`）、
`z_crossarm_check.png`、`z_sinusoid_upper_limit.png`、`aligned_sum.png`。

- 时序形态：I_R 平稳；I_L 早期（t≲60–80 ps）浅缺口后恢复（同原批次的过阻尼
  弛豫形态，更浅）；全程无衰减包络（reservoir-fed），两臂无镜像反相摆动。
- **两臂交叉校验（z 去趋势涨落）**：r = 0.068（rebin 2，置换 p = 0.155）、
  0.111（rebin 4，p = 0.095）——不显著（原批次 0.030/0.095）。
- **锁相正弦振幅 95% 上限（3–200 ps 扫描，合并臂）**：中位 **0.0084**、最差 0.0142
  （原批次 0.0123/0.0193；灵敏度提升受限于补充信号率 ~65%，不及 √3.6）。
  先验周期针对性拟合（`hbt_z_targeted_fits.py`，2026-08-07 起周期按 2 ns 窗口取
  74.4/83.7 ps）：9.84 ps |A| = 0.0023（UL 0.0068）、74.4 ps 0.0062（0.0107）、
  83.7 ps 0.0017（0.0061）。
- **结论：与原批次一致的零结果，上限收紧到 ~1% 强度调制。**

## 6. Fig. 2(c) 式 g⁽²⁾(0)(t) 分析

管线：事件表（73 分钟，**661090 事件**、71971 帧、29 帧 t0 失败；
重建校验：与带级 npz 逐 bin 差 14–39 counts/带 = 带事件的 0.015–0.044%，
与原批次同因（float64 bin 边 ties，原 2 计数容差按事件数 ×2.4 放大后即通过））
→ 计数立方（1/2/4 px）→ g⁽²⁾(0)(t)（三项传播 + plug-in SEM、jackknife pooled）。

定稿图（入 `analysis/26-08-04-derived-hbt-time-reference/`）：
`fig2a_replica_supplement.png`/`.pdf`（intensity 复刻面板）、
`fig2c_replica_supplement_bin4px.png`（g⁽²⁾(0)(t) 4 px 复刻面板）、
`aligned_sum_supplement.png`（71971 帧 t0 对齐总和图，带定案证据）。

产物：`drafts/26-08-06_fig2c_supplement/`：`hbt_event_table.h5`、
`count_cube_bin{1,2,4}px.npy`、`g2_results_bin{1,2,4}px.npz`、`g2_summary.json`、
`fig2c_replica_bin4px.png`、`power_gate_report.json` + 5 图、
`fig2c_science_report.json` + 两宽度判定图、`fig2c_injection_report.json`、
`event_table_sample_overlay.png`。

- **窗口汇总 g⁽²⁾(0)（0–820 ps，jackknife SEM）**：
  - 4 px：L **1.304 ± 0.320**、R **0.990 ± 0.117**（z = +0.95 / −0.08，均与 Poisson 一致）；
  - 2 px：L 1.900 ± 0.476（+1.9σ）、R 0.901 ± 0.163；
  - 1 px：L 1.799 ± 0.662、R 1.029 ± 0.234。
  - 按预登记判据 |z| ≥ 2 才算偏离 Poisson：任何宽度均未达到。
  （原批次 4 px：L 0.31 ± 0.38、R 1.11 ± 0.19。补充 L 臂 2 px 汇总值 +1.9σ 为
  宽度间相关涨落——G6 surrogate 标定下 (a)(b) 两检验均不显著，不构成探测。）
- 逐 bin 显示：与原批次同样全程统计受限（4 px plug-in SEM 中位 L 4.7 / R 1.7，
  原批次 12 / 5.5；提升 ~2× 但仍远超 0.5 展示阈值）。

### G3 功效门槛（补充批次重跑）

`fig2c_power_gate.py --reps 64`（`power_gate_report.json` + 5 图，同目录）：

- 逐 bin plug-in SEM 中位：4 px L 4.67 / R 1.73；2 px L 9.41 / R 3.42；1 px L 19.23 / R 6.83
  —— 逐 bin 曲线任何宽度全程统计受限（同原批次）。
- 窗口汇总 SEM：4 px L 0.314、R 0.117（与实测 jackknife 0.320/0.117 吻合）。
- 锁相正弦振幅 95% 最小可分辨（合并反相）：4 px 中位 **0.52**、9.84 ps 处 0.66
  （原批次 0.78 / 1.01；提升 ~1.5×，受补充信号率 ~65% 限制，不及 √3.6）。
- MC 验证（64 rep）：4 px 与 2 px 两对 UL95 解析/MC ∈ [1.000, 1.024] 全部通过；
  1 px L 不通过（重尾，同原批次）；pooled jackknife/MC 0.91–1.04。
- **判定：borderline**——论文幅度区间 0.2–1.2 仍只盖住上半；盖住 |A| = 0.2 需
  ~(0.52/0.2)² ≈ 6.8× 本批次帧数（≈4.9×10⁵ 帧）。用户已指示"用同样的管线"，
  即沿用原批次 borderline 后的降级路径：交付上限报告并跑完全管线（G6/G7）。

### G6 四项预登记检验（200 surrogates）

判据同 [26-08-04-fig2c.md](26-08-04-fig2c.md) §G6（2026-08-06 预登记）；零假设
surrogate：保持逐帧速率与 t0 结构、循环平移一臂打破同帧符合；窗口 0–820 ps；
4 px 为主（MC 验证通过）、2 px 对照。`fig2c_science_report.json` + 两宽度判定图。

**全部四项检验零结果（与原批次一致）：**

- **(a) bunching 存在性**（逐 bin max|z|）：4 px L p = 0.592、R p = 0.303；
  2 px L p = 0.746、R p = 0.393——无显著 bin。
- **(b) 振荡**（3–200 ps 双正交扫描 + 反相合并，max |A|/σ 置换判定）：
  4 px L p = 0.075、R p = 0.373、合并 p = 0.612；2 px 合并 p = 0.065——均不显著
  （两个 0.07 上下的边缘值在不同宽度/不同臂，不满足 p < 0.05，且互相不印证）。
  95% 上限（|A|+2σ，反相合并）：4 px 中位 **0.50**；9.84 ps 处 0.58、
  74.4 ps 处 0.64、83.7 ps 处 0.54（2 px：0.64；0.56 / 1.10 / 0.62）
  （74.4/83.7 ps 为 2026-08-07 按 2 ns 窗口更正的 k 空间实测周期；重跑后非针对性
  统计量与原报告逐位一致）。
  实测灵敏度与 G3 预测（0.52）吻合。针对性拟合振幅 4 px：9.84 ps A = 0.185 ± 0.200、
  74.4 ps A = 0.339 ± 0.151、83.7 ps A = 0.232 ± 0.152——74.4 ps 处 +2.2σ 的
  单点涨落，经多重性置换（扫描内 max）不显著。
- **(c) 反相**：去趋势涨落零延迟相关 4 px r = +0.033（p = 0.667）、
  2 px r = +0.014（p = 0.622）——无反相证据（方向与符号均不支持）。
- **(d) 包络增长**（晚窗−早窗 pooled g⁽²⁾−1）：4 px L −0.196（p = 0.517）、
  R +0.250（p = 0.154）；2 px L −0.983（p = 0.881）、R −0.122（p = 0.657）——不显著。
- **汇总测量**：见 §6 开头；任何宽度 |z| < 2，与 Poisson 一致。

### G7 合成注入验证

`fig2c_synthetic_injection.py`（立方级注入，实测速率+残余背景，9.84 ps 反相、
增长包络；`fig2c_injection_report.json`）：

- 强注入（峰值 6→8，16 rep）：集合振幅偏差 L +0.2%、R −4.1%（判据 <10%）；
  周期恢复 2% 内 L 15/16、R 16/16；全部 rep 显著。**通过。**
- 论文强度（峰值 1.5→2，16 rep）：集合恢复无偏（L −7.4%、R +4.3%，集合 SEM 内）；
  单 rep 检出率（|z|/σ > 2.45）L 19%、R **81%**（原批次 L 19%、R 50%）——
  G3 "borderline" 的直接实证：论文幅度在本批次 R 臂大多可检出，L 臂仍只有部分。

## 7. 边界（沿用 fig2c 文档 G9，补充批次同适用）

1. 人口振荡零结果在两批次各自独立成立；"与 Josephson 振荡同相"只能以参考脉冲
   相位为基准；
2. g⁽²⁾(0) > 1 ≠ quadrature squeezing（measurement-claims.md §4）；
3. 物理区间差异：非共振 reservoir-fed + Stark kick、1 ns 窗口，复刻的是
   测量管线与观测量形式；
4. 两批次物理状态差异（信号率、L/R 比）未合并分析；如需联合上限，按
   逆方差加权需另行论证状态一致性；
5. 补充事件表与计数立方已保留，二维 g⁽²⁾(t₁,t₂) 与 cross-mesa 关联入口同前。

## 复现

```bash
# worktree: kimi/hbt-supplement-26-08-06
# 1) 转换（生成补充 HDF5；拒绝覆盖，先删旧文件）
python -u analysis/26-08-04-derived-hbt-time-reference/convert_hbt_supplement.py \
  "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-hbt-time-reference/New folder (10)" \
  "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference.h5" \
  "E:/Documents/Codex/joris-homodyne/DATA/26-08-04-derived-hbt-time-reference/26-08-04_hbt_time_reference_supplement.h5"
# 2) fig2a 迹线（~94 min）+ 带重聚合 + 渲染
python -u analysis/26-08-04-derived-hbt-time-reference/extract_hbt_traces.py \
  --h5 .../26-08-04_hbt_time_reference_supplement.h5 --out-dir drafts/26-08-06_fig2a_supplement
python -u analysis/26-08-04-derived-hbt-time-reference/reaggregate_supplement_bands.py
python -u analysis/26-08-04-derived-hbt-time-reference/render_fig2a.py \
  --npz drafts/26-08-06_fig2a_supplement/hbt_traces_bands.npz \
  --bands-l 0,2 --bands-r 1,3 --rebin 8 --t-min 0 --t-max 820 --snr-threshold 0 \
  --out drafts/26-08-06_fig2a_supplement/fig2a_replica_supplement.png
# z 交叉臂 / 正弦上限
python -u tmp/fig2a_z_crossarm_cfd50_supplement.py
python -u tmp/fig2a_upper_limit_cfd50_supplement.py
# 3) fig2c 事件表（~73 min）+ 立方 + g2 + 渲染（带窗覆盖）
#    事件表末尾自动对照 hbt_traces_bands.npz 重建校验（容差按事件数 ×2.4 放宽）
export FIG2C_BAND_WINDOWS="350,440,540,680,770,860,970,1100"
python -u analysis/26-08-04-derived-hbt-time-reference/extract_hbt_event_table.py \
  --h5 .../26-08-04_hbt_time_reference_supplement.h5 \
  --out drafts/26-08-06_fig2c_supplement/hbt_event_table.h5 \
  --npz drafts/26-08-06_fig2a_supplement/hbt_traces_bands.npz
python -u analysis/26-08-04-derived-hbt-time-reference/compute_g2_fig2c.py \
  --table drafts/26-08-06_fig2c_supplement/hbt_event_table.h5 --out-dir drafts/26-08-06_fig2c_supplement
python -u analysis/26-08-04-derived-hbt-time-reference/render_fig2c.py --width 4 --out-dir drafts/26-08-06_fig2c_supplement
# 4) 功效门槛 / 科学判定 / 注入
python -u analysis/26-08-04-derived-hbt-time-reference/fig2c_power_gate.py \
  --npz drafts/26-08-06_fig2a_supplement/hbt_traces_bands.npz --out-dir drafts/26-08-06_fig2c_supplement --reps 64
python -u analysis/26-08-04-derived-hbt-time-reference/fig2c_science_tests.py --out-dir drafts/26-08-06_fig2c_supplement --surrogates 200
python -u analysis/26-08-04-derived-hbt-time-reference/fig2c_synthetic_injection.py --out-dir drafts/26-08-06_fig2c_supplement
```
