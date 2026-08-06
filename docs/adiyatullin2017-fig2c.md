# Adiyatullin 2017 Fig. 2(c)：论文要点存档与复刻目标分解

本文件是论文侧事实的存档，供后续 Fig. 2(c) 相关任务直接引用，不必再回读论文 PDF。
文献标识、样品/激发/探测参数与 Fig. 2 整体布局见 [adiyatullin2017-fig2a.md](adiyatullin2017-fig2a.md)
§1–§3，本文件只补充 Fig. 2(c) 特有内容并给出目标分解。当前状态：**已完成（降级为上限报告，2026-08-06）**——G3 功效门槛判定 borderline（只能覆盖论文幅度区间上半），经用户确认按 §5 G3 的降级路径执行完毕；四项预登记检验全部零结果；管线、判定与复现命令见 [26-08-04-fig2c.md](26-08-04-fig2c.md)。

## 1. Fig. 2(c) 的目标规格（复刻要复现什么）

图注原文："Second-order correlation function g⁽²⁾(0) of the light emission from the left
(blue shaded) and right (red shaded) mesa shows that the light statistics changes from
Poissonian (g⁽²⁾(0) = 1) to bunched (g⁽²⁾(0) > 1) in phase with Josephson oscillations.
Shaded areas represent the error bars calculated as the standard errors of the mean values.
The corresponding results of the theoretical simulations are shown with blue and red lines."

- **内容**：左（蓝色阴影带）、右（红色阴影带）两个 mesa 各自的时间分辨零延迟二阶关联
  g⁽²⁾(0)(t)。阴影带 = ±1 SEM（standard error of the mean，均值标准误）；蓝/红实线 =
  理论模拟（复刻不要求做模拟，但要知道论文图上压着模拟线）。
- **坐标**：纵轴 g⁽²⁾(0)，线性刻度，标出 1 和 2 两档；横轴与 Fig. 2 其余面板共用 0–30 ps。
- **物理特征（判定复刻成败的核心，读自原图与正文）**：
  1. 两臂 g⁽²⁾(0) 从 ≈1（Poisson）出发，随时间出现一个个 bunching 峰（g⁽²⁾(0) > 1）；
  2. 峰的节奏与 Josephson 人口振荡相同（周期约 10 ps，对照 Fig. 2(b) z(t) 峰位
     约 3.5 / 13 / 23 ps）——正文："changes from Poissonian to bunched in phase with
     the Josephson oscillations"；
  3. 两臂**反相**："when the emission of the left mesa is bunched, the right mesa emits
     Poissonian light, and vice versa"；
  4. 峰高随时间**增大**："the magnitude of the bunching gets higher as the polariton
     population decreases"——原图上峰包络从早期 ≈1.2–1.5 升到 25 ps 之后 ≳2；
  5. 误差带随时间变宽（计数衰减，统计变差）。
- **定义与判读阈值**（正文）：g⁽²⁾(t₁,t₂) 为标准正规序定义（论文 Eq. 1）；
  g⁽²⁾(0)(t) ≡ g⁽²⁾(t, t)。g⁽²⁾(0) = 1 为 Poisson，>1 为 super-Poisson（bunched），
  <1 为 antibunched（量子）。
- **配套事实**（解释与界限，来自正文 Discussion 与 SI Note 2）：
  - 大占据数极限下 1 − 1/ᾱ² < g⁽²⁾(0) ≤ 3；论文实验中 t < 30 ps 时上限 3 未达到，
    峰高增长机制是总占据数 n_tot 下降（Eq. 3 分母变小）；g⁽²⁾(0) > 3 的 superbunching
    只在 n_tot ≈ 1 时才可能。
  - t = 30 ps 时系统占据数约 100；可探测的最弱光对应约 50 个 polariton，该处最小
    可分辨 g⁽²⁾(0) ≈ 0.98——**antibunching 在他们的装置上不可宣称**，Fig. 2(c) 全程
    是 Poisson↔bunching 之间的切换。
  - bunching 峰对应 squeezing 相位关系 cos(θ−2φ) < 0（论文 Eq. 3 与 Fig. 2e）；Fig. 2(d)(e)
    的 squeezing 参数曲线是**模型派生**，不是直接测量（见
    [measurement-claims.md](measurement-claims.md) §4：强度 bunching 本身不是 quadrature
    squeezing 的测量）。复刻 Fig. 2(c) = 复刻强度关联测量，不得引申为 squeezing 见证。
- **与 Fig. 3 的关系**：Fig. 3(a) 是左 mesa 的二维 g⁽²⁾(t₁,t₂) 图，Fig. 2(c) 即其对角线
  g⁽²⁾(t,t)。二维图是后续自然延伸，不在本任务范围。

## 2. 论文给出的复刻相关事实（方法，SI Note 1）

- 估计量：g⁽²⁾(t₁,t₂) = ⟨n₁(t₁)·n₂(t₂)⟩ / (⟨n₁(t₁)⟩⟨n₂(t₂)⟩)，其中 n₁、n₂ 是 HBT
  （Hanbury Brown–Twiss，50:50 分束）两臂各自探测到的光子流，平均对激光脉冲数进行。
  对 Fig. 2(c) 的某 mesa，n₁、n₂ = **该 mesa** 落在两臂中的光子数（arm 分辨）；
  Fig. 2(a) 的强度则是臂求和 I_L,R = n_L + n_R。
- 时间分箱：每脉冲每臂 n(t) = 落在 (t−Δt/2, t+Δt/2) 内的光子数；**Δt = 1.296 ps =
  其 streak 屏 2 像素**（其换算 0.648 ps/px；增大 Δt 提高每 bin 计数、降低时间分辨率）。
  30 ps 窗口 ≈ 23 个 bin。
- 误差：g⁽²⁾(0) 的误差由 ⟨n₁⟩、⟨n₂⟩、⟨n₁·n₂⟩ 各自的标准误（SEM）推得
  （具体传播公式未写）。
- 统计规模：>3.9×10⁶ 个激光脉冲（34 h）。Fig. 2(a) 的 "SNR<20" 灰带 = 关联测量不可信的
  低计数区；Fig. 2(c) 的有效数据因此基本止于 ~30 ps。
- 光子计数：自定义算法，CCD 检测阈值"足够高以分辨近距光子对"（数值未公开）。
- 时间抖动：分一小束激发激光做时间参考校正逐帧 jitter（算法未公开）。
- **臂效率无需标定**：n_i = ε_i·m 为同一光子流的独立稀释时，ε₁ε₂ 在
  g⁽²⁾ = ⟨m(m−1)⟩/⟨m⟩² 的比值构造中精确相消——正规序 + 独立臂估计量对两臂效率
  不对称免疫（与 Fig. 2(a) 强度求和的抵消机制不同，这里是比值抵消）。

## 3. 论文未写明、本项目必须自行决定的方法点

| 方法点 | 论文状态 | 本项目的待定项 | 备注 |
|---|---|---|---|
| Δt 在本项目像素尺度 | 其 2 px = 1.296 ps（0.648 ps/px） | 本项目 0.984 ps/px：1.296 ps ≈ 1.32 px 非整数；候选 1/2/4 px，做分箱宽度灵敏度扫描 | 26-07-31 已有多宽度分箱（1/3/6 px）先例 |
| 逐帧 t0 对齐与分箱的先后 | 未写（jitter 校正后计数） | 先用 CFD50 t0 把每帧事件平移到公共时间栅（亚行），再落入 bin | t0 失败帧（8/20000）排除并记录 |
| 事件位置精度 | 未写 | 连通组分质心（亚像素）或整数像素；平移后按 (y−t0) 分箱 | 与 fig2a 管线事件检测一致 |
| 背景处理 | 未写 | 检测前逐帧扣实测 `StreakCamBackG`；是否再做残余背景 β 修正（26-07-31 做法）由数据定 | fig2a G2/G5 |
| SEM 传播公式 | "由三个 SEM 推得"，无公式 | 三项误差传播（图级）+ leave-one-frame-out jackknife（汇总量）对照 | 26-07-31 `compute_g2_matrix.py` 已实现两种 |
| 统计不足区的标示 | SNR<20 灰带（SNR 算法未给） | 以每 bin SEM(g⁽²⁾) 阈值标示/截断（26-07-31 用 MAX_DISPLAY_SEM=0.5） | 语义须写清 |
| 物理区间 | 共振脉冲激发、30 ps、z(0)≈−0.66 | 本项目非共振 reservoir-fed + Stark kick、1 ns 窗口；复刻目标是**测量管线与观测量形式**，不是同一物理状态 | 同 fig2a G12 保留项 |

## 4. 本项目现有资产（复用，勿重造）

- `analysis/26-07-31-derived-single-sweep/compute_g2_matrix.py`：同一 SI Note 1 估计量的
  完整实现——逐帧逐类逐 bin 计数立方、背景 β 修正、三项 SEM 传播 + jackknife、多宽度分箱、
  统计不足 bin 标示；配套单元测试 `tests/test_compute_g2_matrix.py`。
- 26-07-31 先验结果：1000 帧、全窗口汇总 g⁽²⁾(0) 全部与 Poisson 一致（~1.5σ 内），
  same-mesa 值偏高 +0.05–0.08 不构成探测——**统计受限的零结果**。这是 G3 功效门槛
  必须先行的直接动因。
- 26-08-04 管线（fig2a G1–G6 定案）：派生 HDF5（20000 帧主条件）；实测背景扣除；
  CFD50 亚行 t0（偏差 −0.10 px / std 0.33 px）；光子事件检测参数（8-连通、25 counts、
  min 2 px）；四条带 B1 x=380–450、B2 570–680、B3 810–880、B4 1010–1110，臂配对
  {B1,B2}/{B3,B4}，**L 的两条 HBT 臂 = B1、B3；R 的两条臂 = B2、B4**；时间换算
  0.984 ps/px。累计迹线 `drafts/26-08-06_fig2a_cfd50/hbt_traces.npz` 供交叉校验。
- Fig. 2 第 3 页原图（仓内 PDF，300 dpi 渲染件 `tmp/adiyatullin2017_page3-3.png`，
  Git 忽略），fig2a 裁剪脚本 `tmp/paper_crop/crop_paper_fig2a.py` 可照搬出 2(c) 裁剪。

## 5. 复刻目标分解：可完成、可验证的细则

每条给出验证标准。状态（2026-08-06）：全部完成（降级为上限报告）——G3 borderline；
G1/G2/G4/G5/G7/G8 验收通过；G6 四项检验均零结果。逐条结论与复现命令见
[26-08-04-fig2c.md](26-08-04-fig2c.md)。

- **G1 光子事件表提取** —— 对 20000 帧 HBT 主条件重跑检测并**保存逐光子事件表**
  （帧号、亚像素 x/y、幅度），复用 fig2a 的背景扣除与检测参数（fig2a 只存了累计迹线，
  未存事件表）；同时存每帧 CFD50 t0。
  验证：事件表按带汇总重建的 photon 迹线与 `hbt_traces.npz` photon 源逐 bin 一致
  （数值容差内）；报告总事件数与四带分布；抽样帧事件叠加目检。
- **G2 对齐分箱计数立方** —— 事件按 (y−t0) 平移到公共时间栅后分箱，得
  （帧 × 4 带 × bin）计数立方；默认 Δt = 2 px（1.97 ps），另扫 1/4 px。
  验证：立方沿帧求和 = 对应带迹线；无信号区 bin 计数率与背景预期一致；t0 失败帧
  排除记录完整。
- **G3 统计功效门槛（go/no-go，先于全量计算）** —— 用实测计数率推算每 bin 占据数
  ⟨n_c(t)⟩、符合率与 SEM(g⁽²⁾)(t)；给出本数据集可分辨的最小 g⁽²⁾−1 振荡振幅
  （95% 置信，周期 3–200 ps 扫描），与论文特征幅度（峰 1.2→≳2，即 g⁽²⁾−1 ≈ 0.2–1.2）
  对比。判定：若最小可分辨振幅盖不住论文量级，写明所需帧数、把目标降级为上限报告，
  并在此停一下与用户确认；若盖得住，继续。
  验证：解析估算 + 合成 Poisson 数据的 Monte Carlo 复现预测 SEM（偏差 <10%）；
  报告含关键中间数字（每 bin 计数、占据数）。
- **G4 g⁽²⁾(0)(t) 计算** —— same-mesa 臂间关联 g⁽²⁾_LL(0)(t)（B1×B3）与
  g⁽²⁾_RR(0)(t)（B2×B4）+ SEM 带（三项传播，jackknife 对照）；估计器适配
  26-07-31 模块并补单元测试。
  验证：合成 Poisson 数据上 g⁽²⁾≡1 无系统偏置（χ²/ndf≈1，|mean−1| < 2σ）；
  新代码测试通过、既有测试不回归。
- **G5 论文式渲染** —— 蓝=L、红=R 阴影 ±1 SEM 带、线性纵轴、g⁽²⁾=1 参考线；
  统计不足 bin 按 G3 阈值标示或截断（语义写入图注）；论文原图 2(c) 裁剪件与复刻面板
  并排（裁剪脚本照 fig2a 例入 `tmp/paper_crop/`，裁剪件入 `analysis/`）。
  验证：按仓库惯例自审输出图（轴量、刻度、可读性、面板数正确）。
- **G6 科学判定（预先登记判据）** —— 四项检验，每项给定量结论：
  (a) bunching 存在性：是否存在 g⁽²⁾−1 > 0 显著的时间 bin（逐 bin z 分数，以置换检验
  控制多重性）；
  (b) 振荡：与参考脉冲锁相的 g⁽²⁾ 正弦振幅扫描（周期 3–200 ps），显著则报振幅/周期，
  否则给 95% 上限（fig2a G9 同法）；
  (c) 反相：g⁽²⁾_L 与 g⁽²⁾_R 去趋势涨落的互相关符号检验（反相应为显著负相关）；
  (d) 包络增长：峰包络随时间上升的单调性/早晚窗对比。
  验证：判据在跑数据前写入执行文档；结论为阳性（数值+p 值）或零结果（上限）。
- **G7 灵敏度验证** —— 合成注入论文式 g⁽²⁾ 调制（峰值 1.5–2、周期 9.84 ps、两臂反相、
  包络随时间增长、真实背景与 t0 漂移）经同一管线处理。
  验证：恢复周期与注入值差 <2%、调制振幅恢复偏差 <10%；证明若结论为零结果则非管线
  盲区（fig2a G10 的 g⁽²⁾ 版）。
- **G8 定稿与执行文档** —— 对比合成图（论文裁剪 vs 复刻）入
  `analysis/26-08-04-derived-hbt-time-reference/`；新建执行文档
  `docs/26-08-04-fig2c.md` 记录管线、判定、复现命令；中间产物入 `drafts/`（Git 忽略）。
  验证：文档内命令可复跑；`docs/README.md` 与 `MEMORY.md` 索引更新。
- **G9 边界与缺口（写入执行文档，不另行验证）** —
  1. 本数据集人口振荡本身是零结果（fig2a G9）：Fig. 2(c) "in phase with Josephson
     oscillations" 的锁相判定只能以参考脉冲相位为基准，不能对照实测 z(t)；
  2. g⁽²⁾(0) > 1 不等于 quadrature squeezing（measurement-claims.md §4）；论文
     Fig. 2(d)(e) 的 squeezing 曲线为模型派生，不在复刻范围；
  3. 物理区间差异同 fig2a G12 第 3 条（共振 vs 非共振 + Stark kick）；
  4. 二维 g⁽²⁾(t₁,t₂)（论文 Fig. 3）与 cross-mesa 关联（LR/RL）是后续任务入口，
     不在本分解范围。

## 6. 一句话现状

已完成（2026-08-06，上限报告）：G3 门槛判定 borderline（逐 bin 全程统计受限；相干振荡
95% 灵敏度中位 0.78，只盖住论文幅度区间 0.2–1.2 的上半），经用户确认降级执行；
G1 事件表（275218 事件）→ G2 立方 → G4 g2 → G5 渲染 → G6 四项预登记检验全部零结果
（pooled g⁽²⁾(0)：L 0.31±0.38、R 1.11±0.19；振荡上限中位 0.81，9.84 ps 处 1.09）→
G7 注入验证通过 → 文档定稿于 [26-08-04-fig2c.md](26-08-04-fig2c.md)。盖住论文下限
0.2 需 ~3×10⁵ 帧；二维 g⁽²⁾(t₁,t₂) 与 cross-mesa 关联为后续任务入口（事件表与计数
立方已保留）。
