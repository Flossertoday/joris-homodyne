# Adiyatullin 2017 Fig. 2(a)：论文要点存档与复刻目标分解

本文件是论文侧事实的存档，供后续 Fig. 2(a) 相关任务直接引用，不必再回读论文 PDF。
执行侧（本项目管线的实现、运行结果、复现命令）记录在 [26-08-04-fig2a.md](26-08-04-fig2a.md)，
定时方法论证在 [08-04_sync_oscillation.md](08-04_sync_oscillation.md)。本文件只重复其中与目标判定
有关的结论数字，详情以那两份文件为准。

## 1. 文献标识与本地副本

- Adiyatullin, Anderson, Flayac, Portella-Oberli, Jabeen, Ouellet-Plamondon, Sallen, Deveaud,
  "Periodic squeezing in a polariton Josephson junction", *Nature Communications* **8**:1329 (2017),
  DOI 10.1038/s41467-017-01331-8，开放获取（CC BY 4.0）。
- 仓内 PDF：`analysis/Adiyatullin et al. - 2017 - Periodic squeezing in a polariton Josephson junction.pdf`；
  Fig. 2 在第 3 页，正文相关段落在第 2–5 页（Results / Methods），补充信息（SI）为在线附件。
- 文本提取件（Git 忽略，缺失时从 PDF 重新提取）：`tmp/adiyatullin2017.txt`（主文）、
  `tmp/adiyatullin2017_si.txt`（SI，含 Supplementary Notes 1–2 与 Supplementary Figs. 1–3）。
- Fig. 2(a) 原图裁剪件已入库：`analysis/26-08-04-derived-hbt-time-reference/fig2a_paper_original.png`
  （300 dpi 裁自 PDF 第 3 页）。

## 2. Fig. 2(a) 的目标规格（复刻要复现什么）

图注原文："Measured intensity of the light emission from the left (blue) and right (red) mesa.
The gray area indicates the region where signal-to-noise ratio (SNR) is insufficient for confident
correlation measurements."

- **内容**：左右两个 mesa（耦合 polariton 凝聚体）的累积发光强度 I_L(t)（蓝点）、I_R(t)（红点）
  随时间的变化。是对 >3.9×10⁶ 个激光脉冲的光子计数求和后的平均发射动力学。
- **坐标**：纵轴 Intensity (a.u.)，对数刻度，标出 10⁰、10⁻¹、10⁻² 三档；任意单位，数据峰值约
  0.5–0.6（未归一到 1）。横轴与 Fig. 2(b–e) 共用 0–30 ps（刻度标注只在面板 e 底部 "Time (ps)"）。
- **物理特征（判定复刻成败的核心）**：两臂强度**反相交替振荡**——L 的峰对 R 的谷；周期约 10 ps
  （配套面板 b 的失衡曲线峰位约在 3.5 / 13 / 23 ps）；包络指数式衰减，30 ps 内下降约两个数量级，
  与 polariton 寿命 τ = 5.2 ps 一致；早期 I_L > I_R（激光聚焦主要照在一个 mesa 上，约 1.5:1）。
- **底部灰色横带**：低计数区，标注 "SNR<20"，表示该区域的信噪比（SNR）不足以做有信心的**关联**
  测量（针对 g² 而言；强度图本身照常画出）。论文与 SI 均未给出 SNR 的具体算法。
- **配套面板**（复刻 2(a) 不必重做，但物理判定要对照）：(b) 人口失衡 z(t)=(I_L−I_R)/(I_L+I_R)，
  品红测量点振幅约 ±0.4–0.5，橙色为模拟，黑色为模拟与 FWHM=3.4 ps 高斯（实验时间分辨率）的卷积；
  (c) 两臂 g⁽²⁾(0)(t) 随 Josephson 振荡在 Poisson 与 bunching 间切换、两臂反相；
  (d)(e) 模拟的 squeezing 参数演化。

## 3. 论文给出的复刻相关事实

### 3.1 样品与激发

- 两个圆形 mesa，直径 2 μm，中心间距 2.5 μm，隧穿耦合 J = 0.4 meV，形成 bonding/antibonding
  正交模（Fig. 1）。样品为 GaAs λ 腔、单 10 nm In₀.₀₆Ga₀.₉₄As 量子阱、24/20 对 Bragg 镜；
  Rabi 分裂 3.3 meV，激子-光子失谐 −3 meV，polariton 寿命 τ = 5.2 ps（κ = ħ/τ = 125 μeV）。
- Ti:Sapphire 锁模激光 λ = 848 nm，脉宽 ≈300 fs，重复频率 80.5 MHz；电光调制器（EOM，消光比
  >10³）加声光调制器（AOM，~2.5 µs 窗口）做 pulse picking，保证 streak 相机每帧只来一个脉冲；
  单脉冲重复率 64 Hz，受 ORCA-R2 CCD 帧率限制（SI Note 1）。
- 脉冲整形器把激光谱宽压到 0.7 meV（避免激发高能态；脉宽因此变长但仍远小于 streak 分辨率
  3.4 ps FWHM）。脉冲能量 5 pJ，×50 物镜聚成 3 μm 光斑、主要照在一个 mesa 上（以获得高对比度
  Josephson 振荡）。共振激发 bonding/antibonding 态，E = 1.462 eV。样品 5.1 K；全程 34 h，
  每 10⁴ 个脉冲暂停采集做空间稳定（压电调节，光斑漂移 <≈500 nm）。

### 3.2 探测

- 透射几何收集，×50、0.42 NA 物镜；50:50 分束器构成 Hanbury Brown–Twiss（HBT）装置；
  streak 相机以 synchroscan 模式工作，充当 HBT 两臂共用的单光子探测器，同时给出发射的实空间像。
- 实验时间分辨率 3.4 ps FWHM（正文与 Fig. 2b 图注）。
- 时间抖动校正：分一小束激发激光打到 streak 狭缝上作为时间参考（SI Note 1 仅此一句，
  **未写逐帧定位算法**）。
- 累计统计量：>3.9×10⁶ 个激光脉冲（正文 Results 节，指向 SI Note 1）。

### 3.3 数据处理（SI Note 1，Fig. 2(a) 的直接方法依据）

- 自定义光子计数算法处理时间分辨图像；CCD 上的光子检测阈值选得足够高，以分辨相距很近的
  光子对（具体阈值数值未公开）。
- 对每个脉冲，左/右两臂的光子数 n_L,R(t) = 落在时间窗 (t−Δt/2, t+Δt/2) 内的光子数；
  **Δt = 1.296 ps，对应 streak 屏上 2 个像素**（增大 Δt 提高每 bin 计数、降低时间分辨率）。
  （据此其扫描换算约 0.648 ps/px，与本项目 0.984 ps/px 的窗口不同。）
- **强度图（Fig. 2a）算法：I_L,R(t) = n_L + n_R**，即某个 mesa 的强度 = 该 mesa 光子在 HBT
  两臂的计数之和。
- **臂间无归一化**：正文与 SI 均未提任何两臂/两 mesa 间的效率标定或亮度归一化。I = n_L+n_R
  的求和构造使 HBT 分束器两臂的效率不对称在强度中自动抵消；两 mesa 的亮度差是故意为之
  （光斑主要照在一个 mesa 上），直接保留，Fig. 2b 的 z(t) 按原始测量值画。
- g⁽²⁾(t₁,t₂) = ⟨n₁(t₁)·n₂(t₂)⟩ / (⟨n₁(t₁)⟩⟨n₂(t₂)⟩)，平均对脉冲数进行；g⁽²⁾(0) 误差由
  ⟨n₁⟩、⟨n₂⟩、⟨n₁·n₂⟩ 的标准误（SEM）推得。

### 3.4 模拟参数（仅备查；复刻 2(a) 不需要做模拟）

U = 1.4 μeV，J = 0.4 meV，ħΔ_L,R = −0.6 meV，κ = 125 μeV，p_L = 50κ（初始蓝移
μ = UN₀ ≈ 0.7 meV），z(0) = (p_L−p_R)/(p_L+p_R) = −0.66。这些参数用于复现实验中测得的
振荡周期、初始蓝移、激光失谐与 polariton 寿命。

## 4. 论文未写明、本项目必须自行决定的方法点

| 方法点 | 论文状态 | 本项目的决定 | 依据 |
|---|---|---|---|
| 逐帧时间零点算法 | 未公开（仅"use it as a time reference"） | 参考脉冲行廓线前沿 50% 穿越的恒比定时（CFD50），亚行插值 | [08-04_sync_oscillation.md](08-04_sync_oscillation.md) |
| 光子检测阈值 | "足够高以分辨近距光子对"，无数值 | 25 counts、8-连通、min 2 px，由灵敏度扫描拐点确定 | 26-08-04-fig2a.md |
| 背景扣除 | 未写 | 实测 `StreakCamBackG`（≈105 counts/px）逐帧扣除 | 26-08-04-fig2a.md |
| L/R 与探测器像素的对应 | 未写 | 无脉冲实空间强度标定（L/R 比 0.578）配对 | 26-08-04-fig2a.md |
| 行→ps 换算 | 未给（其窗口 0.648 ps/px） | 采集 metadata：1 ns / 1016 px = 0.984 ps/px | `metadata-conventions.toml` |
| 纵轴归一 | a.u.，未写归一方式（峰值 ≈0.5–0.6） | 复刻图取峰值归一；可选 `matched`（R 臂按窗口积分和缩放到 L，消除恒定亮度偏移） | 26-08-04-fig2a.md |
| SNR 算法与灰区画法 | 未给算法；原图为底部横向灰带 | 本项目 SNR 定义见 render_fig2a.py；SNR 中位 >20 时定稿不画灰区 | 26-08-04-fig2a.md |

## 5. 复刻目标分解：可完成、可验证的细则

每条给出验证标准；状态栏为 2026-08-06 全量 CFD50 运行后的现状（证据细节见
26-08-04-fig2a.md 相应条目）。

- **G1 数据就绪** —— 派生 HDF5 含 HBT 主条件（20000 帧 single-sweep、侧面带激光时间参考）、
  实测背景条件与无脉冲强度标定条件；逐像素验证通过。
  验证：`DATA/26-08-04-derived-hbt-time-reference/README.md` 的哈希与逐像素记录。✅ 已完成。
- **G2 背景扣除** —— 每帧减实测 `StreakCamBackG`。
  验证：未曝光区域扣后均值 ≈0（管线诊断图）。✅ 已完成。
- **G3 逐帧时间零点** —— 参考 ROI（行 110–210 × 列 170–215）行廓线上 CFD50 穿越行为 t0，
  峰未高出噪声 5σ 的帧丢弃。
  验证标准：t0 失败帧占比 <0.1%（实测 8/20000）；合成注入实验中 t0 偏差 |bias| ≤0.2 px、
  逐帧 std ≤0.5 px（真实形状参考点实测 −0.10 px / 0.33 px）。✅ 已完成。
- **G4 时间轴换算** —— 0.984 ps/px；样品脉冲与参考点延迟实测 0.32±0.20 ps ≈ 0，不设额外偏移。
  保留：行→ps 换算与 row 0 约定来自采集 metadata，未经数据内标定（开放缺口，见 G12）。
  验证：延迟标定条件的测量记录。✅ 已完成（附保留）。
- **G5 光子事件检测** —— 8-连通组分、阈值 25 counts、min 2 px。
  验证：灵敏度扫描图（threshold_sensitivity）存在明确拐点；analog（扣背景计数）与 photon
  （事件计数）两种迹线定性一致。✅ 已完成。
- **G6 空间带检测与 L/R 配对** —— 四条带 B1–B4（x=380–450、570–680、810–880、1010–1110）；
  臂配对 {B1,B2}/{B3,B4}（臂内间距 210/215 px ≈ 双 mesa 标定 204 px）；L={B1,B3}、R={B2,B4}。
  验证：两臂 L/R 强度比 0.586/0.579 与无脉冲标定 0.578 一致；交叉配对给出 0.59/1.51，可排除。
  ✅ 已完成。
- **G7 累积迹线** —— 各帧按 t0 亚行插值对齐到公共时间栅后累计 mean/SEM，photon 与 analog 双源并存。
  验证：npz 产物含双源迹线与 SEM。✅ 已完成。
- **G8 论文式渲染** —— 蓝=L、红=R 散点、对数强度轴、峰值归一、rebin 8（≈7.9 ps/bin）、
  0–820 ps 窗口；SNR<20 灰区可选（定稿 SNR 中位 L 22 / R 29 >20，未画）。
  验证：渲染输出与论文原图并排目检；矩阵类附图带色标（仓库约定）。✅ 已完成。
- **G9 科学判定：反相交替是否存在** —— 预设判据：z=(I_L−I_R)/(I_L+I_R) 去趋势涨落的两臂
  相关显著（置换检验 p<0.05），或与参考脉冲锁相的正弦振幅显著非零。
  实测：r=0.030（rebin 2，p=0.525）/ 0.095（rebin 4，p=0.145），均不显著；3–200 ps 周期
  z 正弦振幅 95% 上限中位 0.0123、最差 0.0193；k 空间实测周期 37.2/41.8 ps 处 |A|≤0.004、
  上限 ≤0.012；t≲60–150 ps 的 L 侧缺口是单调恢复（过阻尼弛豫），非振荡；加性偏移 +
  线性放大直接目检（`fig2a_offset_zoom.png`，详见 26-08-04-fig2a.md）进一步给出平稳段
  两臂去均值涨落弱共模**正**相关（r=+0.174，置换 p=0.023），无镜像反相结构。
  ✅ 已完成，**结论为零结果（上限约 1–2% 强度调制）**。
- **G10 灵敏度验证** —— 合成帧（真实实测参考脉冲形状+呼吸、真实背景、双态 t0 漂移）注入
  9.84 ps 周期、z 振幅 0.5 的反相振荡，经同一管线处理。
  验证标准：恢复周期与注入值差 <2%（实测 9.87 vs 9.84 ps）；证明零结果不是机器灵敏度问题
  （论文 Fig. 2b 的实测振幅 |z|≈0.4–0.5 恰在注入量级）。✅ 已完成。
- **G11 定稿对比图** —— 论文原图裁剪件（左）与复刻面板（右，同一 npz 原生重绘）合成，
  附英文图注注明 CFD50 时间校准。
  验证：`analysis/26-08-04-derived-hbt-time-reference/fig2a_replica.png` 已按仓库惯例自审。
  ✅ 已完成。
- **G12 未闭合缺口（后续任务入口）** —
  1. 绝对时间轴的数据内标定未做（G4 的保留项）；
  2. 平均迹线之外的逐帧/层级分析未做（[experiment-plan.md](experiment-plan.md) Gate 1 第二问：
     单次 Stark 事件轨迹内是否有可分辨振荡）；
  3. 物理区间差异：论文是共振脉冲激发 bonding/antibonding、30 ps 窗口、z(0)≈−0.66；本项目是
     非共振 reservoir-fed 双阱 + Stark kick、1 ns 窗口。当前复刻的目标是**测量管线与图像形式**，
     不是同一物理状态；同区间对比需要共振激发数据，当前平台不采集（experiment-plan.md「当前平台」）。
  4. k 空间同步扫描图记录到 37.2 ps 的相位缠绕条纹（相位动力学），与实空间人口零结果并存、
     不矛盾；Fig. 2(a) 的观测量（人口）对相对相位不敏感。

## 6. 一句话现状

管线（G1–G8、G10–G11）全部完成并验证；核心判定 G9 为零结果——平均 single-sweep 实空间
数据中无论文式反相人口交替，95% 上限约 1–2% 强度调制；灵敏度验证表明该上限由数据本身决定，
而非管线盲区。细节与复现命令见 [26-08-04-fig2a.md](26-08-04-fig2a.md)。
