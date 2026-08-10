# 26-08-11 交叉 g2(condensate 1×2,LR/RL 对)

补项目 g2 矩阵的最后一块:$g^{(2)}_{LL}$、$g^{(2)}_{RR}$(同 mesa 自 g2,L=B1×B3 /
R=B2×B4)已在 Fig. 2(c) 复刻中完成并均为零结果;**交叉 g2($g^{(2)}_{LR}$ /
$g^{(2)}_{RL}$)** 此前是刻意留出的后续入口(`adiyatullin2017-fig2c.md` G9 第 4 条、
项目报告 §8.3、MEMORY.md 的 full-matrix 目标),本文档记录其首次计算。

## 配对定义与依据

- 四带几何:2 mesa × 2 HBT 臂,臂 A = {B1, B2},臂 B = {B3, B4}(臂间距 430 px,
  mesa 间距 210 px,证据见 `26-08-04-fig2a.md` / 报告 §2.2)。
- 交叉对取**跨臂**:**LR = B1×B4**(臂 A 的 mesa 1 × 臂 B 的 mesa 2)、
  **RL = B2×B3**(臂 A 的 mesa 2 × 臂 B 的 mesa 1)。二者测同一个
  $g^{(2)}_{LR}(0)(t)$,互为独立测量(带互不重叠)——内部一致性检验。
- 同臂交叉对 B1×B2 / B3×B4 **不采用**:两臂共享同一条 streak 迹线与固定图案
  噪声,通道不独立;跨臂对的归一化 g2 不受分束器 50% 分光影响(分子分母同比
  缩减)。
- 估计器、残余背景 β 修正、plug-in SEM、pooled jackknife 全部复用
  `fig2c_g2_estimator.py` 的同 mesa 实现,仅配对不同
  (`CROSS_MESA_PAIRS = {"LR": (0,3), "RL": (1,2)}`)。零假设下
  $g^{(2)}_{LR}=1$;数守恒的 Josephson 布居振荡使两凝聚体反相关
  ($g^{(2)}_{LR}<1$ 方向)。

## 结果(窗口 0–820 ps pooled,jackknife SEM)

| 批次 | 帧数 | 4 px (3.94 ps) LR | 4 px RL | 合并(逆方差) |
|---|---|---|---|---|
| 主批次 | 19 992 | 0.659 ± 0.265 | 0.875 ± 0.281 | 0.76 ± 0.19 |
| 补充批次 | 71 971 | 1.110 ± 0.190 | 1.065 ± 0.200 | 1.09 ± 0.14 |
| 数据集A 主批次(局部基底 β) | 19 992 | 0.770 ± 0.179 | 0.936 ± 0.203 | 0.84 ± 0.13 |
| 数据集A 补充批次(局部基底 β) | 71 971 | 1.075 ± 0.174 | 1.076 ± 0.183 | 1.08 ± 0.13 |

数据集A 由 `compute_g2_cross_mesa_dataset_a.py` 计算(复用 `dataset_a_g2.py` 的
载入/立方/局部基底 β,并带旧全局 β 变体作系统误差夹逼;全部数值见
`drafts/26-08-11_g2_cross_mesa/dataset_a_{main,supplement}/g2_cross_summary.json`)。

1/2 px 分档数值见各批次 `g2_cross_summary.json`(趋势一致,全部与 1 兼容)。

- **两批次、两配对全部与 1 兼容**(最大偏离 1.3σ,主批次 LR 偏低方向;
  补充批次 LR 偏高 0.6σ)——交叉关联为零结果,分辨率与同 mesa 同量级。
- LR 与 RL 在各自批次内互相一致(同帧不同带,零假设下近似独立;合并值仅供
  参考,跨批次不合并——两批次物理状态不同,见 `26-08-06-hbt-supplement.md`
  边界第 4 条)。
- 与同 mesa 结果并读:主批次 L 0.31±0.38 / R 1.11±0.19,补充批次
  L 1.30±0.32 / R 0.99±0.12——四个自/交叉观测量同为零结果。
- 局限:当前占据率下逐 bin g2 退化(plug-in SEM 中位 2.7–4.1,全部 bin 统计受
  限),只有 pooled 值有判别力;本次未做振荡/反相的预先登记式检验,若要做需
  先补功效门槛(power gate)再定判据。

## 产物与复现

产物(Git 忽略):`drafts/26-08-11_g2_cross_mesa/{main,supplement}/` —
`g2_cross_results_bin{1,2,4}px.npz`、`g2_cross_summary.json`、
`fig2c_cross_mesa_bin4px.png`。

```bash
# worktree: kimi/g2-cross
# 主批次(事件表随 Git 跟踪)
python -u analysis/26-08-04-derived-hbt-time-reference/compute_g2_cross_mesa.py \
  --out-dir drafts/26-08-11_g2_cross_mesa/main
python -u analysis/26-08-04-derived-hbt-time-reference/render_fig2c_cross_mesa.py \
  --out-dir drafts/26-08-11_g2_cross_mesa/main --width 4 --title "26-08-04 HBT main batch"
# 补充批次(事件表未跟踪,在主检出;带窗覆盖同 26-08-06 定案)
export FIG2C_BAND_WINDOWS="350,440,540,680,770,860,970,1100"
python -u analysis/26-08-04-derived-hbt-time-reference/compute_g2_cross_mesa.py \
  --table E:/Documents/Codex/joris-homodyne/analysis/26-08-04-derived-hbt-time-reference/detected-peaks-supplement/26-08-04_hbt_time_reference_supplement_events.h5 \
  --out-dir drafts/26-08-11_g2_cross_mesa/supplement
python -u analysis/26-08-04-derived-hbt-time-reference/render_fig2c_cross_mesa.py \
  --out-dir drafts/26-08-11_g2_cross_mesa/supplement --width 4 --title "26-08-04 HBT supplement batch"
```

验证:配对接线测试(共享 B1×B4 符合只出现在 LR、不进 RL)、CLI 端到端
Poisson→1 测试入 `tests/test_fig2c_cross_mesa.py`;`pair_counts` 交叉对断言入
`tests/test_fig2c_g2_estimator.py`;全量 181 项测试通过。
项目报告整合(2026-08-11,本分支):§3 第 7 步改写为含交叉对、§6.2/§6.3 补交叉
数值、摘要与图 10(C7)面板扩为 L/R/LR/RL 四系列、§8.3 移除该项后续——中英双版
同步并经程序化核验(图链接/编号标题/表格竖线/插入数字计数一致)。注意:
project-report 分支另有在途修订(未提交),合并时报告 md 与 PDF 需解冲突。
