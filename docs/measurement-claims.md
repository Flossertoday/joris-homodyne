# 测量与物理主张

本项目必须持续区分以下五层结论。

## 1. Shot-resolved Josephson dynamics

证据是左右人口的反相轨迹、较弱的 common-mode response、与 coupling/detuning/Stark pulse area 一致的参数变化，以及相对相位控制。它是动力学结论，不是量子见证。

## 2. Inter-mode intensity correlation

\(g^{(2)}_{LR}(t_1,t_2)\) 描述两阱强度涨落的联系，既可能来自量子态，也可能来自 reservoir noise、mode competition 或共同探测器漂移。只测 LR 还不是 CSI；必须同时取得匹配的 LL 和 RR，并校准背景、时间响应和通道串扰。

## 3. Intensity Cauchy–Schwarz inequality

采用

\[
\mathcal C(t_1,t_2)=
\frac{[g^{(2)}_{LR}(t_1,t_2)]^2}
{g^{(2)}_{LL}(t_1,t_1)g^{(2)}_{RR}(t_2,t_2)}.
\]

正 Glauber–Sudarshan \(P\) 表示的经典场满足 \(\mathcal C\le 1\)。显著的 \(\mathcal C>1\) 能排除这一类经典场模型，但不自动证明 quadrature squeezing、two-mode squeezing 或 entanglement。若改用倒数定义，经典不等式和违反方向也必须同时反转；旧方案曾在这里出现符号逻辑错误。

## 4. Single-mode quadrature squeezing

2017 论文模型中的量

\[
r=\frac{|\langle\delta a^2\rangle|-\langle\delta a^\dagger\delta a\rangle}{2}
=\frac14-\min_\phi\operatorname{Var}X_\phi
\]

确实对应某一单模正交量低于 vacuum variance。直接检验需要与 condensate 相位相关、时空模式匹配的 local oscillator、相位扫描、vacuum/shot-noise 标定和平衡探测。强度 \(g^{(2)}>1\) 本身不是该测量。

## 5. Two-mode squeezing or entanglement

这需要联合的 phase-sensitive quadratures 或适用的 inseparability criterion。强度互相关或 CSI violation 与它并不等价。

## 对 2017 论文的证据表述

论文实验上测得 Josephson population dynamics、局域同模 \(g^{(2)}(t,t)\) 和局域 two-time \(g^{(2)}(t_1,t_2)\)。squeezing 曲线来自选定量子模型，而不是直接 homodyne 数据。原始稀疏事件在原理上很可能容许回算左右 cross-correlation，但论文没有充分记录四个空间窗、跨模式串扰标定或现存原始文件，因此这仍是推断。

