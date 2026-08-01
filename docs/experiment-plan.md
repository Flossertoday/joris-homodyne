# 实验路线

## 当前平台

- 非共振泵浦维持 condensate dimer；近期不切换到共振激发。
- 独立脉冲激光施加 optical-Stark kick，提供可重复的扰动和时间原点。
- 重复/积累 streak-camera 测量已经看到左右 condensate 的强度振荡。
- 现阶段没有相位锁定的 local oscillator，也没有 shot-noise-limited balanced-homodyne 接收链。

因此，当前平台最自然的研究对象是 Stark-triggered、driven-dissipative Josephson dynamics，而不是直接复现 2017 年论文的 resonant-pulsed quadrature-squeezing 测量。

## Gate 1：single-sweep 轨迹

先回答两个不同问题：

1. 许多 single-sweep 帧的平均是否在匹配探测器响应后复现既有积累曲线？
2. 单次 Stark 事件的轨迹中是否存在统计上可分辨的振荡，而不是只在平均后出现？

每帧优先提取

\[
I_L(t),\quad I_R(t),\quad S(t)=I_L(t)+I_R(t),\quad D(t)=I_L(t)-I_R(t).
\]

比较振荡相位、幅度、频率、阻尼和初始不平衡的分布。单帧不必肉眼可见完整振荡；低光子数时应采用逐帧似然或层级模型，而不是只筛选“好看”的帧。

在称其为 Josephson exchange 前，需要看到：

- 左右强度以反相分量为主；
- 总强度的振荡明显弱于差分强度；
- 频率和相位随 coupling、detuning 与 Stark pulse area 合理变化；
- 时间基准、sweep response、MCP gain、背景和空间串扰已经标定；
- detuning、polarization 和 fluence controls 能排除 reservoir excitation 或普通检测瞬态。

## Gate 2：相对相位

若 Gate 1 成立，再测 time-resolved relative phase，并检验人口—相位耦合是否与同一组参数一致。Stark kick 能提供时间原点并造成相位印记，但不自动提供绝对光学相位参考，也不能代替 homodyne local oscillator。

## Gate 3：完整强度关联

随后测量

\[
g^{(2)}_{ij}(t_1,t_2),\qquad i,j\in\{L,R\},
\]

同时保留 LL、RR、LR、RL，而不是只输出一个 CSI 比值。正常序 \(g^{(2)}\) 需要许多稀疏 photon-counting realization、同一 Stark 事件的身份以及独立 HBT arms；一张 bright analog single-sweep 图不是一个关联矩阵。

CSI 或 relative-number noise 可以从完整数据中派生，但只有在样品特定模型预测的效应大于统计和系统误差时，才适合作为量子见证目标。

## 暂不开展

- 在没有可检验预测前围绕 CSI violation 重建光路。
- 把强度 bunching 或不违反 CSI 解读为 quadrature squeezing 的证实或否定。
- 在没有相位相关 LO 和平衡探测链时声称直接测得单模 quadrature squeezing。

