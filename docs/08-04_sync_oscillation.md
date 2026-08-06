目标: 复刻论文 "E:\Documents\Codex\joris-homodyne\analysis\Adiyatullin et al. - 2017 - Periodic squeezing in a polariton Josephson junction.pdf" 的 Fig. 2(a), 实空间上的强度震荡. 

复刻 Fig. 2(a) 的步骤是:

1. 获得一组有脉冲校正 timing jitter 的数据组(已有, 在 "\26-08-04" 的 "\01_HBT").
2. 定位脉冲时间. 逐帧做图像平移，把每帧按参考脉冲的位置对齐到统一时间原点
3. 把各帧按脉冲时间相加.
4. 当然还需要本底/暗噪声扣除和时间轴定标 (像素 → ps) 等必要步骤。