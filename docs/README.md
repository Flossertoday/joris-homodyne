# 文档导航

本目录只保留需要长期维护、用户能够直接审阅的文档：

- [experiment-plan.md](experiment-plan.md)：当前装置下的分阶段实验路线和判定门槛。
- [measurement-claims.md](measurement-claims.md)：Josephson 动力学、强度关联、CSI 与 squeezing 各自能支持什么结论。
- [literature.md](literature.md)：核心论文、相邻先例及其对本项目的实际约束。
- [data.md](data.md)：原始 MAT、分析 HDF5、环境和可重复性约定。
- [streak-time-axis.md](streak-time-axis.md)：streak 纵轴时间窗的 metadata 调查、1 ns 窗口依据与轴约定。
- [adiyatullin2017-fig2a.md](adiyatullin2017-fig2a.md)：Adiyatullin 2017 论文侧事实存档（免去回读 PDF）与 Fig. 2(a) 复刻目标的可验证分解（G1–G12）；对 08-04 数据集的复现已收尾、未成功，留待 26-07-31 数据集重试。
- [adiyatullin2017-fig2c.md](adiyatullin2017-fig2c.md)：Fig. 2(c)（g⁽²⁾(0) 时间演化）论文侧事实存档与复刻目标分解（G1–G9）；G3 统计功效门槛已执行（borderline），经确认降级为上限报告。
- [26-08-04-fig2a.md](26-08-04-fig2a.md)：26-08-04 HBT 批次复刻 Adiyatullin 2017 Fig. 2(a) 的管线、配对证据与交替振荡零结果（含灵敏度验证）。
- [26-08-04-fig2c.md](26-08-04-fig2c.md)：26-08-04 HBT 批次复刻 Fig. 2(c)（g⁽²⁾(0)(t)）的执行文档：G3 门槛结论、G6 判据预登记、上限报告结果（四项检验全部零结果）与复现命令。
- [26-08-06-hbt-supplement.md](26-08-06-hbt-supplement.md)：26-08-04 HBT 补充批次（72000 帧）的 g2 与 intensity 同管线分析：采集条件确认、带窗重定案、零结果（上限收紧）、G3 borderline / G6 全零 / G7 通过。
- [26-08-08-dataset-a.md](26-08-08-dataset-a.md)：数据集A——"排除一切不可能是数据的峰"的策展（横向 ROI/缺陷列/触边；形态学不可分的证据）、噪声局部地板模型、只引用数据集A的 Fig. 2(a)/2(c) 重算（零结果不变、上限收紧）与帧级/立方级注入验证（算法对论文强度信号灵敏）。
- [sem-limitations.md](sem-limitations.md)：g2 测量 SEM 的限制因素分解（帧数/信号率/底噪/bin 宽×周期/t0）、底噪问题的证据与量级、帧数需求表。
- [peak-localization.md](peak-localization.md)：单帧峰定位算法（Crocker–Grier/trackpy）、圆心语义与二维高斯拟合对照证据。
- [peak-localization-error.md](peak-localization-error.md)：峰定位误差分析——社区做法调查（Savin–Doyle/trackpy、DAOPHOT artificial star tests、SMLM）、逐峰 `ep` 与人工斑注入验证结果。
- [sandbox-initialization.md](sandbox-initialization.md)：托管 Windows 沙盒的初始化入口、实测问题和剩余权限边界。

精确哈希、文件数量和逐像素验证结果放在对应实验的 derived 目录中，不在这里重复，以免多处副本逐渐不一致。

