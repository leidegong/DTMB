# DTMB reference-signal reconstruction

面向无源雷达参考信号的 MATLAB 解调—再调制实现。此仓库发布当前算法源码，不包含实测 IQ/ADC 数据、训练模型、本地缓存或旧版压缩包。

## 两条处理链路

- **低时延支路**：`run_low_latency_reference_pipeline.m`。同步帧输入，默认采用固定预算 PN 稀疏信道估计/均衡、TPS/CPE 跟踪、64QAM 硬判决、IFFT、PN945 重组和状态化 SRRC。运行路径不使用学习型逆信道模型或 FEC。
- **完整 FEC 支路**：`run_optimized_full_chain.m`。均衡帧输入，包含时间解交织、软解映射、LDPC/BCH、重编码和再调制。当前入口输出符号率重构帧，不包含输出 SRRC 成形。
- **同步前端**：`DTMB/real_DTMB_joint_sync_main.m`。需要自行配置采样数据及同步参数；它是离线前端，不是已验收的因果实时接收机。

两个运行入口目前仍是独立批处理程序。单帧状态化核心可供流式集成，但从原始采样到参考输出的端到端实时吞吐尚未验证。

## 环境与测试

已有验证环境为 MATLAB R2022a，需 Signal Processing Toolbox 和 Communications Toolbox。未声明 GNU Octave 兼容。

在仓库根目录执行无需实测数据的测试：

```matlab
addpath('DTMB');
test_DTMB_demod_remod_core;
report = test_full_chain_throughput_optimization();
```

完整链路测试在本地没有历史重构结果时，会跳过实测波形回归比较，只执行合成输入的软解映射与 BCH 检查。`test_low_latency_demod_remod` 另需本地 `fast_16frame_sync_input.mat`，该数据不随仓库发布。

## 使用自己的数据

同步输入 MAT 需包含 `Frame_equalization_packet`（每行 4725 个复样点）、`PN_svd_denoised_packet` 和 `PN_phase_index_packet`。相位索引沿用 `PN_gen` 的模 200 约定。

```matlab
result = run_low_latency_reference_pipeline( ...
    'my_sync_packet.mat', 'outputs/fast', struct());
```

若输入已经均衡，设置 `input_is_equalized=true`，并提供 `equalized_reference_frames`（或 `Frame_equalization_packet`）及相位索引。

```matlab
cfg = struct('input_is_equalized', true);
result = run_low_latency_reference_pipeline( ...
    'my_equalized_packet.mat', 'outputs/fast', cfg);

% FEC 输入包需含 optimized_frames；同步包需含 PN_phase_index_packet。
[result, report] = run_optimized_full_chain( ...
    'my_groupwise_equalized.mat', 'my_sync_packet.mat', ...
    'outputs/full', struct());
```

默认 QUICKSTART、全长验证和历史对比脚本需要未公开的实测文件，不能直接在全新克隆中复现实测数字。用于新数据时，不要把另一段采集的历史结果作为配对真值。

## 结果解释与限制

- 历史 1099 帧测试中，快速支路在 79 个 FEC 有效帧上相对完整重构的相关系数约 0.9881、NMSE 约 -16.26 dB。这不是所有新数据上的保证，FEC 校验通过也不等于已知发送比特真值。
- 内层计时中位数、批处理墙钟耗时和实际流式吞吐是不同口径。约 1.05 ms 是同步后处理的历史估算，不是原始 IQ 到输出的端到端实测。实际流式实现仍需测量均值、尾延迟及队列积压；历史报告中的上下界不能视为严格数学界限。
- 完整链路历史计算约 5.26–6.37 ms/输入帧，不含前端同步和均衡。旧模块计时分母和范围不同，文档中的约 8–10 倍加速只可视为历史估算，不能当作同口径基准结论。
- 旧比较中 PSLR 从约 -25.42 dB 变为 -16.34 dB，主要由 PN945 在 511 点延迟处的重复结构及帧头/帧体功率比变化造成。该值是单帧自相关指标，不等于目标检测能力损失；不能仅为改善 PSLR 数值而任意削弱 PN。
- 雷达用途应以同步参考/监视双通道的杂波对消、CAF 目标峰、旁瓣和长 CPI 相干性为验收依据。

## 进一步说明

- [低时延说明](README_Low_Latency_Demod_Remod_Final.md)
- [全长验证记录](DTMB_low_latency_full_length_validation/README.md)
- [完整链路优化记录](FULL_CHAIN_OPTIMIZATION.md)

历史说明中的前端名称、计时口径如与当前实现不一致，以当前代码及本页限制说明为准。PN ROM 是算法查找表，随源码公开；采样数据、参考文献 PDF 与历史实验模型不在发布范围内。
