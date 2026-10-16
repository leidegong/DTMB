# DTMB 完整解调—再调制重构

## 入口与数据范围

主入口为 `real_DTMB_joint_sync_main.m`。默认读取
`482MHz_3.mat/482MHz_3.mat` 中的 `raw_data`，从第 `302400` 个采样点开始，
处理 1100 个信号帧长度的数据：

```matlab
addpath('DTMB');
out = real_DTMB_joint_sync_main();
```

入口接受标量结构体覆盖参数。例如短帧调试：

```matlab
out = real_DTMB_joint_sync_main(struct( ...
    'target_frames', 80, ...
    'plot_diagnostics', false, ...
    'save_diagnostic_figures', false));
```

## 实现链路

程序依次执行 PN945/CFO/定时同步、TC-AOMP 信道估计与均衡、3780 载波
频域解交织、36 符号系统信息判决、52 支路卷积时域解交织、64QAM 软
解映射、GB 20600-2006 三种准循环 LDPC 译码、缩短 BCH(762,752) 校正、
重新编码和穿孔、64QAM 映射、时域/频域重交织、帧体 IFFT 与 PN945 拼接。

LDPC 默认采用归一化最小和译码（`norm-min-sum`，系数 `0.75`），噪声
可靠度随 64QAM 符号一起通过时域解交织。该配置由当前数据的全码块搜索
得到；相较默认 BP 译码，它对均衡后软信息更稳健。

由于卷积交织，模式 1 和模式 2 的单向时延分别为 170 帧和 510 帧。
完整“解交织—译码—重交织”至少需要 341 帧或 1021 帧；边界区域始终生成
判决导向重构。

## 输出与有效性

最终文件为 `DTMB_Full_Demod_Remod_Reconstruction.mat`，变量名 `result`。
主要字段：

- `reconstructed_frames`：复增益匹配后的 4725 符号重构帧；
- `canonical_reconstructed_frames`：标准归一化重构帧；
- `reconstruction_source`：1 为判决导向，2 为通过 FEC 校验的译码再编码；
- `fec_valid_frame_mask`：可以安全使用完整 FEC 重构的帧；
- `system_info`：调制、码率、交织模式及 TPS 相关置信度；
- `fec_stats.parity_ok`：逐 LDPC 码块校验结果；
- `decision_evm_rms`：逐帧最近星座点 EVM。

程序不会用未通过 LDPC 校验的码块覆盖判决导向结果。失败码块仍保存在
诊断字段中，便于后续改进均衡器，但不能当作有效载荷或纯化参考信号使用。

## 本数据实测结果

从第 302400 点开始的 1100 帧读取量得到 1099 个同步有效帧。TPS 判决为
64QAM、0.6 码率、交织模式 2，相关系数约 0.990。全 52 支路相位、方向和
7488 个码块边界穷举均确认标准对齐（方向 `+1`、相位 `0`、比特偏移 `0`）。

优化分三步完成：公共判决导向信道迭代把 LDPC 硬综合比例从 `0.4778`
降到 `0.4174`；32 帧一组的时变信道跟踪进一步降到约 `0.4105`；随后
LDPC 算法搜索发现默认 BP 与当前 LLR 分布不匹配，改用归一化最小和译码
（系数 `0.75`）后，全量 `1767/1767` 个码块通过，BCH 综合全部为零。
模式 2 去除双侧 510 帧交织暖机区后，第 511–589 帧共 79 帧已使用严格的
FEC 译码—重编码结果，其余边界帧继续使用判决导向重构。

最终数据包为
`DTMB_groupwise_equalization/DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat`。
效果测试的中位 EVM 从 TC-AOMP 的 `0.1458` 降到 `0.0203`，改善约
`17.13 dB`；重构信号中位相关系数为 `0.9179`。

各阶段可单独复现：

```matlab
optimize_DTMB_FEC_equalization
refine_DTMB_groupwise_equalization
optimize_DTMB_LDPC_soft_decoding
```

核心自检可运行：

```matlab
test_DTMB_demod_remod_core
```
