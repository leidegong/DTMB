# DTMB 低时延解调—再调制最终版（2026-09-02）

当前正式方案固定为纯算法链路，不使用神经网络、TCN、残差学习器或其他学习型逆信道模型。主入口为 `run_low_latency_reference_pipeline.m`，单帧实时核心为 `DTMB/DTMB_Low_Latency_Demod_Remod_Frame.m`。

## 最终链路

```text
同步帧
  → 常规 TC-AOMP 信道估计/均衡（默认仅前 2 帧；若上游已均衡则旁路）
  → 单帧 3780 点 FFT 与频率解交织
  → 首帧 TPS/CPE 捕获，后续判决导向相位跟踪
  → 仍处于时间交织顺序的 64QAM 最近点硬判决
  → 频率交织、3780 点 IFFT、PN945 组帧
  → 单帧复增益匹配
  → 跨帧保留状态的连续 SRRC FIR
```

明确不包含：学习型逆信道模型、时间解交织、软解映射、LDPC、BCH、业务比特恢复、重编码和时间重交织。历史 `apply_482mhz_*` 与训练脚本仅保留为实验档案，不是最终入口的依赖。

## 运行

```matlab
run('QUICKSTART_LOW_LATENCY_FINAL.m')
```

使用自己的同步数据：

```matlab
cfg = struct('tc_channel_estimation_frames',2, ...
    'save_shaped_waveform',true);
result = run_low_latency_reference_pipeline( ...
    'D:/data/DTMB_Sync_Packet.mat','D:/output',cfg);
```

若上游已经完成均衡，可取消本入口的 TC-AOMP 启动等待：

```matlab
cfg = struct('input_is_equalized',true);
result = run_low_latency_reference_pipeline(inputFile,outputDir,cfg);
```

输入 MAT 文件至少需要 `PN_phase_index_packet` 和帧矩阵。同步后输入还需要 `PN_svd_denoised_packet`；已均衡输入优先读取 `equalized_reference_frames`，不存在时读取 `Frame_equalization_packet`。

## 时延口径

符号率为 7.56 MHz，一帧 4725 符号对应 0.625 ms。重建帧体必须等当前帧接收完整，所以核心不可消除的输入缓存是 1 帧。默认 2 帧 TC-AOMP 公共信道启动缓存对应 1.25 ms；上游已均衡时总启动缓存回到 1 帧。跨度 8、每符号 2 点的 SRRC 群时延为 8 个输出采样点，约 0.000529 ms。

低时延实现使用以下措施：PN ROM 使用预生成 MAT 加载并持久缓存（文本 ROM 作为回退）；单精度 FFT/IFFT；无距离矩阵的 64QAM 算术量化；系统信息首帧锁定；首帧全范围 CPE 捕获，后续单次判决导向校正，超过 8° 时自动重捕获；SRRC 通过 `filter` 的 `zi/zf` 跨帧连续处理。

`result.timing` 同时给出信号缓存时延、首帧冷启动计算时延和跳过前 5 帧初始化后的稳态计算时延。计算时延受 MATLAB 版本、首次 JIT/FFT 初始化及硬件影响；启动响应查看 `first_frame_compute_ms`。

**不要用 `steady_compute_ms_per_frame_*` 判断实时吞吐。** 它是内层计时器的中位数，既不含信道前端也不含批处理开销，只能作下界。1099 帧全长验证（`validate_low_latency_full_length.m`）给出的实时因子区间为 0.891--1.646，判定为 `inconclusive_harness_bounded`：中位数口径下界 0.891，整链墙钟上界 1.646，本批处理测试台无法收窄。全长还显示 p95 为 0.7832 ms、单帧最大 6.3722 ms、16.5% 的帧超出 0.625 ms 单帧预算，这些在 16 帧样本里都不可见。要给出实时结论需要逐帧喂入、不预分配整段输出的流式基准。详见 `DTMB_low_latency_full_length_validation/README.md`。

## 验证

```matlab
report = test_low_latency_demod_remod();
```

两个运行入口分工：`QUICKSTART_LOW_LATENCY_FINAL.m` 是 16 帧冒烟测试（输出到 `DTMB_low_latency_final_output/`）；`validate_low_latency_full_length.m` 是 1099 帧全长验证（输出到 `DTMB_low_latency_full_length_validation/`），额外给出实时余量区间、超预算帧计数、EVM 的 p95、与完整 FEC 重构的 NMSE/相关，以及 3 次计时重复的区间。两者目录分开，不互相覆盖。

全长验证结果：79 个 FEC 有效帧对完整 FEC 重构 NMSE -16.258 dB、相关 0.988093；全部 1099 帧相关 0.987809；CPE 锁定后重捕获 0 次；逐帧计算时间无漂移。

测试覆盖 64QAM 量化等价性、跨帧 SRRC 与整段卷积的一致性、相对旧全网格 CPE 的 EVM 门限，以及对逆信道模型、FEC 和时间交织依赖的静态检查。
