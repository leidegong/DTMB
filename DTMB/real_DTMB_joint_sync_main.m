function pipeline_result = real_DTMB_joint_sync_main(cfg_override)
% real_DTMB_joint_sync_main
% DTMB 参考信号纯化主流程。
%
% 主链路：
%   初始同步
%   -> 第一次 PN/帧体提取
%   -> 多帧 PN 相干积累
%   -> 残留 CFO 精估计
%   -> 全局残留 CFO 补偿
%   -> 闭环二次/多次同步与矩阵重构
%   -> 逐帧公共相位补偿
%   -> SVD 去噪
%   -> 可选 TC-AOMP 信道估计与均衡

if nargin < 1
    cfg_override = struct();
end
close all; clc;
total_timer = tic;

%% 0. 参数配置
cfg = struct();

% 数据读取参数：指定真实接收数据文件、变量名、起始采样点和目标处理帧数。
cfg.filename = fullfile('482MHz_3.mat', '482MHz_3.mat'); % 原始 DTMB 接收数据文件
cfg.data_var_name = 'raw_data';            % MAT 文件中的接收数据变量名
cfg.start_sample = 302400;                 % 截取数据的起始采样点
cfg.target_frames = 1100;                  % 模式2完整FEC重构需至少1021个有效帧

% DTMB 帧结构与采样率参数：后续同步、矩阵构造和频偏换算都依赖这些值。
cfg.pn_length = 945;                       % PN945 帧头长度，单位：符号
cfg.body_length = 3780;                    % 帧体长度，单位：符号
cfg.frame_symbol_len = cfg.pn_length + cfg.body_length; % 单帧总符号数
cfg.sps = 2;                               % 每个符号对应的采样点数
cfg.Rs = 7.56e6;                           % 符号率，单位：Hz
cfg.Fs = cfg.Rs * cfg.sps;                 % 采样率，单位：Hz
cfg.use_srrc_matched_filter = false;       % 是否启用 SRRC 匹配滤波；关闭时原始 IQ 直接进入同步


% 总体流程开关：控制是否运行调试图、联合 CFO 闭环、TC-AOMP 和诊断图保存。
cfg.debug_mode = false;                    % 粗同步等模块的调试输出开关
cfg.debug_plot_each_frame = false;         % 逐帧 PN 精同步调试图开关
cfg.use_joint_pn_cfo_refine = true;        % 是否启用多帧 PN 残留 CFO 闭环精估计
cfg.plot_diagnostics = true;               % 是否绘制必要诊断图
cfg.run_tc_aomp = true;                    % 是否继续运行 TC-AOMP 信道估计与均衡
cfg.run_full_demod_remod = true;           % 是否执行完整解调、FEC译码和再调制重构
cfg.full_recon_output_packet = 'DTMB_Full_Demod_Remod_Reconstruction.mat';
cfg.ldpc_max_iterations = 30;              % LDPC 最大迭代次数
cfg.ldpc_batch_size = 16;                  % 分批译码以限制峰值内存
cfg.ldpc_algorithm = 'norm-min-sum';       % 实测比默认BP更适合当前DTMB软信息
cfg.ldpc_min_sum_scaling_factor = 0.75;    % 归一化最小和系数（全码块搜索最优）
cfg.ldpc_min_sum_offset = 0.5;             % 仅offset-min-sum模式使用
cfg.ldpc_noise_model = 'deinterleaved';    % 噪声可靠度随符号一起通过时域解交织
cfg.ldpc_llr_scale = 0.25;                 % 最小和译码对全局尺度不敏感，保留搜索首选值
cfg.fec_alignment_search_frames = 12;      % 用少量码块搜索换向/共轭/象限歧义
cfg.fec_search_all_commutator_phases = false; % 标准规定每个数据块首符号同步支路0
cfg.fec_commutator_direction = 1;
cfg.fec_commutator_phase = 0;
cfg.pn_to_body_power_ratio = 2;            % 标准PN帧头与帧体的平均功率比
cfg.run_ambiguity_sidelobe_suppression = false; % 是否在参考信号提纯后执行论文第四章副峰消除参考端处理
cfg.sidelobe_num_frames = [];               % 副峰消除处理帧数，空值表示使用 TC-AOMP 已处理帧
cfg.sidelobe_energy_equalization_mode = 'measured'; % PN 能量均衡方式：paper 使用 1/sqrt(2)，measured 使用实测功率比
cfg.sidelobe_output_packet = 'DTMB_Ambiguity_Sidelobe_Suppression_Packet.mat'; % 副峰消除输出数据包
cfg.sidelobe_use_reference_as_echo = true; % 暂无监测通道时，用提纯参考信号代替回波信号验证完整互模糊流程
cfg.sidelobe_caf_doppler_fft = 256;       % 副峰消除互模糊验证的慢时间多普勒 FFT 点数
cfg.sidelobe_caf_max_lag_symbols = 2 * cfg.frame_symbol_len; % 副峰消除互模糊验证的最大距离延迟，单位：符号；覆盖正负 4 帧
cfg.sidelobe_caf_3d_enable = true;        % 是否绘制连续展开原始 CAF 和三分支联合 CAFm 的三维互模糊图
cfg.sidelobe_caf_3d_doppler_bins = cfg.sidelobe_caf_doppler_fft; % 三维 CAF 多普勒网格点数
cfg.sidelobe_caf_3d_max_doppler_hz = 1000; % 三维 CAF 多普勒显示范围，单位：Hz
cfg.sidelobe_caf_3d_lag_step_symbols = 1; % 三维 CAF 距离延迟绘图步长；1 表示每一个距离延迟点都绘制
cfg.sidelobe_caf_3d_max_lag_symbols = cfg.sidelobe_caf_max_lag_symbols; % 三维 CAF 显示的最大距离延迟
cfg.sidelobe_plot_db_floor = -100;       % 副峰消除零多普勒图和三维图的显示下限，单位：dB
cfg.use_pn945_initial_cfo_sync = true;   % 是否使用 PN945 重复段小数 CFO + 频域差分整数 CFO 估计
cfg.fractional_cfo_num_frames = 200;     % PN945 小数 CFO 多帧累加帧数
cfg.integer_cfo_num_frames = 200;        % PN945 整数 CFO 多帧累加帧数
cfg.integer_cfo_search_min_bin = -20;    % 整数 CFO 搜索最小 bin
cfg.integer_cfo_search_max_bin = 20;     % 整数 CFO 搜索最大 bin

% 单帧和多帧频域帧体星座图参数：控制输出哪一帧、叠加多少帧以及最多绘制多少点。
cfg.plot_frame_idx = 5;                    % 指定绘制第几帧的单帧星座图
cfg.tc_num_frames = [];                     % TC-AOMP 最多处理的帧数，空值表示处理全部有效帧
cfg.tc_channel_estimation_frames = 50;      % 长数据仅用前若干帧估计公共信道，再应用到全部帧
cfg.multi_frame_constellation_frames = 50; % 多帧频域星座图叠加的帧数
cfg.multi_frame_constellation_max_points = 250000; % 多帧星座图最多显示的散点数

% 图像保存与均衡参数：控制诊断图输出目录和频域均衡正则项。
cfg.save_diagnostic_figures = true;        % 是否将所有 figure 保存为 PNG
cfg.figure_output_dir = 'joint_sync_figures'; % 诊断图保存目录
cfg.body_eq_lambda = 1e-3;                 % 帧体频域均衡正则项

% 多帧 PN 残留 CFO 闭环参数：控制闭环迭代次数、停止门限和每轮补偿力度。
cfg.sync_loop_max_iter = 3;                % 残留 CFO 闭环最大迭代次数
cfg.sync_loop_cfo_tol_hz = 0.05;           % 残留 CFO 收敛门限，单位：Hz
cfg.sync_loop_min_improve_hz = 0.01;       % CFO 改善过小时停止迭代的阈值，单位：Hz
cfg.sync_loop_damping = 1.0;               % 每轮 CFO 补偿阻尼系数，1 表示全量补偿
cfg.cfo_estimator = 'fft_phasefit';        % 残留 CFO 估计方式：'phasefit' 或 'fft_phasefit'
cfg.cfo_fft_oversample = 128;              % 帧间 FFT 零填充倍率，用于提高谱峰定位密度
cfg.cfo_fft_window = 'hann';               % 帧间相参 FFT 加窗方式

% PN 去噪和数据包输出参数：控制 SVD 保留秩以及保存的数据包文件名。
cfg.use_pn_svd_denoise = true;           % 是否启用 PN 矩阵 SVD 去噪；关闭时 PN 直接送入后续信道均衡
cfg.rank_keep = 100;                       % PN 矩阵 SVD 去噪保留秩
cfg.use_full_frame_svd_denoise = false;     % 是否启用实验性整帧矩阵 SVD 去噪
cfg.full_frame_rank_keep = 100;            % 整帧矩阵 SVD 去噪保留秩，最大不超过有效帧数
cfg.output_packet = 'DTMB_PN_SVD_Denoised_Packet_JointSync.mat'; % 输出数据包文件名


% NLMS + 频域结构约束直达波净化与信道补偿参数：不使用解调/再调制，只用 PN 训练滤波器并作用到整帧。
cfg.run_nlms_structured_equalizer = false;   % 是否启用 NLMS + 频域结构约束方法
cfg.nlms_filter_len = 20;                  % NLMS FIR 滤波器抽头数
cfg.nlms_mu = 0.005;                       % NLMS 归一化步长
cfg.nlms_delta = 1e-6;                     % NLMS 防除零正则项
cfg.nlms_train_passes = 1;                 % 每帧 PN 重复训练轮数，过多轮数容易过拟合 PN
cfg.nlms_decision_delay = 4;               % NLMS 输出对齐延迟
cfg.use_freq_structure_constraint = false; % 是否启用频域结构约束，当前优化结果显示关闭更好
cfg.freq_constraint_nfft = 1024;           % 频域结构约束 FFT 点数
cfg.freq_constraint_smooth_win = 81;       % 幅度/相位响应平滑窗口
cfg.freq_constraint_gain_floor = 0.05;     % 频域响应幅度下限，相对中位幅度
cfg.freq_constraint_gain_ceiling = 20;     % 频域响应幅度上限，相对中位幅度
cfg = apply_cfg_override(cfg, cfg_override);
fprintf('======================================================\n');
fprintf('DTMB 联合时频同步主流程\n');
fprintf('采样率 Fs = %.6f MHz，每符号采样点 sps = %d，每帧符号数 = %d\n', ...
    cfg.Fs / 1e6, cfg.sps, cfg.frame_symbol_len);
fprintf('======================================================\n');

    %% 1. 读取数据与可选匹配滤波
    fprintf('\n[1] 读取真实 DTMB 数据并进行可选匹配滤波。\n');

    rx_raw = load_rx_segment(cfg);
    if cfg.use_srrc_matched_filter
        fprintf('启用 SRRC 匹配滤波。\n');
        h_srrc = rcosdesign(0.05, 8, cfg.sps, 'sqrt');
        rx_matched0 = conv(rx_raw, h_srrc, 'same');
    else
        fprintf('已关闭 SRRC 匹配滤波，原始 IQ 直接进入同步。\n');
        rx_matched0 = rx_raw;
    end

    %% 2. 初始粗同步与论文二维联合时频同步
    fprintf('\n[2] 初始粗同步与论文二维联合时频同步。\n');

    [detected_mode, coarse_frame_starts, ~, phase_offset_per_samp] = ...
        Frame_Sync(rx_matched0, cfg.sps, cfg.debug_mode);

    if isempty(coarse_frame_starts)
        error('未检测到有效粗帧头位置。');
    end

    fprintf('检测到的帧头模式：%s\n', detected_mode);
    fprintf('初始粗同步有效帧数：%d\n', length(coarse_frame_starts));

    used_initial_cfo_sync = false;
    initial_cfo_info = [];
    if isfield(cfg, 'use_pn945_initial_cfo_sync') && cfg.use_pn945_initial_cfo_sync
        fprintf('启用 PN945 重复段小数 CFO + 频域差分整数 CFO 估计。\n');
        [phase_index_seed, exact_frame_start_seed] = ...
            Symbol_Sync(rx_matched0, coarse_frame_starts, cfg.sps, detected_mode);
        initial_cfo_info = PN945_Fractional_Integer_CFO_Sync(rx_matched0, exact_frame_start_seed, ...
            phase_index_seed, detected_mode, cfg);
        rx_matched_uncomp = rx_matched0;
        rx_matched = initial_cfo_info.rx_comp;
        exact_frame_start = exact_frame_start_seed;
        phase_index = phase_index_seed;
        used_initial_cfo_sync = true;
        fprintf('PN945 初始 CFO 估计：小数 %.9f Hz，整数 bin %+d (%.9f Hz)，合计 %.9f Hz。\n', ...
            initial_cfo_info.f_frac, initial_cfo_info.m_hat, initial_cfo_info.f_int, initial_cfo_info.total_initial_cfo);
    else
        [rx_matched, rx_matched_uncomp, ~, ~] = ...
            Fractional_CFO_Compensation(rx_matched0, coarse_frame_starts, phase_offset_per_samp);
    end

    %% 3. 符号精同步与定时相位精校正
    fprintf('\n[3] 符号级精同步与插值定时相位精校正。\n');

    if ~used_initial_cfo_sync
        [phase_index, exact_frame_start] = ...
            Symbol_Sync(rx_matched, coarse_frame_starts, cfg.sps, detected_mode);
    end

    [rx_matched, ~, rx_16x, ~, ~, ~, exact_frame_start, timing_info] = ...
        Timing_Phase_Refine_8x(rx_matched, rx_matched_uncomp, exact_frame_start, ...
                               cfg.sps, detected_mode, phase_index);

    print_sync_summary(detected_mode, coarse_frame_starts, exact_frame_start, phase_index);

    %% 4. 第一次逐帧 PN945 精同步与矩阵构造
    fprintf('\n[4] 第一次逐帧 PN 精同步，并构造 PN/帧体矩阵。\n');

    num_frames_est = length(coarse_frame_starts);
    first_pass = refine_and_build_matrices(rx_matched, rx_16x, exact_frame_start, phase_index, ...
        cfg.sps, detected_mode, num_frames_est, timing_info, cfg.debug_plot_each_frame);

    fprintf('第一次 PN 矩阵尺寸：%d x %d\n', size(first_pass.PN, 1), size(first_pass.PN, 2));
    fprintf('第一次帧体矩阵尺寸：%d x %d\n', size(first_pass.Frame, 1), size(first_pass.Frame, 2));

    %% 5. 闭环联合时频同步
    sync_pass = first_pass;
    cfo_loop_info = init_cfo_loop_info();

    if cfg.use_joint_pn_cfo_refine
        fprintf('\n[5] 基于多帧 PN 的残留 CFO 闭环精同步。\n');

        [rx_matched, rx_matched_uncomp, rx_16x, exact_frame_start, phase_index, timing_info, ...
            sync_pass, cfo_loop_info] = run_joint_sync_closed_loop( ...
            rx_matched, rx_matched_uncomp, rx_16x, exact_frame_start, phase_index, timing_info, ...
            coarse_frame_starts, detected_mode, num_frames_est, cfg, first_pass);

        fprintf('闭环停止原因：%s\n', translate_stop_reason(cfo_loop_info.stop_reason));
        fprintf('累计施加的残留 CFO 补偿量：%.9f Hz\n', cfo_loop_info.total_applied_cfo_hz);
        fprintf('最终残留 CFO 估计值：%.9f Hz\n', cfo_loop_info.final_residual_cfo.res_cfo_hz);
        fprintf('最终 PN 矩阵尺寸：%d x %d\n', size(sync_pass.PN, 1), size(sync_pass.PN, 2));
        fprintf('最终帧体矩阵尺寸：%d x %d\n', size(sync_pass.Frame, 1), size(sync_pass.Frame, 2));
    else
        fprintf('\n[5] 已关闭多帧 PN 残留 CFO 精估计。\n');
    end

    if cfg.plot_diagnostics && exist('initial_cfo_info', 'var') && isstruct(initial_cfo_info) && isfield(initial_cfo_info, 'rx_frac')
        Plot_PN945_CFO_Stage_Constellations(rx_matched0, initial_cfo_info, sync_pass.PN, cfg);
    end

%% 6. 逐帧公共相位补偿
fprintf('\n[6] 对 PN 与帧体矩阵执行逐帧公共相位补偿。\n');

assert_matrix_consistency(sync_pass);

[PN_sym_matrix_align, Frame_sym_matrix_sync, frame_cpe_info] = ...
    align_pn_frame_common_phase(sync_pass.PN, sync_pass.Frame, sync_pass.phase_index);

fprintf('公共相位补偿有效帧数：%d\n', size(PN_sym_matrix_align, 1));

%% 7. 可选 SVD 去噪与数据包保存
fprintf('\n[7] 可选 SVD 去噪并保存数据包。\n');

if cfg.use_pn_svd_denoise
    [PN_sym_denoised, svd_info] = SVD_Denoise_PN_Matrix(PN_sym_matrix_align, cfg.rank_keep);
    svd_info.enabled = true;
    fprintf('已启用 PN SVD 去噪。\n');
    fprintf('SVD 去噪后 PN 矩阵尺寸：%d x %d\n', size(PN_sym_denoised, 1), size(PN_sym_denoised, 2));
    fprintf('SVD 保留秩：%d\n', svd_info.rank_keep);
else
    PN_sym_denoised = PN_sym_matrix_align;
    svd_info = struct();
    svd_info.enabled = false;
    svd_info.rank_keep = [];
    svd_info.energy_ratio = [];
    svd_info.note = 'PN SVD 去噪已关闭，后续信道均衡直接使用逐帧公共相位补偿后的 PN 矩阵。';
    fprintf('已关闭 PN SVD 去噪，PN 矩阵直接送入后续信道均衡。\n');
    fprintf('未去噪 PN 矩阵尺寸：%d x %d\n', size(PN_sym_denoised, 1), size(PN_sym_denoised, 2));
end

Frame_sym_matrix_for_equalization = Frame_sym_matrix_sync;
full_frame_svd_info = [];
if cfg.use_pn_svd_denoise && cfg.use_full_frame_svd_denoise
    fprintf('启用实验性整帧矩阵 SVD 去噪。\n');
    [Frame_sym_matrix_for_equalization, full_frame_svd_info] = ...
        SVD_Denoise_PN_Matrix(Frame_sym_matrix_sync, cfg.full_frame_rank_keep);
    fprintf('SVD 去噪后整帧矩阵尺寸：%d x %d\n', ...
        size(Frame_sym_matrix_for_equalization, 1), size(Frame_sym_matrix_for_equalization, 2));
    fprintf('整帧 SVD 保留秩：%d\n', full_frame_svd_info.rank_keep);
elseif cfg.use_full_frame_svd_denoise
    fprintf('PN SVD 去噪已关闭，同时跳过实验性整帧矩阵 SVD 去噪，帧体均衡使用公共相位补偿后的整帧矩阵。\n');
else
    fprintf('未启用整帧矩阵 SVD 去噪，帧体均衡继续使用公共相位补偿后的整帧矩阵。\n');
end

[cfg.shared_qam_rms, shared_qam_scale_info] = Compute_DTMB_Shared_QAM_Scale(Frame_sym_matrix_for_equalization, cfg);
cfg.shared_qam_scale_reference = shared_qam_scale_info.reference;
fprintf('统一 QAM shared scale：%.9g，来源：同步/SVD 后 Frame 帧体频域有效数据子载波全局 RMS，排除 TPS %d 个。\n', ...
    cfg.shared_qam_rms, shared_qam_scale_info.num_excluded_tps_bins);

if cfg.plot_diagnostics
    plot_required_sync_diagnostics(cfo_loop_info, first_pass.PN, sync_pass.PN, ...
        PN_sym_matrix_align, PN_sym_denoised, svd_info);
    save_diagnostic_figures(cfg);
end

PN_svd_denoised_packet = single(PN_sym_denoised);
Frame_equalization_packet = single(Frame_sym_matrix_for_equalization);
Frame_svd_denoised_packet = Frame_equalization_packet;
PN_phase_index_packet = sync_pass.phase_index(:);
PN_frame_heads_packet = sync_pass.frame_heads(:);
PN_corr_peak_packet = sync_pass.corr_peak(:);

sync_diagnostics = struct();
sync_diagnostics.cfg = cfg;
sync_diagnostics.detected_mode = detected_mode;
sync_diagnostics.initial_valid_frames = length(coarse_frame_starts);
sync_diagnostics.first_pass = strip_large_fields(first_pass);
sync_diagnostics.final_pass = strip_large_fields(sync_pass);
sync_diagnostics.cfo_loop_info = cfo_loop_info;
sync_diagnostics.used_pn945_initial_cfo_sync = used_initial_cfo_sync;
sync_diagnostics.initial_cfo_info = strip_large_fields(rmfield_if_present(initial_cfo_info, {'rx_comp','rx_frac'}));
sync_diagnostics.final_signal_state = struct( ...
    'rx_matched_length', numel(rx_matched), ...
    'rx_matched_uncomp_length', numel(rx_matched_uncomp), ...
    'rx_16x_length', numel(rx_16x), ...
    'exact_frame_start', exact_frame_start, ...
    'phase_index', phase_index, ...
    'timing_info', strip_large_fields(timing_info));
sync_diagnostics.frame_cpe_info = frame_cpe_info;
sync_diagnostics.svd_info = rmfield(svd_info, intersect(fieldnames(svd_info), {'U','S','V'}));
sync_diagnostics.shared_qam_scale_info = shared_qam_scale_info;
if ~isempty(full_frame_svd_info)
    sync_diagnostics.full_frame_svd_info = rmfield(full_frame_svd_info, intersect(fieldnames(full_frame_svd_info), {'U','S','V'}));
else
    sync_diagnostics.full_frame_svd_info = [];
end

PN_packet_info = struct();
PN_packet_info.detected_mode = detected_mode;
PN_packet_info.pn_length = size(PN_svd_denoised_packet, 2);
PN_packet_info.num_frames = size(PN_svd_denoised_packet, 1);
PN_packet_info.use_pn_svd_denoise = cfg.use_pn_svd_denoise;
PN_packet_info.rank_keep = cfg.rank_keep;
PN_packet_info.use_full_frame_svd_denoise = cfg.use_full_frame_svd_denoise;
PN_packet_info.full_frame_rank_keep = cfg.full_frame_rank_keep;
PN_packet_info.use_joint_pn_cfo_refine = cfg.use_joint_pn_cfo_refine;
PN_packet_info.total_residual_cfo_comp_hz = cfo_loop_info.total_applied_cfo_hz;
PN_packet_info.final_residual_cfo_hz = get_field_or_nan(cfo_loop_info.final_residual_cfo, 'res_cfo_hz');
PN_packet_info.shared_qam_rms = cfg.shared_qam_rms;
PN_packet_info.shared_qam_scale_reference = cfg.shared_qam_scale_reference;
if cfg.use_pn_svd_denoise
    PN_packet_info.description = ...
        '经过多帧 PN 残留 CFO 闭环精同步和 PN SVD 去噪后的 DTMB PN945 数据包';
else
    PN_packet_info.description = ...
        '经过多帧 PN 残留 CFO 闭环精同步但未执行 PN SVD 去噪的 DTMB PN945 数据包';
end

save(cfg.output_packet, ...
    'PN_svd_denoised_packet', ...
    'Frame_equalization_packet', ...
    'Frame_svd_denoised_packet', ...
    'PN_phase_index_packet', ...
    'PN_frame_heads_packet', ...
    'PN_corr_peak_packet', ...
    'PN_packet_info', ...
    'sync_diagnostics', ...
    '-v7.3');

fprintf('数据包已保存：%s\n', cfg.output_packet);

%% 8. 可选 TC-AOMP
if cfg.run_tc_aomp
    fprintf('\n[8] 执行 TC-AOMP 稀疏信道估计与帧体均衡。\n');
    section_timer = tic;
    tc_result = run_tc_aomp_section(PN_sym_denoised, Frame_sym_matrix_for_equalization, sync_pass.phase_index, cfg);
    fprintf('[耗时] TC-AOMP 与 NLMS 均衡总耗时：%.2f 秒\n', toc(section_timer));
    if cfg.plot_diagnostics
        save_timer = tic;
        save_diagnostic_figures(cfg);
        fprintf('[耗时] 最终诊断图保存耗时：%.2f 秒\n', toc(save_timer));
    end
else
    fprintf('\n[8] 已关闭 TC-AOMP。\n');
    tc_result = [];
end

%% 9. 完整解调—再调制重构
if cfg.run_full_demod_remod
    fprintf('\n[9] 执行频域/时域解交织、FEC译码与再调制重构。\n');
    section_timer = tic;
    if isempty(tc_result) || ~isfield(tc_result, 'equalized_reference_frames') || isempty(tc_result.equalized_reference_frames)
        error('完整解调—再调制需要 TC-AOMP 均衡后的参考信号，请先开启 cfg.run_tc_aomp。');
    end
    full_recon_result = run_full_demod_remod_section( ...
        tc_result.equalized_reference_frames, tc_result.phase_index_list, cfg);
    fprintf('[耗时] 完整解调—再调制总耗时：%.2f 秒\n', toc(section_timer));
else
    fprintf('\n[9] 已关闭完整解调—再调制重构。\n');
    full_recon_result = [];
end

%% 10. 可选副峰消除
if cfg.run_ambiguity_sidelobe_suppression
    fprintf('\n[10] 执行论文第四章副峰消除，输入为 TC-AOMP 均衡后的参考信号。\n');
    section_timer = tic;
    if isempty(tc_result) || ~isfield(tc_result, 'equalized_reference_frames') || isempty(tc_result.equalized_reference_frames)
        error('副峰消除需要 TC-AOMP 均衡后的参考信号，请先开启 cfg.run_tc_aomp。');
    end
    DTMB_Ambiguity_Sidelobe_Suppression(tc_result.equalized_reference_frames, cfg);
    fprintf('[耗时] 副峰消除总耗时：%.2f 秒\n', toc(section_timer));
    if cfg.plot_diagnostics
        save_timer = tic;
        save_diagnostic_figures(cfg);
        fprintf('[耗时] 保存副峰消除诊断图耗时：%.2f 秒\n', toc(save_timer));
    end
else
    fprintf('\n[10] 已关闭副峰消除。\n');
end
fprintf('[耗时] 主流程总耗时：%.2f 秒\n', toc(total_timer));
fprintf('\n完成：联合时频同步主流程结束。\n');

pipeline_result = struct('cfg', cfg, 'tc_result', tc_result, ...
    'full_recon_result', full_recon_result);

end

%% 局部辅助函数

function cfg = apply_cfg_override(cfg, override)
    if isempty(override)
        return;
    end
    if ~isstruct(override) || ~isscalar(override)
        error('cfg_override 必须是标量 struct。');
    end
    names = fieldnames(override);
    for idx = 1:numel(names)
        cfg.(names{idx}) = override.(names{idx});
    end
end







