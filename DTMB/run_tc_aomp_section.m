function tc_result = run_tc_aomp_section(PN_sym_denoised, Frame_sym_matrix_sync, phase_index_use, cfg)
    num_frames_total = size(PN_sym_denoised, 1);
    phase_index_list = phase_index_use(:);

    if length(phase_index_list) ~= num_frames_total
        warning('PN 相位号数量 (%d) 与 PN 矩阵行数 (%d) 不一致，将截断到共同长度。', ...
            length(phase_index_list), num_frames_total);
        num_frames_total = min(num_frames_total, length(phase_index_list));
        PN_sym_denoised = PN_sym_denoised(1:num_frames_total, :);
        Frame_sym_matrix_sync = Frame_sym_matrix_sync(1:num_frames_total, :);
        phase_index_list = phase_index_list(1:num_frames_total);
    end

    fprintf('前若干个 PN 相位号：\n');
    disp(phase_index_list(1:min(10,end)).');

    Ltap = 160;

    opts_tc = struct();
    opts_tc.Kmax = 30;
    opts_tc.main_half_width = 2;
    opts_tc.lambda_eq = cfg.body_eq_lambda;
    opts_tc.min_delta_islr = 0.2;
    opts_tc.patience = 5;
    opts_tc.min_iter = 1;

    opts_tc.islr_opts = struct();
    opts_tc.islr_opts.main_half_width = 3;
    opts_tc.islr_opts.roi_half_width = 40;
    if isfield(cfg, 'tc_use_average_channel')
        opts_tc.use_average_channel = cfg.tc_use_average_channel;
    else
        opts_tc.use_average_channel = true;
    end

    if isfield(cfg, 'tc_num_frames') && ~isempty(cfg.tc_num_frames)
        num_frames_output = min(cfg.tc_num_frames, num_frames_total);
    else
        num_frames_output = num_frames_total;
    end
    if isfield(cfg, 'tc_channel_estimation_frames') && ~isempty(cfg.tc_channel_estimation_frames)
        num_frames_proc = min(cfg.tc_channel_estimation_frames, num_frames_output);
    else
        num_frames_proc = num_frames_output;
    end

    h_hat_all = zeros(num_frames_proc, Ltap);
    best_islr_all = zeros(num_frames_proc, 1);
    req_cell = cell(num_frames_output, 1);

    eps_islr_all = zeros(num_frames_proc, 1);
    stop_iter_all = zeros(num_frames_proc, 1);
    best_iter_all = zeros(num_frames_proc, 1);
    stop_reason_all = strings(num_frames_proc, 1);
    tap_count_all = zeros(num_frames_proc, 1);
    support_cell = cell(num_frames_proc, 1);

    tc_loop_timer = tic;
    for m = 1:num_frames_proc
        fprintf('  TC-AOMP 处理第 %d / %d 帧\n', m, num_frames_proc);

        pn_svd_m = PN_sym_denoised(m, :).';
        local_pn_m = PN_gen(phase_index_list(m)).';
        local_pn_m = local_pn_m(:);

        y_m = Extract_IBI_Free_Observation(pn_svd_m, Ltap);
        Theta_m = Build_Sensing_Matrix_TC(local_pn_m, Ltap);
        rcs_frame = Frame_sym_matrix_sync(m, :).';

        opts_tc.eps_islr = Estimate_Epsilon_From_TC_NoiseFloor(rcs_frame, opts_tc.islr_opts);
        eps_islr_all(m) = opts_tc.eps_islr;

        out_m = TC_AOMP_Channel_Estimate(y_m, Theta_m, rcs_frame, opts_tc);

        stop_iter_all(m) = out_m.stop_iter;
        best_iter_all(m) = out_m.best_iter;
        stop_reason_all(m) = out_m.stop_reason;
        h_hat_all(m, :) = out_m.h_hat(:).';
        best_islr_all(m) = out_m.best_islr;
        tap_count_all(m) = numel(out_m.support);
        support_cell{m} = out_m.support(:).';
        req_cell{m} = out_m.req;
    end
    fprintf('[耗时] TC-AOMP 信道估计循环耗时：%.2f 秒\n', toc(tc_loop_timer));

    h_avg = [];
    avg_channel_tap_count = NaN;
    robust_average_info = struct('enabled',false,'method',"arithmetic_mean", ...
        'inlier_mask',false(num_frames_proc,1),'num_inliers',0);
    if opts_tc.use_average_channel
        valid_h = any(abs(h_hat_all) > 0, 2);
        if any(valid_h)
            use_robust_average = isfield(cfg,'tc_robust_average_channel') && ...
                cfg.tc_robust_average_channel;
            if use_robust_average
                [h_avg,robust_average_info] = robust_average_channel(h_hat_all,valid_h,cfg);
            else
                h_avg = mean(h_hat_all(valid_h, :), 1).';
                robust_average_info.inlier_mask = valid_h;
                robust_average_info.num_inliers = nnz(valid_h);
            end
            avg_channel_tap_count = nnz(abs(h_avg) > max(abs(h_avg)) * 1e-3);
            for m = 1:num_frames_output
                rcs_frame = Frame_sym_matrix_sync(m, :).';
                req_cell{m} = Equalize_By_Channel(rcs_frame, h_avg, opts_tc.lambda_eq);
            end
            fprintf('  公共平均信道均衡       ：已启用，%d 帧信道参与平均\n', sum(valid_h));
            fprintf('  公共平均信道应用帧数   ：%d\n', num_frames_output);
            fprintf('  公共平均信道有效抽头数 ：%d\n', avg_channel_tap_count);
            if robust_average_info.enabled
                fprintf('  鲁棒平均信道内点数     ：%d / %d\n', ...
                    robust_average_info.num_inliers,nnz(valid_h));
            end
        else
            h_avg = zeros(Ltap, 1);
            avg_channel_tap_count = 0;
            warning('TC-AOMP 未得到有效逐帧信道估计，无法启用公共平均信道均衡。');
        end
    else
        fprintf('  公共平均信道均衡       ：未启用，继续使用逐帧信道均衡输出\n');
    end

    fprintf('TC-AOMP 处理完成。\n');
    fprintf('  epsilon 停止帧数     ：%d / %d\n', sum(stop_reason_all == "epsilon"), num_frames_proc);
    fprintf('  no_improve 停止帧数  ：%d / %d\n', sum(stop_reason_all == "no_improve"), num_frames_proc);
    fprintf('  Kmax 停止帧数        ：%d / %d\n', sum(stop_reason_all == "Kmax"), num_frames_proc);
    fprintf('  平均停止迭代次数     ：%.2f\n', mean(stop_iter_all));
    fprintf('  平均有效抽头数       ：%.2f\n', mean(tap_count_all));
    fprintf('  有效抽头数范围       ：%d 到 %d\n', min(tap_count_all), max(tap_count_all));
    fprintf('  最佳 ISLR 范围       ：%.2f 到 %.2f dB\n', min(best_islr_all), max(best_islr_all));

    tc_summary_table = table((1:num_frames_proc).', eps_islr_all, best_islr_all, ...
        best_iter_all, stop_iter_all, tap_count_all, stop_reason_all, ...
        'VariableNames', {'frame','eps_islr','best_islr','best_iter','stop_iter','tap_count','stop_reason'});
    disp(tc_summary_table);

    if cfg.plot_diagnostics
        plot_tc_aomp_required_diagnostics(h_hat_all, Frame_sym_matrix_sync, req_cell, tap_count_all, support_cell, cfg);
    end


    if isfield(cfg, 'run_nlms_structured_equalizer') && cfg.run_nlms_structured_equalizer
        NLMS_Structure_Constrained_Equalize_DTMB(PN_sym_denoised, Frame_sym_matrix_sync, ...
            phase_index_list, h_hat_all, req_cell, opts_tc.islr_opts, cfg);
    end
    tc_result = struct();
    tc_result.equalized_reference_frames = req_cell;
    tc_result.h_hat_all = h_hat_all;
    tc_result.h_avg = h_avg;
    tc_result.use_average_channel = opts_tc.use_average_channel;
    tc_result.avg_channel_tap_count = avg_channel_tap_count;
    tc_result.best_islr_all = best_islr_all;
    tc_result.eps_islr_all = eps_islr_all;
    tc_result.best_iter_all = best_iter_all;
    tc_result.stop_iter_all = stop_iter_all;
    tc_result.tap_count_all = tap_count_all;
    tc_result.support_cell = support_cell;
    tc_result.stop_reason_all = stop_reason_all;
    tc_result.phase_index_list = phase_index_list(1:num_frames_output);
    tc_result.robust_average_info = robust_average_info;
end

function [h_avg,info] = robust_average_channel(h_all,valid_h,cfg)
% Reject support/delay outliers, then use a component-wise complex median.
idx = find(valid_h);
H = h_all(idx,:);
n = size(H,1);
shapeCorr = eye(n);
for a = 1:n
    for b = a+1:n
        shapeCorr(a,b) = abs(H(a,:)*H(b,:)')/(norm(H(a,:))*norm(H(b,:))+eps);
        shapeCorr(b,a) = shapeCorr(a,b);
    end
end
[~,medoid] = max(median(shapeCorr,2));
corrThreshold = get_cfg_local(cfg,'tc_robust_corr_threshold',0.99);
tapTolerance = get_cfg_local(cfg,'tc_robust_tap_tolerance',2);
[~,dominantTap] = max(abs(H),[],2);
referenceTap = dominantTap(medoid);
corrToMedoid = shapeCorr(:,medoid);
inlierLocal = corrToMedoid>=corrThreshold & ...
    abs(dominantTap-referenceTap)<=tapTolerance;
minInliers = max(1,ceil(get_cfg_local(cfg,'tc_robust_min_inlier_fraction',0.5)*n));
if nnz(inlierLocal)<minInliers
    [~,order] = sort(corrToMedoid,'descend');
    inlierLocal(order(1:minInliers)) = true;
end
Hin = H(inlierLocal,:);
h_avg = complex(median(real(Hin),1),median(imag(Hin),1)).';
if ~any(abs(h_avg)>0)
    h_avg = mean(Hin,1).';
end
mask = false(size(valid_h));
mask(idx(inlierLocal)) = true;
info = struct('enabled',true,'method',"complex_component_median", ...
    'inlier_mask',mask,'num_inliers',nnz(mask),'medoid_valid_index',medoid, ...
    'corr_to_medoid',corrToMedoid,'dominant_tap',dominantTap-1, ...
    'corr_threshold',corrThreshold,'tap_tolerance',tapTolerance);
end

function value = get_cfg_local(cfg,name,defaultValue)
if isfield(cfg,name) && ~isempty(cfg.(name)), value=cfg.(name); else, value=defaultValue; end
end

function plot_tc_aomp_required_diagnostics(h_hat_all, Frame_sym_matrix_sync, req_cell, tap_count_all, support_cell, cfg)
    num_frames_proc = size(h_hat_all, 1);
    if num_frames_proc == 0
        return;
    end

    % TC-AOMP 必要诊断图：指定帧信道估计、单帧频域星座、多帧频域星座。
    plot_frame_idx = cfg.plot_frame_idx;
    plot_frame_idx = max(1, min(plot_frame_idx, num_frames_proc));

    plot_specified_frame_channel_estimate(h_hat_all, plot_frame_idx, tap_count_all(plot_frame_idx));
    save_all_frame_channel_estimates(h_hat_all, tap_count_all, support_cell, cfg);
    plot_tc_reference_autocorrelation(Frame_sym_matrix_sync, req_cell, plot_frame_idx, cfg);
    plot_frame_body_freq_constellation_equalization(Frame_sym_matrix_sync, req_cell, plot_frame_idx, cfg);
    plot_multi_frame_body_freq_constellation_equalization(Frame_sym_matrix_sync, req_cell, cfg);
end

function plot_specified_frame_channel_estimate(h_hat_all, plot_frame_idx, tap_count)
    % 对应输出图：指定帧信道估计结果图。
    figure('Name', sprintf('第 %d 帧信道估计结果', plot_frame_idx), 'Color', 'w');
    stem(abs(h_hat_all(plot_frame_idx, :)), 'filled');
    grid on;
    xlabel('信道抽头序号');
    ylabel('|h|');
    title(sprintf('第 %d 帧 TC-AOMP 稀疏信道估计结果，有效抽头数 %d', plot_frame_idx, tap_count));
end

function save_all_frame_channel_estimates(h_hat_all, tap_count_all, support_cell, cfg)
    if isfield(cfg, 'plot_tc_each_frame_channel_estimates') && ~cfg.plot_tc_each_frame_channel_estimates
        return;
    end

    if isfield(cfg, 'figure_output_dir') && ~isempty(cfg.figure_output_dir)
        base_dir = cfg.figure_output_dir;
    else
        base_dir = fullfile(pwd, 'joint_sync_figures');
    end
    out_dir = fullfile(base_dir, 'tc_aomp_channel_estimates');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    num_frames_proc = size(h_hat_all, 1);
    for frame_idx = 1:num_frames_proc
        h_abs = abs(h_hat_all(frame_idx, :));
        fig = figure('Visible', 'off', 'Color', 'w', ...
            'Name', sprintf('TC-AOMP 第 %d 帧信道估计', frame_idx));
        stem(h_abs, 'filled');
        grid on;
        xlabel('信道抽头序号');
        ylabel('|h|');
        title(sprintf('第 %d 帧 TC-AOMP 信道估计，有效抽头数 %d', frame_idx, tap_count_all(frame_idx)));
        if ~isempty(support_cell{frame_idx})
            hold on;
            support_idx = support_cell{frame_idx};
            plot(support_idx, h_abs(support_idx), 'ro', 'MarkerSize', 5, 'LineWidth', 1.0);
            hold off;
        end
        saveas(fig, fullfile(out_dir, sprintf('tc_aomp_channel_frame_%03d.png', frame_idx)));
        close(fig);
    end
    fprintf('  每帧 TC-AOMP 信道估计图已保存：%s\n', out_dir);
end

function plot_tc_reference_autocorrelation(Frame_sym_matrix_sync, req_cell, plot_frame_idx, cfg)
    if plot_frame_idx > size(Frame_sym_matrix_sync, 1) || plot_frame_idx > numel(req_cell) || isempty(req_cell{plot_frame_idx})
        warning('第 %d 帧均衡结果不存在，无法绘制 TC-AOMP 参考信号自相关图。', plot_frame_idx);
        return;
    end

    if isfield(cfg, 'tc_autocorr_lag_half_width') && ~isempty(cfg.tc_autocorr_lag_half_width)
        lag_half_width = cfg.tc_autocorr_lag_half_width;
    else
        lag_half_width = 511;
    end

    if isfield(cfg, 'figure_output_dir') && ~isempty(cfg.figure_output_dir)
        base_dir = cfg.figure_output_dir;
    else
        base_dir = fullfile(pwd, 'joint_sync_figures');
    end
    if ~exist(base_dir, 'dir')
        mkdir(base_dir);
    end

    ref_before = Frame_sym_matrix_sync(plot_frame_idx, :).';
    ref_after = req_cell{plot_frame_idx}(:);

    [lag_before, tc_before_db, peak_before, aux_before] = calc_autocorr_db_for_plot(ref_before, cfg);
    [lag_after, tc_after_db, peak_after, aux_after] = calc_autocorr_db_for_plot(ref_after, cfg);

    fig = figure('Name', sprintf('第 %d 帧 TC-AOMP 参考信号自相关', plot_frame_idx), 'Color', 'w');

    subplot(2, 1, 1);
    plot_autocorr_panel(lag_before, tc_before_db, peak_before, aux_before, lag_half_width, ...
        sprintf('第 %d 帧均衡前参考信号自相关', plot_frame_idx));

    subplot(2, 1, 2);
    plot_autocorr_panel(lag_after, tc_after_db, peak_after, aux_after, lag_half_width, ...
        sprintf('第 %d 帧 TC-AOMP 均衡后参考信号自相关', plot_frame_idx));

    out_path = fullfile(base_dir, sprintf('tc_aomp_reference_autocorrelation_frame_%03d_lag_%d.png', ...
        plot_frame_idx, lag_half_width));
    saveas(fig, out_path);
    fprintf('  TC-AOMP 参考信号自相关图已保存：%s\n', out_path);
end

function [lag_rel, tc_db, peak_idx, aux] = calc_autocorr_db_for_plot(ref_signal, cfg)
    opts = struct();
    opts.main_half_width = 3;
    opts.roi_half_width = 40;
    if isstruct(cfg) && isfield(cfg, 'tc_autocorr_main_half_width') && ~isempty(cfg.tc_autocorr_main_half_width)
        opts.main_half_width = cfg.tc_autocorr_main_half_width;
    end
    if isstruct(cfg) && isfield(cfg, 'tc_autocorr_roi_half_width') && ~isempty(cfg.tc_autocorr_roi_half_width)
        opts.roi_half_width = cfg.tc_autocorr_roi_half_width;
    end

    [~, tc_pow, peak_idx, aux] = Calc_ISLR_DTMB(ref_signal, opts);
    tc_db = 10 * log10(tc_pow / (max(tc_pow) + eps) + eps);
    lag_rel = (1:numel(tc_pow)).' - peak_idx;
end

function plot_autocorr_panel(lag_rel, tc_db, peak_idx, aux, lag_half_width, title_text)
    show_mask = lag_rel >= -lag_half_width & lag_rel <= lag_half_width;
    plot(lag_rel(show_mask), tc_db(show_mask), 'b-', 'LineWidth', 1.0);
    grid on;
    xlim([-lag_half_width, lag_half_width]);
    ylim([max(min(tc_db(show_mask)) - 3, -120), 3]);
    xlabel('相对主峰延迟 lag / 符号');
    ylabel('归一化自相关功率 / dB');
    title(title_text);
    hold on;

    main_lag = lag_rel(aux.main_mask);
    roi_lag = lag_rel(aux.roi_mask);
    if ~isempty(roi_lag)
        xline(min(roi_lag), '--k', 'ROI', 'LineWidth', 0.8);
        xline(max(roi_lag), '--k', 'ROI', 'LineWidth', 0.8);
    end
    if ~isempty(main_lag)
        xline(min(main_lag), '--r', '主瓣', 'LineWidth', 0.9);
        xline(max(main_lag), '--r', '主瓣', 'LineWidth', 0.9);
    end
    plot(0, tc_db(peak_idx), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    hold off;
end
function plot_frame_body_freq_constellation_equalization(Frame_sym_matrix_sync, req_cell, plot_frame_idx, cfg)
    if plot_frame_idx > size(Frame_sym_matrix_sync, 1) || plot_frame_idx > numel(req_cell) || isempty(req_cell{plot_frame_idx})
        warning('第 %d 帧均衡结果不存在，无法绘制帧体频域星座图。', plot_frame_idx);
        return;
    end

    frame_after = req_cell{plot_frame_idx};
    if length(frame_after) < 4725
        warning('第 %d 帧长度不足，无法绘制帧体频域星座图。', plot_frame_idx);
        return;
    end

    [metric_before, metric_after] = get_frame_body_qam_metrics(Frame_sym_matrix_sync, req_cell, plot_frame_idx, cfg);

    fprintf('第 %d 帧频域帧体 64QAM 指标：均衡前 EVM = %.2f%%，均衡后 EVM = %.2f%%\n', ...
        plot_frame_idx, 100 * metric_before.evm_rms, 100 * metric_after.evm_rms);
    fprintf('第 %d 帧均衡后 64QAM 校正量：均衡路径 %s，正则 %.1e，公共相位 %.2f 度\n', ...
        plot_frame_idx, metric_after.method_name, metric_after.lambda_eq, ...
        rad2deg(metric_after.phase_rad));
    print_modulation_order_evm(metric_after.sym_corr);

    % 对应输出图：指定帧帧体频域均衡前后星座图。
    figure('Name', sprintf('第 %d 帧帧体频域均衡前后星座图', plot_frame_idx), 'Color', 'w');
    subplot(1,2,1);
    scatter(real(metric_before.sym_corr), imag(metric_before.sym_corr), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title(sprintf('第 %d 帧均衡前帧体频域星座，EVM %.2f%%', ...
        plot_frame_idx, 100 * metric_before.evm_rms));

    subplot(1,2,2);
    scatter(real(metric_after.sym_corr), imag(metric_after.sym_corr), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title(sprintf('第 %d 帧均衡后帧体频域星座，EVM %.2f%%', ...
        plot_frame_idx, 100 * metric_after.evm_rms));
end

function plot_multi_frame_body_freq_constellation_equalization(Frame_sym_matrix_sync, req_cell, cfg)
    num_frames_total = min(size(Frame_sym_matrix_sync, 1), numel(req_cell));
    if num_frames_total == 0
        return;
    end

    if isfield(cfg, 'multi_frame_constellation_frames')
        num_frames_plot = min(cfg.multi_frame_constellation_frames, num_frames_total);
    else
        num_frames_plot = min(10, num_frames_total);
    end

    sym_before_all = [];
    sym_after_all = [];
    evm_before_all = nan(num_frames_plot, 1);
    evm_after_all = nan(num_frames_plot, 1);
    used_frames = false(num_frames_plot, 1);

    for frame_idx = 1:num_frames_plot
        if frame_idx > numel(req_cell) || isempty(req_cell{frame_idx})
            continue;
        end

        try
            [metric_before, metric_after] = get_frame_body_qam_metrics_fast( ...
                Frame_sym_matrix_sync, req_cell, frame_idx, cfg);
        catch ME
            warning('第 %d 帧多帧星座图指标计算失败：%s', frame_idx, ME.message);
            continue;
        end

        sym_before_all = [sym_before_all; metric_before.sym_corr(:)]; %#ok<AGROW>
        sym_after_all = [sym_after_all; metric_after.sym_corr(:)]; %#ok<AGROW>
        evm_before_all(frame_idx) = metric_before.evm_rms;
        evm_after_all(frame_idx) = metric_after.evm_rms;
        used_frames(frame_idx) = true;
    end

    if isempty(sym_before_all) || isempty(sym_after_all)
        warning('没有可用于绘制多帧频域帧体星座图的有效帧。');
        return;
    end

    sym_before_plot = limit_constellation_points(sym_before_all, cfg);
    sym_after_plot = limit_constellation_points(sym_after_all, cfg);
    valid_before_evm = evm_before_all(used_frames);
    valid_after_evm = evm_after_all(used_frames);
    used_frame_idx = find(used_frames);

    fprintf('多帧频域帧体 64QAM 指标：叠加 %d 帧（第 %d 到第 %d 帧），均衡前平均 EVM = %.2f%%，均衡后平均 EVM = %.2f%%\n', ...
        numel(used_frame_idx), used_frame_idx(1), used_frame_idx(end), ...
        100 * mean(valid_before_evm), 100 * mean(valid_after_evm));

    % 对应输出图：多帧帧体频域均衡前后星座图。
    figure('Name', sprintf('前 %d 帧帧体频域均衡前后星座图', numel(used_frame_idx)), 'Color', 'w');
    subplot(1,2,1);
    scatter(real(sym_before_plot), imag(sym_before_plot), 2, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title(sprintf('前 %d 帧均衡前帧体频域星座，平均 EVM %.2f%%', ...
        numel(used_frame_idx), 100 * mean(valid_before_evm)));

    subplot(1,2,2);
    scatter(real(sym_after_plot), imag(sym_after_plot), 2, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title(sprintf('前 %d 帧均衡后帧体频域星座，平均 EVM %.2f%%', ...
        numel(used_frame_idx), 100 * mean(valid_after_evm)));
end

function sym_plot = limit_constellation_points(sym, cfg)
    sym = sym(:);
    valid = isfinite(real(sym)) & isfinite(imag(sym));
    sym = sym(valid);

    if isfield(cfg, 'multi_frame_constellation_max_points')
        max_points = cfg.multi_frame_constellation_max_points;
    else
        max_points = 25000;
    end

    if numel(sym) > max_points
        pick_idx = round(linspace(1, numel(sym), max_points));
        sym_plot = sym(pick_idx);
    else
        sym_plot = sym;
    end
end

function print_modulation_order_evm(sym_corr)
    evm_qpsk = evaluate_square_qam_order_evm(sym_corr, 4);
    evm_16qam = evaluate_square_qam_order_evm(sym_corr, 16);
    evm_64qam = evaluate_square_qam_order_evm(sym_corr, 64);
    fprintf('均衡后调制阶数对比 EVM：QPSK %.2f%%，16QAM %.2f%%，64QAM %.2f%%\n', ...
        100 * evm_qpsk, 100 * evm_16qam, 100 * evm_64qam);
end

function evm_best = evaluate_square_qam_order_evm(sym, M)
    sym = sym(:);
    valid = isfinite(real(sym)) & isfinite(imag(sym));
    sym = sym(valid);

    phi_grid = linspace(-pi/4, pi/4, 91);
    evm_best = inf;
    for phi = phi_grid
        sym_rot = sym * exp(-1j * phi);
        sym_dec = Slice_Square_QAM(sym_rot, M);
        evm = rms(sym_rot - sym_dec) / (rms(sym_dec) + eps);
        evm_best = min(evm_best, evm);
    end
end
function sym_hat = Slice_Square_QAM(sym, M)
    switch M
        case 4
            levels = [-1 1];
        case 16
            levels = [-3 -1 1 3];
        case 64
            levels = [-7 -5 -3 -1 1 3 5 7];
        otherwise
            error('暂不支持 %dQAM 判决。', M);
    end

    avg_energy = mean(abs((levels(:) + 1j * levels(:).')).^2, 'all');
    norm_factor = sqrt(avg_energy);
    sym_u = sym(:) * norm_factor;

    I = real(sym_u);
    Q = imag(sym_u);
    I_hat = zeros(size(I));
    Q_hat = zeros(size(Q));

    for k = 1:numel(sym_u)
        [~, idxI] = min(abs(I(k) - levels));
        [~, idxQ] = min(abs(Q(k) - levels));
        I_hat(k) = levels(idxI);
        Q_hat(k) = levels(idxQ);
    end

    sym_hat = (I_hat + 1j * Q_hat) / norm_factor;
end

function [metric_before, metric_after] = get_frame_body_qam_metrics(Frame_sym_matrix_sync, req_cell, frame_idx, cfg)
    Lpn = 945;
    Ldata = 3780;
    Nfft = 3780;
    body_idx = Lpn + 1 : Lpn + Ldata;

    frame_before = Frame_sym_matrix_sync(frame_idx, :).';
    body_before_td = frame_before(body_idx) - mean(frame_before(body_idx));
    body_before_fd = fftshift(fft(body_before_td, Nfft));

    frame_after = req_cell{frame_idx}(:);
    body_after_td = frame_after(body_idx) - mean(frame_after(body_idx));
    body_after_fd = fftshift(fft(body_after_td, Nfft));

    lambda_eq = get_current_body_eq_lambda(cfg);
    idx_valid = make_dtmb_body_subcarrier_mask(Nfft, cfg);

    metric_before = evaluate_64qam_constellation_fast(body_before_fd(idx_valid), cfg);
    metric_before.method_name = '均衡前';
    metric_before.lambda_eq = lambda_eq;

    metric_after = evaluate_64qam_constellation_fast(body_after_fd(idx_valid), cfg);
    metric_after.method_name = 'TC-AOMP 输出 req';
    metric_after.lambda_eq = lambda_eq;
end
function [metric_before, metric_after] = get_frame_body_qam_metrics_fast(Frame_sym_matrix_sync, req_cell, frame_idx, cfg)
    Lpn = 945;
    Ldata = 3780;
    Nfft = 3780;
    body_idx = Lpn + 1 : Lpn + Ldata;

    frame_before = Frame_sym_matrix_sync(frame_idx, :).';
    body_before_td = frame_before(body_idx) - mean(frame_before(body_idx));
    body_before_fd = fftshift(fft(body_before_td, Nfft));

    lambda_eq = get_current_body_eq_lambda(cfg);
    frame_after = req_cell{frame_idx}(:);
    body_after_td = frame_after(body_idx) - mean(frame_after(body_idx));
    body_after_fd = fftshift(fft(body_after_td, Nfft));

    idx_valid = make_dtmb_body_subcarrier_mask(Nfft, cfg);
    metric_before = evaluate_64qam_constellation_fast(body_before_fd(idx_valid), cfg);
    metric_before.method_name = '均衡前';
    metric_before.lambda_eq = lambda_eq;

    metric_after = evaluate_64qam_constellation_fast(body_after_fd(idx_valid), cfg);
    metric_after.method_name = 'TC-AOMP 输出 req';
    metric_after.lambda_eq = lambda_eq;
end
function metric = evaluate_64qam_constellation_fast(sym, cfg)
    sym = sym(:);
    valid = isfinite(real(sym)) & isfinite(imag(sym));
    sym = sym(valid);
    sym = Normalize_DTMB_QAM_Symbols(sym, cfg);
    phi = Estimate_CPE_64QAM(sym);
    sym_corr = sym * exp(-1j * phi);
    sym_dec = Slice_64QAM(sym_corr);

    metric = struct();
    metric.sym_corr = sym_corr;
    metric.sym_dec = sym_dec;
    metric.phase_rad = phi;
    metric.timing_samp = 0;
    metric.evm_rms = rms(sym_corr - sym_dec) / (rms(sym_dec) + eps);
end
function idx_valid = make_dtmb_body_subcarrier_mask(Nfft, cfg)
    include_tps = isfield(cfg, 'qam_eval_include_tps') && cfg.qam_eval_include_tps;
    idx_valid = DTMB_Body_Frequency_Mask(include_tps, Nfft);
end
function lambda_eq = get_current_body_eq_lambda(cfg)
    lambda_eq = cfg.body_eq_lambda;
end
