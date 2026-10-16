function report = optimize_DTMB_FEC_equalization(sync_packet, reconstruction_packet, output_dir)
% optimize_DTMB_FEC_equalization 以LDPC硬综合为目标优化DTMB帧体均衡。
%
% 比较原4725点整帧均衡与3780点帧体均衡，搜索相对正则系数，随后进行
% 有保护的判决导向公共信道迭代。只有降低LDPC硬综合的更新才会被接受。

    if nargin < 1 || isempty(sync_packet)
        sync_packet = 'DTMB_1100Frame_Sync_Equalized_Packet.mat';
    end
    if nargin < 2 || isempty(reconstruction_packet)
        reconstruction_packet = 'DTMB_Full_Demod_Remod_Reconstruction.mat';
    end
    if nargin < 3 || isempty(output_dir)
        output_dir = 'DTMB_equalization_optimization';
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    score_blocks = 256;
    decision_iterations = 10;
    damping_candidates = [0.25 0.5 0.75 1.0];
    lambda_relative_grid = [1e-8 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1];

    fprintf('[均衡优化 1/6] 读取同步数据并恢复TC-AOMP公共信道。\n');
    S = load(sync_packet, 'PN_svd_denoised_packet', 'Frame_equalization_packet', ...
        'PN_phase_index_packet', 'sync_diagnostics');
    R = load(reconstruction_packet, 'result');
    cfg = S.sync_diagnostics.cfg;
    cfg.plot_diagnostics = false;
    cfg.save_diagnostic_figures = false;
    cfg.tc_num_frames = [];
    cfg.tc_channel_estimation_frames = min(50, size(S.PN_svd_denoised_packet, 1));
    cfg.tc_use_average_channel = true;
    tc_log = evalc('tc = run_tc_aomp_section(S.PN_svd_denoised_packet, S.Frame_equalization_packet, S.PN_phase_index_packet, cfg);'); %#ok<NASGU>

    num_frames = min([size(S.Frame_equalization_packet, 1), ...
        numel(tc.equalized_reference_frames), R.result.num_frames]);
    raw_frames = double(S.Frame_equalization_packet(1:num_frames, :));
    raw_body = raw_frames(:, 946:4725);
    raw_fd = fft(raw_body, 3780, 2);
    system_info = R.result.system_info;
    system_symbols = system_info.symbols;
    check_model = make_check_model(system_info.rate_index);

    fprintf('[均衡优化 2/6] 评价原4725点整帧均衡基线。\n');
    tc_fd = complex(zeros(num_frames, 3780));
    for frame_idx = 1:num_frames
        frame = tc.equalized_reference_frames{frame_idx}(:);
        tc_fd(frame_idx, :) = fft(frame(946:4725), 3780).';
    end
    [baseline_metric, baseline_aux] = score_equalized_fd(tc_fd, system_symbols, ...
        system_info.M, check_model, score_blocks);
    fprintf('  原整帧均衡：硬综合 %.6f，最近点EVM中位数 %.5f。\n', ...
        baseline_metric, baseline_aux.evm_median);

    fprintf('[均衡优化 3/6] 搜索3780点帧体均衡正则系数。\n');
    H_initial = fft(double(tc.h_avg(:)), 3780).';
    channel_power_reference = median(abs(H_initial).^2);
    lambda_absolute_grid = lambda_relative_grid * max(channel_power_reference, eps);
    grid_metric = nan(numel(lambda_absolute_grid), 1);
    grid_evm = nan(numel(lambda_absolute_grid), 1);
    for idx = 1:numel(lambda_absolute_grid)
        eq_fd = mmse_equalize(raw_fd, H_initial, lambda_absolute_grid(idx));
        [grid_metric(idx), aux] = score_equalized_fd(eq_fd, system_symbols, ...
            system_info.M, check_model, score_blocks);
        grid_evm(idx) = aux.evm_median;
        fprintf('  lambda/median|H|^2=%-8.1e：硬综合 %.6f，EVM %.5f。\n', ...
            lambda_relative_grid(idx), grid_metric(idx), grid_evm(idx));
    end
    [best_grid_metric, best_grid_idx] = min(grid_metric);
    best_lambda = lambda_absolute_grid(best_grid_idx);
    best_lambda_relative = lambda_relative_grid(best_grid_idx);
    H_best = H_initial;
    best_metric = best_grid_metric;

    fprintf('[均衡优化 4/6] 判决导向公共信道迭代。\n');
    iteration_number = 0;
    accepted_alpha = 0;
    iteration_metric = best_metric;
    iteration_evm = grid_evm(best_grid_idx);
    for iteration = 1:decision_iterations
        current_fd = mmse_equalize(raw_fd, H_best, best_lambda);
        [~, current_aux] = score_equalized_fd(current_fd, system_symbols, ...
            system_info.M, check_model, score_blocks);
        predicted_fd = build_predicted_fd(current_aux.data_normalized, system_symbols);
        H_estimated = estimate_common_channel(raw_fd, predicted_fd, H_best, 160);

        candidate_metric = nan(numel(damping_candidates), 1);
        candidate_evm = nan(numel(damping_candidates), 1);
        for alpha_idx = 1:numel(damping_candidates)
            alpha = damping_candidates(alpha_idx);
            H_candidate = (1-alpha)*H_best + alpha*H_estimated;
            candidate_fd = mmse_equalize(raw_fd, H_candidate, best_lambda);
            [candidate_metric(alpha_idx), candidate_aux] = score_equalized_fd( ...
                candidate_fd, system_symbols, system_info.M, check_model, score_blocks);
            candidate_evm(alpha_idx) = candidate_aux.evm_median;
        end
        [local_metric, local_idx] = min(candidate_metric);
        local_alpha = damping_candidates(local_idx);
        fprintf('  迭代 %d：最佳阻尼 %.2f，硬综合 %.6f（当前 %.6f）。\n', ...
            iteration, local_alpha, local_metric, best_metric);
        if local_metric < best_metric - 1e-5
            H_best = (1-local_alpha)*H_best + local_alpha*H_estimated;
            best_metric = local_metric;
            iteration_number(end+1, 1) = iteration; %#ok<AGROW>
            accepted_alpha(end+1, 1) = local_alpha; %#ok<AGROW>
            iteration_metric(end+1, 1) = local_metric; %#ok<AGROW>
            iteration_evm(end+1, 1) = candidate_evm(local_idx); %#ok<AGROW>
        else
            fprintf('  更新没有降低LDPC硬综合，停止迭代。\n');
            break;
        end
    end

    fprintf('[均衡优化 5/6] 生成最优3780点均衡帧并运行完整FEC。\n');
    optimized_fd = mmse_equalize(raw_fd, H_best, best_lambda);
    [optimized_metric, best_aux] = score_equalized_fd(optimized_fd, system_symbols, ...
        system_info.M, check_model, score_blocks);
    optimized_frames = raw_frames;
    optimized_frames(:, 946:4725) = ifft(optimized_fd, 3780, 2);
    optimized_cell = mat2cell(optimized_frames, ones(num_frames, 1), 4725);

    cfg.full_recon_output_packet = fullfile(output_dir, ...
        'DTMB_Optimized_Demod_Remod_Reconstruction.mat');
    cfg.fec_search_all_commutator_phases = false;
    cfg.fec_commutator_direction = 1;
    cfg.fec_commutator_phase = 0;
    cfg.ldpc_algorithm = 'norm-min-sum';
    cfg.ldpc_min_sum_scaling_factor = 0.75;
    cfg.ldpc_noise_model = 'deinterleaved';
    cfg.ldpc_llr_scale = 0.25;
    optimized_result = run_full_demod_remod_section(optimized_cell, ...
        S.PN_phase_index_packet(1:num_frames), cfg);

    fprintf('[均衡优化 6/6] 保存报告与诊断图。\n');
    grid_table = table(lambda_relative_grid(:), lambda_absolute_grid(:), ...
        grid_metric, grid_evm, 'VariableNames', ...
        {'lambda_relative','lambda_absolute','syndrome_metric','evm_median'});
    iteration_table = table(iteration_number, accepted_alpha, iteration_metric, ...
        iteration_evm, 'VariableNames', {'iteration','accepted_alpha', ...
        'syndrome_metric','evm_median'});

    report = struct();
    report.num_frames = num_frames;
    report.score_blocks = score_blocks;
    report.baseline_metric = baseline_metric;
    report.baseline_evm_median = baseline_aux.evm_median;
    report.grid_table = grid_table;
    report.best_lambda_relative = best_lambda_relative;
    report.best_lambda_absolute = best_lambda;
    report.iteration_table = iteration_table;
    report.optimized_metric = optimized_metric;
    report.optimized_evm_median = best_aux.evm_median;
    report.metric_improvement = baseline_metric - optimized_metric;
    report.fec_successful_blocks = nnz(optimized_result.fec_stats.parity_ok);
    report.fec_total_blocks = optimized_result.fec_stats.num_blocks;
    report.fec_valid_frames = nnz(optimized_result.fec_valid_frame_mask);
    report.H_initial = H_initial;
    report.H_optimized = H_best;
    report.description = ['3780点帧体MMSE均衡与判决导向公共信道优化；' ...
        '全部候选用LDPC硬综合评价，只有改善的迭代被接受。'];

    writetable(grid_table, fullfile(output_dir, 'lambda_grid.csv'));
    writetable(iteration_table, fullfile(output_dir, 'decision_directed_iterations.csv'));
    save(fullfile(output_dir, 'DTMB_Equalization_Optimization_Report.mat'), ...
        'report', '-v7.3');
    save(fullfile(output_dir, 'DTMB_Optimized_Equalized_Frames.mat'), ...
        'optimized_frames', 'H_best', 'best_lambda', '-v7.3');
    make_optimization_figure(grid_table, iteration_table, baseline_metric, ...
        baseline_aux.evm_median, output_dir);

    fprintf('\n均衡优化结论：\n');
    fprintf('  硬综合：%.6f -> %.6f，下降 %.6f。\n', ...
        baseline_metric, optimized_metric, report.metric_improvement);
    fprintf('  EVM中位数：%.5f -> %.5f。\n', ...
        baseline_aux.evm_median, best_aux.evm_median);
    fprintf('  LDPC通过：%d/%d块；FEC有效帧：%d。\n', ...
        report.fec_successful_blocks, report.fec_total_blocks, report.fec_valid_frames);
    fprintf('  结果目录：%s\n', output_dir);
end

function eq_fd = mmse_equalize(raw_fd, H, lambda)
    eq_fd = raw_fd .* conj(H) ./ (abs(H).^2 + lambda);
end

function [metric, aux] = score_equalized_fd(eq_fd, system_symbols, M, model, score_blocks)
    [z_to_y, tps_idx] = DTMB_Frequency_Permutation();
    z = eq_fd(:, z_to_y);
    data_mask = true(3780, 1);
    data_mask(tps_idx) = false;
    system_rx = z(:, tps_idx);
    data = z(:, data_mask);
    system_symbols = system_symbols(:);
    gain = (system_rx * conj(system_symbols)) / (system_symbols' * system_symbols + eps);
    normalized = data .* exp(-1j*angle(gain));
    scale = sqrt(mean(abs(normalized).^2, 2)) + eps;
    normalized = normalized ./ scale;
    for iteration = 1:2
        decisions = Slice_64QAM(normalized);
        phase = angle(sum(conj(decisions).*normalized, 2));
        normalized = normalized .* exp(-1j*phase);
    end
    decisions = Slice_64QAM(normalized);
    evm = sqrt(mean(abs(normalized-decisions).^2, 2));

    deinterleaved = DTMB_Time_Deinterleave(single(normalized), M, 1, 0);
    use_blocks = min(score_blocks, size(deinterleaved, 1)*3);
    use_symbols = ceil(use_blocks * 7488 / 6);
    symbol_stream = reshape(deinterleaved.', [], 1);
    hard_bits = hard_demap_fast(symbol_stream(1:use_symbols));
    blocks = reshape(hard_bits, 7488, use_blocks);
    metric = syndrome_metric_fast(blocks, model);
    aux = struct('evm_per_frame', evm, 'evm_median', median(evm), ...
        'data_normalized', normalized, 'num_score_blocks', use_blocks, ...
        'num_score_symbols', use_symbols);
end

function predicted_fd = build_predicted_fd(data_normalized, system_symbols)
    num_frames = size(data_normalized, 1);
    [z_to_y, tps_idx] = DTMB_Frequency_Permutation();
    data_mask = true(3780, 1);
    data_mask(tps_idx) = false;
    z = complex(zeros(num_frames, 3780));
    z(:, data_mask) = Slice_64QAM(data_normalized);
    z(:, tps_idx) = repmat(system_symbols(:).', num_frames, 1);
    predicted_fd = complex(zeros(size(z)));
    predicted_fd(:, z_to_y) = z;
end

function H_estimated = estimate_common_channel(raw_fd, predicted_fd, H_current, Ltap)
    predicted_received = predicted_fd .* H_current;
    frame_gain = sum(conj(predicted_received).*raw_fd, 2) ./ ...
        (sum(abs(predicted_received).^2, 2) + eps);
    valid = isfinite(frame_gain) & abs(frame_gain) > eps;
    adjusted = raw_fd(valid, :) ./ frame_gain(valid);
    predicted = predicted_fd(valid, :);
    H_raw = sum(conj(predicted).*adjusted, 1) ./ ...
        (sum(abs(predicted).^2, 1) + eps);
    h_raw = ifft(H_raw, 3780);
    h_truncated = complex(zeros(1, 3780));
    h_truncated(1:Ltap) = h_raw(1:Ltap);
    H_estimated = fft(h_truncated, 3780);
    gain = (H_estimated * H_current') / (H_estimated * H_estimated' + eps);
    H_estimated = gain * H_estimated;
end

function model = make_check_model(rate_index)
    [H, meta] = DTMB_LDPC_Parity_Check_Matrix(rate_index);
    puncture_patterns = logical(dec2bin(0:31, 5) - '0').';
    puncture_syndromes = mod(double(H(:, 1:5)) * double(puncture_patterns), 2) ~= 0;
    affected = any(H(:, 1:5), 2);
    model.H_transmitted = H(:, 6:end);
    model.puncture_syndromes_affected = puncture_syndromes(affected, :);
    model.affected = affected;
    model.num_checks = meta.parity_length;
end

function metric = syndrome_metric_fast(transmitted_bits, model)
    syndrome = mod(double(model.H_transmitted) * double(transmitted_bits), 2) ~= 0;
    fixed_weight = sum(syndrome(~model.affected, :), 1);
    affected_syndrome = syndrome(model.affected, :);
    best_affected = inf(1, size(transmitted_bits, 2));
    for pattern_idx = 1:size(model.puncture_syndromes_affected, 2)
        candidate_weight = sum(xor(affected_syndrome, ...
            model.puncture_syndromes_affected(:, pattern_idx)), 1);
        best_affected = min(best_affected, candidate_weight);
    end
    metric = mean(fixed_weight + best_affected) / model.num_checks;
end

function bits = hard_demap_fast(symbols)
    symbols = symbols(:) * sqrt(42);
    levels = [-7 -5 -3 -1 1 3 5 7];
    labels = logical([0 0 0; 0 0 1; 0 1 1; 0 1 0; ...
                      1 1 0; 1 1 1; 1 0 1; 1 0 0]);
    [~, idx_i] = min(abs(real(symbols) - levels), [], 2);
    [~, idx_q] = min(abs(imag(symbols) - levels), [], 2);
    i_label = labels(idx_i, :);
    q_label = labels(idx_q, :);
    per_symbol = [i_label(:, 3:-1:1), q_label(:, 3:-1:1)];
    bits = reshape(per_symbol.', [], 1);
end

function make_optimization_figure(grid_table, iteration_table, baseline_metric, baseline_evm, output_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 800]);
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    semilogx(grid_table.lambda_relative, grid_table.syndrome_metric, '-o', ...
        'LineWidth', 1.2); hold on;
    yline(baseline_metric, '--r', '原4725点整帧均衡'); hold off; grid on;
    xlabel('lambda / median|H|^2'); ylabel('LDPC硬综合比例');
    title('3780点帧体均衡正则参数搜索');
    nexttile;
    yyaxis left;
    plot(iteration_table.iteration, iteration_table.syndrome_metric, '-o', 'LineWidth', 1.2);
    ylabel('LDPC硬综合比例');
    yyaxis right;
    plot(iteration_table.iteration, iteration_table.evm_median, '-s', 'LineWidth', 1.2);
    yline(baseline_evm, '--', '原均衡EVM'); ylabel('最近点EVM中位数');
    xlabel('已接受的判决导向迭代'); grid on; title('受保护的信道迭代');
    exportgraphics(fig, fullfile(output_dir, 'equalization_optimization.png'), 'Resolution', 180);
    close(fig);
end
