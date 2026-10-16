function report = refine_DTMB_groupwise_equalization(sync_packet, reconstruction_packet, ...
        common_report_packet, output_dir)
% refine_DTMB_groupwise_equalization 对公共信道优化结果进行分组时变跟踪。
%
% 每个候选均用LDPC校验矩阵的硬综合比例评价。分组更新和后续迭代只有在
% 该指标下降时才接受，从而抑制判决反馈造成的发散。

    if nargin < 1 || isempty(sync_packet)
        sync_packet = 'DTMB_1100Frame_Sync_Equalized_Packet.mat';
    end
    if nargin < 2 || isempty(reconstruction_packet)
        reconstruction_packet = 'DTMB_Full_Demod_Remod_Reconstruction.mat';
    end
    if nargin < 3 || isempty(common_report_packet)
        common_report_packet = fullfile('DTMB_equalization_optimization', ...
            'DTMB_Equalization_Optimization_Report.mat');
    end
    if nargin < 4 || isempty(output_dir)
        output_dir = 'DTMB_groupwise_equalization';
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    score_blocks = 512;
    group_sizes = [8 16 32 64 128];
    damping_candidates = [0.25 0.5 0.75 1.0];
    group_iterations = 4;
    channel_taps = 160;

    fprintf('[分组优化 1/6] 读取公共信道结果与原始帧体。\n');
    S = load(sync_packet, 'Frame_equalization_packet', ...
        'PN_phase_index_packet', 'sync_diagnostics');
    R = load(reconstruction_packet, 'result');
    C = load(common_report_packet, 'report');
    num_frames = min([size(S.Frame_equalization_packet, 1), ...
        numel(S.PN_phase_index_packet), R.result.num_frames, C.report.num_frames]);
    raw_frames = double(S.Frame_equalization_packet(1:num_frames, :));
    raw_fd = fft(raw_frames(:, 946:4725), 3780, 2);
    system_info = R.result.system_info;
    system_symbols = system_info.symbols;
    H_common = C.report.H_optimized(:).';
    lambda = C.report.best_lambda_absolute;
    check_model = make_check_model(system_info.rate_index);

    fprintf('[分组优化 2/6] 复核公共信道基线并生成初始判决。\n');
    common_fd = mmse_equalize(raw_fd, H_common, lambda);
    [common_metric, common_aux] = score_equalized_fd(common_fd, system_symbols, ...
        system_info.M, check_model, score_blocks);
    predicted_fd = build_predicted_fd(common_aux.data_normalized, system_symbols);
    fprintf('  公共信道：硬综合 %.6f，最近点EVM中位数 %.5f（%d块评价）。\n', ...
        common_metric, common_aux.evm_median, common_aux.num_score_blocks);

    fprintf('[分组优化 3/6] 搜索分组长度与更新阻尼。\n');
    grid_group_size = repelem(group_sizes(:), numel(damping_candidates));
    grid_alpha = repmat(damping_candidates(:), numel(group_sizes), 1);
    grid_metric = nan(size(grid_alpha));
    grid_evm = nan(size(grid_alpha));
    best_metric = common_metric;
    best_evm = common_aux.evm_median;
    best_group_size = num_frames;
    best_alpha = 0;
    best_H_bank = {H_common};
    best_starts = 1;
    best_ends = num_frames;
    row = 0;
    for size_idx = 1:numel(group_sizes)
        group_size = group_sizes(size_idx);
        [starts, ends] = make_groups(num_frames, group_size);
        H_reference = repmat({H_common}, numel(starts), 1);
        H_estimated = estimate_group_bank(raw_fd, predicted_fd, H_reference, ...
            starts, ends, channel_taps);
        for alpha_idx = 1:numel(damping_candidates)
            row = row + 1;
            alpha = damping_candidates(alpha_idx);
            H_candidate = blend_channel_banks(H_reference, H_estimated, alpha);
            candidate_fd = equalize_channel_bank(raw_fd, H_candidate, ...
                starts, ends, lambda);
            [grid_metric(row), aux] = score_equalized_fd(candidate_fd, ...
                system_symbols, system_info.M, check_model, score_blocks);
            grid_evm(row) = aux.evm_median;
            fprintf('  %3d帧/组，阻尼 %.2f：硬综合 %.6f，EVM %.5f。\n', ...
                group_size, alpha, grid_metric(row), grid_evm(row));
            if grid_metric(row) < best_metric - 1e-5
                best_metric = grid_metric(row);
                best_evm = grid_evm(row);
                best_group_size = group_size;
                best_alpha = alpha;
                best_H_bank = H_candidate;
                best_starts = starts;
                best_ends = ends;
            end
        end
    end

    fprintf('[分组优化 4/6] 对最优分组继续受保护迭代。\n');
    iteration_number = 0;
    iteration_alpha = best_alpha;
    iteration_metric = best_metric;
    iteration_evm = best_evm;
    if best_alpha == 0
        fprintf('  所有分组候选均未改善公共信道基线，跳过分组迭代。\n');
    else
        for iteration = 1:group_iterations
            current_fd = equalize_channel_bank(raw_fd, best_H_bank, ...
                best_starts, best_ends, lambda);
            [~, current_aux] = score_equalized_fd(current_fd, system_symbols, ...
                system_info.M, check_model, score_blocks);
            predicted_fd = build_predicted_fd(current_aux.data_normalized, system_symbols);
            H_estimated = estimate_group_bank(raw_fd, predicted_fd, best_H_bank, ...
                best_starts, best_ends, channel_taps);

            candidate_metric = nan(numel(damping_candidates), 1);
            candidate_evm = nan(numel(damping_candidates), 1);
            candidate_banks = cell(numel(damping_candidates), 1);
            for alpha_idx = 1:numel(damping_candidates)
                alpha = damping_candidates(alpha_idx);
                candidate_banks{alpha_idx} = blend_channel_banks( ...
                    best_H_bank, H_estimated, alpha);
                candidate_fd = equalize_channel_bank(raw_fd, ...
                    candidate_banks{alpha_idx}, best_starts, best_ends, lambda);
                [candidate_metric(alpha_idx), aux] = score_equalized_fd( ...
                    candidate_fd, system_symbols, system_info.M, ...
                    check_model, score_blocks);
                candidate_evm(alpha_idx) = aux.evm_median;
            end
            [local_metric, local_idx] = min(candidate_metric);
            local_alpha = damping_candidates(local_idx);
            fprintf('  迭代 %d：最佳阻尼 %.2f，硬综合 %.6f（当前 %.6f）。\n', ...
                iteration, local_alpha, local_metric, best_metric);
            if local_metric < best_metric - 1e-5
                best_H_bank = candidate_banks{local_idx};
                best_metric = local_metric;
                best_evm = candidate_evm(local_idx);
                iteration_number(end+1, 1) = iteration; %#ok<AGROW>
                iteration_alpha(end+1, 1) = local_alpha; %#ok<AGROW>
                iteration_metric(end+1, 1) = local_metric; %#ok<AGROW>
                iteration_evm(end+1, 1) = best_evm; %#ok<AGROW>
            else
                fprintf('  更新未降低LDPC硬综合，停止分组迭代。\n');
                break;
            end
        end
    end

    fprintf('[分组优化 5/6] 运行全量FEC和再调制重构。\n');
    optimized_fd = equalize_channel_bank(raw_fd, best_H_bank, ...
        best_starts, best_ends, lambda);
    [optimized_metric, optimized_aux] = score_equalized_fd(optimized_fd, ...
        system_symbols, system_info.M, check_model, score_blocks);
    optimized_frames = raw_frames;
    optimized_frames(:, 946:4725) = ifft(optimized_fd, 3780, 2);
    cfg = S.sync_diagnostics.cfg;
    cfg.plot_diagnostics = false;
    cfg.save_diagnostic_figures = false;
    cfg.full_recon_output_packet = fullfile(output_dir, ...
        'DTMB_Groupwise_Demod_Remod_Reconstruction.mat');
    cfg.fec_search_all_commutator_phases = false;
    cfg.fec_commutator_direction = 1;
    cfg.fec_commutator_phase = 0;
    cfg.ldpc_algorithm = 'norm-min-sum';
    cfg.ldpc_min_sum_scaling_factor = 0.75;
    cfg.ldpc_noise_model = 'deinterleaved';
    cfg.ldpc_llr_scale = 0.25;
    optimized_result = run_full_demod_remod_section(optimized_frames, ...
        S.PN_phase_index_packet(1:num_frames), cfg);

    fprintf('[分组优化 6/6] 保存报告、均衡帧和诊断图。\n');
    grid_table = table(grid_group_size, grid_alpha, grid_metric, grid_evm, ...
        'VariableNames', {'group_size','alpha','syndrome_metric','evm_median'});
    iteration_table = table(iteration_number, iteration_alpha, iteration_metric, ...
        iteration_evm, 'VariableNames', {'iteration','accepted_alpha', ...
        'syndrome_metric','evm_median'});
    report = struct();
    report.num_frames = num_frames;
    report.score_blocks = score_blocks;
    report.common_metric = common_metric;
    report.common_evm_median = common_aux.evm_median;
    report.grid_table = grid_table;
    report.best_group_size = best_group_size;
    report.initial_group_alpha = best_alpha;
    report.iteration_table = iteration_table;
    report.optimized_metric = optimized_metric;
    report.optimized_evm_median = optimized_aux.evm_median;
    report.metric_improvement = common_metric - optimized_metric;
    report.fec_successful_blocks = nnz(optimized_result.fec_stats.parity_ok);
    report.fec_total_blocks = optimized_result.fec_stats.num_blocks;
    report.fec_valid_frames = nnz(optimized_result.fec_valid_frame_mask);
    report.group_starts = best_starts;
    report.group_ends = best_ends;
    report.H_group = best_H_bank;
    report.description = ['公共判决导向信道之后的分组时变信道跟踪；' ...
        '网格搜索和迭代均由LDPC硬综合保护。'];
    writetable(grid_table, fullfile(output_dir, 'groupwise_grid.csv'));
    writetable(iteration_table, fullfile(output_dir, 'groupwise_iterations.csv'));
    save(fullfile(output_dir, 'DTMB_Groupwise_Equalization_Report.mat'), ...
        'report', '-v7.3');
    save(fullfile(output_dir, 'DTMB_Groupwise_Equalized_Frames.mat'), ...
        'optimized_frames', 'best_H_bank', 'best_starts', 'best_ends', ...
        'lambda', '-v7.3');
    make_groupwise_figure(grid_table, iteration_table, common_metric, ...
        common_aux.evm_median, output_dir);

    fprintf('\n分组优化结论：\n');
    fprintf('  最优分组：%d帧/组，初始阻尼 %.2f。\n', ...
        best_group_size, best_alpha);
    fprintf('  硬综合：%.6f -> %.6f，下降 %.6f。\n', ...
        common_metric, optimized_metric, report.metric_improvement);
    fprintf('  EVM中位数：%.5f -> %.5f。\n', ...
        common_aux.evm_median, optimized_aux.evm_median);
    fprintf('  LDPC通过：%d/%d块；FEC有效帧：%d。\n', ...
        report.fec_successful_blocks, report.fec_total_blocks, ...
        report.fec_valid_frames);
    fprintf('  结果目录：%s\n', output_dir);
end

function [starts, ends] = make_groups(num_frames, group_size)
    starts = (1:group_size:num_frames).';
    ends = min(starts + group_size - 1, num_frames);
end

function H_bank = estimate_group_bank(raw_fd, predicted_fd, H_reference, ...
        starts, ends, Ltap)
    H_bank = cell(numel(starts), 1);
    for group_idx = 1:numel(starts)
        rows = starts(group_idx):ends(group_idx);
        H_bank{group_idx} = estimate_group_channel(raw_fd(rows, :), ...
            predicted_fd(rows, :), H_reference{group_idx}, Ltap);
    end
end

function H_estimated = estimate_group_channel(raw_fd, predicted_fd, H_reference, Ltap)
    predicted_received = predicted_fd .* H_reference;
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
    gain = (H_estimated * H_reference') / ...
        (H_estimated * H_estimated' + eps);
    H_estimated = gain * H_estimated;
end

function H_bank = blend_channel_banks(H_old, H_new, alpha)
    H_bank = cell(size(H_old));
    for group_idx = 1:numel(H_old)
        H_bank{group_idx} = (1-alpha)*H_old{group_idx} + ...
            alpha*H_new{group_idx};
    end
end

function eq_fd = equalize_channel_bank(raw_fd, H_bank, starts, ends, lambda)
    eq_fd = complex(zeros(size(raw_fd)));
    for group_idx = 1:numel(starts)
        rows = starts(group_idx):ends(group_idx);
        eq_fd(rows, :) = mmse_equalize(raw_fd(rows, :), ...
            H_bank{group_idx}, lambda);
    end
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
    gain = (system_rx * conj(system_symbols)) / ...
        (system_symbols' * system_symbols + eps);
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

function model = make_check_model(rate_index)
    [H, meta] = DTMB_LDPC_Parity_Check_Matrix(rate_index);
    puncture_patterns = logical(dec2bin(0:31, 5) - '0').';
    puncture_syndromes = mod(double(H(:, 1:5)) * ...
        double(puncture_patterns), 2) ~= 0;
    affected = any(H(:, 1:5), 2);
    model.H_transmitted = H(:, 6:end);
    model.puncture_syndromes_affected = puncture_syndromes(affected, :);
    model.affected = affected;
    model.num_checks = meta.parity_length;
end

function metric = syndrome_metric_fast(transmitted_bits, model)
    syndrome = mod(double(model.H_transmitted) * ...
        double(transmitted_bits), 2) ~= 0;
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

function make_groupwise_figure(grid_table, iteration_table, common_metric, ...
        common_evm, output_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 800]);
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    groups = unique(grid_table.group_size, 'stable');
    for idx = 1:numel(groups)
        mask = grid_table.group_size == groups(idx);
        plot(grid_table.alpha(mask), grid_table.syndrome_metric(mask), ...
            '-o', 'LineWidth', 1.1, 'DisplayName', ...
            sprintf('%d帧/组', groups(idx))); hold on;
    end
    yline(common_metric, '--k', '公共信道基线'); hold off; grid on;
    xlabel('分组信道更新阻尼'); ylabel('LDPC硬综合比例');
    title('分组长度与阻尼搜索'); legend('Location', 'best');
    nexttile;
    yyaxis left;
    plot(iteration_table.iteration, iteration_table.syndrome_metric, ...
        '-o', 'LineWidth', 1.2); ylabel('LDPC硬综合比例');
    yyaxis right;
    plot(iteration_table.iteration, iteration_table.evm_median, ...
        '-s', 'LineWidth', 1.2); hold on;
    yline(common_evm, '--', '公共信道EVM'); hold off;
    ylabel('最近点EVM中位数'); xlabel('已接受的分组迭代');
    grid on; title('受保护的分组判决导向迭代');
    exportgraphics(fig, fullfile(output_dir, 'groupwise_optimization.png'), ...
        'Resolution', 180);
    close(fig);
end
