function report = diagnose_DTMB_FEC_alignment(sync_packet, reconstruction_packet, output_dir)
% diagnose_DTMB_FEC_alignment 穷举DTMB时域解交织相位和LDPC码块边界。
%
% 诊断分两级进行：
%   1) 两种换向方向 x 52个支路相位，使用多个完整帧评价LDPC硬综合；
%   2) 对最佳解交织输出穷举全部7488个比特边界，先用24块筛选，
%      再用256块复核最优候选。

    if nargin < 1 || isempty(sync_packet)
        sync_packet = 'DTMB_1100Frame_Sync_Equalized_Packet.mat';
    end
    if nargin < 2 || isempty(reconstruction_packet)
        reconstruction_packet = 'DTMB_Full_Demod_Remod_Reconstruction.mat';
    end
    if nargin < 3 || isempty(output_dir)
        output_dir = 'DTMB_alignment_diagnostic';
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    phase_search_frames = 12;
    boundary_screen_blocks = 24;
    boundary_refine_blocks = 256;
    top_boundary_count = 32;
    boundary_batch_size = 64;

    fprintf('[对齐诊断 1/5] 读取数据并重建TC-AOMP均衡输出。\n');
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

    num_frames = min([numel(tc.equalized_reference_frames), R.result.num_frames]);
    data_normalized = complex(zeros(num_frames, 3744, 'single'));
    fprintf('[对齐诊断 2/5] 逐帧恢复3744个归一化64QAM数据符号。\n');
    for frame_idx = 1:num_frames
        frame = tc.equalized_reference_frames{frame_idx}(:);
        body_fd = fftshift(fft(frame(946:4725), 3780));
        data_rx = DTMB_Frequency_Deinterleave(body_fd);
        phase = R.result.frame_phase_rad(frame_idx);
        scale = R.result.frame_qam_scale(frame_idx) + eps;
        data_normalized(frame_idx, :) = single((data_rx .* exp(-1j*phase) / scale).');
    end

    rate_index = R.result.system_info.rate_index;
    M = R.result.system_info.M;
    delay_frames = round(52 * 51 * M / 3744);
    check_model = make_check_model(rate_index);

    fprintf('[对齐诊断 3/5] 穷举2种方向 x 52个换向器相位。\n');
    search_input_frames = min(num_frames, delay_frames + phase_search_frames + 1);
    search_input = data_normalized(1:search_input_frames, :);
    directions = [1, -1];
    num_phase_candidates = numel(directions) * 52;
    direction_column = zeros(num_phase_candidates, 1);
    phase_column = zeros(num_phase_candidates, 1);
    phase_metric = nan(num_phase_candidates, 1);
    phase_mean_unsatisfied = nan(num_phase_candidates, 1);
    row = 0;
    for direction = directions
        for phase = 0:51
            row = row + 1;
            candidate = DTMB_Time_Deinterleave(search_input, M, direction, phase);
            use_frames = min(phase_search_frames, size(candidate, 1));
            hard_bits = hard_demap_fast(reshape(candidate(1:use_frames, :).', [], 1));
            blocks = reshape(hard_bits, 7488, []);
            [phase_metric(row), weights] = syndrome_metric_fast(blocks, check_model);
            direction_column(row) = direction;
            phase_column(row) = phase;
            phase_mean_unsatisfied(row) = mean(weights);
        end
        fprintf('  已完成方向 %+d 的52个相位。\n', direction);
    end
    phase_table = table(direction_column, phase_column, phase_metric, ...
        phase_mean_unsatisfied, 'VariableNames', ...
        {'direction','commutator_phase','syndrome_metric','mean_unsatisfied_checks'});
    phase_table = sortrows(phase_table, 'syndrome_metric');
    best_direction = phase_table.direction(1);
    best_phase = phase_table.commutator_phase(1);

    fprintf('[对齐诊断 4/5] 最佳方向 %+d、相位 %d；穷举7488个码块边界。\n', ...
        best_direction, best_phase);
    best_deinterleaved = DTMB_Time_Deinterleave(data_normalized, M, ...
        best_direction, best_phase);
    symbol_stream = reshape(best_deinterleaved.', [], 1);
    hard_stream = hard_demap_fast(symbol_stream);
    boundary_offsets = (0:7487).';
    screen_metric = nan(7488, 1);
    for batch_start = 1:boundary_batch_size:7488
        take = batch_start:min(batch_start + boundary_batch_size - 1, 7488);
        offsets = boundary_offsets(take).';
        blocks = extract_boundary_blocks(hard_stream, offsets, boundary_screen_blocks);
        [~, weights] = syndrome_metric_fast(blocks, check_model);
        weights = reshape(weights, numel(offsets), boundary_screen_blocks);
        screen_metric(take) = mean(weights, 2) / check_model.num_checks;
        if mod(batch_start - 1, 1024) == 0
            fprintf('  边界扫描进度：%d/7488。\n', min(take(end), 7488));
        end
    end
    boundary_screen_table = table(boundary_offsets, floor(boundary_offsets/6), ...
        mod(boundary_offsets, 6), screen_metric, 'VariableNames', ...
        {'bit_offset','symbol_offset','in_symbol_bit_offset','screen_syndrome_metric'});
    boundary_screen_table = sortrows(boundary_screen_table, 'screen_syndrome_metric');

    top_count = min(top_boundary_count, height(boundary_screen_table));
    % 标称0边界必须进入长块复核，避免3块初筛的统计波动将它漏掉。
    refine_offsets = unique([0; boundary_screen_table.bit_offset(1:top_count)], 'stable').';
    refine_count = numel(refine_offsets);
    refine_blocks = extract_boundary_blocks(hard_stream, refine_offsets, boundary_refine_blocks);
    [~, refine_weights] = syndrome_metric_fast(refine_blocks, check_model);
    refine_weights = reshape(refine_weights, refine_count, boundary_refine_blocks);
    refine_metric = mean(refine_weights, 2) / check_model.num_checks;
    boundary_refine_table = table(refine_offsets.', floor(refine_offsets.'/6), ...
        mod(refine_offsets.', 6), screen_metric(refine_offsets.' + 1), ...
        refine_metric, 'VariableNames', {'bit_offset','symbol_offset', ...
        'in_symbol_bit_offset','screen_syndrome_metric','refined_syndrome_metric'});
    boundary_refine_table = sortrows(boundary_refine_table, 'refined_syndrome_metric');

    nominal_blocks = extract_boundary_blocks(hard_stream, 0, boundary_refine_blocks);
    nominal_metric = syndrome_metric_fast(nominal_blocks, check_model);
    best_boundary = boundary_refine_table.bit_offset(1);
    best_boundary_metric = boundary_refine_table.refined_syndrome_metric(1);

    fprintf('[对齐诊断 5/5] 保存表格、图和MAT报告。\n');
    report = struct();
    report.num_input_frames = num_frames;
    report.rate_index = rate_index;
    report.M = M;
    report.delay_frames = delay_frames;
    report.phase_search_frames = phase_search_frames;
    report.phase_table = phase_table;
    report.best_direction = best_direction;
    report.best_commutator_phase = best_phase;
    report.best_phase_metric = phase_table.syndrome_metric(1);
    report.second_phase_metric = phase_table.syndrome_metric(2);
    report.phase_separation = report.second_phase_metric - report.best_phase_metric;
    report.boundary_screen_blocks = boundary_screen_blocks;
    report.boundary_refine_blocks = boundary_refine_blocks;
    report.boundary_screen_table = boundary_screen_table;
    report.boundary_refine_table = boundary_refine_table;
    report.best_bit_offset = best_boundary;
    report.best_symbol_offset = floor(best_boundary/6);
    report.best_in_symbol_bit_offset = mod(best_boundary, 6);
    report.best_boundary_metric = best_boundary_metric;
    report.nominal_boundary_metric = nominal_metric;
    report.boundary_improvement = nominal_metric - best_boundary_metric;
    report.standard_alignment_selected = best_direction == 1 && best_phase == 0 && best_boundary == 0;
    report.interleaver_misalignment_likely = report.boundary_improvement > 0.01 || ...
        best_direction ~= 1 || best_phase ~= 0;

    writetable(phase_table, fullfile(output_dir, 'phase_search.csv'));
    writetable(boundary_screen_table, fullfile(output_dir, 'boundary_search_all.csv'));
    writetable(boundary_refine_table, fullfile(output_dir, 'boundary_top_refined.csv'));
    make_diagnostic_figure(phase_table, boundary_offsets, screen_metric, ...
        boundary_refine_table, nominal_metric, boundary_screen_blocks, ...
        boundary_refine_blocks, output_dir);
    save(fullfile(output_dir, 'DTMB_Alignment_Boundary_Diagnostic.mat'), 'report', '-v7.3');

    fprintf('\n诊断结论：\n');
    fprintf('  最佳换向：方向 %+d，相位 %d，硬综合 %.6f（次优 %.6f）。\n', ...
        best_direction, best_phase, report.best_phase_metric, report.second_phase_metric);
    fprintf('  最佳边界：偏移 %d bit = %d symbol + %d bit，%d块综合 %.6f。\n', ...
        best_boundary, report.best_symbol_offset, report.best_in_symbol_bit_offset, ...
        boundary_refine_blocks, best_boundary_metric);
    fprintf('  标称边界：偏移 0 bit，%d块综合 %.6f；相对最佳差值 %.6f。\n', ...
        boundary_refine_blocks, nominal_metric, nominal_metric - best_boundary_metric);
    fprintf('  是否疑似交织/边界错位：%d。\n', report.interleaver_misalignment_likely);
    fprintf('  结果目录：%s\n', output_dir);
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

function [metric, weights] = syndrome_metric_fast(transmitted_bits, model)
    syndrome = mod(double(model.H_transmitted) * double(transmitted_bits), 2) ~= 0;
    fixed_weight = sum(syndrome(~model.affected, :), 1);
    affected_syndrome = syndrome(model.affected, :);
    best_affected = inf(1, size(transmitted_bits, 2));
    for pattern_idx = 1:size(model.puncture_syndromes_affected, 2)
        candidate_weight = sum(xor(affected_syndrome, ...
            model.puncture_syndromes_affected(:, pattern_idx)), 1);
        best_affected = min(best_affected, candidate_weight);
    end
    weights = fixed_weight + best_affected;
    metric = mean(weights) / model.num_checks;
end

function blocks = extract_boundary_blocks(stream, offsets, num_blocks)
    offsets = double(offsets(:).');
    bit_index = reshape((1:7488).', 7488, 1, 1) + ...
        reshape(offsets, 1, [], 1) + reshape((0:num_blocks-1)*7488, 1, 1, []);
    if max(bit_index, [], 'all') > numel(stream)
        error('用于边界复核的解交织比特流长度不足。');
    end
    blocks = reshape(stream(bit_index), 7488, []);
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

function make_diagnostic_figure(phase_table, offsets, boundary_metric, refined, nominal, screen_blocks, refine_blocks, output_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 850]);
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    hold on;
    for direction = [1, -1]
        rows = phase_table.direction == direction;
        [phase, order] = sort(phase_table.commutator_phase(rows));
        values = phase_table.syndrome_metric(rows);
        plot(phase, values(order), '-o', 'MarkerSize', 3, 'DisplayName', ...
            sprintf('方向 %+d', direction));
    end
    hold off; grid on; xlabel('换向器相位'); ylabel('LDPC硬综合比例');
    title('52支路相位穷举'); legend('Location', 'best');

    nexttile;
    plot(offsets, boundary_metric, 'Color', [0.2 0.45 0.75]); hold on;
    scatter(refined.bit_offset, refined.refined_syndrome_metric, 28, 'r', 'filled');
    yline(nominal, '--k', sprintf('标称边界%d块指标', refine_blocks));
    hold off; grid on; xlabel('码块起点偏移 / bit'); ylabel('LDPC硬综合比例');
    title(sprintf('全部7488个码块边界（蓝：%d块穷举；红：%d块复核）', ...
        screen_blocks, refine_blocks));
    exportgraphics(fig, fullfile(output_dir, 'alignment_boundary_diagnostic.png'), 'Resolution', 180);
    close(fig);
end
