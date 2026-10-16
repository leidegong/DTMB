function report = optimize_DTMB_LDPC_soft_decoding(groupwise_frames_packet, ...
        groupwise_reconstruction_packet, output_dir)
% optimize_DTMB_LDPC_soft_decoding 搜索DTMB LDPC软信息和译码算法参数。
%
% 使用固定的分组均衡符号，先比较四类MATLAB LDPC译码器，再搜索LLR
% 全局尺度与逐符号噪声跟踪。小样本搜索后在全部1767个码块上复核。

    if nargin < 1 || isempty(groupwise_frames_packet)
        groupwise_frames_packet = fullfile('DTMB_groupwise_equalization', ...
            'DTMB_Groupwise_Equalized_Frames.mat');
    end
    if nargin < 2 || isempty(groupwise_reconstruction_packet)
        groupwise_reconstruction_packet = fullfile('DTMB_groupwise_equalization', ...
            'DTMB_Groupwise_Demod_Remod_Reconstruction.mat');
    end
    if nargin < 3 || isempty(output_dir)
        output_dir = 'DTMB_LDPC_soft_optimization';
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    evaluation_blocks = 256;
    search_iterations = 30;
    validation_iterations = 60;
    llr_scales = [0.25 0.5 0.75 1 1.5 2 3 4];

    fprintf('[软译码优化 1/5] 读取分组均衡符号和归一化参数。\n');
    G = load(groupwise_frames_packet, 'optimized_frames');
    Q = load(groupwise_reconstruction_packet, 'result');
    result = Q.result;
    num_frames = min(size(G.optimized_frames, 1), result.num_frames);
    frames = double(G.optimized_frames(1:num_frames, :));
    raw_fd = fft(frames(:, 946:4725), 3780, 2);
    [z_to_y, tps_idx] = DTMB_Frequency_Permutation();
    z = raw_fd(:, z_to_y);
    data_mask = true(3780, 1);
    data_mask(tps_idx) = false;
    data = z(:, data_mask);
    normalized = data .* exp(-1j*result.frame_phase_rad(1:num_frames));
    normalized = normalized ./ result.frame_qam_scale(1:num_frames);

    fprintf('[软译码优化 2/5] 精确传播解交织前的逐帧噪声可靠度。\n');
    M = result.system_info.M;
    direction = result.fec_alignment_info.direction;
    phase = result.fec_alignment_info.commutator_phase;
    [deinterleaved, time_info] = DTMB_Time_Deinterleave(single(normalized), ...
        M, direction, phase);
    source_noise = repmat(single(result.noise_variance(1:num_frames)), 1, 3744);
    deinterleaved_noise = DTMB_Time_Deinterleave(source_noise, M, direction, phase);
    num_pre_frames = size(deinterleaved, 1);
    num_blocks = num_pre_frames * 3;

    distance_llr = zeros(7488, num_blocks, 'single');
    legacy_noise = zeros(7488, num_blocks, 'single');
    tracked_noise = zeros(7488, num_blocks, 'single');
    delay_frames = time_info.delay_frames;
    for pre_idx = 1:num_pre_frames
        symbols = deinterleaved(pre_idx, :).';
        base_llr = DTMB_64QAM_Soft_Demap(symbols, 1);
        block_columns = (pre_idx - 1)*3 + (1:3);
        distance_llr(:, block_columns) = single(reshape(base_llr, 7488, 3));

        local_noise = median(result.noise_variance(max(1, pre_idx): ...
            min(num_frames, pre_idx + delay_frames)), 'omitnan');
        legacy_noise(:, block_columns) = single(local_noise);
        symbol_noise = double(deinterleaved_noise(pre_idx, :).');
        bit_noise = repelem(max(symbol_noise, 1e-4), 6);
        tracked_noise(:, block_columns) = single(reshape(bit_noise, 7488, 3));
    end
    legacy_llr = distance_llr ./ max(legacy_noise, 1e-4);
    tracked_llr = distance_llr ./ max(tracked_noise, 1e-4);
    use_blocks = min(evaluation_blocks, num_blocks);
    fprintf('  已生成 %d 个码块；参数搜索使用前 %d 块。\n', num_blocks, use_blocks);

    fprintf('[软译码优化 3/5] 比较BP、分层BP、归一化及偏移最小和译码。\n');
    algorithm = ["bp"; "layered-bp"; repmat("norm-min-sum", 5, 1); ...
        repmat("offset-min-sum", 5, 1)];
    algorithm_parameter = [NaN; NaN; 0.50; 0.65; 0.75; 0.85; 0.95; ...
        0.10; 0.25; 0.50; 0.75; 1.00];
    algorithm_passes = zeros(size(algorithm));
    algorithm_mean_iterations = zeros(size(algorithm_parameter));
    for idx = 1:numel(algorithm)
        [ok, iterations] = decode_parity_only(legacy_llr(:, 1:use_blocks), ...
            result.system_info.rate_index, search_iterations, algorithm(idx), ...
            algorithm_parameter(idx));
        algorithm_passes(idx) = nnz(ok);
        algorithm_mean_iterations(idx) = mean(iterations);
        fprintf('  %-14s 参数 %-4s：通过 %3d/%d，平均迭代 %.1f。\n', ...
            char(algorithm(idx)), parameter_text(algorithm_parameter(idx)), ...
            algorithm_passes(idx), use_blocks, algorithm_mean_iterations(idx));
    end
    algorithm_table = table(algorithm, algorithm_parameter, algorithm_passes, ...
        algorithm_passes/use_blocks, algorithm_mean_iterations, ...
        'VariableNames', {'algorithm','parameter','passed_blocks', ...
        'pass_rate','mean_iterations'});
    [~, algorithm_order] = sortrows([-algorithm_passes, algorithm_mean_iterations], [1 2]);
    best_algorithm_idx = algorithm_order(1);
    best_algorithm = algorithm(best_algorithm_idx);
    best_parameter = algorithm_parameter(best_algorithm_idx);

    fprintf('[软译码优化 4/5] 搜索噪声传播模型和LLR尺度。\n');
    noise_model = [repmat("legacy-window", numel(llr_scales), 1); ...
        repmat("deinterleaved", numel(llr_scales), 1)];
    llr_scale = repmat(llr_scales(:), 2, 1);
    scale_passes = zeros(size(llr_scale));
    scale_mean_iterations = zeros(size(llr_scale));
    for idx = 1:numel(llr_scale)
        if noise_model(idx) == "legacy-window"
            base = legacy_llr(:, 1:use_blocks);
        else
            base = tracked_llr(:, 1:use_blocks);
        end
        [ok, iterations] = decode_parity_only(llr_scale(idx)*base, ...
            result.system_info.rate_index, search_iterations, best_algorithm, ...
            best_parameter);
        scale_passes(idx) = nnz(ok);
        scale_mean_iterations(idx) = mean(iterations);
        fprintf('  %-13s，LLR×%-4.2g：通过 %3d/%d，平均迭代 %.1f。\n', ...
            char(noise_model(idx)), llr_scale(idx), scale_passes(idx), ...
            use_blocks, scale_mean_iterations(idx));
    end
    scale_table = table(noise_model, llr_scale, scale_passes, ...
        scale_passes/use_blocks, scale_mean_iterations, ...
        'VariableNames', {'noise_model','llr_scale','passed_blocks', ...
        'pass_rate','mean_iterations'});
    [~, scale_order] = sortrows([-scale_passes, scale_mean_iterations], [1 2]);
    best_scale_idx = scale_order(1);
    best_noise_model = noise_model(best_scale_idx);
    best_llr_scale = llr_scale(best_scale_idx);

    fprintf('[软译码优化 5/5] 在全部码块上复核最优配置。\n');
    if best_noise_model == "legacy-window"
        validation_llr = best_llr_scale * legacy_llr;
    else
        validation_llr = best_llr_scale * tracked_llr;
    end
    [validation_ok, validation_iteration] = decode_parity_only(validation_llr, ...
        result.system_info.rate_index, validation_iterations, best_algorithm, ...
        best_parameter);
    validation_passes = nnz(validation_ok);
    fprintf('  全量结果：%d/%d块通过（%.2f%%），平均迭代 %.1f。\n', ...
        validation_passes, num_blocks, 100*validation_passes/num_blocks, ...
        mean(validation_iteration));

    report = struct();
    report.num_frames = num_frames;
    report.num_blocks = num_blocks;
    report.evaluation_blocks = use_blocks;
    report.search_iterations = search_iterations;
    report.validation_iterations = validation_iterations;
    report.algorithm_table = algorithm_table;
    report.scale_table = scale_table;
    report.best_algorithm = best_algorithm;
    report.best_algorithm_parameter = best_parameter;
    report.best_noise_model = best_noise_model;
    report.best_llr_scale = best_llr_scale;
    report.validation_parity_ok = validation_ok;
    report.validation_iterations_per_block = validation_iteration;
    report.validation_passed_blocks = validation_passes;
    report.validation_pass_rate = validation_passes/num_blocks;
    report.previous_passed_blocks = nnz(result.fec_stats.parity_ok);
    report.description = ['固定分组均衡结果，搜索LDPC算法、噪声传播模型与' ...
        'LLR全局尺度；最终结果在全部码块上复核。'];
    writetable(algorithm_table, fullfile(output_dir, 'decoder_algorithm_grid.csv'));
    writetable(scale_table, fullfile(output_dir, 'llr_scale_grid.csv'));
    save(fullfile(output_dir, 'DTMB_LDPC_Soft_Optimization_Report.mat'), ...
        'report', '-v7.3');
    make_soft_decoding_figure(algorithm_table, scale_table, output_dir);

    fprintf('\n软译码优化结论：\n');
    fprintf('  最优算法：%s，参数 %s。\n', char(best_algorithm), ...
        parameter_text(best_parameter));
    fprintf('  最优噪声模型：%s，LLR尺度 %.3g。\n', ...
        char(best_noise_model), best_llr_scale);
    fprintf('  LDPC通过：%d -> %d/%d块。\n', ...
        report.previous_passed_blocks, validation_passes, num_blocks);
    fprintf('  结果目录：%s\n', output_dir);
end

function [parity_ok, iterations] = decode_parity_only(received_llr, rate_index, ...
        max_iterations, algorithm, parameter)
    [H_standard, meta] = DTMB_LDPC_Parity_Check_Matrix(rate_index);
    parity_length = meta.parity_length;
    H_matlab = H_standard(:, [parity_length+1:end, 1:parity_length]);
    decoder_cfg = ldpcDecoderConfig(H_matlab);
    decoder_cfg.Algorithm = char(algorithm);
    num_blocks = size(received_llr, 2);
    parity_ok = false(num_blocks, 1);
    iterations = zeros(num_blocks, 1);
    batch_size = 16;
    for batch_start = 1:batch_size:num_blocks
        take = batch_start:min(batch_start + batch_size - 1, num_blocks);
        count = numel(take);
        llr_standard = [zeros(5, count); double(received_llr(:, take))];
        llr_matlab = [llr_standard(parity_length+1:end, :); ...
            llr_standard(1:parity_length, :)];
        args = {'OutputFormat','whole','DecisionType','hard', ...
            'Termination','early'};
        if algorithm == "norm-min-sum"
            args = [args, {'MinSumScalingFactor', parameter}]; %#ok<AGROW>
        elseif algorithm == "offset-min-sum"
            args = [args, {'MinSumOffset', parameter}]; %#ok<AGROW>
        end
        [~, local_iterations, final_checks] = ldpcDecode(llr_matlab, ...
            decoder_cfg, max_iterations, args{:});
        parity_ok(take) = ~any(final_checks, 1).';
        iterations(take) = local_iterations(:);
    end
end

function value = parameter_text(parameter)
    if isnan(parameter)
        value = '-';
    else
        value = sprintf('%.2f', parameter);
    end
end

function make_soft_decoding_figure(algorithm_table, scale_table, output_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 820]);
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    bar(100*algorithm_table.pass_rate); grid on;
    labels = cellstr(algorithm_table.algorithm);
    for idx = 1:height(algorithm_table)
        if isfinite(algorithm_table.parameter(idx))
            labels{idx} = sprintf('%s %.2g', labels{idx}, algorithm_table.parameter(idx));
        end
    end
    xticks(1:height(algorithm_table)); xticklabels(labels); xtickangle(25);
    ylabel('通过率 (%)'); title('LDPC译码算法搜索');
    nexttile;
    models = unique(scale_table.noise_model, 'stable');
    for idx = 1:numel(models)
        mask = scale_table.noise_model == models(idx);
        semilogx(scale_table.llr_scale(mask), 100*scale_table.pass_rate(mask), ...
            '-o', 'LineWidth', 1.2, 'DisplayName', char(models(idx))); hold on;
    end
    hold off; grid on; xlabel('LLR全局尺度'); ylabel('通过率 (%)');
    title('噪声可靠度传播与LLR尺度搜索'); legend('Location', 'best');
    exportgraphics(fig, fullfile(output_dir, 'ldpc_soft_optimization.png'), ...
        'Resolution', 180);
    close(fig);
end
