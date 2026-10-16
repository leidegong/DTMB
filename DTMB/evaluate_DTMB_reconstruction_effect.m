function report = evaluate_DTMB_reconstruction_effect(sync_packet, reconstruction_packet, output_dir)
% evaluate_DTMB_reconstruction_effect 评估DTMB重构前后的星座与相关性能。

    if nargin < 1 || isempty(sync_packet)
        sync_packet = 'DTMB_1100Frame_Sync_Equalized_Packet.mat';
    end
    if nargin < 2 || isempty(reconstruction_packet)
        reconstruction_packet = 'DTMB_Full_Demod_Remod_Reconstruction.mat';
    end
    if nargin < 3 || isempty(output_dir)
        output_dir = 'DTMB_effect_test';
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('[效果测试] 读取同步包和重构包。\n');
    S = load(sync_packet, 'PN_svd_denoised_packet', 'Frame_equalization_packet', ...
        'PN_phase_index_packet', 'sync_diagnostics');
    R = load(reconstruction_packet, 'result');
    cfg = S.sync_diagnostics.cfg;
    cfg.plot_diagnostics = false;
    cfg.save_diagnostic_figures = false;
    cfg.tc_num_frames = [];
    cfg.tc_channel_estimation_frames = min(50, size(S.PN_svd_denoised_packet, 1));
    cfg.tc_use_average_channel = true;

    fprintf('[效果测试] 重建TC-AOMP公共信道均衡输出。\n');
    tc_log = evalc('tc = run_tc_aomp_section(S.PN_svd_denoised_packet, S.Frame_equalization_packet, S.PN_phase_index_packet, cfg);'); %#ok<NASGU>
    num_frames = min([size(S.Frame_equalization_packet, 1), ...
        numel(tc.equalized_reference_frames), size(R.result.reconstructed_frames, 1)]);
    pre_frames = double(S.Frame_equalization_packet(1:num_frames, :));
    tc_frames = complex(zeros(num_frames, 4725));
    for frame_idx = 1:num_frames
        tc_frames(frame_idx, :) = tc.equalized_reference_frames{frame_idx}(:).';
    end
    reconstructed_frames = double(R.result.reconstructed_frames(1:num_frames, :));

    candidate = DTMB_System_Info_Candidates();
    system_symbols = candidate(R.result.system_info.index).symbols;
    fprintf('[效果测试] 计算逐帧EVM、拟合误差和相关系数。\n');
    [evm_pre, constellation_pre] = frame_evm(pre_frames, system_symbols, round(num_frames/2));
    [evm_tc, constellation_tc] = frame_evm(tc_frames, system_symbols, round(num_frames/2));
    [evm_reconstructed, constellation_reconstructed] = ...
        frame_evm(reconstructed_frames, system_symbols, round(num_frames/2));

    nmse = sum(abs(tc_frames - reconstructed_frames).^2, 2) ./ ...
        (sum(abs(tc_frames).^2, 2) + eps);
    correlation = abs(sum(conj(tc_frames).*reconstructed_frames, 2)) ./ ...
        (sqrt(sum(abs(tc_frames).^2, 2).*sum(abs(reconstructed_frames).^2, 2)) + eps);
    canonical = double(R.result.canonical_reconstructed_frames(1:num_frames, :));
    % PSLR is dominated by the PN945 cyclic-extension sidelobe at lag +-511,
    % so it is a monotone function of the PN:body power ratio.  Measure that
    % ratio on every compared signal; comparing PSLR across signals whose
    % ratios differ is meaningless.  See sweep_pn_body_power_ratio.
    pn_body_power_ratio = pn_body_ratio(canonical);
    pn_body_power_ratio_pre = pn_body_ratio(pre_frames);
    pn_body_power_ratio_tc = pn_body_ratio(tc_frames);

    fprintf('[效果测试] 计算50帧自相关旁瓣指标。\n');
    autocorr_frames = unique(round(linspace(1, num_frames, min(50, num_frames))));
    [islr_pre, pslr_pre] = autocorrelation_metrics(pre_frames(autocorr_frames, :));
    [islr_tc, pslr_tc] = autocorrelation_metrics(tc_frames(autocorr_frames, :));
    [islr_reconstructed, pslr_reconstructed] = ...
        autocorrelation_metrics(reconstructed_frames(autocorr_frames, :));

    representative_frame = round(num_frames/2);
    make_constellation_figure(constellation_pre, constellation_tc, ...
        constellation_reconstructed, representative_frame, output_dir);
    make_autocorrelation_figure(pre_frames(representative_frame, :), ...
        tc_frames(representative_frame, :), reconstructed_frames(representative_frame, :), ...
        representative_frame, output_dir);
    make_metric_figure(evm_pre, evm_tc, evm_reconstructed, nmse, correlation, output_dir);

    report = struct();
    report.num_frames = num_frames;
    report.representative_frame = representative_frame;
    report.evm_pre_median = median(evm_pre, 'omitnan');
    report.evm_tc_median = median(evm_tc, 'omitnan');
    report.evm_reconstructed_median = median(evm_reconstructed, 'omitnan');
    report.evm_improvement_tc_to_reconstructed_db = 20*log10( ...
        report.evm_tc_median/(report.evm_reconstructed_median + eps));
    report.nmse_db_median = median(10*log10(nmse + eps), 'omitnan');
    report.correlation_median = median(correlation, 'omitnan');
    report.pn_body_power_ratio_median = median(pn_body_power_ratio, 'omitnan');
    report.pn_body_power_ratio_pre_median = median(pn_body_power_ratio_pre, 'omitnan');
    report.pn_body_power_ratio_tc_median = median(pn_body_power_ratio_tc, 'omitnan');
    report.pslr_comparison_ratio_matched = ...
        max([report.pn_body_power_ratio_pre_median, ...
        report.pn_body_power_ratio_tc_median, ...
        report.pn_body_power_ratio_median]) / ...
        min([report.pn_body_power_ratio_pre_median, ...
        report.pn_body_power_ratio_tc_median, ...
        report.pn_body_power_ratio_median]) < 1.25;
    report.islr_pre_median_db = median(islr_pre, 'omitnan');
    report.islr_tc_median_db = median(islr_tc, 'omitnan');
    report.islr_reconstructed_median_db = median(islr_reconstructed, 'omitnan');
    report.pslr_pre_median_db = median(pslr_pre, 'omitnan');
    report.pslr_tc_median_db = median(pslr_tc, 'omitnan');
    report.pslr_reconstructed_median_db = median(pslr_reconstructed, 'omitnan');
    report.fec_successful_blocks = nnz(R.result.fec_stats.parity_ok);
    report.fec_total_blocks = R.result.fec_stats.num_blocks;
    report.fec_valid_frames = nnz(R.result.fec_valid_frame_mask);
    report.reconstruction_source_counts = [nnz(R.result.reconstruction_source == 1), ...
        nnz(R.result.reconstruction_source == 2)];
    report.evm_pre_per_frame = evm_pre;
    report.evm_tc_per_frame = evm_tc;
    report.evm_reconstructed_per_frame = evm_reconstructed;
    report.nmse_per_frame = nmse;
    report.correlation_per_frame = correlation;
    report.islr_pre_db = islr_pre;
    report.islr_tc_db = islr_tc;
    report.islr_reconstructed_db = islr_reconstructed;
    report.pslr_pre_db = pslr_pre;
    report.pslr_tc_db = pslr_tc;
    report.pslr_reconstructed_db = pslr_reconstructed;

    metric = ["EVM_pre_median"; "EVM_TC_median"; "EVM_reconstructed_median"; ...
        "EVM_improvement_TC_to_reconstructed"; "NMSE_median"; "correlation_median"; ...
        "PN_to_body_power_ratio_pre"; "PN_to_body_power_ratio_TC"; ...
        "PN_to_body_power_ratio_reconstructed"; "PSLR_comparison_ratio_matched"; ...
        "ISLR_pre_median"; "ISLR_TC_median"; ...
        "ISLR_reconstructed_median"; "PSLR_pre_median"; "PSLR_TC_median"; ...
        "PSLR_reconstructed_median"; "FEC_successful_blocks"; "FEC_valid_frames"];
    value = [report.evm_pre_median; report.evm_tc_median; report.evm_reconstructed_median; ...
        report.evm_improvement_tc_to_reconstructed_db; report.nmse_db_median; ...
        report.correlation_median; report.pn_body_power_ratio_pre_median; ...
        report.pn_body_power_ratio_tc_median; report.pn_body_power_ratio_median; ...
        double(report.pslr_comparison_ratio_matched); ...
        report.islr_pre_median_db; report.islr_tc_median_db; ...
        report.islr_reconstructed_median_db; report.pslr_pre_median_db; ...
        report.pslr_tc_median_db; report.pslr_reconstructed_median_db; ...
        report.fec_successful_blocks; report.fec_valid_frames];
    unit = ["linear"; "linear"; "linear"; "dB"; "dB"; "linear"; "linear"; ...
        "linear"; "linear"; "flag"; "dB"; "dB"; "dB"; "dB"; "dB"; "dB"; ...
        "blocks"; "frames"];
    summary_table = table(metric, value, unit);
    report.summary_table = summary_table;
    writetable(summary_table, fullfile(output_dir, 'effect_summary.csv'));
    save(fullfile(output_dir, 'DTMB_Effect_Test_Report.mat'), 'report', '-v7.3');
    disp(summary_table);
    fprintf('[效果测试] 结果已保存到：%s\n', output_dir);
end

function [evm, selected_constellation] = frame_evm(frames, system_symbols, selected_frame)
    num_frames = size(frames, 1);
    evm = nan(num_frames, 1);
    selected_constellation = [];
    for frame_idx = 1:num_frames
        body_fd = fftshift(fft(frames(frame_idx, 946:4725).', 3780));
        [data, system_rx] = DTMB_Frequency_Deinterleave(body_fd);
        phase0 = angle((system_symbols' * system_rx)/(system_symbols' * system_symbols + eps));
        phase_corrected = data*exp(-1j*phase0);
        scale = rms(phase_corrected) + eps;
        normalized = phase_corrected/scale;
        residual_phase = Estimate_CPE_64QAM(normalized);
        normalized = data*exp(-1j*(phase0 + residual_phase))/(scale + eps);
        sliced = Slice_64QAM(normalized);
        evm(frame_idx) = sqrt(mean(abs(normalized - sliced).^2, 'omitnan'));
        if frame_idx == selected_frame
            selected_constellation = normalized;
        end
    end
end

function ratio = pn_body_ratio(frames)
    ratio = mean(abs(frames(:, 1:945)).^2, 2) ./ ...
        (mean(abs(frames(:, 946:end)).^2, 2) + eps);
end

function [islr_db, pslr_db] = autocorrelation_metrics(frames)
    num_frames = size(frames, 1);
    islr_db = nan(num_frames, 1);
    pslr_db = nan(num_frames, 1);
    main_half_width = 3;
    roi_half_width = 511;
    for frame_idx = 1:num_frames
        power_corr = abs(xcorr(frames(frame_idx, :).')).^2;
        [peak_power, peak_idx] = max(power_corr);
        lag = (1:numel(power_corr)).' - peak_idx;
        roi = abs(lag) <= roi_half_width;
        main = abs(lag) <= main_half_width;
        side = roi & ~main;
        islr_db(frame_idx) = 10*log10((mean(power_corr(side)) + eps) / ...
            (sum(power_corr(main)) + eps));
        pslr_db(frame_idx) = 10*log10((max(power_corr(side)) + eps) / (peak_power + eps));
    end
end

function make_constellation_figure(pre, tc, reconstructed, frame_idx, output_dir)
    max_points = 3744;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 480]);
    tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_constellation_panel(pre(1:min(end,max_points)), 'TC-AOMP前');
    plot_constellation_panel(tc(1:min(end,max_points)), 'TC-AOMP后');
    plot_constellation_panel(reconstructed(1:min(end,max_points)), '判决重构后');
    sgtitle(sprintf('DTMB 第%d帧 64QAM星座对比', frame_idx));
    exportgraphics(fig, fullfile(output_dir, 'constellation_comparison.png'), 'Resolution', 180);
    close(fig);
end

function plot_constellation_panel(symbols, title_text)
    nexttile;
    scatter(real(symbols), imag(symbols), 5, '.', 'MarkerEdgeAlpha', 0.35);
    axis equal; grid on; xlim([-1.4 1.4]); ylim([-1.4 1.4]);
    xlabel('I'); ylabel('Q'); title(title_text);
end

function make_autocorrelation_figure(pre, tc, reconstructed, frame_idx, output_dir)
    [lag, pre_db] = normalized_autocorrelation(pre);
    [~, tc_db] = normalized_autocorrelation(tc);
    [~, recon_db] = normalized_autocorrelation(reconstructed);
    show = abs(lag) <= 600;
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 650]);
    plot(lag(show), pre_db(show), 'Color', [0.55 0.55 0.55], 'LineWidth', 0.9); hold on;
    plot(lag(show), tc_db(show), 'b-', 'LineWidth', 1.0);
    plot(lag(show), recon_db(show), 'r-', 'LineWidth', 1.0); hold off;
    grid on; ylim([-90 3]); xlim([-600 600]);
    xlabel('延迟 / 符号'); ylabel('归一化自相关功率 / dB');
    title(sprintf('DTMB 第%d帧自相关对比', frame_idx));
    legend('TC-AOMP前', 'TC-AOMP后', '重构后', 'Location', 'best');
    exportgraphics(fig, fullfile(output_dir, 'autocorrelation_comparison.png'), 'Resolution', 180);
    close(fig);
end

function [lag, corr_db] = normalized_autocorrelation(frame)
    power_corr = abs(xcorr(frame(:))).^2;
    [peak, peak_idx] = max(power_corr);
    lag = (1:numel(power_corr)).' - peak_idx;
    corr_db = 10*log10(power_corr/(peak + eps) + eps);
end

function make_metric_figure(evm_pre, evm_tc, evm_recon, nmse, correlation, output_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 850]);
    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    plot(evm_pre, 'Color', [0.55 0.55 0.55]); hold on;
    plot(evm_tc, 'b'); plot(evm_recon, 'r'); hold off; grid on;
    ylabel('RMS EVM'); legend('TC-AOMP前', 'TC-AOMP后', '重构后');
    title('逐帧64QAM最近点EVM');
    nexttile; plot(10*log10(nmse + eps), 'k'); grid on;
    ylabel('NMSE / dB'); title('重构信号相对TC-AOMP均衡信号的拟合误差');
    nexttile; plot(correlation, 'Color', [0.1 0.55 0.2]); grid on;
    xlabel('帧序号'); ylabel('相关系数'); ylim([0 1.02]); title('逐帧归一化相关系数');
    exportgraphics(fig, fullfile(output_dir, 'per_frame_metrics.png'), 'Resolution', 180);
    close(fig);
end
