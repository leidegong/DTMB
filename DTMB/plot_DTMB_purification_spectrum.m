function report = plot_DTMB_purification_spectrum(sync_packet, final_packet, output_dir)
% plot_DTMB_purification_spectrum 绘制DTMB参考信号纯化前后的频谱对比。
%
% 只使用最终fec_valid_frame_mask标记的同一批帧，以单位RMS归一化消除
% 纯增益差异；同时显示连续帧Welch PSD和3780载波平均功率谱。

    if nargin < 1 || isempty(sync_packet)
        sync_packet = 'DTMB_1100Frame_Sync_Equalized_Packet.mat';
    end
    if nargin < 2 || isempty(final_packet)
        final_packet = fullfile('DTMB_groupwise_equalization', ...
            'DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat');
    end
    if nargin < 3 || isempty(output_dir)
        output_dir = 'DTMB_spectrum_comparison';
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    fprintf('[频谱对比 1/4] 读取纯化前与最终重构数据。\n');
    S = load(sync_packet, 'PN_svd_denoised_packet', 'Frame_equalization_packet', ...
        'PN_phase_index_packet', 'sync_diagnostics');
    R = load(final_packet, 'result');
    cfg = S.sync_diagnostics.cfg;
    cfg.plot_diagnostics = false;
    cfg.save_diagnostic_figures = false;
    cfg.tc_num_frames = [];
    cfg.tc_channel_estimation_frames = min(50, size(S.PN_svd_denoised_packet, 1));
    cfg.tc_use_average_channel = true;

    fprintf('[频谱对比 2/4] 重建TC-AOMP中间结果。\n');
    tc_log = evalc(['tc = run_tc_aomp_section(S.PN_svd_denoised_packet, ' ...
        'S.Frame_equalization_packet, S.PN_phase_index_packet, cfg);']); %#ok<NASGU>
    num_frames = min([size(S.Frame_equalization_packet, 1), ...
        numel(tc.equalized_reference_frames), size(R.result.reconstructed_frames, 1)]);
    valid_mask = logical(R.result.fec_valid_frame_mask(1:num_frames));
    selected_frames = find(valid_mask);
    if numel(selected_frames) < 8
        count = min(79, num_frames);
        first = floor((num_frames-count)/2) + 1;
        selected_frames = (first:first+count-1).';
        selection_source = '中心帧（无足够FEC有效标记）';
    else
        selection_source = 'fec_valid_frame_mask';
    end

    before_frames = double(S.Frame_equalization_packet(selected_frames, 1:4725));
    tc_frames = complex(zeros(numel(selected_frames), 4725));
    for idx = 1:numel(selected_frames)
        frame = tc.equalized_reference_frames{selected_frames(idx)};
        tc_frames(idx, :) = double(frame(1:4725)).';
    end
    after_frames = double(R.result.reconstructed_frames(selected_frames, 1:4725));

    fprintf('[频谱对比 3/4] 计算Welch PSD与帧体载波平均谱。\n');
    Fs = 7.56e6;
    nfft = 16384;
    window_length = 8192;
    overlap = window_length/2;
    [psd_before, frequency_hz] = centered_welch(before_frames, Fs, ...
        window_length, overlap, nfft);
    psd_tc = centered_welch(tc_frames, Fs, window_length, overlap, nfft);
    psd_after = centered_welch(after_frames, Fs, window_length, overlap, nfft);
    common_reference = max(psd_before);
    psd_before_db = 10*log10(psd_before/common_reference + eps);
    psd_tc_db = 10*log10(psd_tc/common_reference + eps);
    psd_after_db = 10*log10(psd_after/common_reference + eps);

    [carrier_before, carrier_frequency_hz, carrier_map_before] = ...
        average_body_carrier_spectrum(before_frames, Fs);
    [carrier_tc, ~, ~] = average_body_carrier_spectrum(tc_frames, Fs);
    [carrier_after, ~, carrier_map_after] = average_body_carrier_spectrum(after_frames, Fs);
    carrier_before_db = 10*log10(carrier_before + eps);
    carrier_tc_db = 10*log10(carrier_tc + eps);
    carrier_after_db = 10*log10(carrier_after + eps);

    edge_guard = round(0.02*numel(carrier_before)) + 1;
    interior = edge_guard:(numel(carrier_before)-edge_guard+1);
    flatness_before = spectral_flatness(carrier_before(interior));
    flatness_tc = spectral_flatness(carrier_tc(interior));
    flatness_after = spectral_flatness(carrier_after(interior));
    ripple_before_db = std(carrier_before_db(interior));
    ripple_tc_db = std(carrier_tc_db(interior));
    ripple_after_db = std(carrier_after_db(interior));

    fprintf('[频谱对比 4/4] 保存频谱图和数值表。\n');
    spectrum_table = table(frequency_hz/1e6, psd_before_db, psd_tc_db, ...
        psd_after_db, 'VariableNames', {'frequency_MHz','before_dB', ...
        'tc_aomp_dB','purified_dB'});
    carrier_table = table(carrier_frequency_hz/1e6, carrier_before_db, ...
        carrier_tc_db, carrier_after_db, 'VariableNames', ...
        {'frequency_MHz','before_dB','tc_aomp_dB','purified_dB'});
    writetable(spectrum_table, fullfile(output_dir, 'welch_spectrum.csv'));
    writetable(carrier_table, fullfile(output_dir, 'body_carrier_spectrum.csv'));

    report = struct();
    report.num_selected_frames = numel(selected_frames);
    report.selected_frames = selected_frames;
    report.selection_source = selection_source;
    report.sample_rate_hz = Fs;
    report.psd_normalization = '各时域序列先单位RMS；三条PSD共同以纯化前峰值为0 dB';
    report.carrier_normalization = '各阶段平均载波功率分别归一化为单位总功率';
    report.spectral_flatness = [flatness_before, flatness_tc, flatness_after];
    report.carrier_ripple_std_db = [ripple_before_db, ripple_tc_db, ripple_after_db];
    report.labels = {'纯化前同步帧','TC-AOMP','最终FEC纯化'};
    save(fullfile(output_dir, 'DTMB_Spectrum_Comparison_Report.mat'), ...
        'report', '-v7.3');
    make_spectrum_figure(frequency_hz, psd_before_db, psd_tc_db, ...
        psd_after_db, carrier_frequency_hz, carrier_before_db, carrier_tc_db, ...
        carrier_after_db, carrier_map_before, carrier_map_after, ...
        selected_frames, output_dir);

    fprintf('  选用帧：%d帧（%s）。\n', numel(selected_frames), selection_source);
    fprintf('  载波谱平坦度：%.4f -> %.4f -> %.4f。\n', ...
        flatness_before, flatness_tc, flatness_after);
    fprintf('  载波功率起伏标准差：%.3f -> %.3f -> %.3f dB。\n', ...
        ripple_before_db, ripple_tc_db, ripple_after_db);
    fprintf('  输出目录：%s\n', output_dir);
end

function [psd_value, frequency_hz] = centered_welch(frames, Fs, ...
        window_length, overlap, nfft)
    stream = reshape(frames.', [], 1);
    stream = stream - mean(stream);
    stream = stream/(rms(stream) + eps);
    window = hann(window_length, 'periodic');
    psd_value = pwelch(stream, window, overlap, nfft, Fs, 'twosided');
    psd_value = fftshift(real(psd_value));
    frequency_hz = ((-nfft/2):(nfft/2-1)).' * (Fs/nfft);
end

function [power_spectrum, frequency_hz, normalized_map_db] = ...
        average_body_carrier_spectrum(frames, Fs)
    body = frames(:, 946:4725);
    body = body/(rms(body(:)) + eps);
    carrier = fftshift(fft(body, 3780, 2), 2);
    raw_power = abs(carrier).^2;
    power_spectrum = mean(raw_power, 1).';
    power_spectrum = power_spectrum/(sum(power_spectrum) + eps);
    per_frame_reference = median(raw_power, 2) + eps;
    normalized_map_db = 10*log10(raw_power ./ per_frame_reference);
    frequency_hz = ((-1890):1889).' * (Fs/3780);
end

function value = spectral_flatness(power_spectrum)
    power_spectrum = max(real(power_spectrum(:)), eps);
    value = exp(mean(log(power_spectrum)))/mean(power_spectrum);
end

function make_spectrum_figure(frequency_hz, before_psd, tc_psd, after_psd, ...
        carrier_frequency_hz, before_carrier, tc_carrier, after_carrier, ...
        before_map, after_map, selected_frames, output_dir)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 950]);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    plot(frequency_hz/1e6, before_psd, 'Color', [0.55 0.55 0.55], ...
        'LineWidth', 0.9); hold on;
    plot(frequency_hz/1e6, tc_psd, 'b', 'LineWidth', 0.9);
    plot(frequency_hz/1e6, after_psd, 'r', 'LineWidth', 1.0); hold off;
    grid on; xlabel('基带频率 (MHz)'); ylabel('PSD (dB)');
    title('连续帧Welch功率谱（单位RMS、共同参考）');
    legend({'纯化前同步帧','TC-AOMP','最终FEC纯化'}, 'Location', 'best');
    xlim([frequency_hz(1), frequency_hz(end)]/1e6);

    nexttile;
    plot(carrier_frequency_hz/1e6, before_carrier, ...
        'Color', [0.55 0.55 0.55], 'LineWidth', 0.8); hold on;
    plot(carrier_frequency_hz/1e6, tc_carrier, 'b', 'LineWidth', 0.8);
    plot(carrier_frequency_hz/1e6, after_carrier, 'r', 'LineWidth', 0.9); hold off;
    grid on; xlabel('子载波频率 (MHz)'); ylabel('平均功率 (dB)');
    title('3780点帧体平均载波谱（各自单位总功率）');
    legend({'纯化前同步帧','TC-AOMP','最终FEC纯化'}, 'Location', 'best');

    color_limits = percentile_limits([before_map(:); after_map(:)], [2 98]);
    nexttile;
    imagesc(carrier_frequency_hz/1e6, selected_frames, before_map);
    axis xy; caxis(color_limits); colorbar; colormap(turbo);
    xlabel('子载波频率 (MHz)'); ylabel('原始帧号');
    title('纯化前逐帧载波功率（相对帧内中位数，dB）');

    nexttile;
    imagesc(carrier_frequency_hz/1e6, selected_frames, after_map);
    axis xy; caxis(color_limits); colorbar; colormap(turbo);
    xlabel('子载波频率 (MHz)'); ylabel('原始帧号');
    title('最终FEC纯化后逐帧载波功率（dB）');

    exportgraphics(fig, fullfile(output_dir, ...
        'purification_spectrum_comparison.png'), 'Resolution', 180);
    close(fig);
end

function limits = percentile_limits(values, percent)
    values = sort(values(isfinite(values)));
    if isempty(values)
        limits = [-1 1];
        return;
    end
    index = 1 + round((numel(values)-1)*percent(:)/100);
    index = min(max(index, 1), numel(values));
    limits = values(index).';
    if limits(1) == limits(2)
        limits = limits + [-1 1];
    end
end
