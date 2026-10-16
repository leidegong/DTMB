function info = Joint_Time_Freq_Sync_2D_PN(rx_signal, frame_start_seed, phase_index_seed, detected_mode, cfg)
% Joint_Time_Freq_Sync_2D_PN
% 按论文 2.3.4 的思想，在“帧内时延 d × 帧间频偏 v”二维平面上进行 PN 相参匹配。

    if ~strcmpi(detected_mode, 'PN945')
        error('当前二维联合时频同步实现仅支持 PN945，检测到模式：%s', detected_mode);
    end

    rx_signal = rx_signal(:);
    Lpn_sym = cfg.pn_length;
    Lframe_sym = cfg.frame_symbol_len;
    sps = cfg.sps;
    Fs = cfg.Fs;
    Lpn_samp = Lpn_sym * sps;
    Lframe_samp = Lframe_sym * sps;
    total_len = numel(rx_signal);

    max_frames_by_len = floor((total_len - round(frame_start_seed) - Lpn_samp + 1) / Lframe_samp) + 1;
    num_frames = min([cfg.joint_2d_num_frames, cfg.target_frames, max_frames_by_len]);
    if num_frames < 8
        error('二维联合时频同步可用帧数过少：%d', num_frames);
    end

    search_half = cfg.joint_2d_time_search_samp;
    search_step = cfg.joint_2d_time_step_samp;
    delay_offsets = (-search_half:search_step:search_half).';
    num_delay = numel(delay_offsets);

    h_srrc = rcosdesign(0.05, 8, sps, 'sqrt');
    local_pn_cache = cell(200, 1);
    for ph = 0:199
        local_pn = PN_gen(ph).';
        local_up = zeros(numel(local_pn) * sps, 1);
        local_up(1:sps:end) = local_pn(:);
        local_pn_cache{ph + 1} = conv(local_up, h_srrc, 'same');
    end

    xi = complex(zeros(num_delay, num_frames));
    frame_start_seed = round(frame_start_seed);
    fprintf('二维联合同步搜索：%d 个时延候选 × %d 帧，时延范围 [%d, %d] 采样点。\n', ...
        num_delay, num_frames, delay_offsets(1), delay_offsets(end));

    for id = 1:num_delay
        d = delay_offsets(id);
        for m = 1:num_frames
            idx = frame_start_seed + d + (m - 1) * Lframe_samp;
            stop_idx = idx + Lpn_samp - 1;
            if idx < 1 || stop_idx > total_len
                xi(id, m) = 0;
                continue;
            end
            phase_m = mod(phase_index_seed + m - 1, 200);
            local_pn_samp = local_pn_cache{phase_m + 1};
            rx_pn = rx_signal(idx:stop_idx);
            xi(id, m) = rx_pn(:)' * conj(local_pn_samp(:));
        end
    end

    switch lower(cfg.joint_2d_window)
        case 'hann'
            win = hann(num_frames).';
        otherwise
            win = ones(1, num_frames);
    end
    win = win / (rms(win) + eps);

    nfft = 2^nextpow2(max(num_frames * cfg.joint_2d_fft_oversample, num_frames));
    spectrum = fftshift(fft(xi .* win, nfft, 2), 2);
    mag = abs(spectrum);
    [peak_value, lin_idx] = max(mag(:));
    [best_delay_idx, best_freq_idx] = ind2sub(size(mag), lin_idx);

    freq_bins = (-nfft/2:nfft/2-1) / nfft;
    cycles_per_frame = freq_bins(best_freq_idx);
    frame_period = Lframe_samp / Fs;
    cfo_hz = cycles_per_frame / frame_period;

    best_delay_offset = delay_offsets(best_delay_idx);
    best_frame_start = frame_start_seed + best_delay_offset;

    info = struct();
    info.frame_start_seed = frame_start_seed;
    info.frame_start = best_frame_start;
    info.phase_index_seed = phase_index_seed;
    info.best_delay_offset_samp = best_delay_offset;
    info.cfo_hz = cfo_hz;
    info.peak_value = peak_value;
    info.num_frames = num_frames;
    info.delay_offsets = delay_offsets;
    info.freq_axis_hz = freq_bins / frame_period;
    info.best_delay_index = best_delay_idx;
    info.best_freq_index = best_freq_idx;
    info.delay_profile = max(mag, [], 2);
    info.freq_profile = mag(best_delay_idx, :).';
    info.xi_best_delay = xi(best_delay_idx, :).';

    if isfield(cfg, 'plot_diagnostics') && cfg.plot_diagnostics
        plot_joint_2d_sync_result(info);
    end
end

function plot_joint_2d_sync_result(info)
    figure('Name', '二维联合时频同步估计结果', 'Color', 'w');
    subplot(1, 2, 1);
    plot(info.delay_offsets, info.delay_profile, 'LineWidth', 1.2);
    grid on;
    xlabel('时延候选 / 采样点');
    ylabel('帧间相参峰值');
    title(sprintf('时域同步峰值，最佳偏移 %+d 点', info.best_delay_offset_samp));

    subplot(1, 2, 2);
    plot(info.freq_axis_hz, info.freq_profile, 'LineWidth', 1.2);
    grid on;
    xlabel('载频偏移 / Hz');
    ylabel('相参 FFT 幅度');
    title(sprintf('帧间频偏峰值，CFO %.6f Hz', info.cfo_hz));
end

