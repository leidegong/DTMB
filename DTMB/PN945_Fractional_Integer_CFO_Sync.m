function info = PN945_Fractional_Integer_CFO_Sync(rx_signal, frame_start_seed, phase_index_seed, detected_mode, cfg)
% PN945_Fractional_Integer_CFO_Sync
% 使用 PN945 内部重复结构估计小数倍频偏，并在小数补偿后用 PN945
% 频域差分匹配估计整数倍频偏。该模块替代原二维联合 delay+CFO 估计；
% 不修改 Symbol_Sync 给出的初始帧起点和相位索引。

    if ~strcmpi(detected_mode, 'PN945')
        error('PN945 新频偏估计方法仅支持 PN945，检测到模式：%s', detected_mode);
    end

    rx_signal = rx_signal(:);
    frame_start = round(frame_start_seed);
    sps = cfg.sps;
    Fs = cfg.Fs;
    Lframe_samp = cfg.frame_symbol_len * sps;

    frac = estimate_fractional_cfo_pn945(rx_signal, frame_start, cfg, Lframe_samp);
    rx_frac = compensate_cfo_hz(rx_signal, Fs, frac.fractional_cfo_hz);

    intcfo = estimate_integer_cfo_pn945(rx_frac, frame_start, phase_index_seed, cfg, Lframe_samp);
    rx_comp = compensate_integer_cfo_bin(rx_frac, intcfo.integer_bin, intcfo.fft_len);

    info = struct();
    info.method = 'pn945_repetition_fractional_plus_frequency_differential_integer_cfo';
    info.frame_start_seed = frame_start;
    info.frame_start = frame_start;
    info.phase_index_seed = phase_index_seed;
    info.rx_comp = rx_comp;
    info.rx_frac = rx_frac;
    info.fractional = frac;
    info.integer = intcfo;
    info.f_frac = frac.fractional_cfo_hz;
    info.m_hat = intcfo.integer_bin;
    info.f_int = intcfo.integer_cfo_hz;
    info.total_initial_cfo = info.f_frac + info.f_int;
    info.best_delay_offset_samp = 0;
    info.cfo_hz = info.total_initial_cfo;
    info.peak_value = frac.correlation_abs;

    fprintf('PN945 CFO：小数 %.9f Hz（有效 %d/%d 帧），整数 bin %+d = %.9f Hz（有效 %d/%d 帧），合计 %.9f Hz。\n', ...
        info.f_frac, frac.num_valid, frac.num_frames, ...
        info.m_hat, info.f_int, intcfo.num_valid, intcfo.num_frames, info.total_initial_cfo);

    if isfield(cfg, 'plot_diagnostics') && cfg.plot_diagnostics
        plot_pn945_cfo_result(info);
    end
end

function info = estimate_fractional_cfo_pn945(rx_signal, frame_start, cfg, Lframe_samp)
    repeat_len = 434 * cfg.sps;
    repeat_delay = 511 * cfg.sps;
    max_frames_by_len = floor((numel(rx_signal) - frame_start - repeat_delay - repeat_len + 1) / Lframe_samp) + 1;
    num_frames = min([get_cfg_scalar(cfg, 'fractional_cfo_num_frames', cfg.target_frames), cfg.target_frames, max_frames_by_len]);
    if num_frames <= 0
        error('PN945 小数 CFO 估计可用样本不足。');
    end

    frame_corr = complex(zeros(num_frames, 1));
    valid = false(num_frames, 1);
    for m = 1:num_frames
        idx1 = frame_start + (m - 1) * Lframe_samp;
        idx2 = idx1 + repeat_delay;
        stop2 = idx2 + repeat_len - 1;
        if idx1 < 1 || stop2 > numel(rx_signal)
            continue;
        end
        first = rx_signal(idx1:idx1 + repeat_len - 1);
        second = rx_signal(idx2:idx2 + repeat_len - 1);
        frame_corr(m) = sum(conj(first) .* second);
        valid(m) = isfinite(real(frame_corr(m))) && isfinite(imag(frame_corr(m)));
    end

    if ~any(valid)
        error('PN945 小数 CFO 没有有效重复段。');
    end

    P = sum(frame_corr(valid));
    phase = angle(P);
    info = struct();
    info.method = 'pn945_repeated_ca_fractional_cfo';
    info.frame_start_seed = frame_start;
    info.frame_start = frame_start;
    info.phase_index_seed = [];
    info.fractional_cfo_hz = phase * cfg.Fs / (2 * pi * repeat_delay);
    info.repeat_delay_samples = repeat_delay;
    info.repeat_len_samples = repeat_len;
    info.num_frames = num_frames;
    info.num_valid = sum(valid);
    info.correlation = P;
    info.correlation_abs = abs(P);
    info.correlation_phase_rad = phase;
    info.frame_correlation = frame_corr;
    info.valid = valid;
end

function info = estimate_integer_cfo_pn945(rx_frac, frame_start, phase_index_seed, cfg, Lframe_samp)
    fft_len = cfg.body_length * cfg.sps;
    pn_len = cfg.pn_length * cfg.sps;
    max_frames_by_len = floor((numel(rx_frac) - frame_start - pn_len + 1) / Lframe_samp) + 1;
    num_frames = min([get_cfg_scalar(cfg, 'integer_cfo_num_frames', cfg.target_frames), cfg.target_frames, max_frames_by_len]);
    if num_frames <= 0
        error('PN945 整数 CFO 估计可用样本不足。');
    end

    min_bin = get_cfg_scalar(cfg, 'integer_cfo_search_min_bin', -20);
    max_bin = get_cfg_scalar(cfg, 'integer_cfo_search_max_bin', 20);
    search_bins = (min_bin:max_bin).';
    metric = zeros(numel(search_bins), 1);
    valid = false(num_frames, 1);

    h_srrc = rcosdesign(0.05, 8, cfg.sps, 'sqrt');
    local_pn_cache = cell(200, 1);
    for ph = 0:199
        local_pn = PN_gen(ph).';
        local_pn_up = zeros(length(local_pn) * cfg.sps, 1);
        local_pn_up(1:cfg.sps:end) = local_pn(:);
        local_pn_cache{ph + 1} = conv(local_pn_up, h_srrc, 'same');
    end

    for m = 1:num_frames
        st = frame_start + (m - 1) * Lframe_samp;
        stop_idx = st + pn_len - 1;
        if st < 1 || stop_idx > numel(rx_frac)
            continue;
        end

        phase_m = mod(phase_index_seed + m - 1, 200);
        r_pn = rms_normalize(rx_frac(st:stop_idx));
        s_pn = rms_normalize(local_pn_cache{phase_m + 1});

        X = fft(r_pn, fft_len);
        P = fft(s_pn, fft_len);
        Y = X .* conj(circshift(X, -1));
        Y_ref = P .* conj(circshift(P, -1));
        if ~(all(isfinite(real(Y))) && all(isfinite(imag(Y))) && ...
                all(isfinite(real(Y_ref))) && all(isfinite(imag(Y_ref))))
            continue;
        end

        for k = 1:numel(search_bins)
            metric(k) = metric(k) + abs(sum(Y .* conj(circshift(Y_ref, search_bins(k)))));
        end
        valid(m) = true;
    end

    if ~any(valid)
        error('PN945 整数 CFO 没有有效帧。');
    end

    [best_metric, best_idx] = max(metric);
    integer_bin = search_bins(best_idx);
    info = struct();
    info.method = 'pn945_frequency_differential_integer_cfo';
    info.frame_start_seed = frame_start;
    info.frame_start = frame_start;
    info.phase_index_seed = phase_index_seed;
    info.integer_bin = integer_bin;
    info.integer_cfo_hz = integer_bin * cfg.Fs / fft_len;
    info.fft_len = fft_len;
    info.pn_len_samples = pn_len;
    info.num_frames = num_frames;
    info.num_valid = sum(valid);
    info.search_bins = search_bins;
    info.metric = metric;
    info.best_metric = best_metric;
    info.valid = valid;
end

function y = compensate_cfo_hz(x, Fs, cfo_hz)
    n = (0:numel(x)-1).';
    y = x(:) .* exp(-1j * 2 * pi * cfo_hz / Fs * n);
end

function y = compensate_integer_cfo_bin(x, integer_bin, fft_len)
    n = (0:numel(x)-1).';
    y = x(:) .* exp(-1j * 2 * pi * integer_bin / fft_len * n);
end

function y = rms_normalize(x)
    x = x(:);
    y = x ./ (rms(x) + eps);
end

function value = get_cfg_scalar(cfg, field_name, default_value)
    if isfield(cfg, field_name) && ~isempty(cfg.(field_name))
        value = cfg.(field_name);
    else
        value = default_value;
    end
end

function plot_pn945_cfo_result(info)
    figure('Name', 'PN945 小数+整数 CFO 估计结果', 'Color', 'w');
    subplot(1, 2, 1);
    stem(find(info.fractional.valid), abs(info.fractional.frame_correlation(info.fractional.valid)), 'filled');
    grid on;
    xlabel('帧序号');
    ylabel('|P_{frame}|');
    title(sprintf('小数 CFO %.6f Hz', info.f_frac));

    subplot(1, 2, 2);
    plot(info.integer.search_bins, info.integer.metric, 'o-', 'LineWidth', 1.1);
    grid on;
    xlabel('整数频偏 bin');
    ylabel('多帧累加 metric');
    title(sprintf('整数 bin %+d，%.3f Hz', info.m_hat, info.f_int));
end
