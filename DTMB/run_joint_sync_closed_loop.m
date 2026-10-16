function [rx_matched, rx_matched_uncomp, rx_16x, exact_frame_start, phase_index, timing_info, ...
          current_pass, loop_info] = run_joint_sync_closed_loop( ...
          rx_matched, rx_matched_uncomp, rx_16x, exact_frame_start, phase_index, timing_info, ...
          coarse_frame_starts, detected_mode, num_frames_est, cfg, current_pass)

    loop_info = init_cfo_loop_info();
    loop_info.enabled = true;
    loop_info.stop_reason = 'max_iter';

    for iter = 1:cfg.sync_loop_max_iter
        fprintf('\n  闭环迭代 %d / %d\n', iter, cfg.sync_loop_max_iter);

        cfo_before = estimate_residual_cfo_from_pn(current_pass.PN, current_pass.phase_index, ...
            current_pass.frame_heads, current_pass.corr_peak, cfg.Fs, cfg);

        fprintf('    补偿前残留 CFO：%.9f Hz（%d 个有效帧）\n', ...
            cfo_before.res_cfo_hz, cfo_before.num_valid);

        entry = struct();
        entry.iteration = iter;
        entry.before_cfo_hz = cfo_before.res_cfo_hz;
        entry.after_cfo_hz = nan;
        entry.applied_cfo_hz = 0;
        entry.total_applied_cfo_hz = loop_info.total_applied_cfo_hz;
        entry.before_valid_frames = cfo_before.num_valid;
        entry.after_valid_frames = 0;
        entry.accepted = false;

        if abs(cfo_before.res_cfo_hz) <= cfg.sync_loop_cfo_tol_hz
            loop_info.stop_reason = 'cfo_tol_before_compensation';
            loop_info.final_residual_cfo = cfo_before;
            loop_info.history(end+1) = entry;
            break;
        end

        applied_cfo_hz = cfg.sync_loop_damping * cfo_before.res_cfo_hz;
        rx_matched_try = apply_cfo_compensation(rx_matched, cfg.Fs, applied_cfo_hz);
        rx_matched_uncomp_try = apply_cfo_compensation(rx_matched_uncomp, cfg.Fs, applied_cfo_hz);

        [phase_index_try, exact_frame_start_try] = ...
            Symbol_Sync(rx_matched_try, coarse_frame_starts, cfg.sps, detected_mode);

        [rx_matched_try, ~, rx_16x_try, ~, ~, ~, exact_frame_start_try, timing_info_try] = ...
            Timing_Phase_Refine_8x(rx_matched_try, rx_matched_uncomp_try, exact_frame_start_try, ...
                                   cfg.sps, detected_mode, phase_index_try);

        pass_try = refine_and_build_matrices(rx_matched_try, rx_16x_try, exact_frame_start_try, ...
            phase_index_try, cfg.sps, detected_mode, num_frames_est, timing_info_try, ...
            cfg.debug_plot_each_frame);

        cfo_after = estimate_residual_cfo_from_pn(pass_try.PN, pass_try.phase_index, ...
            pass_try.frame_heads, pass_try.corr_peak, cfg.Fs, cfg);

        fprintf('    补偿后残留 CFO：%.9f Hz（%d 个有效帧）\n', ...
            cfo_after.res_cfo_hz, cfo_after.num_valid);

        before_abs = abs(cfo_before.res_cfo_hz);
        after_abs = abs(cfo_after.res_cfo_hz);
        improve_hz = before_abs - after_abs;

        entry.after_cfo_hz = cfo_after.res_cfo_hz;
        entry.applied_cfo_hz = applied_cfo_hz;
        entry.total_applied_cfo_hz = loop_info.total_applied_cfo_hz + applied_cfo_hz;
        entry.after_valid_frames = cfo_after.num_valid;

        if after_abs > before_abs + cfg.sync_loop_min_improve_hz
            loop_info.stop_reason = 'rejected_worse_cfo';
            loop_info.final_residual_cfo = cfo_before;
            loop_info.history(end+1) = entry;
            fprintf('    拒绝本轮补偿：残留 CFO 增大 %.9f Hz\n', after_abs - before_abs);
            break;
        end

        rx_matched = rx_matched_try;
        rx_matched_uncomp = rx_matched_uncomp_try;
        rx_16x = rx_16x_try;
        exact_frame_start = exact_frame_start_try;
        phase_index = phase_index_try;
        timing_info = timing_info_try;
        current_pass = pass_try;

        loop_info.total_applied_cfo_hz = loop_info.total_applied_cfo_hz + applied_cfo_hz;
        loop_info.final_residual_cfo = cfo_after;
        entry.accepted = true;
        entry.total_applied_cfo_hz = loop_info.total_applied_cfo_hz;
        loop_info.history(end+1) = entry;

        if after_abs <= cfg.sync_loop_cfo_tol_hz
            loop_info.stop_reason = 'cfo_tol_after_compensation';
            break;
        end

        if improve_hz < cfg.sync_loop_min_improve_hz
            loop_info.stop_reason = 'small_cfo_improvement';
            break;
        end
    end

    if isempty(loop_info.final_residual_cfo)
        loop_info.final_residual_cfo = estimate_residual_cfo_from_pn(current_pass.PN, ...
            current_pass.phase_index, current_pass.frame_heads, current_pass.corr_peak, cfg.Fs, cfg);
    end
end

function cfo = estimate_residual_cfo_from_pn(PN_matrix, phase_index, frame_heads, corr_peak, Fs, cfg)
% 基于多帧 PN 相干积累估计残留 CFO。
% 每帧先用对应本地 PN 去调制并相干求和，得到帧间复观测序列。
% 默认方法为“帧间相参 FFT 粗精估计 + 残差相位斜率细化”。

    [M, Lpn] = size(PN_matrix);
    phase_index = phase_index(:);
    frame_heads = double(frame_heads(:));
    corr_peak = double(corr_peak(:));

    M = min([M, numel(phase_index), numel(frame_heads), numel(corr_peak)]);
    PN_matrix = PN_matrix(1:M, :);
    phase_index = phase_index(1:M);
    frame_heads = frame_heads(1:M);
    corr_peak = corr_peak(1:M);

    z = zeros(M, 1);
    amp = zeros(M, 1);
    for k = 1:M
        local_pn = PN_gen(phase_index(k)).';
        local_pn = local_pn(:).';
        if length(local_pn) ~= Lpn
            error('第 %d 帧本地 PN 长度不一致：本地=%d，矩阵=%d。', ...
                k, length(local_pn), Lpn);
        end
        z(k) = sum(PN_matrix(k, :) .* conj(local_pn));
        amp(k) = abs(z(k));
    end

    amp_med = median(amp(amp > 0));
    if isempty(amp_med) || amp_med <= 0
        error('PN 相干积累幅度异常，无法估计残留 CFO。');
    end

    valid = amp > 0.3 * amp_med;
    peak_med = median(corr_peak(corr_peak > 0));
    if ~isempty(peak_med) && peak_med > 0
        valid = valid & (corr_peak > 0.3 * peak_med);
    end

    if sum(valid) < 5
        warning('用于残留 CFO 估计的有效帧偏少，放宽异常帧剔除门限。');
        valid = amp > 0.1 * amp_med;
        if ~isempty(peak_med) && peak_med > 0
            valid = valid & (corr_peak > 0.1 * peak_med);
        end
    end
    if sum(valid) < 3
        error('有效帧数不足，无法进行残留 CFO 估计。');
    end

    z_valid = z(valid);
    amp_valid = amp(valid);
    peak_valid = corr_peak(valid);
    t_abs = frame_heads(valid) / Fs;
    t = t_abs - t_abs(1);

    amp_w = amp_valid / (max(amp_valid) + eps);
    peak_w = peak_valid / (max(peak_valid) + eps);
    w = amp_w .* max(peak_w, eps);

    method = 'phasefit';
    if isfield(cfg, 'cfo_estimator')
        method = cfg.cfo_estimator;
    end

    fft_info = struct();
    if strcmpi(method, 'fft_phasefit')
        [fft_cfo_hz, fft_info] = estimate_cfo_by_frame_fft(z_valid, t, w, cfg);
        z_res = z_valid .* exp(-1j * 2*pi*fft_cfo_hz*t(:));
        [residual_cfo_hz, ~, intercept, ~, fit_phi] = fit_cfo_by_phase(z_res, t, w);
        res_cfo_hz = fft_cfo_hz + residual_cfo_hz;
        fit_phi = fit_phi + 2*pi*fft_cfo_hz*t(:);
        phi = unwrap(angle(z_valid));
    else
        [res_cfo_hz, ~, intercept, phi, fit_phi] = fit_cfo_by_phase(z_valid, t, w);
        fft_cfo_hz = nan;
        residual_cfo_hz = nan;
    end

    cfo = struct();
    cfo.res_cfo_hz = res_cfo_hz;
    cfo.slope_rad_per_sec = 2*pi*res_cfo_hz;
    cfo.intercept = intercept;
    cfo.z = z;
    cfo.amp = amp;
    cfo.valid = valid;
    cfo.t = t;
    cfo.t_abs = t_abs;
    cfo.phi_raw = angle(z);
    cfo.phi_unwrap = phi;
    cfo.fit_phi = fit_phi;
    cfo.frame_heads_valid = frame_heads(valid);
    cfo.phase_index_valid = phase_index(valid);
    cfo.weight = w;
    cfo.num_valid = sum(valid);
    cfo.estimator = method;
    cfo.fft_cfo_hz = fft_cfo_hz;
    cfo.phasefit_residual_cfo_hz = residual_cfo_hz;
    cfo.fft_info = fft_info;
end

function [cfo_hz, slope, intercept, phi, fit_phi] = fit_cfo_by_phase(z_seq, t, w)
    phi = unwrap(angle(z_seq));
    X = [t(:), ones(numel(t), 1)];
    Xw = bsxfun(@times, X, w(:));
    beta = (X' * Xw) \ (X' * (w(:) .* phi(:)));
    slope = beta(1);
    intercept = beta(2);
    fit_phi = X * beta;
    cfo_hz = slope / (2*pi);
end

function [cfo_hz, fft_info] = estimate_cfo_by_frame_fft(z_seq, t, w, cfg)
    z_seq = z_seq(:);
    t = t(:);
    w = w(:);

    M = numel(z_seq);
    Ts = median(diff(t));
    if ~isfinite(Ts) || Ts <= 0
        error('帧间时间轴异常，无法进行帧间相参 FFT。');
    end

    oversample = 128;
    if isfield(cfg, 'cfo_fft_oversample')
        oversample = cfg.cfo_fft_oversample;
    end

    nfft = 2^nextpow2(max(M * oversample, 1024));
    win = make_frame_fft_window(M, cfg);
    x = z_seq .* sqrt(w + eps) .* win;

    spec = fftshift(fft(x, nfft));
    pow = abs(spec).^2;
    freq_axis = ((0:nfft-1).' - floor(nfft/2)) / (nfft * Ts);

    [~, idx_peak] = max(pow);
    bin_offset = parabolic_peak_offset(pow, idx_peak);
    df = 1 / (nfft * Ts);
    cfo_hz = freq_axis(idx_peak) + bin_offset * df;

    fft_info = struct();
    fft_info.freq_axis = freq_axis;
    fft_info.power = pow;
    fft_info.fft_cfo_hz = cfo_hz;
    fft_info.nfft = nfft;
    fft_info.frame_interval_s = Ts;
    fft_info.native_resolution_hz = 1 / (M * Ts);
    fft_info.fft_bin_spacing_hz = df;
    fft_info.oversample = oversample;
end

function win = make_frame_fft_window(M, cfg)
    win_name = 'hann';
    if isfield(cfg, 'cfo_fft_window')
        win_name = cfg.cfo_fft_window;
    end

    switch lower(win_name)
        case {'hann', 'hanning'}
            n = (0:M-1).';
            win = 0.5 - 0.5*cos(2*pi*n/max(M-1, 1));
        otherwise
            win = ones(M, 1);
    end
end

function delta = parabolic_peak_offset(power_spectrum, idx_peak)
    if idx_peak <= 1 || idx_peak >= numel(power_spectrum)
        delta = 0;
        return;
    end

    y1 = power_spectrum(idx_peak - 1);
    y2 = power_spectrum(idx_peak);
    y3 = power_spectrum(idx_peak + 1);
    denom = y1 - 2*y2 + y3;
    if abs(denom) < eps
        delta = 0;
    else
        delta = 0.5 * (y1 - y3) / denom;
        delta = max(min(delta, 0.5), -0.5);
    end
end

