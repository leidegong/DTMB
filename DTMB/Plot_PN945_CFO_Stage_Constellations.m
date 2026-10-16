function Plot_PN945_CFO_Stage_Constellations(rx_initial, initial_cfo_info, PN_after_loop, cfg)
% Plot_PN945_CFO_Stage_Constellations
% 对比 PN945 初始 CFO 补偿各阶段的帧头 PN 星座：
%   1. 匹配滤波后、CFO 补偿前
%   2. 小数 CFO 补偿后
%   3. 小数 + 整数 CFO 补偿后
%   4. 残留 CFO 闭环精同步后

    if nargin < 4 || isempty(cfg)
        cfg = struct();
    end
    if ~isstruct(initial_cfo_info) || ~isfield(initial_cfo_info, 'rx_frac') || ~isfield(initial_cfo_info, 'rx_comp')
        warning('PN945:CFOStagePlotSkipped', '缺少 rx_frac/rx_comp，跳过 PN945 CFO 阶段星座图。');
        return;
    end

    max_points = get_cfg_scalar(cfg, 'multi_frame_constellation_max_points', 250000);
    num_frames = get_cfg_scalar(cfg, 'pn945_cfo_constellation_frames', 20);

    pn_initial = extract_pn_symbol_samples(rx_initial, initial_cfo_info.frame_start_seed, cfg, num_frames);
    pn_frac = extract_pn_symbol_samples(initial_cfo_info.rx_frac, initial_cfo_info.frame_start_seed, cfg, num_frames);
    pn_comp = extract_pn_symbol_samples(initial_cfo_info.rx_comp, initial_cfo_info.frame_start_seed, cfg, num_frames);
    pn_loop = PN_after_loop;

    initial_pts = limit_complex_points(pn_initial, max_points);
    frac_pts = limit_complex_points(pn_frac, max_points);
    comp_pts = limit_complex_points(pn_comp, max_points);
    loop_pts = limit_complex_points(pn_loop, max_points);

    figure('Name', 'PN945 CFO 补偿阶段帧头星座对比', 'Color', 'w');

    subplot(2,2,1);
    scatter(real(initial_pts), imag(initial_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title('初始匹配滤波后');

    subplot(2,2,2);
    scatter(real(frac_pts), imag(frac_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title(sprintf('小数 CFO 补偿后，%.3f Hz', initial_cfo_info.f_frac));

    subplot(2,2,3);
    scatter(real(comp_pts), imag(comp_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title(sprintf('小数+整数 CFO 后，bin %+d', initial_cfo_info.m_hat));

    subplot(2,2,4);
    scatter(real(loop_pts), imag(loop_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title('残留 CFO 闭环精同步后');
end

function pn_mat = extract_pn_symbol_samples(rx_signal, frame_start_seed, cfg, num_frames_req)
    rx_signal = rx_signal(:);
    frame_start = round(frame_start_seed);
    sps = cfg.sps;
    pn_len = cfg.pn_length;
    frame_len_samp = cfg.frame_symbol_len * sps;

    max_frames_by_len = floor((numel(rx_signal) - frame_start - (pn_len - 1) * sps) / frame_len_samp) + 1;
    num_frames = min([num_frames_req, max_frames_by_len, get_cfg_scalar(cfg, 'target_frames', num_frames_req)]);
    if num_frames <= 0
        pn_mat = complex(zeros(0, pn_len));
        return;
    end

    pn_mat = complex(nan(num_frames, pn_len));
    for m = 1:num_frames
        st = frame_start + (m - 1) * frame_len_samp;
        idx = st:sps:(st + (pn_len - 1) * sps);
        if idx(1) < 1 || idx(end) > numel(rx_signal)
            continue;
        end
        pn_mat(m, :) = rx_signal(idx).';
    end
end

function pts = limit_complex_points(x, max_points)
    pts = x(:);
    valid = isfinite(real(pts)) & isfinite(imag(pts));
    pts = pts(valid);
    if numel(pts) > max_points
        idx = round(linspace(1, numel(pts), max_points));
        pts = pts(idx);
    end
end

function value = get_cfg_scalar(cfg, field_name, default_value)
    if isfield(cfg, field_name) && ~isempty(cfg.(field_name))
        value = cfg.(field_name);
    else
        value = default_value;
    end
end
