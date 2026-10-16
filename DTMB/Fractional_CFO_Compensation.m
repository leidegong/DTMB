function [rx_matched_corrected, rx_matched_uncomp, phase_trajectory, inst_cfo_seq] = ...
    Fractional_CFO_Compensation(rx_matched, coarse_frame_starts, phase_offset_per_samp)
% Fractional_CFO_Compensation
% 功能：
%   根据 Frame_Sync 输出的每帧小数频偏估计值，
%   通过线性插值 + 相位积分，对整段信号进行连续小数频偏补偿
%
% 输入：
%   rx_matched            - 匹配滤波后的接收信号
%   coarse_frame_starts   - 粗同步得到的每帧起始位置
%   phase_offset_per_samp - 每帧对应的“每采样点相位偏移”
%
% 输出：
%   rx_matched_corrected  - 补偿后信号
%   rx_matched_uncomp     - 补偿前信号备份
%   phase_trajectory      - 整段信号的累计相位轨迹
%   inst_cfo_seq          - 插值得到的瞬时 CFO 序列

    total_len = length(rx_matched);

    % 保证列向量
    rx_matched = rx_matched(:);
    x_known = coarse_frame_starts(:);
    y_known = phase_offset_per_samp(:);
    x_query = (1:total_len).';

    if isempty(x_known) || isempty(y_known)
        error('coarse_frame_starts 或 phase_offset_per_samp 为空，无法进行小数频偏补偿。');
    end

    % 1. 构造整段瞬时频偏序列
    if length(x_known) < 2
        inst_cfo_seq = ones(total_len, 1) * y_known(1);
    else
        inst_cfo_seq = interp1(x_known, y_known, x_query, 'linear', 'extrap');
        inst_cfo_seq(x_query < x_known(1)) = y_known(1);
        inst_cfo_seq(x_query > x_known(end)) = y_known(end);
    end

    % 2. 相位积分
    phase_trajectory = cumsum(inst_cfo_seq);

    % 3. 去旋转补偿
    rx_matched_uncomp = rx_matched;
    rot_factor = exp(-1j * phase_trajectory);
    rx_matched_corrected = rx_matched .* rot_factor;

    disp('    -> ✅ 逐帧动态跟踪的小数频偏补偿完成！(已消除长序列相位累积误差)');
end