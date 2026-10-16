function [rx_matched_out, rx_matched_corrected, rx_16x, rx_16x_uncomp, ...
          rx_sym_comp, rx_sym_uncomp, exact_frame_start_out, timing_info] = ...
          Timing_Phase_Refine_8x(rx_matched, rx_matched_uncomp, exact_frame_start, ...
                                 sps, detected_mode, phase_index)
% Timing_Phase_Refine_8x
% 功能：
%   对已做小数频偏补偿的信号做 8 倍插值，
%   在 16 倍采样域搜索最佳采样相位，再重构回 2 倍采样信号
%
% 输入：
%   rx_matched            - CFO补偿后的2倍采样信号
%   rx_matched_uncomp     - CFO补偿前的2倍采样信号（用于对比绘图）
%   exact_frame_start     - Symbol_Sync 输出的精确帧头起点
%   sps                   - 当前过采样率（你的代码里为2）
%   detected_mode         - 'PN945' 或 'PN420'
%   phase_index           - 当前 PN 相位号
%
% 输出：
%   rx_matched_out        - 重构后的2倍采样信号
%   rx_matched_corrected  - 与 rx_matched_out 一致，保留旧变量接口
%   rx_16x                - 8倍插值后的补偿后信号
%   rx_16x_uncomp         - 8倍插值后的补偿前信号
%   rx_sym_comp           - 最佳采样相位抽取得到的1倍符号率补偿后符号
%   rx_sym_uncomp         - 最佳采样相位抽取得到的1倍符号率补偿前符号
%   exact_frame_start_out - 修正后的精确帧头位置（对应新2倍采样序列）
%   timing_info           - 结构体，保存中间结果

    disp('>>> 正在执行最佳采样相位校正 (16倍过采样精细寻优)...');

    interp_factor = 8;
    sps_16x = sps * interp_factor;

    disp('    -> 正在对接收信号进行 8 倍插值 (提升至16倍采样)...');
    rx_16x = interp(double(rx_matched(:)), interp_factor);
    rx_16x_uncomp = interp(double(rx_matched_uncomp(:)), interp_factor);

    % 原始帧头位置映射到16倍采样域
    exact_frame_start_16x = (exact_frame_start - 1) * interp_factor + 1;

    % 模式参数
    if strcmp(detected_mode, 'PN945')
        temp_pn_tot = 945;
        temp_pn_eff = 511;
        temp_pn_start = 218;
    elseif strcmp(detected_mode, 'PN420')
        temp_pn_tot = 420;
        temp_pn_eff = 255;
        temp_pn_start = 165;
    else
        error('未知 detected_mode');
    end

    % 本地参考 PN
    local_pn_full = PN_gen(phase_index);
    local_pn_ref = local_pn_full(temp_pn_start : temp_pn_start + temp_pn_eff - 1).';

    search_len_16x = temp_pn_tot * sps_16x;
    search_range = -sps_16x : sps_16x;

    max_corr_val = 0;
    best_offset_16x = 0;

    % 搜索最佳偏移
    for offset = search_range
        idx_test = exact_frame_start_16x + offset;

        if idx_test < 1 || (idx_test + search_len_16x - 1) > length(rx_16x)
            continue;
        end

        temp_samp_16x = rx_16x(idx_test : idx_test + search_len_16x - 1);
        temp_sym = temp_samp_16x(1 : sps_16x : end);
        temp_pn = temp_sym(temp_pn_start : temp_pn_start + temp_pn_eff - 1);

        corr_val = abs(temp_pn(:)' * local_pn_ref(:));

        if corr_val > max_corr_val
            max_corr_val = corr_val;
            best_offset_16x = offset;
        end
    end

    best_idx_16x_abs = exact_frame_start_16x + best_offset_16x;

    if best_idx_16x_abs < 1
        disp('警告：偏移量过小，超出数组左侧边界，已强制归 1');
        best_idx_16x_abs = 1;
    end

    % 计算2倍采样重构的起始抽头
    start_idx_16x = mod(best_idx_16x_abs - 1, interp_factor) + 1;

    % 重构回2倍采样
    rx_matched_out = rx_16x(start_idx_16x : interp_factor : end);
    rx_matched_corrected = rx_matched_out;

    % 修正后的精确帧起点
    exact_frame_start_out = (best_idx_16x_abs - start_idx_16x) / interp_factor + 1;

    % 1倍符号率序列（补偿前/后）
    rx_sym_uncomp = rx_16x_uncomp(start_idx_16x : interp_factor * sps : end);
    rx_sym_comp   = rx_16x(start_idx_16x : interp_factor * sps : end);

    % 输出调试信息
    timing_info = struct();
    timing_info.interp_factor = interp_factor;
    timing_info.sps_16x = sps_16x;
    timing_info.best_offset_16x = best_offset_16x;
    timing_info.best_idx_16x_abs = best_idx_16x_abs;
    timing_info.start_idx_16x = start_idx_16x;
    timing_info.max_corr_val = max_corr_val;

    disp(['    -> ✅ 当前最优偏移量: ', num2str(best_offset_16x), ' (单位: 1/16 符号周期)']);
    disp('    -> ✅ 已将信号最佳相位抽头重构为 2 倍采样，无缝对接后续模块');
end