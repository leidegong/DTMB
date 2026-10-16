function [detected_mode, coarse_frame_starts, metric,phase_offset_per_samp] = Frame_Sync(rx_signal, sps, debug)
% Frame_Sync: 基于延迟自相关的 DTMB 帧头模式盲检与粗同步
% 功能: 同时测试 PN420 和 PN945 的循环特征，自动判决发射模式

    disp('>>> 开始执行帧头模式盲检 (PN420 vs PN945)...');
    % --- 1. 定义两种模式的参数 (符号数) ---
    % 假设 A: PN420
    N1_sym_420 = 165;    N2_sym_420 = 255;    Frame_sym_420 = 4200;
    % 假设 B: PN945
    N1_sym_945 = 434;    N2_sym_945 = 511;    Frame_sym_945 = 4725;
    
    % 映射到采样点级
    N1_samp_420 = N1_sym_420 * sps;  N2_samp_420 = N2_sym_420 * sps;  Frame_samp_420 = Frame_sym_420 * sps;
    N1_samp_945 = N1_sym_945 * sps;  N2_samp_945 = N2_sym_945 * sps;  Frame_samp_945 = Frame_sym_945 * sps;
    
    % --- 2. 并行计算两种假设的滑动自相关度量 ---
    
    % 【计算 PN420 度量】
    rx_delay_420 = zeros(size(rx_signal));
    rx_delay_420(1:end-N2_samp_420) = rx_signal(N2_samp_420+1:end); 
    num_420 = filter(ones(N1_samp_420, 1), 1, rx_signal .* conj(rx_delay_420));
    
    % 计算两路信号各自的能量并取平均作为分母 (防止噪声除零放大)
    energy_sig_420 = filter(ones(N1_samp_420, 1), 1, abs(rx_signal).^2);
    energy_delay_420 = filter(ones(N1_samp_420, 1), 1, abs(rx_delay_420).^2);
    den_420 = (energy_sig_420 + energy_delay_420) / 2;
    metric_420 = abs(num_420) ./ (den_420 + eps);
    metric_420(1:N1_samp_420) = 0;
    
    % 【计算 PN945 度量】
    rx_delay_945 = zeros(size(rx_signal));
    rx_delay_945(1:end-N2_samp_945) = rx_signal(N2_samp_945+1:end); 
    num_945 = filter(ones(N1_samp_945, 1), 1, rx_signal .* conj(rx_delay_945));
    
    % 计算两路信号各自的能量并取平均作为分母
    energy_sig_945 = filter(ones(N1_samp_945, 1), 1, abs(rx_signal).^2);
    energy_delay_945 = filter(ones(N1_samp_945, 1), 1, abs(rx_delay_945).^2);
    den_945 = (energy_sig_945 + energy_delay_945) / 2;
    metric_945 = abs(num_945) ./ (den_945 + eps);
    metric_945(1:N1_samp_945) = 0;
    
    % --- 3. 模式判决 (比较两种度量的全局最大值) ---
    max_peak_420 = max(metric_420);
    max_peak_945 = max(metric_945);
    
    % 设定一个基础门限，防止纯噪声误判 (通常合法信号的相关峰值在 0.5~1 之间)
    detect_threshold = 0.3; 
    
    if (max_peak_420 > max_peak_945) && (max_peak_420 > detect_threshold)
        detected_mode = 'PN420';
        metric = metric_420;
        active_N1 = N1_samp_420;
        active_Frame = Frame_samp_420;
        disp(['    -> 盲检结果: 信号判定为 【PN420】 模式！(峰值:', num2str(max_peak_420), ')']);
        
    elseif (max_peak_945 > max_peak_420) && (max_peak_945 > detect_threshold)
        detected_mode = 'PN945';
        metric = metric_945;
        active_N1 = N1_samp_945;
        active_Frame = Frame_samp_945;
        disp(['    -> 盲检结果: 信号判定为 【PN945】 模式！(峰值:', num2str(max_peak_945), ')']);
        
    else
        detected_mode = 'Unknown';
        coarse_frame_starts = [];
        metric = max(metric_420, metric_945);
        disp('    -> [警告] 盲检失败: 未检测到明显帧头特征，可能是纯噪声或极低信噪比！');
        return;
    end
    
    % --- 4. 在判决出的度量曲线上寻峰定位 ---
    threshold = max(metric) * 0.4; 
    min_dist = active_Frame * 0.9;
    
    [pks, locs] = findpeaks(metric, 'MinPeakHeight', threshold, 'MinPeakDistance', min_dist);
    
    % 补偿 filter 带来的窗口延迟
    coarse_frame_starts = locs - active_N1 + 1;
    
    % 剔除越界异常值
    valid_idx = (coarse_frame_starts > 0) & (coarse_frame_starts < length(rx_signal));
    coarse_frame_starts = coarse_frame_starts(valid_idx);
    locs = locs(valid_idx); pks = pks(valid_idx);
    
    disp(['    -> 粗同步完成: 共提取出 ', num2str(length(coarse_frame_starts)), ' 帧数据起点。']);

    % =========================================================
    %  修改：计算每帧独立的小数频偏 (Fractional Frequency Offset)
    % =========================================================
    if strcmp(detected_mode, 'PN420')
        delay_samp = N2_samp_420;
        complex_peaks = num_420(locs);
    else
        delay_samp = N2_samp_945;
        complex_peaks = num_945(locs);
    end
    
    % 不再求全局均值，而是直接求出每个帧头峰值的独立角度
    peak_angles = angle(complex_peaks); 
    
    % 计算每帧对应的采样点相位偏移步进（此时是一个与 locs 等长的列向量）
    phase_offset_per_samp = peak_angles / delay_samp;
    
    % 【优化】使用滑动平均稍微平滑一下频偏序列，抵抗单帧深度衰落或突发噪声
    if length(phase_offset_per_samp) > 5
        phase_offset_per_samp = smoothdata(phase_offset_per_samp, 'movmean', 5);
    end
    
    avg_cfo = mean(phase_offset_per_samp);
    disp(['    -> CFO估计: 检测到 ', num2str(length(phase_offset_per_samp)), ' 个动态频偏值']);
    disp(['    -> CFO估计: 平均采样点相位偏移为 ', num2str(avg_cfo), ' rad']);
    % =========================================================
    % --- 5. 调试绘图 ---
    if debug
        figure('Name', '帧头模式盲检验证', 'Color', 'w', 'Position', [150, 150, 1000, 600]);
        
        subplot(2,1,1);
        plot(metric_420, 'b-', 'LineWidth', 1);
        title('假设A: PN420 匹配度量曲线'); ylabel('相关度 M(k)'); grid on;
        if strcmp(detected_mode, 'PN420')
            hold on; plot(locs, pks, 'rv', 'MarkerFaceColor', 'r');
            legend('PN420 M(k)', '有效峰值');
        end
        xlim([0, min(length(metric_420), active_Frame * 4)]); 
        
        subplot(2,1,2);
        plot(metric_945, 'g-', 'LineWidth', 1);
        title('假设B: PN945 匹配度量曲线'); xlabel('采样点索引'); ylabel('相关度 M(k)'); grid on;
        if strcmp(detected_mode, 'PN945')
            hold on; plot(locs, pks, 'rv', 'MarkerFaceColor', 'r');
            legend('PN945 M(k)', '有效峰值');
        end
        xlim([0, min(length(metric_945), active_Frame * 4)]); 
    end
end