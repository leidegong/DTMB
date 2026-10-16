function [phase_index, position] = Symbol_Sync(rx_signal, coarse_starts, sps, detected_mode)
% Symbol_Sync: 兼容 PN420 和 PN945 的符号精同步
    
    if length(coarse_starts) < 1
        error('未检测到帧粗同步位置！');
    end

    search_margin = 100 * sps;      % 容错余量
    idx_start = max(1, coarse_starts(1) - search_margin);

    if strcmp(detected_mode, 'PN420')
        % --- 原有 PN420 逻辑 (利用相位滑动差分) ---
        len_sym = 4200; 
        local_pn_0 = PN_gen(0); % 注意：你需要确保你的 PN_gen 能正确生成 PN420
        local_pn_0_up = zeros(length(local_pn_0) * sps, 1);
        local_pn_0_up(1:sps:end) = local_pn_0; % 间隔插零
        h_srrc = rcosdesign(0.05, 8, sps, 'sqrt'); % 生成与接收端一致的滤波器
        local_pn_0_samp = conv(local_pn_0_up, h_srrc, 'same'); % 滤波并保持长度一致
        L_pn = length(local_pn_0_samp);
        
        % 计算 p1
        idx_end_1 = min(length(rx_signal), coarse_starts(1) + L_pn + search_margin); 
        [c1, lags1] = xcorr(rx_signal(idx_start : idx_end_1), local_pn_0_samp);
        [~, max_idx_1] = max(abs(c1));
        p1 = idx_start + lags1(max_idx_1);
        
        % 计算 p2
        idx_start_2 = max(1, coarse_starts(2) - search_margin);
        idx_end_2 = min(length(rx_signal), coarse_starts(2) + L_pn + search_margin); 
        [c2, lags2] = xcorr(rx_signal(idx_start_2 : idx_end_2), local_pn_0_samp);
        [~, max_idx_2] = max(abs(c2));
        p2 = idx_start_2 + lags2(max_idx_2);
        
        % PN420 独有的差分推导 k_pn
        delta_pos_sym = round((p2 - p1) / sps) - len_sym;
        phase_index = -1;
        for k = 0:224
            if k <= 111, expected_delta = (-1)^(k+1) * (k + 1);
            elseif k <= 223, expected_delta = -((-1)^(223-k+1) * (223-k + 1));
            else, expected_delta = 0; end
            if expected_delta == delta_pos_sym
                phase_index = k; break;
            end
        end
        
        % PN420 计算最终 position 的代码略 (保留你原有的 4~5 步计算 position)
        if phase_index <= 112
            if mod(phase_index, 2) == 0 % 偶数
            delta_p = phase_index / 2;
            else % 奇数
            delta_p = -(phase_index + 1) / 2;
            end
        else
        % k_pn > 112，等效为计算 224 - k_pn
            temp_k = 224 - phase_index;
            if mod(temp_k, 2) == 0
                delta_p = temp_k / 2;
            else
                delta_p = -(temp_k + 1) / 2;
            end
        end
    
    % 5. 精确计算绝对起始位置 (换算回采样点)
        position = p1 - (delta_p * sps);
        
    elseif strcmp(detected_mode, 'PN945')
        len_sym = 4725; 
        local_pn_0 = PN_gen(0); 
        local_pn_0_up = zeros(length(local_pn_0) * sps, 1);
        local_pn_0_up(1:sps:end) = local_pn_0;
        h_srrc = rcosdesign(0.05, 8, sps, 'sqrt');
        local_pn_0_samp = conv(local_pn_0_up, h_srrc, 'same'); 
        L_pn = length(local_pn_0_samp);
        
        % 计算 p1
        idx_end_1 = min(length(rx_signal), coarse_starts(1) + L_pn + search_margin); 
        [c1, lags1] = xcorr(rx_signal(idx_start : idx_end_1), local_pn_0_samp);
        [~, max_idx_1] = max(abs(c1));
        p1 = idx_start + lags1(max_idx_1);
        
        % 计算 p2
        idx_start_2 = max(1, coarse_starts(2) - search_margin);
        idx_end_2 = min(length(rx_signal), coarse_starts(2) + L_pn + search_margin); 
        [c2, lags2] = xcorr(rx_signal(idx_start_2 : idx_end_2), local_pn_0_samp);
        [~, max_idx_2] = max(abs(c2));
        p2 = idx_start_2 + lags2(max_idx_2);
        
        delta_pos_sym = round((p2 - p1) / sps)-len_sym;
        phase_index = -1;
        for k = 0:199
            if k <= 99
                expected_delta = (-1)^(k+1) * (k + 1);
            else
                % 根据表格推导出的最简下半区公式
                expected_delta = (-1)^(k+1) * (200 - k); 
            end
    
            if expected_delta == delta_pos_sym
                phase_index = k;
                break;
            end
        end
       
        if phase_index <= 99
            if mod(phase_index, 2) == 0 % 偶数
                delta_p = phase_index / 2;
            else % 奇数
                delta_p = -(phase_index + 1) / 2;
            end
        else
            temp_k = 200 - phase_index;
            if mod(temp_k, 2) == 0
                delta_p = temp_k / 2;
            else
                delta_p = -(temp_k + 1) / 2;
            end
        end
    
    % 5. 精确计算绝对起始位置 (换算回采样点)
        position = p1 - (delta_p * sps);

    else
        error('未知的帧头模式');
    end
    % =========================================================
    % ---> 新增：计算整段接收信号与 0 相位本地信号的全局互相关 <---
    % =========================================================
    disp('    -> 正在计算全部信号与 0 相位本地序列的全局互相关，请稍候...');
    
    local_pn_0_tmp = PN_gen(0);
    if strcmp(detected_mode, 'PN945')
        % 提取 PN945 的核心 511 符号 (索引 218 到 728)
        local_pn_0_core = local_pn_0_tmp(218 : 728); 
    else
        % 提取 PN420 的核心 255 符号 (假设索引为 165 到 419)
        local_pn_0_core = local_pn_0_tmp(165 : 419); 
    end
    
    % 2. 恢复正确的匹配滤波上采样 (弃用 repelem)
    local_pn_0_up_global = zeros(length(local_pn_0_core) * sps, 1);
    local_pn_0_up_global(1:sps:end) = local_pn_0_core;
    h_srrc = rcosdesign(0.05, 8, sps, 'sqrt');
    local_pn_0_global_samp = conv(local_pn_0_up_global, h_srrc, 'same');
    
    % 3. 全局互相关
    [c_all, lags_all] = xcorr(rx_signal, local_pn_0_global_samp);

    if strcmp(detected_mode, 'PN420')
        sym_per_frame = 4200;
    else
        sym_per_frame = 4725;
    end
    
    % 限制 X 轴显示范围 (例如看前 3 帧)
    zoom_samples = 3 * sym_per_frame * sps + 2000;
    % 生成完整的 0相位本地序列 (带有过采样并进行匹配成型)
    % 只取核心部分 (PN945 的 218~728)
%     local_pn_0_tmp = PN_gen(0);
%     local_pn_0_up_global = repelem(local_pn_0_tmp, sps);
% %     local_pn_0_up_global = zeros(length(local_pn_0_tmp) * sps, 1);
% %     local_pn_0_up_global(1:sps:end) = local_pn_0_tmp;
% %     h_srrc = rcosdesign(0.05, 8, sps, 'sqrt');
% %     local_pn_0_global = conv(local_pn_0_up_global, h_srrc, 'same');
%     
%     % 核心：对整个长度为 255 帧的 rx_signal 进行互相关
%     [c_all, lags_all] = xcorr(rx_signal, local_pn_0_up_global);
%     if strcmp(detected_mode, 'PN420')
%         sym_per_frame = 4200;
%     else
%         sym_per_frame = 4725;
%     end
%     % 10帧的采样点数 (加上 2000 个点的余量以便看清边缘)
%     zoom_samples = 10 * sym_per_frame * sps + 2000;

    % 创建绘图
    figure('Name', ['全局互相关: 所有帧 vs ', detected_mode, ' (k=0)'], 'Color', 'w');
    plot(lags_all, abs(c_all), 'b-', 'LineWidth', 1);
    xlim([-1000, zoom_samples]);

    title(['全部接收信号 (255帧) 与 0 相位本地序列 ', detected_mode, ' 全局互相关']);
    xlabel('滞后采样点数 (Lags)');
    ylabel('互相关幅度 (Magnitude)');
    grid on;
end