function [frame_head_positions_exact, phase_index_list, corr_peak_list, frame_status, debug_info] = ...
    Refine_All_Frame_Heads_PN945(rx_signal, first_frame_start, first_phase_index, sps, num_frames, debug_plot)

    if nargin < 6
        debug_plot = false;
    end

    frame_len_sym = 4725;
    pn_len_sym = 945;
    frame_len_samp = frame_len_sym * sps;
    pn_len_samp = pn_len_sym * sps;

    search_margin_sym = 80;
    search_margin_samp = search_margin_sym * sps;

    h_srrc = rcosdesign(0.05, 8, sps, 'sqrt');
    frame_head_positions_exact = nan(num_frames, 1);
    phase_index_list = nan(num_frames, 1);
    corr_peak_list = nan(num_frames, 1);
    frame_status = strings(num_frames,1);

    debug_info = struct();
    debug_info.pred_start_list = nan(num_frames,1);
    debug_info.best_offset_list = nan(num_frames,1);
    debug_info.actual_frame_index = nan(num_frames,1);
    debug_info.pred_frame_len_samp = nan(num_frames,1);
    debug_info.nominal_frame_len_samp = frame_len_samp;

    disp('------------------------------------------------------');
    disp('>>> 开始逐帧精同步 PN945 帧头（不跳过坏帧）...');
    fprintf('%6s %8s %10s %10s %8s %8s %12s\n', ...
        '帧号', '预测起点', '精确起点', '偏移', '相位号', '峰值', '状态');

    for k = 1:num_frames

        % 第k帧固定对应第k个物理帧，不再跳帧
        actual_idx = k - 1;

        [pred_start, frame_len_pred] = predict_frame_start_with_drift( ...
            first_frame_start, frame_len_samp, frame_head_positions_exact, frame_status, k);
        this_phase = mod(first_phase_index + actual_idx, 200);

        debug_info.pred_start_list(k) = pred_start;
        debug_info.actual_frame_index(k) = k;
        debug_info.pred_frame_len_samp(k) = frame_len_pred;

        % 构造本地PN模板
        local_pn = PN_gen(this_phase).';
        local_pn_up = zeros(length(local_pn)*sps, 1);
        local_pn_up(1:sps:end) = local_pn(:);
        local_pn_samp = conv(local_pn_up, h_srrc, 'same');

        % 单次局部搜索
        [best_start, max_val, ok_flag, corr_axis, corr_abs] = ...
            local_search(rx_signal, local_pn_samp, pred_start, search_margin_samp, pn_len_samp);

        best_offset = best_start - pred_start;

        frame_head_positions_exact(k) = best_start;
        phase_index_list(k) = this_phase;
        corr_peak_list(k) = max_val;
        debug_info.best_offset_list(k) = best_offset;

        if ok_flag
            frame_status(k) = "ok";
        else
            frame_status(k) = "fail";
        end

        fprintf('%6d %12d %12.0f %10.0f %10.0f %12.4f %12s\n', ...
            k, debug_info.pred_start_list(k), frame_head_positions_exact(k), ...
            debug_info.best_offset_list(k), phase_index_list(k), ...
            corr_peak_list(k), frame_status(k));

        if debug_plot && ~isempty(corr_abs)
            figure('Name', sprintf('第%d帧精同步验证', k), 'Color', 'w');
            plot(corr_axis, corr_abs ./ max(corr_abs + eps), 'b', 'LineWidth', 1.2); hold on;
            xline(debug_info.pred_start_list(k), '--k', '预测起点', 'LineWidth', 1.0);
            if ~isnan(frame_head_positions_exact(k))
                xline(frame_head_positions_exact(k), '--r', '精确起点', 'LineWidth', 1.2);
            end
            grid on;
            xlabel('采样点位置');
            ylabel('归一化相关幅度');
            title(sprintf('第%d帧 | phase=%d | status=%s', ...
                k, phase_index_list(k), frame_status(k)));
            hold off;
        end
    end
end

function [pred_start, frame_len_pred] = predict_frame_start_with_drift( ...
    first_frame_start, nominal_frame_len_samp, frame_head_positions_exact, frame_status, k)

    actual_idx = k - 1;
    pred_start_nominal = first_frame_start + actual_idx * nominal_frame_len_samp;
    frame_len_pred = nominal_frame_len_samp;

    if k <= 3
        pred_start = round(pred_start_nominal);
        return;
    end

    prev_idx = (0:k-2).';
    ok = (frame_status(1:k-1) == "ok") & isfinite(frame_head_positions_exact(1:k-1));
    if sum(ok) < 3
        pred_start = round(pred_start_nominal);
        return;
    end

    p = polyfit(prev_idx(ok), double(frame_head_positions_exact(ok)), 1);
    frame_len_pred = p(1);

    max_sfo_ppm = 300;
    max_dev = nominal_frame_len_samp * max_sfo_ppm * 1e-6;
    if abs(frame_len_pred - nominal_frame_len_samp) > max_dev
        frame_len_pred = nominal_frame_len_samp;
        pred_start = round(pred_start_nominal);
        return;
    end

    pred_start = round(polyval(p, actual_idx));
end


function [best_start, max_val, ok_flag, corr_axis, corr_abs] = ...
    local_search(rx_signal, local_pn_samp, pred_start, margin_samp, pn_len_samp)

    total_len = length(rx_signal);
    search_start = max(1, pred_start - margin_samp);
    search_end = min(total_len, pred_start + pn_len_samp - 1 + margin_samp);

    if (search_end - search_start + 1) < length(local_pn_samp)
        best_start = nan;
        max_val = 0;
        ok_flag = false;
        corr_axis = [];
        corr_abs = [];
        return;
    end

    search_seg = rx_signal(search_start:search_end);
    [c, lags] = xcorr(search_seg, local_pn_samp);

    corr_abs = abs(c);
    [max_val, idx] = max(corr_abs);
    best_start = search_start + lags(idx);
    corr_axis = search_start + lags;

    ok_flag = (best_start >= 1) && ((best_start + pn_len_samp - 1) <= total_len);
end
