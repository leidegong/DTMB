function pass = refine_and_build_matrices(rx_matched, rx_16x, exact_frame_start, phase_index, ...
                                          sps, detected_mode, num_frames_est, timing_info, ...
                                          debug_plot_each_frame)
    [frame_head_positions_exact, phase_index_list, corr_peak_list, frame_status, refine_debug_info] = ...
        Refine_All_Frame_Heads_PN945(rx_matched, exact_frame_start, phase_index, ...
                                     sps, num_frames_est, debug_plot_each_frame);

    fprintf('逐帧 PN 精同步成功/失败帧数：%d / %d\n', ...
        sum(frame_status == "ok"), sum(frame_status == "fail"));

    valid_sync_idx = (frame_status == "ok");
    frame_heads0 = frame_head_positions_exact(valid_sync_idx);
    phase_index0 = phase_index_list(valid_sync_idx);
    corr_peak0 = corr_peak_list(valid_sync_idx);

    [PN_sym_matrix, Frame_sym_matrix, timing_track, valid_frame_heads] = ...
        Build_PN_Frame_Matrix_From_InterpStream(rx_16x, frame_heads0, phase_index0, ...
                                               timing_info, sps, detected_mode);

    valid_timing = timing_track.valid_mask(:);
    pass = struct();
    pass.PN = PN_sym_matrix;
    pass.Frame = Frame_sym_matrix;
    pass.frame_heads = valid_frame_heads(:);
    pass.phase_index = phase_index0(valid_timing);
    pass.corr_peak = corr_peak0(valid_timing);
    pass.timing_track = timing_track;
    pass.frame_status = frame_status;
    pass.refine_debug_info = refine_debug_info;

    assert_matrix_consistency(pass);
    print_frame_offset_drift_summary(pass);

    fprintf('前若干帧的帧头位置 / PN 相位号 / 相关峰值：\n');
    for k = 1:min(10, numel(pass.frame_heads))
        fprintf('  %3d: 帧头 = %10.2f，相位号 = %3d，峰值 = %.6g\n', ...
            k, pass.frame_heads(k), pass.phase_index(k), pass.corr_peak(k));
    end
end

function print_frame_offset_drift_summary(pass)
    if ~isfield(pass, 'refine_debug_info') || ~isfield(pass.refine_debug_info, 'best_offset_list')
        return;
    end

    offsets = double(pass.refine_debug_info.best_offset_list(:));
    frame_idx = (0:numel(offsets)-1).';
    ok = (pass.frame_status(:) == "ok") & isfinite(offsets);
    if sum(ok) < 3
        return;
    end

    p = polyfit(frame_idx(ok), offsets(ok), 1);
    offset_slope_samp_per_frame = p(1);

    nominal_frame_len = nan;
    if isfield(pass.refine_debug_info, 'nominal_frame_len_samp')
        nominal_frame_len = pass.refine_debug_info.nominal_frame_len_samp;
    end

    if isfinite(nominal_frame_len) && nominal_frame_len > 0
        sfo_ppm = offset_slope_samp_per_frame / nominal_frame_len * 1e6;
        fprintf('帧头偏移线性斜率：%.6f 采样点/帧，等效采样时钟偏差 %.3f ppm\n', ...
            offset_slope_samp_per_frame, sfo_ppm);
    else
        fprintf('帧头偏移线性斜率：%.6f 采样点/帧\n', offset_slope_samp_per_frame);
    end
end

