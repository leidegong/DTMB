function assert_matrix_consistency(pass)
    n_pn = size(pass.PN, 1);
    n_frame = size(pass.Frame, 1);
    n_head = numel(pass.frame_heads);
    n_phase = numel(pass.phase_index);
    n_peak = numel(pass.corr_peak);

    if ~(n_pn == n_frame && n_pn == n_head && n_pn == n_phase && n_pn == n_peak)
        error(['矩阵与同步元数据行数不一致：PN=%d，帧体=%d，帧头=%d，', ...
               '相位号=%d，相关峰=%d。'], n_pn, n_frame, n_head, n_phase, n_peak);
    end
end


