function info = DTMB_Decode_System_Info(system_symbol_matrix)
% DTMB_Decode_System_Info 用多帧归一化相关判决64QAM系统信息模式。

    if size(system_symbol_matrix, 2) ~= 36
        error('系统信息矩阵必须每行包含 36 个符号。');
    end
    candidates = DTMB_System_Info_Candidates();
    num_frames = size(system_symbol_matrix, 1);
    scores = nan(num_frames, numel(candidates));
    gains = nan(num_frames, numel(candidates));
    for frame_idx = 1:num_frames
        r = system_symbol_matrix(frame_idx, :).';
        valid = isfinite(real(r)) & isfinite(imag(r));
        for candidate_idx = 1:numel(candidates)
            c = candidates(candidate_idx).symbols;
            use = valid & isfinite(real(c));
            if nnz(use) < 16
                continue;
            end
            gain = (c(use)' * r(use)) / (c(use)' * c(use) + eps);
            gains(frame_idx, candidate_idx) = gain;
            scores(frame_idx, candidate_idx) = abs(c(use)' * r(use))^2 / ...
                ((c(use)' * c(use)) * (r(use)' * r(use)) + eps);
        end
    end
    aggregate = median(scores, 1, 'omitnan');
    [sorted_scores, order] = sort(aggregate, 'descend');
    best_idx = order(1);
    best = candidates(best_idx);

    info = best;
    info.candidate_scores = aggregate;
    info.per_frame_scores = scores;
    info.per_frame_complex_gain = gains(:, best_idx);
    info.best_score = sorted_scores(1);
    if numel(sorted_scores) >= 2
        info.second_score = sorted_scores(2);
        info.confidence_ratio = sorted_scores(1) / (sorted_scores(2) + eps);
    else
        info.second_score = NaN;
        info.confidence_ratio = Inf;
    end
end
