function [PN_sym_matrix, Frame_sym_matrix, timing_track, valid_frame_heads] = ...
    Build_PN_Frame_Matrix_From_InterpStream(rx_16x, frame_head_positions, phase_index_list, ...
                                           timing_info, sps, detected_mode)
% Build_PN_Frame_Matrix_From_InterpStream
%   Per-frame timing tracking in the interpolated sampling domain. For each
%   refined 2x-sample frame head, search a small 16x-sample neighborhood and
%   extract PN/full-frame symbols at that frame's best timing phase.

    rx_16x = rx_16x(:);
    frame_head_positions = frame_head_positions(:);
    phase_index_list = phase_index_list(:);

    interp_factor = timing_info.interp_factor;
    start_idx_16x = timing_info.start_idx_16x;
    sps_16x = sps * interp_factor;

    switch upper(detected_mode)
        case 'PN420'
            Lpn = 420;
            core_start = 165;
            core_len = 255;
        case 'PN945'
            Lpn = 945;
            core_start = 218;
            core_len = 511;
        otherwise
            error('Unsupported frame head mode: %s', detected_mode);
    end

    Ldata = 3780;
    Lframe = Lpn + Ldata;

    num_frames = min(length(frame_head_positions), length(phase_index_list));
    PN_sym_matrix = zeros(num_frames, Lpn);
    Frame_sym_matrix = zeros(num_frames, Lframe);
    valid_mask = false(num_frames, 1);

    timing_offset_16x = nan(num_frames, 1);
    timing_corr_peak = nan(num_frames, 1);
    timing_start_16x = nan(num_frames, 1);

    search_range = -sps_16x:sps_16x;
    total_len = length(rx_16x);

    row = 0;
    for k = 1:num_frames
        base_idx_16x = (round(frame_head_positions(k)) - 1) * interp_factor + start_idx_16x;
        local_pn = PN_gen(phase_index_list(k)).';
        local_core = local_pn(core_start:core_start + core_len - 1);

        best_corr = -inf;
        best_offset = 0;
        best_pn = [];

        for offset = search_range
            idx = base_idx_16x + offset;
            pn_stop_idx = idx + (Lpn - 1) * sps_16x;

            if idx < 1 || pn_stop_idx > total_len
                continue;
            end

            pn_sym = rx_16x(idx:sps_16x:pn_stop_idx);
            pn_core = pn_sym(core_start:core_start + core_len - 1);
            corr_val = abs(pn_core(:)' * local_core(:));

            if corr_val > best_corr
                best_corr = corr_val;
                best_offset = offset;
                best_pn = pn_sym;
            end
        end

        best_idx = base_idx_16x + best_offset;
        frame_stop_idx = best_idx + (Lframe - 1) * sps_16x;

        if isempty(best_pn) || best_idx < 1 || frame_stop_idx > total_len
            continue;
        end

        frame_sym = rx_16x(best_idx:sps_16x:frame_stop_idx);
        if length(frame_sym) ~= Lframe
            continue;
        end

        row = row + 1;
        PN_sym_matrix(row, :) = best_pn(:).';
        Frame_sym_matrix(row, :) = frame_sym(:).';
        valid_mask(k) = true;
        timing_offset_16x(k) = best_offset;
        timing_corr_peak(k) = best_corr;
        timing_start_16x(k) = best_idx;
    end

    PN_sym_matrix = PN_sym_matrix(1:row, :);
    Frame_sym_matrix = Frame_sym_matrix(1:row, :);
    valid_frame_heads = frame_head_positions(valid_mask);

    timing_track = struct();
    timing_track.offset_16x = timing_offset_16x(valid_mask);
    timing_track.offset_symbols = timing_track.offset_16x / sps_16x;
    timing_track.corr_peak = timing_corr_peak(valid_mask);
    timing_track.start_16x = timing_start_16x(valid_mask);
    timing_track.valid_mask = valid_mask;

    disp(['    -> 16x per-frame timing extracted PN/frame matrices: ', ...
        num2str(size(PN_sym_matrix,1)), ' x ', num2str(size(PN_sym_matrix,2)), ...
        ' / ', num2str(size(Frame_sym_matrix,1)), ' x ', num2str(size(Frame_sym_matrix,2))]);
    fprintf('    -> 16x timing offset range: %.4f to %.4f symbols\n', ...
        min(timing_track.offset_symbols), max(timing_track.offset_symbols));
end
