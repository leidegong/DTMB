function [output_frames, valid_frame_mask, info] = DTMB_Time_Interleave(input_frames, M, commutator_direction, commutator_phase)
% DTMB_Time_Interleave 52支路卷积交织；暖机区以NaN标记。

    if nargin < 3 || isempty(commutator_direction)
        commutator_direction = 1;
    end
    if nargin < 4 || isempty(commutator_phase)
        commutator_phase = 0;
    end
    if ~ismember(commutator_direction, [-1, 1])
        error('换向器方向只能为 +1 或 -1。');
    end
    [num_frames, symbols_per_frame] = size(input_frames);
    if symbols_per_frame ~= 3744
        error('时域交织要求每帧 3744 个数据符号。');
    end
    B = 52;
    delay_symbols = B * (B - 1) * M;
    delay_frames = round(delay_symbols / symbols_per_frame);
    stream = reshape(input_frames.', [], 1);
    output = nan(size(stream), 'like', stream);
    valid = false(size(stream));
    t0 = (0:(numel(stream)-1)).';
    branch = mod(commutator_direction * t0 + commutator_phase, B);
    source0 = t0 - branch * M * B;
    use = source0 >= 0;
    output(t0(use) + 1) = stream(source0(use) + 1);
    valid(t0(use) + 1) = true;
    output_frames = reshape(output, symbols_per_frame, []).';
    valid_matrix = reshape(valid, symbols_per_frame, []).';
    valid_frame_mask = all(valid_matrix, 2);
    info = struct('B', B, 'M', M, 'delay_symbols', delay_symbols, ...
        'delay_frames', delay_frames, 'commutator_direction', commutator_direction, ...
        'commutator_phase', mod(commutator_phase, B), ...
        'input_frames', num_frames, ...
        'fully_valid_frames', nnz(valid_frame_mask));
end
