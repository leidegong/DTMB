function [output_frames, info] = DTMB_Time_Deinterleave(input_frames, M, commutator_direction, commutator_phase)
% DTMB_Time_Deinterleave 52支路卷积解交织并去掉统一暖机时延。

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
        error('时域解交织要求每帧 3744 个数据符号。');
    end
    B = 52;
    delay_symbols = B * (B - 1) * M;
    delay_frames = delay_symbols / symbols_per_frame;
    if delay_frames ~= round(delay_frames)
        error('交织时延不是整数帧。');
    end
    stream = reshape(input_frames.', [], 1);
    delayed = nan(size(stream), 'like', stream);
    t0 = (delay_symbols:(numel(stream)-1)).';
    branch = mod(commutator_direction * t0 + commutator_phase, B);
    source0 = t0 - (B - 1 - branch) * M * B;
    delayed(t0 + 1) = stream(source0 + 1);
    cropped = delayed(delay_symbols + 1:end);
    output_frames = reshape(cropped, symbols_per_frame, []).';
    info = struct('B', B, 'M', M, 'delay_symbols', delay_symbols, ...
        'delay_frames', round(delay_frames), 'commutator_direction', commutator_direction, ...
        'commutator_phase', mod(commutator_phase, B), ...
        'input_frames', num_frames, ...
        'output_frames', size(output_frames, 1));
end
