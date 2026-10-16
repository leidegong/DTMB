function [data_symbols, system_symbols, z_symbols] = DTMB_Frequency_Deinterleave(body_fd_shifted)
% DTMB_Frequency_Deinterleave 将FFT移位顺序的Y载波恢复为Z/X顺序。

    body_fd_shifted = body_fd_shifted(:);
    if numel(body_fd_shifted) ~= 3780
        error('频域解交织要求 3780 个载波，实际为 %d。', numel(body_fd_shifted));
    end
    [z_to_y, tps_idx] = DTMB_Frequency_Permutation();
    % 工程FFT输出先做fftshift用于观察，标准Y序列对应还原后的FFT索引。
    y_symbols = ifftshift(body_fd_shifted);
    z_symbols = y_symbols(z_to_y);
    data_mask = true(3780, 1);
    data_mask(tps_idx) = false;
    system_symbols = z_symbols(tps_idx);
    data_symbols = z_symbols(data_mask);
end
