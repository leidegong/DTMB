function [body_fd_shifted, y_symbols, z_symbols] = DTMB_Frequency_Interleave(data_symbols, system_symbols)
% DTMB_Frequency_Interleave 将36个系统信息符号和3744个数据符号映射到Y载波。

    data_symbols = data_symbols(:);
    system_symbols = system_symbols(:);
    if numel(data_symbols) ~= 3744 || numel(system_symbols) ~= 36
        error('频域交织输入应为 3744 个数据符号和 36 个系统信息符号。');
    end
    [z_to_y, tps_idx] = DTMB_Frequency_Permutation();
    data_mask = true(3780, 1);
    data_mask(tps_idx) = false;
    z_symbols = zeros(3780, 1, 'like', data_symbols + system_symbols(1));
    z_symbols(tps_idx) = system_symbols;
    z_symbols(data_mask) = data_symbols;
    y_symbols = zeros(3780, 1, 'like', z_symbols);
    y_symbols(z_to_y) = z_symbols;
    body_fd_shifted = fftshift(y_symbols);
end
