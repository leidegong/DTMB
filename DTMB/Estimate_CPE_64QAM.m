function phi_best = Estimate_CPE_64QAM(sym)
% Estimate_CPE_64QAM
% 基于 64QAM 判决辅助的公共相位误差估计
%
% 输入:
%   sym       当前帧频域符号，建议已做 rms 归一化
%
% 输出:
%   phi_best  估计得到的公共相位误差，单位 rad
%
% 使用方法:
%   phi = Estimate_CPE_64QAM(sym);
%   sym_corr = sym * exp(-1j * phi);

    sym = sym(:);

    % 去掉 NaN / Inf
    valid = isfinite(real(sym)) & isfinite(imag(sym));
    sym = sym(valid);

    if numel(sym) < 100
        phi_best = 0;
        return;
    end

    % 剔除幅度过小和过大的异常点，避免边缘频点或坏点影响估计
    amp = abs(sym);
    low_th  = prctile(amp, 5);
    high_th = prctile(amp, 95);

    mask = amp > low_th & amp < high_th;
    sym_use = sym(mask);

    if numel(sym_use) < 100
        phi_best = 0;
        return;
    end

    % 64QAM 具有 pi/2 对称性，先在 [-pi/4, pi/4] 内搜索
    phi_grid = linspace(-pi/4, pi/4, 91);
    % 一次性评估全部候选角度。与逐角度循环完全等价，但全量帧测试
    % 的执行时间可从数十分钟降低到分钟以内。
    sym_rot = sym_use .* exp(-1j * phi_grid);
    sym_dec = Slice_64QAM(sym_rot);
    metric = mean(abs(sym_rot - sym_dec).^2, 1);

    [~, idx_best] = min(metric);
    phi_best = phi_grid(idx_best);
end
