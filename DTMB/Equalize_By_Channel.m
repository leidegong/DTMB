function req = Equalize_By_Channel(rcs, h_hat, lambda)
% Equalize_By_Channel
% 功能：
%   用估计信道 h_hat 对参考信号 rcs 做简单频域均衡
%
% 输入：
%   rcs    - 当前帧参考信号（列向量）
%   h_hat  - 当前估计信道
%   lambda - 正则项
%
% 输出：
%   req    - 均衡后的参考信号

    if nargin < 3
        lambda = 1e-3;
    end

    rcs = rcs(:);
    h_hat = h_hat(:);

    N = length(rcs);
    H = fft(h_hat, N);
    R = fft(rcs, N);

    Req = R .* conj(H) ./ (abs(H).^2 + lambda);
    req = ifft(Req);
end