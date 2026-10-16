function sym_hat = Slice_64QAM(sym)
%SLICE_64QAM Quantize to the nearest normalized square-64QAM point.
%   The arithmetic quantizer is equivalent to an eight-level distance
%   search except at exact decision boundaries, where it selects the lower
%   level.  It avoids allocating an N-by-8 distance matrix and is suitable
%   for frame-rate streaming code.
%
% 输入:
%   sym      待判决复符号
%
% 输出:
%   sym_hat  判决后的归一化 64QAM 符号
%
    normFactor = sqrt(42);
    scaled = sym * normFactor;
    finiteMask = isfinite(real(scaled)) & isfinite(imag(scaled));

    % Boundaries are the even integers.  2*ceil(x/2)-1 maps every interval
    % to its nearest odd level; clipping supplies the two outer levels.
    levelI = min(max(2*ceil(real(scaled)/2)-1,-7),7);
    levelQ = min(max(2*ceil(imag(scaled)/2)-1,-7),7);
    levelI(~finiteMask) = -1;
    levelQ(~finiteMask) = -1;
    sym_hat = complex(levelI,levelQ) / normFactor;
end
