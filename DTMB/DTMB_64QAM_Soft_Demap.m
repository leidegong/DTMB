function [llr, hard_bits, metrics] = DTMB_64QAM_Soft_Demap(symbols, noise_variance)
% DTMB_64QAM_Soft_Demap 最大对数软解映射；正LLR表示比特0更可信。

    symbols = symbols(:);
    if nargin < 2 || isempty(noise_variance)
        sliced = Slice_64QAM(symbols);
        noise_variance = median(abs(symbols - sliced).^2, 'omitnan');
    end
    noise_variance = double(noise_variance(:));
    if isscalar(noise_variance)
        noise_variance = repmat(noise_variance, numel(symbols), 1);
    elseif numel(noise_variance) ~= numel(symbols)
        error('噪声方差必须是标量或与输入符号等长的向量。');
    end
    valid_noise = isfinite(noise_variance) & noise_variance > 0;
    if any(valid_noise)
        fallback_noise = median(noise_variance(valid_noise));
    else
        sliced = Slice_64QAM(symbols);
        fallback_noise = median(abs(symbols - sliced).^2, 'omitnan');
    end
    noise_variance(~valid_noise) = fallback_noise;
    noise_variance = max(noise_variance, 1e-4);
    % Square 64QAM is separable.  For an I-axis bit, minimization over Q is
    % identical for the bit-0 and bit-1 hypotheses and cancels in d1-d0
    % (and vice versa).  Two N-by-8 real distance matrices therefore give
    % exactly the same max-log LLRs as six N-by-32 complex matrices.
    persistent levels axis_labels bit_to_axis_column
    if isempty(levels)
        levels = (-7:2:7)/sqrt(42);
        axis_labels = logical([0 0 0;0 0 1;0 1 1;0 1 0; ...
            1 1 0;1 1 1;1 0 1;1 0 0]);
        % DTMB_64QAM_Map maps serial bits 1:3 through label columns 3:1.
        bit_to_axis_column = [3 2 1];
    end
    distance_i = (real(symbols)-levels).^2;
    distance_q = (imag(symbols)-levels).^2;
    llr_matrix = zeros(numel(symbols),6);
    for axis_bit = 1:3
        label_column = bit_to_axis_column(axis_bit);
        use0 = ~axis_labels(:,label_column).';
        use1 = axis_labels(:,label_column).';
        llr_matrix(:,axis_bit) = ...
            (min(distance_i(:,use1),[],2)-min(distance_i(:,use0),[],2)) ...
            ./ noise_variance;
        llr_matrix(:,axis_bit+3) = ...
            (min(distance_q(:,use1),[],2)-min(distance_q(:,use0),[],2)) ...
            ./ noise_variance;
    end
    llr = reshape(llr_matrix.', [], 1);
    hard_bits = llr < 0;
    metrics = struct('noise_variance', noise_variance, ...
        'noise_variance_median', median(noise_variance), ...
        'evm_rms', sqrt(mean(abs(symbols - Slice_64QAM(symbols)).^2, 'omitnan')));
end
