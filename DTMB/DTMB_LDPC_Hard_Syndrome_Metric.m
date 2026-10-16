function [metric, detail] = DTMB_LDPC_Hard_Syndrome_Metric(transmitted_bits, rate_index)
% DTMB_LDPC_Hard_Syndrome_Metric 对5个穿孔位穷举后的最小校验综合比例。
% transmitted_bits 每列为一个7488比特硬判决码字。

    transmitted_bits = logical(transmitted_bits);
    if size(transmitted_bits, 1) ~= 7488
        error('综合评价输入必须每列为 7488 比特。');
    end
    [H, meta] = DTMB_LDPC_Parity_Check_Matrix(rate_index);
    puncture_patterns = logical(dec2bin(0:31, 5) - '0').';
    puncture_syndromes = mod(double(H(:, 1:5)) * double(puncture_patterns), 2);
    weights = zeros(size(transmitted_bits, 2), 1);
    best_pattern = zeros(size(transmitted_bits, 2), 1);
    for block_idx = 1:size(transmitted_bits, 2)
        received_syndrome = mod(double(H(:, 6:end)) * double(transmitted_bits(:, block_idx)), 2);
        candidate_syndromes = mod(received_syndrome + puncture_syndromes, 2);
        [weights(block_idx), best_pattern(block_idx)] = min(sum(candidate_syndromes, 1));
    end
    metric = mean(weights) / meta.parity_length;
    detail = struct('unsatisfied_checks', weights, 'best_puncture_pattern', best_pattern, ...
        'num_checks', meta.parity_length, 'num_blocks', size(transmitted_bits, 2));
end
