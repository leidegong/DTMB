function [reencoded_bits, payload_bits, stats] = DTMB_LDPC_Decode_Reencode(received_llr, rate_index, max_iterations, batch_size, decoder_options)
% DTMB_LDPC_Decode_Reencode 对7488比特穿孔码字译码、BCH校正并重新编码。
% received_llr 每列一个传输码字，正值表示0，负值表示1。

    if nargin < 3 || isempty(max_iterations)
        max_iterations = 30;
    end
    if nargin < 4 || isempty(batch_size)
        batch_size = 16;
    end
    if nargin < 5 || isempty(decoder_options)
        decoder_options = struct();
    end
    if size(received_llr, 1) ~= 7488
        error('LDPC 输入必须每列为 7488 个传输LLR。');
    end
    [H_standard, meta] = DTMB_LDPC_Parity_Check_Matrix(rate_index);
    parity_length = meta.parity_length;
    information_length = meta.information_length;
    H_matlab = H_standard(:, [parity_length+1:end, 1:parity_length]);
    encoder_cfg = ldpcEncoderConfig(H_matlab);
    decoder_cfg = ldpcDecoderConfig(H_matlab);
    algorithm = get_option(decoder_options, 'algorithm', 'bp');
    decoder_cfg.Algorithm = algorithm;
    min_sum_scaling_factor = get_option(decoder_options, ...
        'min_sum_scaling_factor', 0.75);
    min_sum_offset = get_option(decoder_options, 'min_sum_offset', 0.5);

    num_blocks = size(received_llr, 2);
    reencoded_bits = false(7488, num_blocks);
    bch_blocks_per_ldpc = information_length / 762;
    payload_bits = false(752 * bch_blocks_per_ldpc, num_blocks);
    iterations = zeros(num_blocks, 1);
    parity_ok = false(num_blocks, 1);
    bch_fail_before = zeros(num_blocks, 1);
    bch_corrected = zeros(num_blocks, 1);
    bch_uncorrectable = zeros(num_blocks, 1);
    changed_bits = zeros(num_blocks, 1);
    ldpc_information_all = false(information_length,num_blocks);

    decode_args = {'OutputFormat','whole','DecisionType','hard', ...
        'Termination','early'};
    if strcmp(algorithm,'norm-min-sum')
        decode_args = [decode_args, ...
            {'MinSumScalingFactor',min_sum_scaling_factor}];
    elseif strcmp(algorithm,'offset-min-sum')
        decode_args = [decode_args,{'MinSumOffset',min_sum_offset}];
    end

    for batch_start = 1:batch_size:num_blocks
        take = batch_start:min(batch_start + batch_size - 1, num_blocks);
        batch_count = numel(take);
        llr_standard = [zeros(5, batch_count); double(received_llr(:, take))];
        llr_matlab = [llr_standard(parity_length+1:end, :); llr_standard(1:parity_length, :)];
        [whole_matlab, batch_iterations, final_checks] = ldpcDecode( ...
            llr_matlab, decoder_cfg, max_iterations, decode_args{:});
        iterations(take) = batch_iterations(:);
        parity_ok(take) = ~any(final_checks, 1).';
        ldpc_information_all(:,take) = logical( ...
            whole_matlab(1:information_length,:));
    end

    % BCH arithmetic is vectorized over columns.  Running it once for the
    % complete burst avoids repeating 752 polynomial steps per LDPC batch.
    bch_matrix = reshape(ldpc_information_all,762,[]);
    [bch_corrected_words,bch_payload,bch_stats] = ...
        DTMB_BCH762_Decode_Encode(bch_matrix);
    corrected_information_all = reshape( ...
        bch_corrected_words,information_length,num_blocks);
    payload_bits(:,:) = reshape( ...
        bch_payload,752*bch_blocks_per_ldpc,num_blocks);
    bch_fail_before(:) = sum(reshape( ...
        bch_stats.syndrome_fail_before_mask, ...
        bch_blocks_per_ldpc,num_blocks),1).';
    bch_corrected(:) = sum(reshape( ...
        bch_stats.single_error_corrected_mask, ...
        bch_blocks_per_ldpc,num_blocks),1).';
    bch_uncorrectable(:) = sum(reshape( ...
        bch_stats.uncorrectable_mask, ...
        bch_blocks_per_ldpc,num_blocks),1).';

    for batch_start = 1:batch_size:num_blocks
        take = batch_start:min(batch_start+batch_size-1,num_blocks);
        codeword_matlab = logical(ldpcEncode( ...
            corrected_information_all(:,take),encoder_cfg));
        codeword_standard = [codeword_matlab(information_length+1:end, :); ...
                             codeword_matlab(1:information_length, :)];
        reencoded_bits(:, take) = codeword_standard(6:end, :);
        changed_bits(take) = sum(reencoded_bits(:, take) ~= (received_llr(:, take) < 0), 1).';
    end

    stats = struct();
    stats.rate_index = rate_index;
    stats.rate = meta.rate;
    stats.num_blocks = num_blocks;
    stats.max_iterations = max_iterations;
    stats.batch_size = batch_size;
    stats.algorithm = algorithm;
    stats.min_sum_scaling_factor = min_sum_scaling_factor;
    stats.min_sum_offset = min_sum_offset;
    stats.iterations = iterations;
    stats.parity_ok = parity_ok;
    stats.parity_success_rate = mean(parity_ok);
    stats.bch_syndrome_fail_before = bch_fail_before;
    stats.bch_single_error_corrected = bch_corrected;
    stats.bch_uncorrectable = bch_uncorrectable;
    stats.changed_bits = changed_bits;
    stats.changed_bit_rate = sum(changed_bits) / numel(reencoded_bits);
end

function value = get_option(options, field_name, default_value)
    if isstruct(options) && isfield(options, field_name) && ...
            ~isempty(options.(field_name))
        value = options.(field_name);
    else
        value = default_value;
    end
end
