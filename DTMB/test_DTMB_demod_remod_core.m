function test_DTMB_demod_remod_core()
% test_DTMB_demod_remod_core 核心映射、交织和FEC的确定性自检。

    rng(20260815);
    bits = logical(randi([0 1], 600, 1));
    symbols = DTMB_64QAM_Map(bits);
    [~, hard_bits] = DTMB_64QAM_Soft_Demap(symbols, 1e-4);
    assert(isequal(bits, hard_bits), '64QAM 映射往返失败。');

    candidates = DTMB_System_Info_Candidates();
    data = randn(3744, 1) + 1j*randn(3744, 1);
    body_fd = DTMB_Frequency_Interleave(data, candidates(4).symbols);
    [data_back, system_back] = DTMB_Frequency_Deinterleave(body_fd);
    assert(max(abs(data - data_back)) < 1e-12, '频域交织往返失败。');
    decoded_info = DTMB_Decode_System_Info(repmat(system_back.', 4, 1));
    assert(decoded_info.index == 4, '系统信息判决失败。');

    x = reshape(complex(1:(400*3744)), 400, 3744);
    for direction = [1, -1]
        for phase = [0, 17]
            [y, ~] = DTMB_Time_Interleave(x, 240, direction, phase);
            [z, ~] = DTMB_Time_Deinterleave(y, 240, direction, phase);
            valid = all(isfinite(z), 2);
            assert(max(abs(z(valid, :) - x(1:nnz(valid), :)), [], 'all') == 0, ...
                '时域交织往返失败。');
        end
    end

    for rate_index = 1:3
        [H, meta] = DTMB_LDPC_Parity_Check_Matrix(rate_index);
        bch_count = meta.information_length / 762;
        payload = logical(randi([0 1], 752, bch_count));
        bch_words = DTMB_BCH762_Encode(payload);
        p = meta.parity_length;
        encoder_cfg = ldpcEncoderConfig(H(:, [p+1:end, 1:p]));
        matlab_codeword = logical(ldpcEncode(bch_words(:), encoder_cfg));
        standard_codeword = [matlab_codeword(meta.information_length+1:end); ...
                             matlab_codeword(1:meta.information_length)];
        transmitted = standard_codeword(6:end);
        llr = (1 - 2*double(transmitted))*20;
        [reencoded, payload_back, stats] = ...
            DTMB_LDPC_Decode_Reencode(llr, rate_index, 10, 1);
        assert(isequal(transmitted, reencoded), 'LDPC 再编码不一致。');
        assert(isequal(payload(:), payload_back), 'BCH 载荷不一致。');
        assert(stats.parity_success_rate == 1, '无噪声LDPC校验失败。');
    end
    fprintf('DTMB 解调—再调制核心自检全部通过。\n');
end
