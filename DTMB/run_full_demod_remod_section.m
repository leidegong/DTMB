function result = run_full_demod_remod_section(equalized_reference_frames, phase_index_list, cfg)
% run_full_demod_remod_section DTMB完整解调、FEC译码和再调制重构。
%
% 始终生成判决导向重构；当输入帧数大于两倍时域交织时延时，
% 中心帧再由“解交织-LDPC/BCH-重编码-重交织”的结果替换。

    if nargin < 3
        cfg = struct();
    end
    matrix_input = ~iscell(equalized_reference_frames);
    phase_index_list = phase_index_list(:);
    if matrix_input
        num_input_frames = size(equalized_reference_frames,1);
    else
        num_input_frames = numel(equalized_reference_frames);
    end
    num_frames = min(num_input_frames,numel(phase_index_list));
    if num_frames == 0
        error('没有可供解调重构的均衡帧。');
    end
    if matrix_input
        if size(equalized_reference_frames,2) < 4725
            error('均衡帧矩阵的列数为 %d，小于 4725。', ...
                size(equalized_reference_frames,2));
        end
        equalized_matrix = double(equalized_reference_frames(1:num_frames,1:4725));
    else
        equalized_reference_frames = equalized_reference_frames(1:num_frames);
        equalized_matrix = complex(zeros(num_frames,4725));
        for frame_idx = 1:num_frames
            frame = equalized_reference_frames{frame_idx}(:);
            if numel(frame) < 4725
                error('第 %d 个均衡帧长度为 %d，小于 4725。', ...
                    frame_idx,numel(frame));
            end
            equalized_matrix(frame_idx,:) = double(frame(1:4725)).';
        end
    end
    phase_index_list = phase_index_list(1:num_frames);
    frame_erasure_mask = get_cfg(cfg, 'fec_frame_erasure_mask', false(num_frames, 1));
    frame_erasure_mask = logical(frame_erasure_mask(:));
    if numel(frame_erasure_mask) < num_frames
        frame_erasure_mask(end+1:num_frames, 1) = false;
    else
        frame_erasure_mask = frame_erasure_mask(1:num_frames);
    end
    nfft = 3780;
    frame_length = 4725;
    stage_timing = struct();

    fprintf('  [重构 1/5] 逐帧FFT、频域解交织和系统信息提取。\n');
    t_stage = tic;
    [z_to_y,tps_idx] = DTMB_Frequency_Permutation();
    data_mask = true(nfft,1);
    data_mask(tps_idx) = false;
    % DTMB_Frequency_Deinterleave applies ifftshift immediately after the
    % diagnostic fftshift.  Operating in native FFT order removes both
    % permutations and batches all frames in one FFT call.
    body_y = fft(equalized_matrix(:,946:frame_length),nfft,2);
    body_z = body_y(:,z_to_y);
    data_rx = body_z(:,data_mask);
    system_rx = body_z(:,tps_idx);
    clear body_y body_z
    stage_timing.frequency_extract_s = toc(t_stage);

    system_rx(frame_erasure_mask, :) = complex(NaN, NaN);
    system_info = DTMB_Decode_System_Info(system_rx);
    fprintf('  系统信息判决：64QAM，码率 %.1f，时域交织模式 %d (M=%d)，相关 %.4f，置信比 %.3f。\n', ...
        system_info.rate, system_info.interleaver_mode, system_info.M, ...
        system_info.best_score, system_info.confidence_ratio);

    fprintf('  [重构 2/5] TPS定相、64QAM自动增益和软信息计算。\n');
    t_stage = tic;
    data_normalized = complex(zeros(size(data_rx)));
    hard_interleaved = complex(zeros(size(data_rx)));
    frame_phase = zeros(num_frames, 1);
    frame_tps_phase_raw = zeros(num_frames, 1);
    frame_qam_quarter_turns = zeros(num_frames, 1);
    frame_scale = zeros(num_frames, 1);
    noise_variance = zeros(num_frames, 1);
    decision_evm = zeros(num_frames, 1);
    track_phase_continuity = get_cfg(cfg, ...
        'fec_track_frame_phase_continuity', false);
    cpe_method = lower(string(get_cfg(cfg,'fec_cpe_method', ...
        'decision_directed')));
    previous_total_phase = NaN;
    previous_residual_phase = NaN;
    cpe_reacquired = false(num_frames,1);
    for frame_idx = 1:num_frames
        if frame_erasure_mask(frame_idx)
            data_normalized(frame_idx, :) = 0;
            hard_interleaved(frame_idx, :) = Slice_64QAM(complex(zeros(3744, 1))).';
            frame_phase(frame_idx) = 0;
            frame_scale(frame_idx) = 1;
            noise_variance(frame_idx) = 1e9;
            decision_evm(frame_idx) = NaN;
            previous_residual_phase = NaN;
            continue;
        end
        initial_phase_raw = angle(system_info.per_frame_complex_gain(frame_idx));
        initial_phase = initial_phase_raw;
        % TPS complex gain can have a +/-1 ambiguity.  Resolve it against
        % the preceding frame before estimating the four-fold-symmetric
        % 64-QAM residual phase.
        if track_phase_continuity && isfinite(previous_total_phase)
            initial_phase = closest_periodic_phase(initial_phase_raw, ...
                previous_total_phase, pi);
        end
        phase_corrected = data_rx(frame_idx, :).'*exp(-1j*initial_phase);
        scale = rms(phase_corrected(isfinite(phase_corrected))) + eps;
        initial_normalized = phase_corrected / scale;
        if cpe_method == "decision_directed"
            [residual_phase,cpe_info,fast_sliced] = ...
                Estimate_CPE_64QAM_Low_Latency(initial_normalized, ...
                previous_residual_phase,cfg);
            cpe_reacquired(frame_idx) = cpe_info.reacquired;
        elseif cpe_method == "grid"
            residual_phase = Estimate_CPE_64QAM(initial_normalized);
        else
            error('未知 fec_cpe_method：%s。',cpe_method);
        end
        total_phase = initial_phase + residual_phase;
        if track_phase_continuity && isfinite(previous_total_phase)
            continuous_phase = closest_periodic_phase(total_phase, ...
                previous_total_phase, pi/2);
            frame_qam_quarter_turns(frame_idx) = round( ...
                (continuous_phase-total_phase)/(pi/2));
            total_phase = continuous_phase;
        end
        if cpe_method == "decision_directed" && ...
                frame_qam_quarter_turns(frame_idx) == 0
            % Phase rotation preserves RMS. Reuse the hard decisions already
            % produced by the CPE tracker instead of slicing twice.
            normalized = initial_normalized*exp(-1j*residual_phase);
            sliced = fast_sliced;
        else
            phase_corrected = data_rx(frame_idx,:).'*exp(-1j*total_phase);
            scale = rms(phase_corrected(isfinite(phase_corrected)))+eps;
            normalized = phase_corrected/scale;
            sliced = Slice_64QAM(normalized);
        end
        residual = normalized - sliced;
        nv = median(abs(residual).^2, 'omitnan');
        data_normalized(frame_idx, :) = normalized.';
        hard_interleaved(frame_idx, :) = sliced.';
        frame_phase(frame_idx) = total_phase;
        frame_tps_phase_raw(frame_idx) = initial_phase_raw;
        frame_scale(frame_idx) = scale;
        noise_variance(frame_idx) = max(nv, 1e-4);
        decision_evm(frame_idx) = sqrt(mean(abs(residual).^2, 'omitnan'));
        previous_total_phase = total_phase;
        previous_residual_phase = residual_phase;
    end
    stage_timing.phase_normalize_s = toc(t_stage);

    fprintf('  [重构 3/5] 生成全帧判决导向基线。\n');
    t_stage = tic;
    candidates = DTMB_System_Info_Candidates();
    system_clean = candidates(system_info.index).symbols;
    decision_canonical = build_frames(hard_interleaved, system_clean, phase_index_list, cfg);
    reconstruction_canonical = decision_canonical;
    reconstruction_source = ones(num_frames, 1, 'uint8');
    stage_timing.decision_rebuild_s = toc(t_stage);

    delay_frames = 52 * 51 * system_info.M / 3744;
    fec_available = num_frames > 2 * delay_frames;
    fec_valid_mask = false(num_frames, 1);
    fec_stats = struct();
    decoded_payload_bits = false(0, 0);
    reencoded_blocks = false(7488, 0);
    llr_blocks = zeros(7488, 0, 'single');
    fec_alignment_info = struct();

    if fec_available
        t_stage = tic;
        fprintf('  [重构 4/5] 时域解交织、64QAM软解映射、LDPC/BCH译码与重编码。\n');
        t_fec_part = tic;
        [deinterleaved_data, fec_alignment_info] = DTMB_Select_FEC_Alignment( ...
            data_normalized, system_info.M, system_info.rate_index, cfg);
        stage_timing.fec_alignment_deinterleave_s = toc(t_fec_part);
        time_deinterleave_info = fec_alignment_info.time_deinterleave_info;
        num_pre_frames = size(deinterleaved_data, 1);
        noise_model = get_cfg(cfg, 'ldpc_noise_model', 'frame-window');
        t_fec_part = tic;
        if strcmp(noise_model, 'deinterleaved')
            noise_input = repmat(single(noise_variance), 1, 3744);
            deinterleaved_noise = DTMB_Time_Deinterleave(noise_input, ...
                system_info.M, fec_alignment_info.direction, ...
                fec_alignment_info.commutator_phase);
        else
            deinterleaved_noise = [];
        end
        stage_timing.fec_noise_prepare_s = toc(t_fec_part);
        llr_scale = get_cfg(cfg, 'ldpc_llr_scale', 1);
        llr_blocks = zeros(7488,num_pre_frames*3,'single');
        demap_batch_frames = max(1,round(get_cfg(cfg, ...
            'fec_demap_batch_frames',16)));
        if ~strcmp(noise_model,'deinterleaved')
            local_noise_per_frame = zeros(num_pre_frames,1);
            for pre_idx = 1:num_pre_frames
                local_noise_per_frame(pre_idx) = median(noise_variance( ...
                    pre_idx:min(num_frames,pre_idx+delay_frames)),'omitnan');
            end
        end
        t_fec_part = tic;
        for batch_first = 1:demap_batch_frames:num_pre_frames
            take = batch_first:min(num_pre_frames, ...
                batch_first+demap_batch_frames-1);
            symbols = reshape(deinterleaved_data(take,:).',[],1);
            if strcmp(noise_model,'deinterleaved')
                local_noise = reshape(double(deinterleaved_noise(take,:).'),[],1);
            else
                local_noise = repelem(local_noise_per_frame(take),3744);
            end
            llr = DTMB_64QAM_Soft_Demap(symbols,local_noise);
            block_columns = (take(1)-1)*3+(1:3*numel(take));
            llr_blocks(:,block_columns) = single(llr_scale * ...
                reshape(llr,7488,3*numel(take)));
        end
        stage_timing.fec_soft_demap_s = toc(t_fec_part);
        max_iterations = get_cfg(cfg, 'ldpc_max_iterations', 30);
        batch_size = get_cfg(cfg, 'ldpc_batch_size', 16);
        decoder_options = struct( ...
            'algorithm', get_cfg(cfg, 'ldpc_algorithm', 'bp'), ...
            'min_sum_scaling_factor', get_cfg(cfg, ...
                'ldpc_min_sum_scaling_factor', 0.75), ...
            'min_sum_offset', get_cfg(cfg, 'ldpc_min_sum_offset', 0.5));
        t_fec_part = tic;
        [reencoded_blocks, decoded_payload_bits, fec_stats] = ...
            DTMB_LDPC_Decode_Reencode(llr_blocks, system_info.rate_index, ...
                max_iterations, batch_size, decoder_options);
        stage_timing.fec_ldpc_bch_reencode_s = toc(t_fec_part);
        fec_stats.noise_model = noise_model;
        fec_stats.llr_scale = llr_scale;
        t_fec_part = tic;
        mapped_symbols = DTMB_64QAM_Map(reencoded_blocks(:));
        clean_pre_frames = reshape(mapped_symbols,3744,num_pre_frames).';
        block_ok = reshape(fec_stats.parity_ok,3,num_pre_frames).';
        pre_symbol_fec_ok = repelem(block_ok,1,1248);
        stage_timing.fec_qam_map_s = toc(t_fec_part);
        % 未通过LDPC校验的码块不能作为“纯化”结果参与重构。
        clean_pre_frames(~pre_symbol_fec_ok) = complex(NaN, NaN);
        t_fec_part = tic;
        [fec_interleaved_short, local_valid, time_interleave_info] = ...
            DTMB_Time_Interleave(clean_pre_frames, system_info.M, ...
                fec_alignment_info.direction, fec_alignment_info.commutator_phase);
        local_valid = local_valid & all(isfinite(fec_interleaved_short), 2);
        stage_timing.fec_time_interleave_s = toc(t_fec_part);
        fec_interleaved = complex(nan(num_frames, 3744));
        fec_interleaved(1:num_pre_frames, :) = fec_interleaved_short;
        fec_valid_mask(1:num_pre_frames) = local_valid;
        t_fec_part = tic;
        fec_canonical = build_frames(fec_interleaved(fec_valid_mask, :), system_clean, ...
            phase_index_list(fec_valid_mask), cfg);
        reconstruction_canonical(fec_valid_mask, :) = fec_canonical;
        reconstruction_source(fec_valid_mask) = 2;
        stage_timing.fec_valid_frame_rebuild_s = toc(t_fec_part);
        fec_stats.time_interleave_info = time_interleave_info;
        fec_stats.valid_frame_indices = find(fec_valid_mask);
        fec_stats.pre_symbol_fec_ok = pre_symbol_fec_ok;
        fec_stats.num_successful_blocks = nnz(fec_stats.parity_ok);
        if any(fec_valid_mask)
            fprintf('  FEC完整重构有效帧：%d 到 %d，共 %d 帧；LDPC校验通过率 %.2f%%。\n', ...
                find(fec_valid_mask, 1, 'first'), find(fec_valid_mask, 1, 'last'), ...
                nnz(fec_valid_mask), 100*fec_stats.parity_success_rate);
        else
            fprintf('  LDPC校验通过率 %.2f%%（%d/%d 码块），没有可安全采用的完整FEC重构帧；保留判决导向结果。\n', ...
                100*fec_stats.parity_success_rate, fec_stats.num_successful_blocks, fec_stats.num_blocks);
        end
        stage_timing.fec_decode_reencode_s = toc(t_stage);
    else
        [~, time_deinterleave_info] = DTMB_Time_Deinterleave(data_normalized, system_info.M);
        fprintf('  [重构 4/5] 输入仅 %d 帧，完整FEC重构需要至少 %d 帧；本次保留判决导向结果。\n', ...
            num_frames, 2*delay_frames + 1);
        stage_timing.fec_decode_reencode_s = 0;
        stage_timing.fec_alignment_deinterleave_s = 0;
        stage_timing.fec_noise_prepare_s = 0;
        stage_timing.fec_soft_demap_s = 0;
        stage_timing.fec_ldpc_bch_reencode_s = 0;
        stage_timing.fec_qam_map_s = 0;
        stage_timing.fec_time_interleave_s = 0;
        stage_timing.fec_valid_frame_rebuild_s = 0;
    end

    fprintf('  [重构 5/5] 帧头拼接、复增益匹配和质量统计。\n');
    t_stage = tic;
    match_gain = sum(conj(reconstruction_canonical).*equalized_matrix,2) ./ ...
        (sum(abs(reconstruction_canonical).^2,2)+eps);
    reconstructed_frames = single(reconstruction_canonical.*match_gain);
    stage_timing.gain_match_s = toc(t_stage);

    result = struct();
    result.reconstructed_frames = reconstructed_frames;
    result.canonical_reconstructed_frames = single(reconstruction_canonical);
    result.reconstruction_source = reconstruction_source;
    result.source_labels = {'decision_directed', 'fec_decode_reencode'};
    result.fec_valid_frame_mask = fec_valid_mask;
    result.frame_erasure_mask = frame_erasure_mask;
    result.fec_available = fec_available;
    result.system_info = system_info;
    result.time_deinterleave_info = time_deinterleave_info;
    result.fec_alignment_info = fec_alignment_info;
    result.fec_stats = fec_stats;
    result.decoded_payload_bits = decoded_payload_bits;
    result.reencoded_ldpc_blocks = reencoded_blocks;
    if get_cfg(cfg, 'save_ldpc_llr_blocks', false)
        result.ldpc_llr_blocks = llr_blocks;
    end
    result.frame_phase_rad = frame_phase;
    result.frame_tps_phase_raw_rad = frame_tps_phase_raw;
    result.frame_qam_quarter_turns = frame_qam_quarter_turns;
    result.phase_continuity_tracking_enabled = track_phase_continuity;
    result.cpe_method = cpe_method;
    result.cpe_reacquired = cpe_reacquired;
    result.frame_qam_scale = frame_scale;
    result.noise_variance = noise_variance;
    result.decision_evm_rms = decision_evm;
    result.complex_match_gain = match_gain;
    result.start_sample = get_cfg(cfg, 'start_sample', NaN);
    result.num_frames = num_frames;
    stage_timing.total_measured_s = stage_timing.frequency_extract_s + ...
        stage_timing.phase_normalize_s + stage_timing.decision_rebuild_s + ...
        stage_timing.fec_decode_reencode_s + stage_timing.gain_match_s;
    stage_timing.per_input_frame_ms = ...
        1e3*stage_timing.total_measured_s/num_frames;
    result.timing = stage_timing;
    result.description = ['DTMB解调-再调制重构；source=2为完整FEC重构，' ...
        'source=1为交织边界或帧数不足时的判决导向重构。'];

    output_packet = get_cfg(cfg, 'full_recon_output_packet', 'DTMB_Full_Demod_Remod_Reconstruction.mat');
    save(output_packet, 'result', '-v7.3');
    fprintf('  完整解调—再调制数据包已保存：%s\n', output_packet);
end

function phase_out = closest_periodic_phase(phase_in, reference, period)
% Choose the periodic representative nearest to the preceding frame.
    phase_out = phase_in + period*round((reference-phase_in)/period);
end

function frames = build_frames(interleaved_data, system_symbols, phase_index_list, cfg)
    num_frames = size(interleaved_data, 1);
    pn_to_body_power_ratio = get_cfg(cfg, 'pn_to_body_power_ratio', 2);
    body_ifft_gain = sqrt(3780 / pn_to_body_power_ratio);
    [z_to_y,tps_idx] = DTMB_Frequency_Permutation();
    data_mask = true(3780,1);
    data_mask(tps_idx) = false;
    z = complex(zeros(num_frames,3780));
    z(:,tps_idx) = repmat(system_symbols(:).',num_frames,1);
    z(:,data_mask) = interleaved_data;
    y = complex(zeros(num_frames,3780));
    y(:,z_to_y) = z;
    body_td = body_ifft_gain*ifft(y,3780,2);
    pn = double(PN_gen(phase_index_list(:)));
    frames = [pn,body_td];
end

function value = get_cfg(cfg, field_name, default_value)
    if isstruct(cfg) && isfield(cfg, field_name) && ~isempty(cfg.(field_name))
        value = cfg.(field_name);
    else
        value = default_value;
    end
end
