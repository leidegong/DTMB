function [deinterleaved_data, info] = DTMB_Select_FEC_Alignment(data_normalized, M, rate_index, cfg)
% DTMB_Select_FEC_Alignment 用LDPC硬综合自动消除换向、共轭和象限歧义。

    if nargin < 4
        cfg = struct();
    end
    search_frames = get_cfg(cfg, 'fec_alignment_search_frames', 12);
    if get_cfg(cfg, 'fec_search_all_commutator_phases', false)
        directions = [1, -1];
        phases = 0:51;
    else
        directions = get_cfg(cfg, 'fec_commutator_direction', 1);
        phases = get_cfg(cfg, 'fec_commutator_phase', 0);
    end
    conjugations = [false, true];
    quadrants = 0:3;
    branch_records = repmat(struct('direction', 0, 'commutator_phase', 0, ...
        'syndrome_metric', Inf), numel(directions)*numel(phases), 1);
    branch_idx = 0;
    branch_best_metric = Inf;
    search_input_frames = min(size(data_normalized, 1), ...
        round(52*51*M/3744) + search_frames + 1);
    search_input = data_normalized(1:search_input_frames, :);
    for direction = directions
        for phase = phases
            [candidate_data, ~] = DTMB_Time_Deinterleave(search_input, M, direction, phase);
            metric = evaluate_candidate(candidate_data, rate_index, search_frames);
            branch_idx = branch_idx + 1;
            branch_records(branch_idx).direction = direction;
            branch_records(branch_idx).commutator_phase = phase;
            branch_records(branch_idx).syndrome_metric = metric;
            if metric < branch_best_metric
                branch_best_metric = metric;
                best_direction = direction;
                best_phase = phase;
            end
        end
    end

    % Deinterleaving is a pure sample permutation, so conjugation and
    % quadrant rotation commute with it.  Deinterleave the full matrix once
    % instead of repeating the same large allocation for all 8 hypotheses.
    [base_deinterleaved,best_time_info] = DTMB_Time_Deinterleave( ...
        data_normalized,M,best_direction,best_phase);
    records = repmat(struct('direction', best_direction, 'commutator_phase', best_phase, ...
        'conjugated', false, 'quadrant', 0, 'syndrome_metric', Inf), 8, 1);
    record_idx = 0;
    best_metric = Inf;
    for conjugated = conjugations
            use_frames = min(search_frames,size(base_deinterleaved,1));
            base_search = base_deinterleaved(1:use_frames,:);
            if conjugated, base_search = conj(base_search); end
            for quadrant = quadrants
                candidate_search = base_search*exp(-1j*quadrant*pi/2);
                metric = evaluate_candidate(candidate_search,rate_index,search_frames);
                record_idx = record_idx + 1;
                records(record_idx).direction = best_direction;
                records(record_idx).commutator_phase = best_phase;
                records(record_idx).conjugated = conjugated;
                records(record_idx).quadrant = quadrant;
                records(record_idx).syndrome_metric = metric;
                if metric < best_metric
                    best_metric = metric;
                    best_record_idx = record_idx;
                end
            end
    end
    [sorted_metric, order] = sort([records.syndrome_metric]);
    best = records(best_record_idx);
    deinterleaved_data = base_deinterleaved;
    if best.conjugated, deinterleaved_data = conj(deinterleaved_data); end
    deinterleaved_data = deinterleaved_data*exp(-1j*best.quadrant*pi/2);
    info = best;
    info.time_deinterleave_info = best_time_info;
    info.candidates = records;
    info.branch_candidates = branch_records;
    info.second_metric = sorted_metric(2);
    info.separation_ratio = sorted_metric(2) / (sorted_metric(1) + eps);
    info.candidate_order = order;
    fprintf('  FEC对齐搜索：方向 %+d，支路相位 %d，共轭 %d，象限 %d；硬综合 %.4f，次优 %.4f。\n', ...
        info.direction, info.commutator_phase, info.conjugated, info.quadrant, ...
        info.syndrome_metric, info.second_metric);
end

function metric = evaluate_candidate(candidate_data, rate_index, search_frames)
    use_frames = min(search_frames, size(candidate_data, 1));
    symbols = reshape(candidate_data(1:use_frames,:).',[],1);
    [~,hard_bits] = DTMB_64QAM_Soft_Demap(symbols,1);
    hard_blocks = reshape(hard_bits,7488,use_frames*3);
    metric = DTMB_LDPC_Hard_Syndrome_Metric(hard_blocks, rate_index);
end

function value = get_cfg(cfg, name, default_value)
    if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
        value = cfg.(name);
    else
        value = default_value;
    end
end
