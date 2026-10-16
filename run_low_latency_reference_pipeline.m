function result = run_low_latency_reference_pipeline(syncPacket,outputDir,cfgOverride)
%RUN_LOW_LATENCY_REFERENCE_PIPELINE Final algorithmic demod/remod pipeline.
%   The default synchronized-input path emits frame one immediately using
%   a fixed-budget PN channel estimate, then refines the cached equalizer as
%   frame two arrives. Set CFG.INPUT_IS_EQUALIZED=true when equalization is
%   supplied upstream. No learned inverse-channel model is loaded or called.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
if isfolder(fullfile(fileparts(root),'example'))
    artifactRoot = fileparts(root);
else
    artifactRoot = root;
end
if nargin < 1 || isempty(syncPacket)
    if strcmp(artifactRoot,root)
        syncPacket = fullfile(root,'fast_16frame_sync_input.mat');
    else
        syncPacket = fullfile(artifactRoot,'example','fast_16frame_sync_input.mat');
    end
end
if nargin < 2 || isempty(outputDir)
    if strcmp(artifactRoot,root)
        outputDir = fullfile(root,'DTMB_low_latency_final_output');
    else
        outputDir = fullfile(artifactRoot,'outputs');
    end
end
if nargin < 3, cfgOverride = struct(); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
assert(isfile(syncPacket),'Input packet not found: %s',syncPacket);

cfg = defaults();
cfg = merge_cfg(cfg,cfgOverride);
S = load(syncPacket);
assert(isfield(S,'PN_phase_index_packet'), ...
    'Input packet is missing PN_phase_index_packet.');

prewarmTiming = 0;
if cfg.prewarm_runtime
    prewarmInfo = DTMB_Low_Latency_Prewarm(cfg);
    prewarmTiming = prewarmInfo.elapsed_s;
end
tPipeline = tic;
if cfg.input_is_equalized
    if isfield(S,'equalized_reference_frames')
        equalizedFrames = S.equalized_reference_frames;
        inputVariable = 'equalized_reference_frames';
    else
        assert(isfield(S,'Frame_equalization_packet'), ...
            'Equalized input variable was not found.');
        equalizedFrames = S.Frame_equalization_packet;
        inputVariable = 'Frame_equalization_packet';
    end
    channelTiming = 0;
    channelStartupFrames = 0;
    channelEstimateCount = 0;
    channelComputeMs = zeros(0,1);
    frontend = struct('mode',"upstream_equalized",'input_variable',inputVariable, ...
        'channel_estimation_frames',0,'uses_learned_inverse_channel_model',false);
elseif cfg.streaming_channel_bootstrap
    required = {'PN_svd_denoised_packet','Frame_equalization_packet'};
    for k = 1:numel(required)
        assert(isfield(S,required{k}),'Input packet is missing %s.',required{k});
    end
    numFrames = min([size(S.PN_svd_denoised_packet,1), ...
        size(S.Frame_equalization_packet,1),numel(S.PN_phase_index_packet)]);
    if ~isempty(cfg.tc_num_frames), numFrames = min(numFrames,cfg.tc_num_frames); end
    equalizedFrames = complex(zeros(numFrames,4725,'single'));
    channelComputeMs = zeros(numFrames,1);
    channelState = [];
    fprintf('[low-latency 1/2] Immediate-output cached PN/sparse equalization.\n');
    tChannel = tic;
    for k = 1:numFrames
        [equalized,channelState,channelInfo] = ...
            DTMB_Low_Latency_Channel_Equalize_Frame( ...
            S.PN_svd_denoised_packet(k,:),S.Frame_equalization_packet(k,:), ...
            S.PN_phase_index_packet(k),channelState,cfg);
        equalizedFrames(k,:) = equalized.';
        channelComputeMs(k) = 1e3*channelInfo.total_s;
    end
    channelTiming = toc(tChannel);
    channelStartupFrames = 1;
    channelEstimateCount = channelState.num_estimates + ...
        channelState.num_periodic_updates;
    tapThreshold = max(abs(channelState.h_average))*1e-3;
    frontend = struct('mode',"fast_pn_sparse_streaming", ...
        'input_variable','Frame_equalization_packet', ...
        'channel_estimation_frames',channelEstimateCount, ...
        'startup_buffer_frames',channelStartupFrames, ...
        'channel_estimation_compute_s',channelTiming, ...
        'average_channel_tap_count',nnz(abs(channelState.h_average)>tapThreshold), ...
        'fast_channel_num_taps',cfg.fast_channel_num_taps, ...
        'fast_channel_method',string(cfg.fast_channel_method), ...
        'channel_update_interval_frames',cfg.channel_update_interval_frames, ...
        'uses_learned_inverse_channel_model',false);
else
    required = {'PN_svd_denoised_packet','Frame_equalization_packet'};
    for k = 1:numel(required)
        assert(isfield(S,required{k}),'Input packet is missing %s.',required{k});
    end
    tcCfg = cfg;
    tcCfg.plot_diagnostics = false;
    tcCfg.save_diagnostic_figures = false;
    tcCfg.run_nlms_structured_equalizer = false;
    fprintf('[low-latency 1/2] Two-frame conventional channel estimate/equalization.\n');
    tChannel = tic;
    tc = run_tc_aomp_section(S.PN_svd_denoised_packet, ...
        S.Frame_equalization_packet,S.PN_phase_index_packet,tcCfg);
    channelTiming = toc(tChannel);
    equalizedFrames = tc.equalized_reference_frames;
    channelStartupFrames = size(tc.h_hat_all,1);
    channelEstimateCount = channelStartupFrames;
    channelComputeMs = zeros(0,1);
    frontend = struct('mode',"tc_aomp_algorithmic", ...
        'input_variable','Frame_equalization_packet', ...
        'channel_estimation_frames',channelStartupFrames, ...
        'channel_estimation_compute_s',channelTiming, ...
        'average_channel_tap_count',tc.avg_channel_tap_count, ...
        'uses_learned_inverse_channel_model',false);
end

fprintf('[low-latency 2/2] Stateful demodulation/remodulation and continuous SRRC.\n');
result = run_low_latency_demod_remod( ...
    equalizedFrames,S.PN_phase_index_packet,cfg);
result.frontend = frontend;
result.timing.channel_frontend_compute_s = channelTiming;
result.timing.runtime_prewarm_s = prewarmTiming;
if isempty(channelComputeMs)
    result.timing.channel_compute_ms_per_frame_median = NaN;
    result.timing.channel_compute_ms_per_frame_p95 = NaN;
    result.timing.channel_first_frame_compute_ms = 0;
else
    steadyChannel = channelComputeMs(min(numel(channelComputeMs), ...
        cfg.tc_channel_estimation_frames)+1:end);
    if isempty(steadyChannel), steadyChannel = channelComputeMs; end
    sortedChannel = sort(steadyChannel);
    result.timing.channel_compute_ms_per_frame_median = median(steadyChannel);
    result.timing.channel_compute_ms_per_frame_p95 = ...
        sortedChannel(max(1,ceil(0.95*numel(sortedChannel))));
    result.timing.channel_first_frame_compute_ms = channelComputeMs(1);
end
result.timing.steady_end_to_end_ms_per_frame_estimate = ...
    result.timing.steady_compute_ms_per_frame_median + ...
    result.timing.channel_compute_ms_per_frame_median;
result.timing.first_frame_end_to_end_compute_ms = ...
    result.timing.channel_first_frame_compute_ms + ...
    result.timing.first_frame_compute_ms;
result.timing.pipeline_wall_s = toc(tPipeline);
result.timing.channel_startup_buffer_frames = channelStartupFrames;
result.timing.total_startup_buffer_frames = max(1,channelStartupFrames);
frameDurationMs = result.timing.minimum_input_buffer_ms;
result.timing.total_startup_buffer_ms = ...
    result.timing.total_startup_buffer_frames*frameDurationMs;
result.timing.minimum_signal_latency_ms = ...
    result.timing.total_startup_buffer_ms + result.timing.srrc_group_delay_ms;
result.timing.steady_response_latency_ms_estimate = ...
    result.timing.minimum_signal_latency_ms + ...
    result.timing.steady_end_to_end_ms_per_frame_estimate;
result.algorithm_version = "passive_radar_low_latency_v2_20260902";
result.runtime_configuration = struct( ...
    'fast_channel_method',string(cfg.fast_channel_method), ...
    'fast_channel_num_taps',cfg.fast_channel_num_taps, ...
    'channel_bootstrap_estimates',cfg.tc_channel_estimation_frames, ...
    'channel_update_interval_frames',cfg.channel_update_interval_frames, ...
    'cpe_update_interval_frames',cfg.cpe_update_interval_frames, ...
    'uses_fec',false,'uses_inverse_channel_model',false);

summary = table(result.num_frames,string(frontend.mode), ...
    result.uses_inverse_channel_model,channelEstimateCount, ...
    cfg.cpe_update_interval_frames, ...
    result.timing.total_startup_buffer_frames, ...
    result.timing.minimum_signal_latency_ms, ...
    result.timing.first_frame_end_to_end_compute_ms, ...
    result.timing.steady_compute_ms_per_frame_median, ...
    result.timing.steady_compute_ms_per_frame_p95, ...
    result.timing.steady_end_to_end_ms_per_frame_estimate, ...
    result.timing.steady_response_latency_ms_estimate, ...
    result.timing.time_interleaver_buffer_frames, ...
    median(result.decision_evm_rms,'omitnan'),nnz(result.cpe_reacquired(2:end)), ...
    'VariableNames',{'num_frames','frontend','uses_inverse_channel_model', ...
    'channel_estimation_frames','cpe_update_interval_frames', ...
    'startup_buffer_frames', ...
    'minimum_signal_latency_ms','first_frame_compute_ms', ...
    'steady_core_compute_ms_median','steady_core_compute_ms_p95', ...
    'steady_end_to_end_ms_estimate', ...
    'steady_response_latency_ms_estimate', ...
    'time_interleaver_buffer_frames','decision_evm_median', ...
    'cpe_reacquisitions_after_lock'});
result.summary = summary;
writetable(summary,fullfile(outputDir,'low_latency_final_summary.csv'));
if cfg.save_result
    save(fullfile(outputDir,'DTMB_Low_Latency_Demod_Remod_Final.mat'), ...
        'result','-v7.3');
end
disp(summary);
end

function cfg = defaults()
cfg = struct();
cfg.input_is_equalized = false;
cfg.prewarm_runtime = true;
cfg.streaming_channel_bootstrap = true;
cfg.body_eq_lambda = 1e-3;
cfg.tc_num_frames = [];
cfg.tc_channel_estimation_frames = 2;
cfg.fast_channel_num_taps = 32;
cfg.fast_channel_method = 'threshold_ls';
cfg.channel_update_interval_frames = 0;
cfg.channel_update_alpha = 0.25;
cfg.tc_use_average_channel = true;
cfg.tc_robust_average_channel = false;
cfg.phase_continuity_tracking = true;
cfg.phase_rate_prediction = false;
cfg.cpe_tracking_method = 'decision_directed';
cfg.cpe_dd_max_correction_deg = 8;
cfg.cpe_update_interval_frames = 4;
cfg.cpe_initial_update_frames = 2;
cfg.pn_to_body_power_ratio = 2;
cfg.shaping_sps = 2;
cfg.shaping_rolloff = 0.05;
cfg.shaping_span = 8;
cfg.symbol_rate_hz = 7.56e6;
cfg.save_shaped_waveform = true;
cfg.flush_srrc_tail = true;
cfg.save_result = true;
end

function cfg = merge_cfg(cfg,override)
names = fieldnames(override);
for k = 1:numel(names), cfg.(names{k}) = override.(names{k}); end
end
