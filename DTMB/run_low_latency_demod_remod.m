function result = run_low_latency_demod_remod(equalizedFrames,phaseIndices,cfg)
%RUN_LOW_LATENCY_DEMOD_REMOD Run the stateful frame core over one burst.
%   This convenience wrapper retains batch outputs for analysis while using
%   the exact same one-frame API and persistent FIR state as deployment.

if nargin < 3, cfg = struct(); end
if iscell(equalizedFrames)
    numInput = numel(equalizedFrames);
    getFrame = @(k) equalizedFrames{k};
else
    numInput = size(equalizedFrames,1);
    getFrame = @(k) equalizedFrames(k,:).';
end
phaseIndices = phaseIndices(:);
numFrames = min(numInput,numel(phaseIndices));
assert(numFrames > 0,'No equalized frames were supplied.');

frameLength = 4725;
dataSymbols = 3744;
sps = get_cfg(cfg,'shaping_sps',2);
collectShaped = get_cfg(cfg,'save_shaped_waveform',true);
reconstructed = complex(zeros(numFrames,frameLength,'single'));
canonical = complex(zeros(numFrames,frameLength,'single'));
hard = complex(zeros(numFrames,dataSymbols,'single'));
if collectShaped
    shaped = complex(zeros(numFrames*frameLength*sps,1,'single'));
else
    shaped = complex(zeros(0,1,'single'));
end
framePhase = zeros(numFrames,1);
tpsPhase = zeros(numFrames,1);
cpePhase = zeros(numFrames,1);
qamScale = zeros(numFrames,1);
decisionEvm = zeros(numFrames,1);
matchGain = complex(zeros(numFrames,1));
computeMs = zeros(numFrames,1);
cpeCandidates = zeros(numFrames,1);
cpeReacquired = false(numFrames,1);
cpeSearchMode = strings(numFrames,1);

state = [];
tTotal = tic;
for k = 1:numFrames
    [frameOut,state] = DTMB_Low_Latency_Demod_Remod_Frame( ...
        getFrame(k),phaseIndices(k),state,cfg);
    reconstructed(k,:) = frameOut.reconstructed_frame;
    canonical(k,:) = frameOut.canonical_reconstructed_frame;
    hard(k,:) = frameOut.hard_interleaved_symbols;
    if collectShaped
        first = (k-1)*frameLength*sps+1;
        shaped(first:first+frameLength*sps-1) = frameOut.shaped_reference;
    end
    framePhase(k) = frameOut.total_phase_rad;
    tpsPhase(k) = frameOut.tracked_tps_phase_rad;
    cpePhase(k) = frameOut.tracked_residual_cpe_rad;
    qamScale(k) = frameOut.qam_scale;
    decisionEvm(k) = frameOut.decision_evm_rms;
    matchGain(k) = frameOut.complex_match_gain;
    computeMs(k) = 1e3*frameOut.timing.total_s;
    cpeCandidates(k) = frameOut.cpe.candidate_count;
    cpeReacquired(k) = frameOut.cpe.reacquired;
    cpeSearchMode(k) = frameOut.cpe.search_mode;
end
if collectShaped && get_cfg(cfg,'flush_srrc_tail',true)
    [tail,state] = DTMB_Low_Latency_Flush(state);
    shaped = [shaped;tail];
end
wallTime = toc(tTotal);

symbolRate = get_cfg(cfg,'symbol_rate_hz',7.56e6);
filterGroupDelaySamples = (numel(state.shaping_filter)-1)/2;
timing = struct();
timing.total_wall_s = wallTime;
timing.total_core_compute_s = state.total_compute_s;
timing.compute_ms_per_frame_median = median(computeMs);
sortedComputeMs = sort(computeMs);
timing.compute_ms_per_frame_p95 = ...
    sortedComputeMs(max(1,ceil(0.95*numel(sortedComputeMs))));
warmupFrames = min(get_cfg(cfg,'performance_warmup_frames',5), ...
    max(0,numFrames-1));
steadyComputeMs = computeMs(warmupFrames+1:end);
sortedSteadyMs = sort(steadyComputeMs);
timing.first_frame_compute_ms = computeMs(1);
timing.performance_warmup_frames = warmupFrames;
timing.steady_compute_ms_per_frame_median = median(steadyComputeMs);
timing.steady_compute_ms_per_frame_p95 = ...
    sortedSteadyMs(max(1,ceil(0.95*numel(sortedSteadyMs))));
timing.compute_ms_per_frame = computeMs;
timing.minimum_input_buffer_frames = 1;
timing.minimum_input_buffer_ms = 1e3*frameLength/symbolRate;
timing.time_interleaver_buffer_frames = 0;
timing.time_interleaver_buffer_ms = 0;
timing.srrc_group_delay_samples = filterGroupDelaySamples;
timing.srrc_group_delay_ms = 1e3*filterGroupDelaySamples/(symbolRate*sps);

result = struct();
result.reconstructed_frames = reconstructed;
result.canonical_reconstructed_frames = canonical;
result.shaped_reference = shaped;
result.shaping_filter = state.shaping_filter;
result.shaping_sps = state.sps;
result.hard_interleaved_symbols = hard;
result.system_info_index = state.system_info_index;
result.frame_phase_rad = framePhase;
result.tracked_tps_phase_rad = tpsPhase;
result.tracked_residual_cpe_rad = cpePhase;
result.frame_qam_scale = qamScale;
result.decision_evm_rms = decisionEvm;
result.complex_match_gain = matchGain;
result.cpe_candidate_count = cpeCandidates;
result.cpe_reacquired = cpeReacquired;
result.cpe_search_mode = cpeSearchMode;
result.phase_index_list = phaseIndices(1:numFrames);
result.num_frames = numFrames;
result.final_state = state;
result.timing = timing;
result.uses_inverse_channel_model = false;
result.description = ['Low-latency algorithmic DTMB demodulation/remodulation: ' ...
    'one-frame FFT, TPS/tracked CPE, nearest-64QAM projection, IFFT/PN ' ...
    'regeneration, and stateful continuous SRRC. No learned inverse-channel ' ...
    'model, FEC, or time interleaving.'];
end

function value = get_cfg(cfg,name,defaultValue)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
