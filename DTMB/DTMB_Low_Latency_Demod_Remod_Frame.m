function [out,state] = DTMB_Low_Latency_Demod_Remod_Frame(equalizedFrame,phaseIndex,state,cfg)
%DTMB_LOW_LATENCY_DEMOD_REMOD_FRAME Process one equalized DTMB frame.
%   [OUT,STATE] = ... preserves TPS/CPE tracking and SRRC FIR state between
%   calls.  The path is strictly algorithmic: FFT, carrier permutation,
%   TPS/CPE estimation, nearest-64QAM projection, IFFT, PN assembly and
%   stateful SRRC shaping.  It has no FEC, time interleaver, neural network,
%   learned residual, or inverse-channel-model dependency.

if nargin < 3 || isempty(state), state = struct(); end
if nargin < 4, cfg = struct(); end
frameLength = 4725;
nfft = 3780;
frame = equalizedFrame(:);
assert(numel(frame) >= frameLength,'Equalized frame must contain 4725 symbols.');
frame = single(frame(1:frameLength));
if ~isfield(state,'initialized') || ~state.initialized
    state = initialize_state(cfg);
end

tTotal = tic;
t = tic;
[zToY,tpsIndex] = DTMB_Frequency_Permutation();
dataMask = true(nfft,1);
dataMask(tpsIndex) = false;
ySymbols = fft(frame(946:4725),nfft);
zSymbols = ySymbols(zToY);
data = zSymbols(dataMask);
systemRx = zSymbols(tpsIndex);
timing.frequency_frontend_s = toc(t);

t = tic;
if isempty(state.system_symbols)
    candidates = DTMB_System_Info_Candidates();
    if isempty(state.system_info_index)
        decoded = DTMB_Decode_System_Info(double(systemRx.'));
        state.system_info_index = decoded.index;
    else
        decoded = candidates(state.system_info_index);
    end
    state.system_info = candidates(state.system_info_index);
    state.system_symbols = single(state.system_info.symbols);
else
    decoded = state.system_info;
end
systemClean = state.system_symbols;
systemGain = (systemClean'*systemRx)/(systemClean'*systemClean+eps);
systemScore = abs(systemClean'*systemRx)^2 / ...
    ((systemClean'*systemClean)*(systemRx'*systemRx)+eps);
rawTpsPhase = angle(systemGain);
trackedTpsPhase = rawTpsPhase;
tpsSignFlip = false;
if state.frame_count > 0 && state.track_phase
    targetTps = predict_next(state.tps_history,state.cfg);
    trackedTpsPhase = closest_periodic(rawTpsPhase,targetTps,pi);
    tpsSignFlip = mod(round((trackedTpsPhase-rawTpsPhase)/pi),2) ~= 0;
end

phaseCorrected = data*exp(-1j*trackedTpsPhase);
initialScale = sqrt(mean(abs(phaseCorrected(isfinite(phaseCorrected))).^2))+eps;
normalizedInitial = phaseCorrected/initialScale;
if state.frame_count > 0 && state.track_phase
    targetCpe = predict_next(state.cpe_history,state.cfg);
else
    targetCpe = NaN;
end
frameNumber = state.frame_count+1;
cpeInterval = max(1,round(get_cfg(state.cfg,'cpe_update_interval_frames',1)));
cpeInitial = max(1,round(get_cfg(state.cfg,'cpe_initial_update_frames',2)));
updateCpe = frameNumber <= cpeInitial || ...
    mod(frameNumber-cpeInitial,cpeInterval) == 0;
if updateCpe || ~isfinite(targetCpe)
    [trackedCpe,cpeInfo,hardSymbols] = Estimate_CPE_64QAM_Low_Latency( ...
        normalizedInitial,targetCpe,state.cfg);
else
    trackedCpe = targetCpe;
    corrected = normalizedInitial*exp(-1j*trackedCpe);
    hardSymbols = Slice_64QAM(corrected);
    fullMetric = mean(abs(corrected-hardSymbols).^2);
    cpeInfo = struct('reacquired',false,'candidate_count',0, ...
        'metric',NaN,'full_metric',fullMetric,'search_mode',"hold");
end
totalPhase = trackedTpsPhase + trackedCpe;
qamScale = initialScale;
decisionEvm = sqrt(cpeInfo.full_metric);
timing.phase_projection_s = toc(t);

t = tic;
bodyIfftGain = sqrt(nfft/state.pn_to_body_power_ratio);
cleanZ = complex(zeros(nfft,1,'like',data));
cleanZ(tpsIndex) = systemClean;
cleanZ(dataMask) = hardSymbols;
cleanY = complex(zeros(nfft,1,'like',data));
cleanY(zToY) = cleanZ;
body = bodyIfftGain*ifft(cleanY,nfft);
pn = single(PN_gen(phaseIndex));
canonicalFrame = [pn(:);body(:)];
matchGain = (canonicalFrame'*frame)/(canonicalFrame'*canonicalFrame+eps);
reconstructedFrame = matchGain*canonicalFrame;
timing.frame_regeneration_s = toc(t);

t = tic;
upsampled = complex(zeros(frameLength*state.sps,1,'single'));
upsampled(1:state.sps:end) = single(reconstructedFrame);
[shaped,state.srrc_zi] = filter(state.shaping_filter,1,upsampled,state.srrc_zi);
timing.srrc_s = toc(t);
timing.total_s = toc(tTotal);

state.frame_count = state.frame_count + 1;
state.tps_history = append_history(state.tps_history,trackedTpsPhase);
state.cpe_history = append_history(state.cpe_history,trackedCpe);
state.total_compute_s = state.total_compute_s + timing.total_s;
state.reacquisition_count = state.reacquisition_count + double(cpeInfo.reacquired);

decoded.per_frame_complex_gain = systemGain;
decoded.per_frame_score = systemScore;
out = struct();
out.reconstructed_frame = single(reconstructedFrame.');
out.canonical_reconstructed_frame = single(canonicalFrame.');
out.shaped_reference = shaped;
out.hard_interleaved_symbols = single(hardSymbols.');
out.system_info = decoded;
out.raw_tps_phase_rad = rawTpsPhase;
out.tracked_tps_phase_rad = trackedTpsPhase;
out.tracked_residual_cpe_rad = trackedCpe;
out.total_phase_rad = totalPhase;
out.tps_sign_flip = tpsSignFlip;
out.qam_scale = qamScale;
out.decision_evm_rms = decisionEvm;
out.complex_match_gain = matchGain;
out.cpe = cpeInfo;
out.timing = timing;
end

function state = initialize_state(cfg)
state = struct();
state.initialized = true;
state.cfg = cfg;
state.sps = get_cfg(cfg,'shaping_sps',2);
rolloff = get_cfg(cfg,'shaping_rolloff',0.05);
span = get_cfg(cfg,'shaping_span',8);
state.shaping_filter = single(rcosdesign(rolloff,span,state.sps,'sqrt').');
state.srrc_zi = complex(zeros(numel(state.shaping_filter)-1,1,'single'));
state.pn_to_body_power_ratio = get_cfg(cfg,'pn_to_body_power_ratio',2);
state.track_phase = get_cfg(cfg,'phase_continuity_tracking',true);
state.system_info_index = get_cfg(cfg,'system_info_index',[]);
state.system_info = struct();
state.system_symbols = [];
state.tps_history = zeros(0,1);
state.cpe_history = zeros(0,1);
state.frame_count = 0;
state.total_compute_s = 0;
state.reacquisition_count = 0;
end

function history = append_history(history,value)
history = [history(:);value];
if numel(history) > 2, history = history(end-1:end); end
end

function target = predict_next(history,cfg)
useRate = get_cfg(cfg,'phase_rate_prediction',false);
if useRate && numel(history) >= 2 && all(isfinite(history(end-1:end)))
    target = history(end) + history(end)-history(end-1);
elseif ~isempty(history)
    target = history(end);
else
    target = NaN;
end
end

function value = closest_periodic(value,target,period)
value = value + round((target-value)/period)*period;
end

function value = get_cfg(cfg,name,defaultValue)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
