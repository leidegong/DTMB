function result = run_fast_decision_reconstruction(equalized_reference_frames, phase_index_list, cfg)
%RUN_FAST_DECISION_RECONSTRUCTION Low-latency DTMB reference reconstruction.
%   This branch deliberately omits time deinterleaving, soft demapping,
%   LDPC/BCH, re-encoding and time re-interleaving.  It projects the
%   equalized, still-time-interleaved QAM symbols directly to the nearest
%   constellation points and rebuilds the corresponding transmitted frame.

if nargin < 3, cfg = struct(); end
if iscell(equalized_reference_frames)
    numInput = numel(equalized_reference_frames);
    getFrame = @(k) equalized_reference_frames{k}(:);
else
    numInput = size(equalized_reference_frames,1);
    getFrame = @(k) equalized_reference_frames(k,:).';
end
phase_index_list = phase_index_list(:);
numFrames = min(numInput,numel(phase_index_list));
assert(numFrames > 0,'No equalized frames were supplied.');

nfft = 3780;
frameLength = 4725;
dataRx = complex(zeros(numFrames,3744,'single'));
systemRx = complex(zeros(numFrames,36,'single'));
frameErasure = logical(get_cfg(cfg,'frame_erasure_mask',false(numFrames,1)));
frameErasure = frameErasure(:);
frameErasure(end+1:numFrames,1) = false;
frameErasure = frameErasure(1:numFrames);

timing = struct();
tTotal = tic;
t = tic;
for k = 1:numFrames
    frame = double(getFrame(k));
    assert(numel(frame) >= frameLength,'Frame %d is shorter than 4725 samples.',k);
    fd = fftshift(fft(frame(946:4725),nfft));
    [data,sys] = DTMB_Frequency_Deinterleave(fd);
    dataRx(k,:) = single(data(:).');
    systemRx(k,:) = single(sys(:).');
end
timing.fft_frequency_deinterleave_s = toc(t);

systemForDecode = double(systemRx);
systemForDecode(frameErasure,:) = complex(NaN,NaN);
systemInfo = DTMB_Decode_System_Info(systemForDecode);
candidates = DTMB_System_Info_Candidates();
systemClean = candidates(systemInfo.index).symbols;

t = tic;
hardSymbols = complex(zeros(numFrames,3744,'single'));
framePhase = zeros(numFrames,1);
frameScale = zeros(numFrames,1);
decisionEvm = nan(numFrames,1);
rawInitialPhase = nan(numFrames,1);
trackedInitialPhase = nan(numFrames,1);
rawResidualPhase = nan(numFrames,1);
trackedResidualPhase = nan(numFrames,1);
tpsSignFlip = false(numFrames,1);
qamQuarterTurns = zeros(numFrames,1);
trackPhase = get_cfg(cfg,'phase_continuity_tracking',false);
trackTpsSign = trackPhase && get_cfg(cfg,'tps_sign_disambiguation',true);
trackQamModulo = trackPhase && get_cfg(cfg,'qam_phase_modulo_unwrap',true);
for k = 1:numFrames
    if frameErasure(k)
        hardSymbols(k,:) = single(Slice_64QAM(complex(zeros(3744,1))).');
        frameScale(k) = 1;
        continue
    end
    data = double(dataRx(k,:).');
    rawInitialPhase(k) = angle(systemInfo.per_frame_complex_gain(k));
    initialPhase = rawInitialPhase(k);
    if trackTpsSign && k>1
        target = phase_prediction(trackedInitialPhase,k);
        initialPhase = closest_periodic_phase(initialPhase,target,pi);
        tpsSignFlip(k) = mod(round((initialPhase-rawInitialPhase(k))/pi),2)~=0;
    end
    trackedInitialPhase(k) = initialPhase;
    phaseCorrected = data*exp(-1j*initialPhase);
    initialScale = rms(phaseCorrected(isfinite(phaseCorrected)))+eps;
    rawResidualPhase(k) = Estimate_CPE_64QAM(phaseCorrected/initialScale);
    residualPhase = rawResidualPhase(k);
    if trackQamModulo && k>1
        target = phase_prediction(trackedResidualPhase,k);
        residualPhase = closest_periodic_phase(residualPhase,target,pi/2);
        qamQuarterTurns(k) = round((residualPhase-rawResidualPhase(k))/(pi/2));
    end
    trackedResidualPhase(k) = residualPhase;
    totalPhase = initialPhase+residualPhase;
    phaseCorrected = data*exp(-1j*totalPhase);
    scale = rms(phaseCorrected(isfinite(phaseCorrected)))+eps;
    normalized = phaseCorrected/scale;
    sliced = Slice_64QAM(normalized);
    hardSymbols(k,:) = single(sliced(:).');
    framePhase(k) = totalPhase;
    frameScale(k) = scale;
    decisionEvm(k) = sqrt(mean(abs(normalized-sliced).^2,'omitnan'));
end
timing.cpe_projection_s = toc(t);

pnToBodyPowerRatio = get_cfg(cfg,'pn_to_body_power_ratio',2);
bodyIfftGain = sqrt(nfft/pnToBodyPowerRatio);
canonicalFrames = complex(zeros(numFrames,frameLength,'single'));
t = tic;
for k = 1:numFrames
    bodyFd = DTMB_Frequency_Interleave(double(hardSymbols(k,:).'),systemClean);
    bodyTd = bodyIfftGain*ifft(ifftshift(bodyFd),nfft);
    pn = PN_gen(phase_index_list(k));
    canonicalFrames(k,:) = single([pn(:);bodyTd(:)].');
end
timing.frame_body_ifft_assembly_s = toc(t);

% Match one complex scalar per frame to the equalized observation.  This
% preserves the clean constellation while keeping the local replica in the
% observed phase/amplitude convention.
t = tic;
reconstructedFrames = complex(zeros(size(canonicalFrames),'single'));
matchGain = complex(zeros(numFrames,1));
for k = 1:numFrames
    clean = double(canonicalFrames(k,:).');
    received = double(getFrame(k));
    received = received(1:frameLength);
    gain = (clean'*received)/(clean'*clean+eps);
    matchGain(k) = gain;
    reconstructedFrames(k,:) = single((gain*clean).');
end
timing.complex_gain_match_s = toc(t);

sps = get_cfg(cfg,'shaping_sps',2);
rolloff = get_cfg(cfg,'shaping_rolloff',0.05);
span = get_cfg(cfg,'shaping_span',8);
saveShaped = get_cfg(cfg,'save_shaped_waveform',true);
t = tic;
shapingFilter = rcosdesign(rolloff,span,sps,'sqrt');
if saveShaped
    serialFrames = reshape(reconstructedFrames.',[],1);
    upsampled = zeros(numel(serialFrames)*sps,1,'single');
    upsampled(1:sps:end) = serialFrames;
    shapedReference = conv(upsampled,single(shapingFilter(:)),'full');
else
    shapedReference = complex(zeros(0,1,'single'));
end
timing.srrc_shaping_s = toc(t);
timing.total_s = toc(tTotal);
timing.per_frame_ms = 1e3*timing.total_s/numFrames;
timing.interleaver_buffer_frames = 0;
timing.interleaver_buffer_ms = 0;

result = struct();
result.reconstructed_frames = reconstructedFrames;
result.canonical_reconstructed_frames = canonicalFrames;
result.shaped_reference = shapedReference;
result.shaping_filter = shapingFilter;
result.shaping_sps = sps;
result.hard_interleaved_symbols = hardSymbols;
result.system_info = systemInfo;
result.frame_phase_rad = framePhase;
result.raw_tps_phase_rad = rawInitialPhase;
result.tracked_tps_phase_rad = trackedInitialPhase;
result.raw_residual_cpe_rad = rawResidualPhase;
result.tracked_residual_cpe_rad = trackedResidualPhase;
result.tps_sign_flip = tpsSignFlip;
result.qam_quarter_turns = qamQuarterTurns;
result.phase_continuity_tracking = trackPhase;
result.frame_qam_scale = frameScale;
result.decision_evm_rms = decisionEvm;
result.complex_match_gain = matchGain;
result.frame_erasure_mask = frameErasure;
result.phase_index_list = phase_index_list(1:numFrames);
result.num_frames = numFrames;
result.timing = timing;
result.description = ['Fast DTMB decision reconstruction: FFT/frequency deinterleave, ' ...
    'TPS/CPE correction, nearest-64QAM projection, frequency interleave/IFFT, ' ...
    'PN frame assembly and continuous SRRC shaping. No FEC or time interleaver.'];

outputPacket = get_cfg(cfg,'output_packet','');
if ~isempty(outputPacket)
    save(outputPacket,'result','-v7.3');
end
end

function value = get_cfg(cfg,name,defaultValue)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end

function target = phase_prediction(history,k)
% Constant-rate one-step prediction within the current burst.
if k>=3 && isfinite(history(k-1)) && isfinite(history(k-2))
    target = history(k-1)+(history(k-1)-history(k-2));
else
    target = history(k-1);
end
end

function y = closest_periodic_phase(x,target,period)
% Choose the equivalent x+n*period closest to the continuity prediction.
y = x+round((target-x)/period)*period;
end
