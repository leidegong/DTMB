function [metrics, detail] = evaluate_DTMB_radar_reference_metrics(rawReference, cleanReference, surveillance, fs, opts)
%EVALUATE_DTMB_RADAR_REFERENCE_METRICS Radar-oriented reference quality.
%   [METRICS,DETAIL] = evaluate_DTMB_radar_reference_metrics(RAW,CLEAN,SURV,FS,OPTS)
%   reports waveform NMSE/EVM/correlation, delay-tap ECA cancellation (CA),
%   and CAF peak-to-background SNR. CLEAN and SURV must be time aligned.
%
%   SURV may be []: waveform metrics are still returned and ECA/CAF fields
%   are NaN with an explanatory status.  CAF target SNR is only labelled as
%   target SNR when both opts.target_delay_samples and
%   opts.target_doppler_hz are supplied. Otherwise it is a detected-peak
%   CAF SNR and must not be presented as a known-target result.

if nargin < 5 || isempty(opts), opts = struct(); end
if nargin < 4 || isempty(fs), fs = NaN; end
opts = defaults(opts, struct( ...
    'start_sample',1, 'max_samples',100000, ...
    'eca_delays',0:31, 'eca_doppler_hz',0, ...
    'samples_per_pulse',4725, 'max_pulses',128, ...
    'caf_max_delay_samples',255, 'caf_doppler_fft',128, ...
    'target_delay_samples',[], 'target_doppler_hz',[], ...
    'target_guard_delay',2, 'target_guard_doppler',2, ...
    'clutter_guard_delay',2, 'clutter_guard_doppler',1));

rawReference = double(rawReference(:));
cleanReference = double(cleanReference(:));
n = min(numel(rawReference),numel(cleanReference));
assert(n > 0, 'Reference inputs must not be empty.');
rawReference = rawReference(1:n) - mean(rawReference(1:n));
cleanReference = cleanReference(1:n) - mean(cleanReference(1:n));

% Fit the unavoidable complex gain before comparing waveform quality.
gain = (cleanReference' * rawReference) / (cleanReference' * cleanReference + eps);
fittedClean = gain * cleanReference;
err = rawReference - fittedClean;
metrics = struct();
metrics.waveform_nmse_db = 10*log10(sum(abs(err).^2)/(sum(abs(rawReference).^2)+eps));
metrics.evm_rms = sqrt(sum(abs(err).^2)/(sum(abs(fittedClean).^2)+eps));
metrics.reference_correlation = abs(rawReference'*fittedClean) / ...
    (norm(rawReference)*norm(fittedClean)+eps);
metrics.complex_gain = gain;
metrics.eca_ca_db = NaN;
metrics.caf_target_snr_db = NaN;
metrics.caf_detected_peak_snr_db = NaN;
metrics.status = "waveform-only: surveillance input not supplied";
detail = struct('eca_residual',[],'eca_clutter_estimate',[], ...
    'caf_power',[],'caf_delay_samples',[],'caf_doppler_hz',[]);

if nargin < 3 || isempty(surveillance)
    return
end
assert(isfinite(fs) && fs > 0, 'A positive sampling rate fs is required for ECA/CAF.');
surveillance = double(surveillance(:));
nAll = min([numel(cleanReference),numel(surveillance)]);
first = max(1,round(opts.start_sample));
last = min(nAll,first+round(opts.max_samples)-1);
assert(last >= first, 'Selected evaluation interval is empty.');
x = cleanReference(first:last);
y = surveillance(first:last);
x = x-mean(x); y = y-mean(y);

% Delay/Doppler dictionary ECA. CA is input surveillance power divided by
% residual power over exactly the same valid interval.
delays = unique(round(opts.eca_delays(:).'));
delays = delays(delays >= 0);
maxDelay = max(delays);
assert(numel(x) > maxDelay+1, 'ECA interval is shorter than its delay span.');
valid = (maxDelay+1:numel(x)).';
dopplers = opts.eca_doppler_hz(:).';
X = complex(zeros(numel(valid),numel(delays)*numel(dopplers)));
t = (valid-1)/fs;
col = 0;
for f = dopplers
    tone = exp(1j*2*pi*f*t);
    for d = delays
        col = col+1;
        X(:,col) = x(valid-d).*tone;
    end
end
yValid = y(valid);
beta = X\yValid;
clutter = X*beta;
residual = yValid-clutter;
metrics.eca_ca_db = 10*log10(sum(abs(yValid).^2)/(sum(abs(residual).^2)+eps));
detail.eca_residual = residual;
detail.eca_clutter_estimate = clutter;

% Frame/pulse-based cross ambiguity. Each row is a delay bin and the slow
% time FFT across columns provides Doppler integration.
nFast = round(opts.samples_per_pulse);
nPulse = min(floor(numel(residual)/nFast),round(opts.max_pulses));
if nPulse < 4
    metrics.status = "ECA available; interval too short for CAF";
    return
end
nUse = nFast*nPulse;
xCaf = reshape(x(valid(1:nUse)),nFast,nPulse);
yCaf = reshape(residual(1:nUse),nFast,nPulse);
maxCafDelay = min(round(opts.caf_max_delay_samples),nFast-2);
slow = complex(zeros(maxCafDelay+1,nPulse));
for d = 0:maxCafDelay
    slow(d+1,:) = sum(yCaf(d+1:end,:).*conj(xCaf(1:end-d,:)),1);
end
nDoppler = max(nPulse,round(opts.caf_doppler_fft));
caf = fftshift(fft(slow,nDoppler,2),2);
cafPower = abs(caf).^2;
dopplerAxis = (-floor(nDoppler/2):ceil(nDoppler/2)-1) * ...
    (fs/nFast/nDoppler);
delayAxis = (0:maxCafDelay).';
detail.caf_power = cafPower;
detail.caf_delay_samples = delayAxis;
detail.caf_doppler_hz = dopplerAxis;

knownTarget = ~isempty(opts.target_delay_samples) && ~isempty(opts.target_doppler_hz);
if knownTarget
    [~,r0] = min(abs(delayAxis-opts.target_delay_samples));
    [~,c0] = min(abs(dopplerAxis-opts.target_doppler_hz));
else
    searchMask = true(size(cafPower));
    searchMask(1:min(end,opts.clutter_guard_delay+1),:) = false;
    dc = floor(nDoppler/2)+1;
    searchMask(:,max(1,dc-opts.clutter_guard_doppler): ...
        min(nDoppler,dc+opts.clutter_guard_doppler)) = false;
    candidate = cafPower;
    candidate(~searchMask) = -Inf;
    [~,linearIndex] = max(candidate(:));
    [r0,c0] = ind2sub(size(candidate),linearIndex);
end

noiseMask = true(size(cafPower));
noiseMask(max(1,r0-opts.target_guard_delay):min(end,r0+opts.target_guard_delay), ...
    max(1,c0-opts.target_guard_doppler):min(nDoppler,c0+opts.target_guard_doppler)) = false;
noiseFloor = median(cafPower(noiseMask),'omitnan');
peakSnr = 10*log10(cafPower(r0,c0)/(noiseFloor+eps));
metrics.caf_detected_delay_samples = delayAxis(r0);
metrics.caf_detected_doppler_hz = dopplerAxis(c0);
metrics.caf_detected_peak_snr_db = peakSnr;
if knownTarget
    metrics.caf_target_snr_db = peakSnr;
    metrics.status = "waveform, ECA and known-target CAF metrics available";
else
    metrics.status = "waveform and ECA available; CAF value is detected-peak, not target truth";
end
end

function out = defaults(in, d)
out = in;
names = fieldnames(d);
for k = 1:numel(names)
    if ~isfield(out,names{k}) || isempty(out.(names{k}))
        out.(names{k}) = d.(names{k});
    end
end
end
