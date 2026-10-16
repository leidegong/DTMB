function report = benchmark_low_latency_passive_radar()
%BENCHMARK_LOW_LATENCY_PASSIVE_RADAR Tune latency without using FEC output.
%   Sweeps the fixed-budget PN channel estimator and residual-CPE update
%   interval against the existing complete-link reconstruction.  The full
%   link is used only as offline truth for this benchmark, never at runtime.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
sync = load(fullfile(root,'DTMB_1100Frame_Sync_Equalized_Packet.mat'));
gold = load(fullfile(root,'DTMB_groupwise_equalization', ...
    'DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat'),'result');
gold = gold.result;
n = min([size(sync.PN_svd_denoised_packet,1), ...
    size(sync.Frame_equalization_packet,1), ...
    numel(sync.PN_phase_index_packet),size(gold.canonical_reconstructed_frames,1)]);
phase = sync.PN_phase_index_packet(1:n);
mask = logical(gold.fec_valid_frame_mask(1:n));
assert(any(mask),'The complete-link file contains no FEC-valid frames.');

tapCandidates = [12 16 20 24 32 40 48];
tapRows = repmat(empty_row(),numel(tapCandidates),1);
for i = 1:numel(tapCandidates)
    cfg = base_cfg();
    cfg.fast_channel_num_taps = tapCandidates(i);
    [equalized,channelMedianMs,firstChannelMs,estimatorMedianMs] = ...
        equalize_all(sync,phase,n,cfg);
    low = run_low_latency_demod_remod(equalized,phase,cfg);
    tapRows(i) = compare_row("channel_taps",tapCandidates(i),low,gold,mask, ...
        channelMedianMs,firstChannelMs,estimatorMedianMs);
end
tapSweep = struct2table(tapRows);
bestCorrelation = max([tapRows.correlation_to_full]);
bestNmse = min([tapRows.nmse_to_full_db]);
tapEligible = find([tapRows.correlation_to_full] >= bestCorrelation-5e-4 & ...
    [tapRows.nmse_to_full_db] <= bestNmse+0.20);
bestTap = tapCandidates(tapEligible(1));

cpeCandidates = [1 2 4 8];
cpeRows = repmat(empty_row(),numel(cpeCandidates),1);
cfg = base_cfg();
cfg.fast_channel_num_taps = bestTap;
[bestEqualized,channelMedianMs,firstChannelMs,estimatorMedianMs] = ...
    equalize_all(sync,phase,n,cfg);
lowCandidates = cell(numel(cpeCandidates),1);
for i = 1:numel(cpeCandidates)
    cfg.cpe_update_interval_frames = cpeCandidates(i);
    lowCandidates{i} = run_low_latency_demod_remod(bestEqualized,phase,cfg);
    cpeRows(i) = compare_row("cpe_interval",cpeCandidates(i), ...
        lowCandidates{i},gold,mask,channelMedianMs,firstChannelMs, ...
        estimatorMedianMs);
end
phaseReference = lowCandidates{1}.tracked_residual_cpe_rad;
for i = 1:numel(cpeCandidates)
    phaseDelta = angle(exp(1j*(lowCandidates{i}.tracked_residual_cpe_rad- ...
        phaseReference)))*180/pi;
    cpeRows(i).phase_delta_rms_deg = sqrt(mean(phaseDelta.^2));
    sortedDelta = sort(abs(phaseDelta));
    cpeRows(i).phase_delta_p95_deg = ...
        sortedDelta(max(1,ceil(0.95*numel(sortedDelta))));
end
cpeSweep = struct2table(cpeRows);

% Preserve phase coherence: accept an interval only when its correlation
% loss is negligible relative to per-frame CPE updates.
corrFloor = cpeRows(1).correlation_to_full-1e-4;
nmseCeiling = cpeRows(1).nmse_to_full_db+0.10;
eligible = find([cpeRows.correlation_to_full] >= corrFloor & ...
    [cpeRows.nmse_to_full_db] <= nmseCeiling & ...
    [cpeRows.phase_delta_p95_deg] <= 1.0);
selectedIndex = eligible(end);
selectedCpeInterval = cpeCandidates(selectedIndex);
selected = lowCandidates{selectedIndex};
[selectedIslr,selectedPslr] = ambiguity_sidelobes( ...
    selected.canonical_reconstructed_frames(mask,:));
[fullIslr,fullPslr] = ambiguity_sidelobes( ...
    gold.canonical_reconstructed_frames(mask,:));

report = struct();
report.algorithm_version = "passive_radar_low_latency_v2_20260902";
report.num_frames = n;
report.num_fec_valid_truth_frames = nnz(mask);
report.selected_channel_taps = bestTap;
report.selected_cpe_update_interval_frames = selectedCpeInterval;
report.tap_sweep = tapSweep;
report.cpe_sweep = cpeSweep;
report.selected_islr_db = selectedIslr;
report.selected_pslr_db = selectedPslr;
report.full_link_islr_db = fullIslr;
report.full_link_pslr_db = fullPslr;
report.note = ['FEC is used only to define offline truth frames. The selected ' ...
    'runtime path contains no FEC or inverse-channel model.'];

outputDir = fullfile(root,'DTMB_low_latency_final_output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
writetable(tapSweep,fullfile(outputDir,'passive_radar_channel_tap_sweep.csv'));
writetable(cpeSweep,fullfile(outputDir,'passive_radar_cpe_interval_sweep.csv'));
save(fullfile(outputDir,'passive_radar_latency_benchmark.mat'),'report');
disp(tapSweep);
disp(cpeSweep);
fprintf(['Selected: %d channel taps, CPE every %d frame(s), ' ...
    'ISLR %.4f dB, PSLR %.4f dB.\n'],bestTap,selectedCpeInterval, ...
    selectedIslr,selectedPslr);
end

function cfg = base_cfg()
cfg = struct('fast_channel_num_taps',4, ...
    'fast_channel_method','threshold_ls', ...
    'tc_channel_estimation_frames',2, ...
    'channel_update_interval_frames',0, ...
    'cpe_update_interval_frames',1, ...
    'cpe_initial_update_frames',2, ...
    'cpe_tracking_method','decision_directed', ...
    'save_shaped_waveform',false,'flush_srrc_tail',false, ...
    'performance_warmup_frames',8);
end

function [equalized,steadyMedianMs,firstMs,estimatorMedianMs] = ...
        equalize_all(sync,phase,n,cfg)
equalized = complex(zeros(n,4725,'single'));
computeMs = zeros(n,1);
state = [];
% Warm MATLAB dispatch once, then report the repeatable estimate cost.
DTMB_Fast_PN_Channel_Estimate(sync.PN_svd_denoised_packet(1,:),phase(1),cfg);
estimatorMs = zeros(7,1);
for i = 1:numel(estimatorMs)
    t = tic;
    DTMB_Fast_PN_Channel_Estimate( ...
        sync.PN_svd_denoised_packet(1,:),phase(1),cfg);
    estimatorMs(i) = 1e3*toc(t);
end
estimatorMedianMs = median(estimatorMs);
for k = 1:n
    [frame,state,info] = DTMB_Low_Latency_Channel_Equalize_Frame( ...
        sync.PN_svd_denoised_packet(k,:), ...
        sync.Frame_equalization_packet(k,:),phase(k),state,cfg);
    equalized(k,:) = frame.';
    computeMs(k) = 1e3*info.total_s;
end
firstMs = computeMs(1);
steadyMedianMs = median(computeMs(min(n,2)+1:end));
end

function row = compare_row(kind,value,low,gold,mask,channelMedianMs, ...
        firstChannelMs,estimatorMedianMs)
a = double(low.canonical_reconstructed_frames(mask,:));
b = double(gold.canonical_reconstructed_frames(mask,:));
a = a(:); b = b(:);
gain = (a'*b)/(a'*a+eps);
error = gain*a-b;
row = empty_row();
row.sweep = string(kind);
row.value = value;
row.decision_evm_median = median(low.decision_evm_rms(mask),'omitnan');
row.nmse_to_full_db = 10*log10(sum(abs(error).^2)/(sum(abs(b).^2)+eps));
row.correlation_to_full = abs(a'*b)/(norm(a)*norm(b)+eps);
row.channel_median_ms = channelMedianMs;
row.core_median_ms = low.timing.steady_compute_ms_per_frame_median;
row.end_to_end_median_ms = channelMedianMs+row.core_median_ms;
row.first_channel_ms = firstChannelMs;
row.channel_estimator_median_ms = estimatorMedianMs;
end

function row = empty_row()
row = struct('sweep',"",'value',NaN,'decision_evm_median',NaN, ...
    'nmse_to_full_db',NaN,'correlation_to_full',NaN, ...
    'channel_median_ms',NaN,'core_median_ms',NaN, ...
    'end_to_end_median_ms',NaN,'first_channel_ms',NaN, ...
    'channel_estimator_median_ms',NaN, ...
    'phase_delta_rms_deg',NaN,'phase_delta_p95_deg',NaN);
end

function [islrDb,pslrDb] = ambiguity_sidelobes(frames)
numUse = min(size(frames,1),64);
indices = unique(round(linspace(1,size(frames,1),numUse)));
powerSum = zeros(2*size(frames,2)-1,1);
for k = indices
    x = double(frames(k,:)).';
    x = x-mean(x);
    powerSum = powerSum+abs(xcorr(x,x)).^2;
end
powerMean = powerSum/numel(indices);
[mainPower,peak] = max(powerMean);
side = powerMean;
side(peak) = [];
islrDb = 10*log10(mean(side)/(mainPower+eps));
pslrDb = 10*log10(max(side)/(mainPower+eps));
end
