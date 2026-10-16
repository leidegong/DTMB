function report = run_fast_reconstruction_experiment()
%RUN_FAST_RECONSTRUCTION_EXPERIMENT Build and test the independent fast path.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
outDir = fullfile(root,'DTMB_fast_reconstruction');
if ~exist(outDir,'dir'), mkdir(outDir); end

groupFile = fullfile(root,'DTMB_groupwise_equalization', ...
    'DTMB_Groupwise_Equalized_Frames.mat');
syncFile = fullfile(root,'DTMB_1100Frame_Sync_Equalized_Packet.mat');
fullFile = fullfile(root,'DTMB_groupwise_equalization', ...
    'DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat');
profileFile = fullfile(root,'DTMB_module_profile','DTMB_Module_Profile_Report.mat');

fprintf('[1/5] Loading equalized frames and full-FEC baseline.\n');
G = load(groupFile,'optimized_frames');
S = load(syncFile,'PN_phase_index_packet');
F = load(fullFile,'result');
numFrames = min([size(G.optimized_frames,1),numel(S.PN_phase_index_packet), ...
    F.result.num_frames]);
inputFrames = G.optimized_frames(1:numFrames,:);
phaseIndex = S.PN_phase_index_packet(1:numFrames);

cfg = struct('save_shaped_waveform',true,'shaping_sps',2, ...
    'shaping_rolloff',0.05,'shaping_span',8,'output_packet','');

fprintf('[2/5] Memory-profile pass (also warms MATLAB JIT).\n');
profile clear;
profile on -memory;
profiledResult = run_fast_decision_reconstruction(inputFrames,phaseIndex,cfg); %#ok<NASGU>
profile off;
p = profile('info');
peakBytes = 0;
if ~isempty(p.FunctionTable), peakBytes = max([p.FunctionTable.PeakMem]); end
w = whos('profiledResult');
peakBytes = max(peakBytes,w.bytes);
clear profiledResult

fprintf('[3/5] Timed fast reconstruction pass.\n');
cfg.output_packet = fullfile(outDir,'DTMB_Fast_Reconstruction.mat');
fastResult = run_fast_decision_reconstruction(inputFrames,phaseIndex,cfg);

fprintf('[4/5] Computing matched-frame quality metrics.\n');
fullMask = logical(F.result.fec_valid_frame_mask(1:numFrames));
assert(any(fullMask),'The full baseline contains no FEC-valid frames.');
fullFrames = F.result.reconstructed_frames(1:numFrames,:);
fastFrames = fastResult.reconstructed_frames;
systemSymbols = fastResult.system_info.symbols;

fastAll = frame_quality(double(inputFrames),double(fastFrames),systemSymbols);
inputValid = double(inputFrames(fullMask,:));
fastValidFrames = double(fastFrames(fullMask,:));
fullValidFrames = double(fullFrames(fullMask,:));
fastValid = frame_quality(inputValid,fastValidFrames,systemSymbols);
fullValid = frame_quality(inputValid,fullValidFrames,systemSymbols);

selected = unique(round(linspace(1,nnz(fullMask),min(50,nnz(fullMask)))));
[islrInput,pslrInput] = autocorrelation_quality(inputValid(selected,:));
[islrFast,pslrFast] = autocorrelation_quality(fastValidFrames(selected,:));
[islrFull,pslrFull] = autocorrelation_quality(fullValidFrames(selected,:));

[directNmseDb,directCorrelation] = waveform_match(fastValidFrames,fullValidFrames);
symbolAgreement = compare_frame_symbols( ...
    double(fastResult.canonical_reconstructed_frames(fullMask,:)), ...
    double(F.result.canonical_reconstructed_frames(fullMask,:)));

method = ["Equalized input";"Fast decision";"Full FEC"];
num_frames = repmat(nnz(fullMask),3,1);
evm_median = [median(fastValid.reference_evm,'omitnan'); ...
    median(fastValid.test_evm,'omitnan');median(fullValid.test_evm,'omitnan')];
waveform_nmse_db_median = [NaN;median(fastValid.nmse_db,'omitnan'); ...
    median(fullValid.nmse_db,'omitnan')];
correlation_median = [1;median(fastValid.correlation,'omitnan'); ...
    median(fullValid.correlation,'omitnan')];
islr_db_median = [median(islrInput,'omitnan');median(islrFast,'omitnan'); ...
    median(islrFull,'omitnan')];
pslr_db_median = [median(pslrInput,'omitnan');median(pslrFast,'omitnan'); ...
    median(pslrFull,'omitnan')];
% PSLR is set by the PN945 cyclic-extension sidelobe at lag +-511, so it
% tracks the PN:body power ratio and says nothing about body reconstruction
% quality unless the ratios match.  The equalized input carries a different
% ratio than the remodulated frames, so the ratio is reported alongside and
% the PSLR column must not be read as a quality ranking on its own.
% See sweep_pn_body_power_ratio.
pn_body_power_ratio = [median(pn_body_ratio(inputValid),'omitnan'); ...
    median(pn_body_ratio(fastValidFrames),'omitnan'); ...
    median(pn_body_ratio(fullValidFrames),'omitnan')];
comparisonTable = table(method,num_frames,evm_median,waveform_nmse_db_median, ...
    correlation_median,pn_body_power_ratio,islr_db_median,pslr_db_median);

fullComputeMs = NaN;
sharedFrontMs = NaN;
fullCacheMs = 638.125;
if exist(profileFile,'file')
    P = load(profileFile,'report');
    fullComputeMs = sum(P.report.module_table.per_frame_ms);
    moduleNames = string(P.report.module_table.module);
    sharedFrontMs = sum(P.report.module_table.per_frame_ms( ...
        moduleNames=="TC-AOMP" | moduleNames=="Equalization"));
end
fastReconstructionMs = fastResult.timing.per_frame_ms;
fastEndToEndMs = sharedFrontMs+fastReconstructionMs;

report = struct();
report.generated_at = datetime('now');
report.num_frames = numFrames;
report.fec_valid_frames = nnz(fullMask);
report.fec_valid_frame_indices = find(fullMask);
report.fast_timing = fastResult.timing;
report.fast_peak_additional_mib = peakBytes/2^20;
report.full_profile_compute_ms_per_frame = fullComputeMs;
report.shared_tc_aomp_equalization_ms_per_frame = sharedFrontMs;
report.fast_reconstruction_only_ms_per_frame = fastReconstructionMs;
report.fast_end_to_end_compute_ms_per_frame = fastEndToEndMs;
report.full_interleaver_cache_ms = fullCacheMs;
report.fast_interleaver_cache_ms = 0;
report.compute_speedup_vs_profile = fullComputeMs/fastEndToEndMs;
report.compute_reduction_percent = 100*(1-fastEndToEndMs/fullComputeMs);
report.fast_all_frames_nmse_db_median = median(fastAll.nmse_db,'omitnan');
report.fast_all_frames_evm_median = median(fastAll.test_evm,'omitnan');
report.fast_all_frames_correlation_median = median(fastAll.correlation,'omitnan');
report.fast_vs_full_nmse_db_median = median(directNmseDb,'omitnan');
report.fast_vs_full_correlation_median = median(directCorrelation,'omitnan');
report.fast_vs_full_qam_symbol_agreement = symbolAgreement;
report.comparison_table = comparisonTable;
report.fast_per_frame = fastAll;
report.fast_valid_per_frame = fastValid;
report.full_valid_per_frame = fullValid;
report.radar_metrics_status = ['NMSE/EVM/correlation/autocorrelation are available. ' ...
    'ECA-CA and known-target CAF-SNR remain pending until a time-aligned ' ...
    'surveillance channel and target ROI are supplied.'];

writetable(comparisonTable,fullfile(outDir,'fast_vs_full_quality.csv'));
save(fullfile(outDir,'DTMB_Fast_Reconstruction_Report.mat'),'report','-v7.3');
make_figures(inputValid,fastValidFrames,fullValidFrames,fastAll,fullMask, ...
    comparisonTable,report,outDir);

fprintf('[5/5] Complete.\n');
disp(comparisonTable);
fprintf(['Fast reconstruction after equalization: %.3f ms/frame; ' ...
    'fast end-to-end estimate: %.3f ms/frame; full profile: %.3f ms/frame.\n'], ...
    fastReconstructionMs,fastEndToEndMs,fullComputeMs);
fprintf('Fast/full symbol agreement on %d FEC-valid frames: %.6f.\n', ...
    nnz(fullMask),symbolAgreement);
end

function q = frame_quality(reference,test,systemSymbols)
n = min(size(reference,1),size(test,1));
reference = reference(1:n,:); test = test(1:n,:);
err = reference-test;
q.nmse_db = 10*log10(sum(abs(err).^2,2)./(sum(abs(reference).^2,2)+eps));
q.correlation = abs(sum(conj(reference).*test,2))./ ...
    (sqrt(sum(abs(reference).^2,2).*sum(abs(test).^2,2))+eps);
q.reference_evm = frame_evm(reference,systemSymbols);
q.test_evm = frame_evm(test,systemSymbols);
end

function evm = frame_evm(frames,systemSymbols)
n = size(frames,1);
evm = nan(n,1);
for k = 1:n
    fd = fftshift(fft(frames(k,946:4725).',3780));
    [data,sys] = DTMB_Frequency_Deinterleave(fd);
    phase0 = angle((systemSymbols'*sys)/(systemSymbols'*systemSymbols+eps));
    z = data*exp(-1j*phase0);
    scale = rms(z)+eps;
    residualPhase = Estimate_CPE_64QAM(z/scale);
    z = data*exp(-1j*(phase0+residualPhase));
    z = z/(rms(z)+eps);
    evm(k) = sqrt(mean(abs(z-Slice_64QAM(z)).^2,'omitnan'));
end
end

function [nmseDb,corrValue] = waveform_match(a,b)
err = a-b;
nmseDb = 10*log10(sum(abs(err).^2,2)./(sum(abs(b).^2,2)+eps));
corrValue = abs(sum(conj(a).*b,2))./ ...
    (sqrt(sum(abs(a).^2,2).*sum(abs(b).^2,2))+eps);
end

function agreement = compare_frame_symbols(a,b)
assert(size(a,1)==size(b,1));
same = 0; total = 0;
for k = 1:size(a,1)
    afd = fftshift(fft(a(k,946:4725).',3780));
    bfd = fftshift(fft(b(k,946:4725).',3780));
    ad = Slice_64QAM(DTMB_Frequency_Deinterleave(afd));
    bd = Slice_64QAM(DTMB_Frequency_Deinterleave(bfd));
    same = same+nnz(abs(ad-bd)<1e-7);
    total = total+numel(ad);
end
agreement = same/total;
end

function ratio = pn_body_ratio(frames)
ratio = mean(abs(frames(:,1:945)).^2,2)./(mean(abs(frames(:,946:end)).^2,2)+eps);
end

function [islrDb,pslrDb] = autocorrelation_quality(frames)
n = size(frames,1);
islrDb = nan(n,1); pslrDb = nan(n,1);
for k = 1:n
    p = abs(xcorr(frames(k,:).')).^2;
    [peak,idx] = max(p);
    lag = (1:numel(p)).'-idx;
    roi = abs(lag)<=511; main = abs(lag)<=3; side = roi & ~main;
    islrDb(k) = 10*log10((mean(p(side))+eps)/(sum(p(main))+eps));
    pslrDb(k) = 10*log10((max(p(side))+eps)/(peak+eps));
end
end

function make_figures(inputFrames,fastFrames,fullFrames,fastAll,fullMask,T,R,outDir)
fig = figure('Visible','off','Color','w','Position',[80 80 1450 780]);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
nexttile;
bar(1:3,T.evm_median); grid on; ylabel('RMS EVM');
xticks(1:3); xticklabels(T.method); title('Matched FEC-valid frames');
nexttile;
yyaxis left; bar(1:2,T.waveform_nmse_db_median(2:3)); ylabel('NMSE / dB');
yyaxis right; plot(1:2,T.correlation_median(2:3),'ko-','LineWidth',1.5); ylim([0 1.02]);
ylabel('correlation'); grid on; xticks(1:2); xticklabels(T.method(2:3));
title('Waveform fit to equalized input');
nexttile;
plot(fastAll.nmse_db,'Color',[0.1 0.45 0.85]); grid on; hold on;
indices = find(fullMask); scatter(indices,fastAll.nmse_db(fullMask),10,'r','filled'); hold off;
xlabel('frame'); ylabel('Fast NMSE / dB'); title('Fast branch: all 1099 frames');
nexttile;
bar([R.full_profile_compute_ms_per_frame,R.fast_end_to_end_compute_ms_per_frame; ...
    R.full_interleaver_cache_ms,R.fast_interleaver_cache_ms].');
grid on; ylabel('ms');
xticks(1:2); xticklabels({'Full FEC','Fast'}); legend('compute/frame','buffer wait');
title('Computation and mandatory buffering');
exportgraphics(fig,fullfile(outDir,'fast_vs_full_summary.png'),'Resolution',180);
close(fig);

k = ceil(size(inputFrames,1)/2);
[a,~] = normalized_constellation(inputFrames(k,:));
[b,~] = normalized_constellation(fastFrames(k,:));
[c,~] = normalized_constellation(fullFrames(k,:));
fig = figure('Visible','off','Color','w','Position',[80 80 1450 450]);
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
constellation_panel(a,'Equalized input');
constellation_panel(b,'Fast decision reconstruction');
constellation_panel(c,'Full FEC reconstruction');
exportgraphics(fig,fullfile(outDir,'fast_vs_full_constellation.png'),'Resolution',180);
close(fig);
end

function [z,sliced] = normalized_constellation(frame)
fd = fftshift(fft(frame(946:4725).',3780));
z = DTMB_Frequency_Deinterleave(fd);
z = z/(rms(z)+eps);
phase = Estimate_CPE_64QAM(z);
z = z*exp(-1j*phase); z = z/(rms(z)+eps);
sliced = Slice_64QAM(z);
end

function constellation_panel(z,titleText)
nexttile; scatter(real(z),imag(z),5,'.'); axis equal; grid on;
xlim([-1.4 1.4]); ylim([-1.4 1.4]); xlabel('I'); ylabel('Q'); title(titleText);
end
