function report = sweep_pn_body_power_ratio(fullChainPacket,equalizedPacket, ...
        syncPacket,outputDir)
%SWEEP_PN_BODY_POWER_RATIO Explain PSLR as a function of the PN:body power ratio.
%   The frame body of the existing full-FEC reconstruction is held fixed; only
%   the PN945 header amplitude is rescaled in memory.  EVM, NMSE and every
%   symbol decision are therefore untouched, so any PSLR/ISLR movement is
%   attributable to the header scaling alone.
%
%   The sweep also decomposes the autocorrelation into PN-only and body-only
%   contributions, which identifies where the peak sidelobe actually comes from.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
if nargin < 1 || isempty(fullChainPacket)
    fullChainPacket = fullfile(root,'DTMB_full_chain_optimization', ...
        'DTMB_Full_Throughput_Optimized.mat');
end
if nargin < 2 || isempty(equalizedPacket)
    equalizedPacket = fullfile(root,'DTMB_groupwise_equalization', ...
        'DTMB_Groupwise_Equalized_Frames.mat');
end
if nargin < 3 || isempty(syncPacket)
    syncPacket = fullfile(root,'DTMB_1100Frame_Sync_Equalized_Packet.mat');
end
if nargin < 4 || isempty(outputDir)
    outputDir = fullfile(root,'DTMB_pn_power_ratio_calibration');
end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

F = load(fullChainPacket,'result');
G = load(equalizedPacket,'optimized_frames');
S = load(syncPacket,'Frame_equalization_packet');
validMask = logical(F.result.fec_valid_frame_mask);
frameIdx = find(validMask);
assert(~isempty(frameIdx),'The full-chain packet contains no FEC-valid frames.');

reconstructed = double(F.result.canonical_reconstructed_frames(frameIdx,:));
groupwiseInput = double(G.optimized_frames(frameIdx,:));
syncInput = double(S.Frame_equalization_packet(frameIdx,:));

% ---- Reference points: the three signals exactly as they are stored -------
referenceNames = ["full_fec_reconstruction";"groupwise_equalized_input"; ...
    "sync_stage_equalized_input"];
referenceData = {reconstructed,groupwiseInput,syncInput};
refRatio = zeros(numel(referenceData),1);
refIslr = zeros(numel(referenceData),1);
refPslr = zeros(numel(referenceData),1);
refLag = zeros(numel(referenceData),1);
for k = 1:numel(referenceData)
    refRatio(k) = median(pn_body_ratio(referenceData{k}));
    [islr,pslr,lag] = autocorrelation_metrics(referenceData{k});
    refIslr(k) = median(islr,'omitnan');
    refPslr(k) = median(pslr,'omitnan');
    refLag(k) = median(lag,'omitnan');
end
referenceTable = table(referenceNames,refRatio,10*log10(refRatio), ...
    refIslr,refPslr,refLag,'VariableNames',{'signal','pn_body_power_ratio', ...
    'pn_body_power_ratio_db','islr_db_median','pslr_db_median', ...
    'peak_sidelobe_lag_symbols'});

% ---- Header / body decomposition of the reconstructed frame --------------
headerOnly = [reconstructed(:,1:945),zeros(size(reconstructed,1),3780)];
bodyOnly = [zeros(size(reconstructed,1),945),reconstructed(:,946:end)];
[islrPn,pslrPn,lagPn] = autocorrelation_metrics(headerOnly);
[islrBody,pslrBody,lagBody] = autocorrelation_metrics(bodyOnly);
decompositionTable = table(["pn945_header_only";"frame_body_only"], ...
    [median(islrPn,'omitnan');median(islrBody,'omitnan')], ...
    [median(pslrPn,'omitnan');median(pslrBody,'omitnan')], ...
    [median(lagPn,'omitnan');median(lagBody,'omitnan')], ...
    'VariableNames',{'component','islr_db_median','pslr_db_median', ...
    'peak_sidelobe_lag_symbols'});

% ---- Sweep: body fixed, header rescaled ----------------------------------
sweepRatios = unique([logspace(log10(0.05),log10(8),33),refRatio(:).',2]);
numRatios = numel(sweepRatios);
islrMedian = zeros(numRatios,1);
pslrMedian = zeros(numRatios,1);
pslrP05 = zeros(numRatios,1);
pslrP95 = zeros(numRatios,1);
lagMedian = zeros(numRatios,1);
header = reconstructed(:,1:945);
body = reconstructed(:,946:end);
headerPower = mean(abs(header).^2,2);
bodyPower = mean(abs(body).^2,2);
for k = 1:numRatios
    alpha = sqrt(sweepRatios(k)*bodyPower./(headerPower+eps));
    [islr,pslr,lag] = autocorrelation_metrics([alpha.*header,body]);
    islrMedian(k) = median(islr,'omitnan');
    pslrMedian(k) = median(pslr,'omitnan');
    pslrP05(k) = prctile(pslr,5);
    pslrP95(k) = prctile(pslr,95);
    lagMedian(k) = median(lag,'omitnan');
end
sweepTable = table(sweepRatios(:),10*log10(sweepRatios(:)),islrMedian, ...
    pslrMedian,pslrP05,pslrP95,lagMedian, ...
    'VariableNames',{'pn_body_power_ratio','pn_body_power_ratio_db', ...
    'islr_db_median','pslr_db_median','pslr_db_p05','pslr_db_p95', ...
    'peak_sidelobe_lag_symbols'});

% ---- The comparison the historical fast-vs-full table got wrong ----------
matchedRatio = refRatio(strcmp(referenceNames,'groupwise_equalized_input'));
alpha = sqrt(matchedRatio*bodyPower./(headerPower+eps));
[islrMatched,pslrMatched] = autocorrelation_metrics([alpha.*header,body]);

report = struct();
report.analysis_version = "pn_body_power_ratio_calibration_v1_20260908";
report.num_frames = numel(frameIdx);
report.frame_index_first = frameIdx(1);
report.frame_index_last = frameIdx(end);
report.reference_table = referenceTable;
report.decomposition_table = decompositionTable;
report.sweep_table = sweepTable;
report.groupwise_input_pn_body_ratio = matchedRatio;
report.reconstruction_pslr_at_groupwise_ratio_db = median(pslrMatched,'omitnan');
report.reconstruction_islr_at_groupwise_ratio_db = median(islrMatched,'omitnan');
report.reconstruction_minus_sync_input_pslr_db = ...
    refPslr(strcmp(referenceNames,'full_fec_reconstruction')) - ...
    refPslr(strcmp(referenceNames,'sync_stage_equalized_input'));

writetable(referenceTable,fullfile(outputDir,'pn_ratio_reference_points.csv'));
writetable(decompositionTable,fullfile(outputDir,'pn_ratio_decomposition.csv'));
writetable(sweepTable,fullfile(outputDir,'pn_ratio_sweep.csv'));
make_sweep_figure(sweepTable,referenceTable, ...
    median(pslrBody,'omitnan'),outputDir);
save(fullfile(outputDir,'PN_Power_Ratio_Calibration_Report.mat'),'report');

fprintf('\n[PN功率比标定] FEC有效帧 %d 帧（第 %d--%d 帧）\n', ...
    report.num_frames,report.frame_index_first,report.frame_index_last);
disp(referenceTable);
disp(decompositionTable);
fprintf('[PN功率比标定] 帧体不变、PN比值改到均衡输入水平 %.4f 后 PSLR=%.3f dB。\n', ...
    matchedRatio,report.reconstruction_pslr_at_groupwise_ratio_db);
fprintf('[PN功率比标定] 同比值下重构相对同步级均衡输入的 PSLR 差为 %.3f dB。\n', ...
    report.reconstruction_minus_sync_input_pslr_db);
end

function ratio = pn_body_ratio(frames)
ratio = mean(abs(frames(:,1:945)).^2,2)./(mean(abs(frames(:,946:end)).^2,2)+eps);
end

function [islrDb,pslrDb,lagAtMax] = autocorrelation_metrics(frames)
% Identical ROI and main-lobe definition to
% evaluate_DTMB_reconstruction_effect and run_fast_reconstruction_experiment,
% with the lag of the peak sidelobe added.
numFrames = size(frames,1);
islrDb = nan(numFrames,1);
pslrDb = nan(numFrames,1);
lagAtMax = nan(numFrames,1);
mainHalfWidth = 3;
roiHalfWidth = 511;
for frameIdx = 1:numFrames
    powerCorr = abs(xcorr(frames(frameIdx,:).')).^2;
    [peakPower,peakIdx] = max(powerCorr);
    lag = (1:numel(powerCorr)).'-peakIdx;
    roi = abs(lag) <= roiHalfWidth;
    main = abs(lag) <= mainHalfWidth;
    side = roi & ~main;
    islrDb(frameIdx) = 10*log10((mean(powerCorr(side))+eps)/ ...
        (sum(powerCorr(main))+eps));
    [sidePeak,sideIdx] = max(powerCorr(side));
    pslrDb(frameIdx) = 10*log10((sidePeak+eps)/(peakPower+eps));
    sideLag = lag(side);
    lagAtMax(frameIdx) = sideLag(sideIdx);
end
end

function make_sweep_figure(sweepTable,referenceTable,bodyOnlyPslrDb,outputDir)
fig = figure('Visible','off','Position',[100 100 1000 400]);
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
plot(sweepTable.pn_body_power_ratio_db,sweepTable.pslr_db_median,'-o', ...
    'LineWidth',1.4,'MarkerSize',4);
yline(bodyOnlyPslrDb,'k--','body-only PSLR floor','FontSize',7, ...
    'LabelHorizontalAlignment','left');
plot(referenceTable.pn_body_power_ratio_db,referenceTable.pslr_db_median, ...
    'rs','MarkerSize',9,'LineWidth',1.4);
labelOffset = [2.6 -2.0 -3.6];
for k = 1:height(referenceTable)
    text(referenceTable.pn_body_power_ratio_db(k)-0.4, ...
        referenceTable.pslr_db_median(k)+labelOffset(min(k,end)), ...
        strrep(char(referenceTable.signal(k)),'_','\_'), ...
        'FontSize',7,'HorizontalAlignment','right');
end
xlim([-15 12]);
xlabel('PN:body power ratio (dB)'); ylabel('PSLR (dB)');
title('PSLR vs PN header scaling (body fixed)');

nexttile; hold on; grid on;
plot(sweepTable.pn_body_power_ratio_db,sweepTable.islr_db_median,'-o', ...
    'LineWidth',1.4,'MarkerSize',4);
xlabel('PN:body power ratio (dB)'); ylabel('ISLR (dB)');
title('ISLR vs PN header scaling (body fixed)');

print(fig,fullfile(outputDir,'pn_power_ratio_sweep.png'),'-dpng','-r150');
close(fig);
end
