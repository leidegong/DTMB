function report = validate_low_latency_full_length(syncPacket,goldPacket, ...
        outputDir,numTimingRepeats)
%VALIDATE_LOW_LATENCY_FULL_LENGTH Run the final low-latency entry over the
%   whole 1099-frame capture instead of the bundled 16-frame sample.
%   QUICKSTART_LOW_LATENCY_FINAL stays the fast smoke test; this entry is the
%   full-length validation.  It reports steady-state timing over a long run,
%   the real-time margin against the 0.625 ms frame duration, and waveform
%   agreement with the complete-FEC reconstruction.  The complete link is used
%   only as offline truth here, never at runtime.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
if nargin < 1 || isempty(syncPacket)
    syncPacket = fullfile(root,'DTMB_1100Frame_Sync_Equalized_Packet.mat');
end
if nargin < 2 || isempty(goldPacket)
    goldPacket = fullfile(root,'DTMB_groupwise_equalization', ...
        'DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat');
end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(root,'DTMB_low_latency_full_length_validation');
end
if nargin < 4 || isempty(numTimingRepeats), numTimingRepeats = 3; end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

% Same runtime configuration as QUICKSTART_LOW_LATENCY_FINAL.
cfg = struct('tc_channel_estimation_frames',2, ...
    'streaming_channel_bootstrap',true, ...
    'fast_channel_num_taps',32, ...
    'fast_channel_method','threshold_ls', ...
    'cpe_update_interval_frames',4, ...
    'cpe_tracking_method','decision_directed', ...
    'save_shaped_waveform',true,'save_result',true);
result = run_low_latency_reference_pipeline(syncPacket,outputDir,cfg);

% Dependency and streaming audit on the full-length result.
dependencySummary = verify_low_latency_demod_remod_final(result);

frameDurationMs = result.timing.minimum_input_buffer_ms;
coreMs = result.timing.compute_ms_per_frame(:);
warmupFrames = result.timing.performance_warmup_frames;
steadyCore = coreMs(min(warmupFrames+1,numel(coreMs)):end);

report = struct();
report.validation_version = "low_latency_full_length_v1_20260908";
report.num_frames = result.num_frames;
report.frame_duration_ms = frameDurationMs;
report.frontend_mode = string(result.frontend.mode);
report.uses_inverse_channel_model = result.uses_inverse_channel_model;
report.dependency_summary = dependencySummary;

report.steady_core_ms_median = median(steadyCore);
report.steady_core_ms_mean = mean(steadyCore);
report.steady_core_ms_p95 = prctile(steadyCore,95);
report.steady_core_ms_max = max(steadyCore);
report.channel_ms_median = result.timing.channel_compute_ms_per_frame_median;
report.channel_ms_p95 = result.timing.channel_compute_ms_per_frame_p95;
report.channel_ms_mean = ...
    1e3*result.timing.channel_frontend_compute_s/result.num_frames;
report.steady_end_to_end_ms_median = ...
    report.steady_core_ms_median + report.channel_ms_median;
% Throughput is decided by the mean, not the median: a stream keeps up only
% if the average per-frame cost stays under the frame duration.  The median
% can sit under budget while the tail pushes the mean over it.
%
% The per-frame core/channel figures above come from inner timers that
% exclude loop and allocation overhead, and the two are not measured the same
% way, so they must not be summed into a throughput verdict.  The real-time
% factor is taken from the pipeline wall clock, which covers the channel
% front-end and the core on one clock and charges every overhead once.
report.pipeline_wall_ms_per_frame = ...
    1e3*result.timing.pipeline_wall_s/result.num_frames;
report.realtime_factor_wall = ...
    report.pipeline_wall_ms_per_frame/frameDurationMs;
report.keeps_up_on_average = report.realtime_factor_wall < 1;
% The two clocks bracket the true streaming cost and neither is it.  The inner
% timers exclude the harness (full-capture preallocation, writing every frame
% into a 1099-frame array, profiling scaffolding) that a streaming deployment
% would not pay, so they are a lower bound.  The pipeline wall clock charges
% all of it, so it is an upper bound.  A verdict on real-time operation needs
% a streaming benchmark, not this batch harness.
report.realtime_factor_lower_bound = ...
    report.steady_end_to_end_ms_median/frameDurationMs;
report.realtime_factor_upper_bound = report.realtime_factor_wall;
if report.realtime_factor_upper_bound < 1
    report.realtime_verdict = "keeps_up";
elseif report.realtime_factor_lower_bound >= 1
    report.realtime_verdict = "does_not_keep_up";
else
    report.realtime_verdict = "inconclusive_harness_bounded";
end
report.first_frame_end_to_end_ms = ...
    result.timing.first_frame_end_to_end_compute_ms;
report.minimum_signal_latency_ms = result.timing.minimum_signal_latency_ms;
report.time_interleaver_buffer_frames = ...
    result.timing.time_interleaver_buffer_frames;

% Real-time margin: how often a single frame misses its own frame budget.
report.core_over_budget_frames = nnz(steadyCore > frameDurationMs);
report.core_over_budget_fraction = ...
    report.core_over_budget_frames/numel(steadyCore);
report.core_realtime_margin_percent = ...
    100*(1-report.steady_core_ms_median/frameDurationMs);
report.end_to_end_realtime_margin_percent = ...
    100*(1-report.steady_end_to_end_ms_median/frameDurationMs);

% Repeat the timed run: single-run MATLAB timings on a loaded desktop move by
% more than 10%, so the real-time verdict must be read from the spread across
% runs rather than from one measurement.  Repeats skip the large saves.
repeatCfg = cfg;
repeatCfg.save_result = false;   % disk write only; it is outside pipeline_wall_s
repeatCore = zeros(numTimingRepeats,1);
repeatWall = zeros(numTimingRepeats,1);
repeatCore(1) = report.steady_core_ms_mean;
repeatWall(1) = report.pipeline_wall_ms_per_frame;
for k = 2:numTimingRepeats
    fprintf('[full-length] Timing repeat %d/%d.\n',k,numTimingRepeats);
    repeatResult = run_low_latency_reference_pipeline( ...
        syncPacket,outputDir,repeatCfg);
    repeatWarmup = repeatResult.timing.performance_warmup_frames;
    repeatSteady = repeatResult.timing.compute_ms_per_frame( ...
        min(repeatWarmup+1,repeatResult.num_frames):end);
    repeatCore(k) = mean(repeatSteady);
    repeatWall(k) = ...
        1e3*repeatResult.timing.pipeline_wall_s/repeatResult.num_frames;
end
report.timing_repeats = numTimingRepeats;
report.repeat_core_ms_mean = repeatCore;
report.repeat_wall_ms_per_frame = repeatWall;
report.repeat_wall_ms_min = min(repeatWall);
report.repeat_wall_ms_max = max(repeatWall);
report.repeat_realtime_factor_min = min(repeatWall)/frameDurationMs;
report.repeat_realtime_factor_max = max(repeatWall)/frameDurationMs;
report.keeps_up_in_every_repeat = report.repeat_realtime_factor_max < 1;

report.decision_evm_median = median(result.decision_evm_rms,'omitnan');
report.decision_evm_p95 = prctile(result.decision_evm_rms,95);
report.cpe_reacquisitions_after_lock = nnz(result.cpe_reacquired(2:end));

% Waveform agreement with the complete-FEC reconstruction (offline truth).
report.gold_comparison_available = false;
if isfile(goldPacket)
    G = load(goldPacket,'result');
    n = min(result.num_frames,G.result.num_frames);
    mask = false(n,1);
    mask(logical(G.result.fec_valid_frame_mask(1:n))) = true;
    report.gold_comparison_available = true;
    report.gold_fec_valid_frames = nnz(mask);
    [report.nmse_to_full_fec_db,report.correlation_to_full_fec] = ...
        match_metrics(result.canonical_reconstructed_frames(mask,:), ...
        G.result.canonical_reconstructed_frames(mask,:));
    [report.nmse_all_frames_db,report.correlation_all_frames] = ...
        match_metrics(result.canonical_reconstructed_frames(1:n,:), ...
        G.result.canonical_reconstructed_frames(1:n,:));
end

% Compare against the bundled 16-frame smoke run, when it is present.
smokeCsv = fullfile(root,'DTMB_low_latency_final_output', ...
    'low_latency_final_summary.csv');
report.smoke_comparison_available = false;
if isfile(smokeCsv)
    smoke = readtable(smokeCsv);
    report.smoke_comparison_available = true;
    report.smoke_num_frames = smoke.num_frames(1);
    report.smoke_steady_core_ms_median = smoke.steady_core_compute_ms_median(1);
    report.smoke_decision_evm_median = smoke.decision_evm_median(1);
end

summary = table(report.num_frames,report.frame_duration_ms, ...
    report.steady_core_ms_median,report.steady_core_ms_mean, ...
    report.steady_core_ms_p95,report.steady_core_ms_max, ...
    report.steady_end_to_end_ms_median,report.pipeline_wall_ms_per_frame, ...
    report.realtime_factor_lower_bound,report.realtime_factor_upper_bound, ...
    report.repeat_realtime_factor_min,report.repeat_realtime_factor_max, ...
    report.realtime_verdict, ...
    report.core_over_budget_frames,report.core_realtime_margin_percent, ...
    report.end_to_end_realtime_margin_percent, ...
    report.minimum_signal_latency_ms,report.decision_evm_median, ...
    report.decision_evm_p95,report.cpe_reacquisitions_after_lock, ...
    'VariableNames',{'num_frames','frame_duration_ms', ...
    'steady_core_ms_median','steady_core_ms_mean', ...
    'steady_core_ms_p95','steady_core_ms_max', ...
    'steady_end_to_end_ms_median','pipeline_wall_ms_per_frame', ...
    'realtime_factor_lower_bound','realtime_factor_upper_bound', ...
    'repeat_realtime_factor_min','repeat_realtime_factor_max', ...
    'realtime_verdict', ...
    'core_over_budget_frames', ...
    'core_realtime_margin_percent','end_to_end_realtime_margin_percent', ...
    'minimum_signal_latency_ms','decision_evm_median','decision_evm_p95', ...
    'cpe_reacquisitions_after_lock'});
report.summary = summary;
writetable(summary,fullfile(outputDir,'low_latency_full_length_summary.csv'));
save(fullfile(outputDir,'Low_Latency_Full_Length_Report.mat'),'report');
make_timing_figure(coreMs,warmupFrames,frameDurationMs,outputDir);

disp(summary);
if report.gold_comparison_available
    fprintf(['[全长验证] 对完整FEC重构：%d 个FEC有效帧 NMSE=%.3f dB，' ...
        '相关=%.6f；全部 %d 帧相关=%.6f。\n'],report.gold_fec_valid_frames, ...
        report.nmse_to_full_fec_db,report.correlation_to_full_fec, ...
        result.num_frames,report.correlation_all_frames);
end
fprintf('[全长验证] 超出单帧预算(%.4f ms)的稳态帧数：%d / %d。\n', ...
    frameDurationMs,report.core_over_budget_frames,numel(steadyCore));
fprintf(['[全长验证] %d 次重复的整链墙钟 %.4f--%.4f ms/帧，' ...
    '实时因子 %.4f--%.4f，%s次重复都跟得上。\n'],numTimingRepeats, ...
    report.repeat_wall_ms_min,report.repeat_wall_ms_max, ...
    report.repeat_realtime_factor_min,report.repeat_realtime_factor_max, ...
    ternary(report.keeps_up_in_every_repeat,'每','并非每'));
fprintf('[全长验证] 实时判据 %.4f--%.4f（内层计时下界--墙钟上界），结论：%s。\n', ...
    report.realtime_factor_lower_bound,report.realtime_factor_upper_bound, ...
    report.realtime_verdict);
end

function out = ternary(condition,trueValue,falseValue)
if condition, out = trueValue; else, out = falseValue; end
end

function [nmseDb,correlation] = match_metrics(testFrames,goldFrames)
a = double(testFrames(:));
b = double(goldFrames(:));
gain = (a'*b)/(a'*a+eps);
residual = gain*a-b;
nmseDb = 10*log10(sum(abs(residual).^2)/(sum(abs(b).^2)+eps));
correlation = abs(a'*b)/(norm(a)*norm(b)+eps);
end

function make_timing_figure(coreMs,warmupFrames,frameDurationMs,outputDir)
fig = figure('Visible','off','Position',[100 100 1000 380]);
tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
plot(coreMs,'.','MarkerSize',4);
yline(frameDurationMs,'r--',sprintf('frame budget %.4f ms',frameDurationMs), ...
    'FontSize',8);
xlabel('frame index'); ylabel('core compute (ms)');
title('Per-frame core compute over the full capture');

nexttile; hold on; grid on;
histogram(coreMs(min(warmupFrames+1,numel(coreMs)):end),60);
xline(frameDurationMs,'r--','FontSize',8);
xlabel('core compute (ms)'); ylabel('frames');
title('Steady-state distribution (warm-up excluded)');

print(fig,fullfile(outputDir,'low_latency_full_length_timing.png'),'-dpng','-r150');
close(fig);
end
