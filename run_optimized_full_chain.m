function [result,report] = run_optimized_full_chain( ...
        equalizedPacket,syncPacket,outputDir,cfgOverride)
%RUN_OPTIMIZED_FULL_CHAIN Reproduce the throughput-optimized DTMB FEC path.
%   The input is the existing channel-equalized frame packet.  Protocol
%   buffering is unchanged; this entry optimizes computation only.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
if nargin < 1 || isempty(equalizedPacket)
    equalizedPacket = fullfile(root,'DTMB_groupwise_equalization', ...
        'DTMB_Groupwise_Equalized_Frames.mat');
end
if nargin < 2 || isempty(syncPacket)
    syncPacket = fullfile(root,'DTMB_1100Frame_Sync_Equalized_Packet.mat');
end
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(root,'DTMB_full_chain_optimization');
end
if nargin < 4, cfgOverride = struct(); end
if ~exist(outputDir,'dir'), mkdir(outputDir); end

G = load(equalizedPacket,'optimized_frames');
S = load(syncPacket,'PN_phase_index_packet','sync_diagnostics');
assert(isfield(G,'optimized_frames'), ...
    'Equalized packet must contain optimized_frames.');
assert(isfield(S,'PN_phase_index_packet'), ...
    'Synchronized packet must contain PN_phase_index_packet.');
if isfield(S,'sync_diagnostics') && isfield(S.sync_diagnostics,'cfg')
    cfg = S.sync_diagnostics.cfg;
else
    cfg = struct();
end
cfg.plot_diagnostics = false;
cfg.save_diagnostic_figures = false;
cfg.fec_search_all_commutator_phases = false;
cfg.fec_commutator_direction = 1;
cfg.fec_commutator_phase = 0;
cfg.ldpc_algorithm = 'norm-min-sum';
cfg.ldpc_min_sum_scaling_factor = 0.75;
cfg.ldpc_noise_model = 'deinterleaved';
cfg.ldpc_llr_scale = 0.25;
cfg.fec_cpe_method = 'decision_directed';
cfg.fec_demap_batch_frames = 16;
cfg.save_ldpc_llr_blocks = false;
cfg = merge_cfg(cfg,cfgOverride);
cfg.full_recon_output_packet = fullfile(outputDir, ...
    'DTMB_Full_Throughput_Optimized.mat');

tWall = tic;
result = run_full_demod_remod_section( ...
    G.optimized_frames,S.PN_phase_index_packet,cfg);
wallSeconds = toc(tWall);

profileCsv = fullfile(root,'DTMB_module_profile','module_timing_memory.csv');
baselineFecMs = NaN;
baselineAllMs = NaN;
if isfile(profileCsv)
    baseline = readtable(profileCsv);
    baselineAllMs = sum(baseline.per_frame_ms);
    exclude = ismember(string(baseline.module),["TC-AOMP","Equalization"]);
    baselineFecMs = sum(baseline.per_frame_ms(~exclude));
end

report = struct();
report.algorithm_version = "full_fec_throughput_v2_20260903";
report.num_input_frames = result.num_frames;
report.num_fec_valid_frames = nnz(result.fec_valid_frame_mask);
report.ldpc_parity_success_rate = result.fec_stats.parity_success_rate;
report.mean_ldpc_iterations = mean(result.fec_stats.iterations);
report.optimized_compute_ms_per_input_frame = ...
    result.timing.per_input_frame_ms;
report.optimized_wall_ms_per_input_frame = ...
    1e3*wallSeconds/result.num_frames;
report.baseline_equalized_input_ms_per_frame = baselineFecMs;
report.baseline_all_profiled_modules_ms_per_frame = baselineAllMs;
report.compute_speedup_vs_equalized_baseline = ...
    baselineFecMs/result.timing.per_input_frame_ms;
report.wall_seconds = wallSeconds;
report.stage_timing = result.timing;
report.time_interleaver_mode = result.system_info.interleaver_mode;
report.time_interleaver_M = result.system_info.M;
delayFrames = 52*51*result.system_info.M/3744;
report.mandatory_decode_reencode_wait_frames = 2*delayFrames+1;
report.mandatory_decode_reencode_wait_ms = ...
    0.625*report.mandatory_decode_reencode_wait_frames;

baselineFile = fullfile(root,'DTMB_groupwise_equalization', ...
    'DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat');
report.baseline_comparison_available = false;
if isfile(baselineFile)
    B = load(baselineFile,'result');
    n = min(result.num_frames,B.result.num_frames);
    common = result.fec_valid_frame_mask(1:n) & ...
        B.result.fec_valid_frame_mask(1:n);
    report.baseline_comparison_available = true;
    report.baseline_fec_valid_frames = nnz(B.result.fec_valid_frame_mask(1:n));
    report.common_fec_valid_frames = nnz(common);
    report.fec_valid_mask_exact = isequal( ...
        result.fec_valid_frame_mask(1:n),B.result.fec_valid_frame_mask(1:n));
    a = result.canonical_reconstructed_frames(common,:);
    b = B.result.canonical_reconstructed_frames(common,:);
    report.fec_valid_reconstruction_exact = isequal(a,b);
    report.fec_valid_reconstruction_correlation = ...
        abs(double(a(:))'*double(b(:))) / ...
        (norm(double(a(:)))*norm(double(b(:)))+eps);
    allA = double(result.canonical_reconstructed_frames(1:n,:));
    allB = double(B.result.canonical_reconstructed_frames(1:n,:));
    report.all_frame_reconstruction_exact = isequal(allA,allB);
    report.all_frame_reconstruction_correlation = ...
        abs(allA(:)'*allB(:))/(norm(allA(:))*norm(allB(:))+eps);
end

summary = table(report.num_input_frames,report.num_fec_valid_frames, ...
    report.ldpc_parity_success_rate, ...
    report.optimized_compute_ms_per_input_frame, ...
    report.optimized_wall_ms_per_input_frame, ...
    report.baseline_equalized_input_ms_per_frame, ...
    report.compute_speedup_vs_equalized_baseline, ...
    report.mandatory_decode_reencode_wait_ms, ...
    'VariableNames',{'input_frames','fec_valid_frames', ...
    'ldpc_parity_success_rate','optimized_compute_ms_per_frame', ...
    'optimized_wall_ms_per_frame','baseline_compute_ms_per_frame', ...
    'compute_speedup','mandatory_protocol_wait_ms'});
report.summary = summary;
writetable(summary,fullfile(outputDir,'full_chain_optimization_summary.csv'));
save(fullfile(outputDir,'Full_Chain_Optimization_Report.mat'),'report');
disp(summary);
end

function cfg = merge_cfg(cfg,override)
names = fieldnames(override);
for k = 1:numel(names), cfg.(names{k}) = override.(names{k}); end
end
