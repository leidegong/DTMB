function report = test_low_latency_demod_remod()
%TEST_LOW_LATENCY_DEMOD_REMOD Numerical, continuity and dependency checks.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
rng(20260902);

% Arithmetic slicer must match the original nearest-level definition away
% from exact boundaries.
x = 2*(randn(10000,1)+1j*randn(10000,1));
levels = [-7 -5 -3 -1 1 3 5 7];
u = x*sqrt(42);
[~,i] = min(abs(real(u)-levels),[],2);
[~,q] = min(abs(imag(u)-levels),[],2);
reference = (levels(i).'+1j*levels(q).')/sqrt(42);
assert(isequal(Slice_64QAM(x),reference),'Fast 64QAM slicer mismatch.');

% Validate the stateful path on the bundled synchronized-frame sample.
inputFile = fullfile(root,'fast_16frame_sync_input.mat');
S = load(inputFile,'Frame_equalization_packet','PN_phase_index_packet');
cfg = struct('save_shaped_waveform',true,'flush_srrc_tail',true);
result = run_low_latency_demod_remod( ...
    S.Frame_equalization_packet,S.PN_phase_index_packet,cfg);
serial = reshape(result.reconstructed_frames.',[],1);
upsampled = complex(zeros(numel(serial)*result.shaping_sps,1,'single'));
upsampled(1:result.shaping_sps:end) = serial;
expected = conv(upsampled,result.shaping_filter,'full');
srrcMaxError = max(abs(expected-result.shaped_reference));
assert(srrcMaxError < 5e-6,'Stateful SRRC is not continuous across frames.');

% Compare decision quality with the validated full-grid implementation.
legacyCfg = struct('phase_continuity_tracking',true, ...
    'tps_sign_disambiguation',true,'qam_phase_modulo_unwrap',true, ...
    'save_shaped_waveform',false);
legacy = run_fast_decision_reconstruction( ...
    S.Frame_equalization_packet,S.PN_phase_index_packet,legacyCfg);
newEvm = median(result.decision_evm_rms,'omitnan');
legacyEvm = median(legacy.decision_evm_rms,'omitnan');
assert(newEvm <= legacyEvm+0.005, ...
    'Low-latency CPE tracking degraded median EVM unexpectedly.');
dependencySummary = verify_low_latency_demod_remod_final(result);

report = struct('passed',true,'num_frames',result.num_frames, ...
    'srrc_max_abs_error',srrcMaxError,'new_evm_median',newEvm, ...
    'legacy_evm_median',legacyEvm, ...
    'core_compute_ms_median',result.timing.compute_ms_per_frame_median, ...
    'core_compute_ms_p95',result.timing.compute_ms_per_frame_p95, ...
    'steady_compute_ms_median',result.timing.steady_compute_ms_per_frame_median, ...
    'steady_compute_ms_p95',result.timing.steady_compute_ms_per_frame_p95, ...
    'dependency_summary',dependencySummary);
fprintf(['Low-latency demod/remod tests passed: SRRC error %.3g, ' ...
    'EVM %.6f (legacy %.6f).\n'],srrcMaxError,newEvm,legacyEvm);
end
