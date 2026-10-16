function summary = verify_low_latency_demod_remod_final(resultInput)
%VERIFY_LOW_LATENCY_DEMOD_REMOD_FINAL Audit dependencies and stream output.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
if isfolder(fullfile(fileparts(root),'example'))
    artifactRoot = fileparts(root);
else
    artifactRoot = root;
end
if nargin < 1 || isempty(resultInput)
    if strcmp(artifactRoot,root)
        resultInput = fullfile(root,'DTMB_low_latency_final_output', ...
            'DTMB_Low_Latency_Demod_Remod_Final.mat');
    else
        resultInput = fullfile(artifactRoot,'outputs', ...
            'DTMB_Low_Latency_Demod_Remod_Final.mat');
    end
end
if isstruct(resultInput)
    result = resultInput;
else
    S = load(resultInput,'result');
    result = S.result;
end

entryFiles = {fullfile(root,'run_low_latency_reference_pipeline.m'), ...
    fullfile(root,'DTMB','DTMB_Low_Latency_Demod_Remod_Frame.m')};
[files,~] = matlab.codetools.requiredFilesAndProducts(entryFiles);
names = lower(string(files));
forbidden = ["apply_482mhz","purifier","neural","tcn","residual", ...
    "ldpc","bch","time_deinterleave","time_interleave", ...
    "run_full_demod_remod"];
hit = false(size(names));
for k = 1:numel(forbidden), hit = hit | contains(names,forbidden(k)); end
assert(~any(hit),'Forbidden model/FEC dependency detected: %s', ...
    strjoin(names(hit),', '));
assert(~result.uses_inverse_channel_model, ...
    'Final result unexpectedly declares an inverse channel model.');
assert(result.timing.time_interleaver_buffer_frames == 0);
assert(result.timing.minimum_input_buffer_frames == 1);
if isfield(result.timing,'total_startup_buffer_frames')
    assert(result.timing.total_startup_buffer_frames == 1, ...
        'Final streaming frontend must emit the first received frame.');
end
if isfield(result,'runtime_configuration')
    assert(~result.runtime_configuration.uses_fec);
end

if ~isempty(result.shaped_reference)
    baseLength = result.num_frames*4725*result.shaping_sps;
    validLengths = [baseLength,baseLength+numel(result.shaping_filter)-1];
    assert(any(numel(result.shaped_reference) == validLengths), ...
        'Stateful SRRC output length is inconsistent.');
end
summary = table(result.num_frames,numel(files), ...
    result.uses_inverse_channel_model, ...
    result.timing.time_interleaver_buffer_frames, ...
    result.timing.minimum_input_buffer_frames, ...
    result.timing.steady_compute_ms_per_frame_median, ...
    'VariableNames',{'num_frames','dependency_files', ...
    'uses_inverse_channel_model','time_interleaver_buffer_frames', ...
    'minimum_input_buffer_frames','steady_core_compute_ms_median'});
disp(summary);
end
