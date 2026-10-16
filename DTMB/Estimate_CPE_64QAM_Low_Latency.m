function [phase,info,hardSymbols] = Estimate_CPE_64QAM_Low_Latency(sym,targetPhase,cfg)
%ESTIMATE_CPE_64QAM_LOW_LATENCY Decision-directed CPE search for streaming.
%   The first frame uses a global pi/2-ambiguity search.  When TARGETPHASE
%   is finite, a single decision-directed correction is normally enough.
%   An unexpectedly large correction automatically triggers a global
%   reacquisition.  A local-grid tracking mode remains available for
%   diagnostics.

if nargin < 2 || isempty(targetPhase), targetPhase = NaN; end
if nargin < 3, cfg = struct(); end
inputSym = sym(:);
valid = isfinite(real(inputSym)) & isfinite(imag(inputSym));
sym = inputSym(valid);
if numel(sym) < 100
    phase = 0;
    info = struct('reacquired',false,'candidate_count',0, ...
        'metric',NaN,'full_metric',NaN,'search_mode',"insufficient_data");
    hardSymbols = Slice_64QAM(inputSym);
    return
end

globalPoints = get_cfg(cfg,'cpe_global_points',91);
localStep = deg2rad(get_cfg(cfg,'cpe_local_step_deg',1));
localRadius = get_cfg(cfg,'cpe_local_radius_steps',6);
reacquired = ~isfinite(targetPhase);
if reacquired
    sym = robust_trim(sym);
    candidates = linspace(-pi/4,pi/4,globalPoints);
    searchMode = "global";
    [phase,~,bestMetric] = search_candidates(sym,candidates);
    candidateCount = numel(candidates);
else
    trackingMethod = lower(string(get_cfg(cfg,'cpe_tracking_method','decision_directed')));
    if trackingMethod == "local_grid"
        localSym = robust_trim(sym);
        candidates = targetPhase + (-localRadius:localRadius)*localStep;
        [phase,bestIndex,bestMetric] = search_candidates(localSym,candidates);
        needsReacquire = get_cfg(cfg,'cpe_reacquire_on_edge',true) && ...
            (bestIndex == 1 || bestIndex == numel(candidates));
        candidateCount = numel(candidates);
        searchMode = "local_grid";
    else
        rotated = sym*exp(-1j*targetPhase);
        decided = Slice_64QAM(rotated);
        correction = angle(decided'*rotated);
        phase = targetPhase + correction;
        bestMetric = NaN;
        maxCorrection = deg2rad(get_cfg(cfg,'cpe_dd_max_correction_deg',8));
        needsReacquire = abs(correction) > maxCorrection;
        candidateCount = 1;
        searchMode = "decision_directed";
    end
    if needsReacquire
        globalCandidates = linspace(-pi/4,pi/4,globalPoints);
        [rawGlobal,~,bestMetric] = search_candidates(robust_trim(sym),globalCandidates);
        phase = rawGlobal + round((targetPhase-rawGlobal)/(pi/2))*(pi/2);
        reacquired = true;
        searchMode = "global_reacquire";
        candidateCount = candidateCount + numel(globalCandidates);
    end
end

corrected = inputSym*exp(-1j*phase);
hardSymbols = Slice_64QAM(corrected);
fullMetric = mean(abs(corrected(valid)-hardSymbols(valid)).^2);
info = struct('reacquired',reacquired,'candidate_count',candidateCount, ...
    'metric',bestMetric,'full_metric',fullMetric,'search_mode',searchMode);
end

function sym = robust_trim(sym)
% Match the validated 5--95 percent amplitude trimming without prctile.
amplitude = abs(sym);
sortedAmplitude = sort(amplitude);
n = numel(sortedAmplitude);
low = sortedAmplitude(max(1,ceil(0.05*n)));
high = sortedAmplitude(min(n,floor(0.95*n)));
use = amplitude > low & amplitude < high;
if nnz(use) >= 100, sym = sym(use); end
end

function [phase,bestIndex,bestMetric] = search_candidates(sym,candidates)
rotated = sym .* exp(-1j*candidates);
decided = Slice_64QAM(rotated);
metric = mean(abs(rotated-decided).^2,1);
[bestMetric,bestIndex] = min(metric);
phase = candidates(bestIndex);
end

function value = get_cfg(cfg,name,defaultValue)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
