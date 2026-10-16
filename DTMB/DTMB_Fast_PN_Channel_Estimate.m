function [h,info] = DTMB_Fast_PN_Channel_Estimate(pnObserved,phaseIndex,cfg)
%DTMB_FAST_PN_CHANNEL_ESTIMATE Fixed-budget PN945 sparse channel estimate.
%   This low-latency OMP variant omits the full-frame equalization and
%   autocorrelation calculation previously performed at every AOMP
%   iteration.  A small fixed tap budget gives deterministic execution.

if nargin < 3, cfg = struct(); end
tapLength = get_cfg(cfg,'channel_impulse_response_length',160);
numTaps = get_cfg(cfg,'fast_channel_num_taps',4);
numTaps = min(max(1,round(numTaps)),tapLength);
method = lower(string(get_cfg(cfg,'fast_channel_method','threshold_ls')));

observed = single(pnObserved(:));
assert(numel(observed) >= 945,'PN observation must contain 945 symbols.');
localPn = single(PN_gen(phaseIndex)).';
y = observed(tapLength:945);
theta = single(Build_Sensing_Matrix_TC(localPn,tapLength));

inputEnergy = sum(abs(y).^2)+eps('single');
if method == "threshold_ls"
    % PN columns are close to orthogonal. One matched-filter pass followed
    % by a small LS fit avoids K full projection passes in classical OMP.
    projection = abs(theta'*y);
    [~,ordered] = sort(projection,'descend');
    support = ordered(1:numTaps).';
    omega = theta(:,support);
    coefficients = omega\y;
    residual = y-omega*coefficients;
    residualRatio = sum(abs(residual).^2)/inputEnergy;
elseif method == "omp"
    residual = y;
    support = zeros(1,numTaps);
    selected = false(tapLength,1);
    coefficients = complex(zeros(0,1,'single'));
    residualRatio = zeros(numTaps,1);
    for k = 1:numTaps
        projection = abs(theta'*residual);
        projection(selected) = -Inf;
        [~,support(k)] = max(projection);
        selected(support(k)) = true;
        omega = theta(:,support(1:k));
        coefficients = omega\y;
        residual = y-omega*coefficients;
        residualRatio(k) = sum(abs(residual).^2)/inputEnergy;
    end
else
    error('Unknown fast_channel_method: %s',method);
end

h = complex(zeros(tapLength,1,'single'));
h(support) = coefficients;
info = struct('support',support,'num_taps',numTaps,'method',method, ...
    'residual_energy_ratio',residualRatio(end), ...
    'residual_energy_history',residualRatio);
end

function value = get_cfg(cfg,name,defaultValue)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
