function info = DTMB_Low_Latency_Prewarm(cfg)
%DTMB_LOW_LATENCY_PREWARM Initialize JIT/FFT/FIR paths before acquisition.
%   The generated states are deliberately discarded. Calling this during
%   receiver initialization moves MATLAB one-time dispatch costs outside
%   the first live DTMB frame without adding signal-buffer latency.

if nargin < 1, cfg = struct(); end
t = tic;
phaseIndex = 0;
pn = single(PN_gen(phaseIndex));
[zToY,tpsIndex] = DTMB_Frequency_Permutation();
dataMask = true(3780,1);
dataMask(tpsIndex) = false;
candidates = DTMB_System_Info_Candidates();
z = complex(zeros(3780,1,'single'));
z(tpsIndex) = single(candidates(1).symbols);
qamPattern = single((repmat([-7 -3 1 5],1,936)+ ...
    1j*repelem([-5 -1 3 7],936)).');
z(dataMask) = qamPattern;
y = complex(zeros(3780,1,'single'));
y(zToY) = z;
body = sqrt(single(3780/get_cfg(cfg,'pn_to_body_power_ratio',2))) * ...
    ifft(y,3780);
dummyFrame = [pn(:);body(:)];
DTMB_Fast_PN_Channel_Estimate(pn,phaseIndex,cfg);
[equalized,~,~] = DTMB_Low_Latency_Channel_Equalize_Frame( ...
    pn,dummyFrame,phaseIndex,[],cfg);
coreState = [];
for k = 1:3
    [~,coreState] = DTMB_Low_Latency_Demod_Remod_Frame( ...
        equalized,phaseIndex,coreState,cfg);
end
info = struct('elapsed_s',toc(t),'states_discarded',true);
end

function value = get_cfg(cfg,name,defaultValue)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
