function [tail,state] = DTMB_Low_Latency_Flush(state)
%DTMB_LOW_LATENCY_FLUSH Emit the finite SRRC tail at end of a burst.
assert(isstruct(state) && isfield(state,'srrc_zi'), ...
    'A state returned by DTMB_Low_Latency_Demod_Remod_Frame is required.');
n = numel(state.shaping_filter)-1;
[tail,state.srrc_zi] = filter(state.shaping_filter,1, ...
    complex(zeros(n,1,'single')),state.srrc_zi);
end
