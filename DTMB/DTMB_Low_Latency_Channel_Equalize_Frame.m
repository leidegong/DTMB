function [equalized,state,info] = DTMB_Low_Latency_Channel_Equalize_Frame( ...
    pnObserved,receivedFrame,phaseIndex,state,cfg)
%DTMB_LOW_LATENCY_CHANNEL_EQUALIZE_FRAME Immediate-output channel frontend.
%   Frame one is equalized with its own PN estimate and emitted immediately.
%   The estimate is refined as the next bootstrap frames arrive, without
%   buffering or revising earlier output.  The frequency equalizer is cached.

if nargin < 4 || isempty(state), state = struct(); end
if nargin < 5, cfg = struct(); end
if ~isfield(state,'initialized') || ~state.initialized
    state = initialize_state(cfg);
end
frame = single(receivedFrame(:));
assert(numel(frame) >= state.frame_length,'Received frame must contain 4725 symbols.');
frame = frame(1:state.frame_length);

frameNumber = state.frame_count+1;
bootstrap = state.num_estimates < state.bootstrap_estimates;
periodic = state.update_interval > 0 && ~bootstrap && ...
    mod(frameNumber-state.bootstrap_estimates,state.update_interval) == 0;
shouldEstimate = bootstrap || periodic;
tEstimate = tic;
estimateInfo = struct('support',zeros(1,0),'num_taps',0, ...
    'residual_energy_ratio',NaN,'residual_energy_history',zeros(0,1));
if shouldEstimate
    [hNew,estimateInfo] = DTMB_Fast_PN_Channel_Estimate( ...
        pnObserved,phaseIndex,state.cfg);
    if bootstrap
        state.num_estimates = state.num_estimates+1;
        if state.num_estimates == 1
            state.h_average = hNew;
        else
            state.h_average = state.h_average + ...
                (hNew-state.h_average)/state.num_estimates;
        end
    else
        alpha = state.update_alpha;
        state.h_average = (1-alpha)*state.h_average+alpha*hNew;
        state.num_periodic_updates = state.num_periodic_updates+1;
    end
    hFrequency = fft(state.h_average,state.frame_length);
    state.equalizer = conj(hFrequency)./(abs(hFrequency).^2+state.lambda);
end
estimateTime = toc(tEstimate);
assert(~isempty(state.equalizer),'Channel equalizer has not been initialized.');

tEqualize = tic;
equalized = ifft(fft(frame,state.frame_length).*state.equalizer);
equalizeTime = toc(tEqualize);
state.frame_count = frameNumber;
state.total_compute_s = state.total_compute_s+estimateTime+equalizeTime;

info = struct('frame_number',frameNumber,'channel_updated',shouldEstimate, ...
    'bootstrap_complete',state.num_estimates >= state.bootstrap_estimates, ...
    'num_channel_estimates',state.num_estimates, ...
    'support',estimateInfo.support,'residual_energy_ratio', ...
    estimateInfo.residual_energy_ratio,'estimation_s',estimateTime, ...
    'equalization_s',equalizeTime,'total_s',estimateTime+equalizeTime);
end

function state = initialize_state(cfg)
cirLength = max(1,round(get_cfg(cfg,'channel_impulse_response_length',160)));
state = struct('initialized',true,'cfg',cfg,'frame_length',4725, ...
    'lambda',single(get_cfg(cfg,'body_eq_lambda',1e-3)), ...
    'bootstrap_estimates',max(1,round(get_cfg(cfg, ...
        'tc_channel_estimation_frames',2))), ...
    'update_interval',max(0,round(get_cfg(cfg, ...
        'channel_update_interval_frames',0))), ...
    'update_alpha',single(get_cfg(cfg,'channel_update_alpha',0.25)), ...
    'frame_count',0,'num_estimates',0,'num_periodic_updates',0, ...
    'h_average',complex(zeros(cirLength,1,'single')),'equalizer',[], ...
    'total_compute_s',0);
end

function value = get_cfg(cfg,name,defaultValue)
if isstruct(cfg) && isfield(cfg,name) && ~isempty(cfg.(name))
    value = cfg.(name);
else
    value = defaultValue;
end
end
