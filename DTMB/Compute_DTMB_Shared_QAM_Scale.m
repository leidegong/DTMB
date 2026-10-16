function [shared_rms, info] = Compute_DTMB_Shared_QAM_Scale(Frame_matrix, cfg)
% Compute_DTMB_Shared_QAM_Scale
% Use the post-sync/SVD frame-body frequency-domain RMS as one shared QAM
% scale for constellation plots, EVM evaluation, and downstream compensation.

    Lpn = get_cfg_scalar(cfg, 'pn_length', 945);
    Ldata = get_cfg_scalar(cfg, 'body_length', 3780);
    Nfft = Ldata;

    info = struct();
    info.reference = 'post-sync-or-full-frame-SVD Frame body frequency-domain global RMS, excluding TPS';
    info.pn_length = Lpn;
    info.body_length = Ldata;
    info.matrix_sps = 1;
    info.nfft = Nfft;
    info.include_tps = false;
    info.num_frames = 0;
    info.num_valid_symbols = 0;
    info.num_excluded_tps_bins = 0;
    info.fallback = false;

    if isempty(Frame_matrix) || size(Frame_matrix, 2) < Lpn + Ldata
        warning('SharedQAMScale:InvalidFrame', ...
            'Frame matrix is empty or too short; using shared RMS = 1.');
        shared_rms = 1;
        info.shared_rms = shared_rms;
        info.fallback = true;
        return;
    end

    body_idx = Lpn + (1:Ldata);
    body = Frame_matrix(:, body_idx);
    body = body - mean(body, 2);
    body_fd = fftshift(fft(body, Nfft, 2), 2);
    mask = DTMB_Body_Frequency_Mask(false, Nfft);
    symbols = body_fd(:, mask);
    valid = isfinite(real(symbols)) & isfinite(imag(symbols));
    values = symbols(valid);

    info.num_frames = size(Frame_matrix, 1);
    info.num_valid_symbols = numel(values);
    info.num_excluded_tps_bins = Nfft - nnz(mask);

    if isempty(values)
        warning('SharedQAMScale:NoValidSymbol', ...
            'No valid frequency-domain body symbols; using shared RMS = 1.');
        shared_rms = 1;
        info.shared_rms = shared_rms;
        info.fallback = true;
        return;
    end

    shared_rms = sqrt(mean(abs(values).^2));
    if ~isfinite(shared_rms) || shared_rms <= 0
        warning('SharedQAMScale:BadScale', ...
            'Invalid shared RMS; using shared RMS = 1.');
        shared_rms = 1;
        info.fallback = true;
    end
    info.shared_rms = shared_rms;
end

function value = get_cfg_scalar(cfg, name, default_value)
    if isstruct(cfg) && isfield(cfg, name) && ~isempty(cfg.(name))
        value = cfg.(name);
    else
        value = default_value;
    end
end
