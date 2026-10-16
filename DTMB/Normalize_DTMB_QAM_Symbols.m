function sym_scaled = Normalize_DTMB_QAM_Symbols(sym, cfg)
% Normalize_DTMB_QAM_Symbols 使用 cfg.shared_qam_rms 统一 QAM 频域符号尺度。

    scale = [];
    if isstruct(cfg) && isfield(cfg, 'shared_qam_rms')
        scale = cfg.shared_qam_rms;
    end
    if isempty(scale) || ~isfinite(scale) || scale <= 0
        scale = rms(sym(:)) + eps;
    end
    sym_scaled = sym ./ (scale + eps);
end