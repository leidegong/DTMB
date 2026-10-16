function eps_islr = Estimate_Epsilon_From_TC_NoiseFloor(rcs_frame, islr_opts)
% 根据接收参考信号 TC 结果的噪底平台估计 epsilon

    [~, tc_pow, ~, aux] = Calc_ISLR_DTMB(rcs_frame, islr_opts);

    tc_db = 10 * log10(tc_pow / (max(tc_pow) + eps) + eps);

    side_db = tc_db(aux.side_mask);

    % 用“平台值”而不是中位数
    eps_islr = prctile(side_db, 75);
end