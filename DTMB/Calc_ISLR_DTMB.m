function [islr_db, tc_pow, peak_idx, aux] = Calc_ISLR_DTMB(req, opts)
% 基于论文式(2-22)的 ISLR 计算
% req  : 均衡后参考信号 req-m
% opts : 
%   opts.main_half_width  - 主瓣半宽，建议 0 或 1
%   opts.roi_idx          - ROI 样点索引（优先）
%   opts.roi_half_width   - 若未给 roi_idx，则用主峰两侧对称窗口定义 ROI

    if nargin < 2
        opts = struct();
    end
    if ~isfield(opts, 'main_half_width')
        opts.main_half_width = 0;
    end
    if ~isfield(opts, 'roi_idx')
        opts.roi_idx = [];
    end
    if ~isfield(opts, 'roi_half_width')
        opts.roi_half_width = [];
    end

    req = req(:);

    % 论文表2.1第7步：对均衡后参考信号做时域自相关
    xi = xcorr(req, req);
    tc_pow = abs(xi).^2;

    [~, peak_idx] = max(tc_pow);
    N = numel(tc_pow);

    % ROI 定义
    if ~isempty(opts.roi_idx)
        roi_mask = false(N,1);
        roi_mask(opts.roi_idx) = true;
    else
        if isempty(opts.roi_half_width)
            error('请设置 opts.roi_idx 或 opts.roi_half_width 来定义 ROI');
        end
        roi_l = max(1, peak_idx - opts.roi_half_width);
        roi_r = min(N, peak_idx + opts.roi_half_width);
        roi_mask = false(N,1);
        roi_mask(roi_l:roi_r) = true;
    end

    % 主瓣区域
    main_l = max(1, peak_idx - opts.main_half_width);
    main_r = min(N, peak_idx + opts.main_half_width);
    main_mask = false(N,1);
    main_mask(main_l:main_r) = true;

    % ROI 内旁瓣
    side_mask = roi_mask & ~main_mask;
    if ~any(side_mask)
        error('ROI 内没有旁瓣样点，无法计算 ISLR');
    end

    E_main = sum(tc_pow(main_mask));
    E_side_avg = mean(tc_pow(side_mask));   % 论文要求：旁瓣平均功率

    islr_db = 10 * log10((E_side_avg + eps) / (E_main + eps));

    aux = struct();
    aux.roi_mask  = roi_mask;
    aux.main_mask = main_mask;
    aux.side_mask = side_mask;
end