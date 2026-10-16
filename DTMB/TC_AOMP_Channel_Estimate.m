function out = TC_AOMP_Channel_Estimate(y, Theta, rcs_frame, opts)
% TC_AOMP_Channel_Estimate
% 功能：
%   基于 TC-AOMP 的稀疏信道估计
%
% 输入：
%   y         - IBI-free 区域观测向量, B×1
%   Theta     - 感知矩阵, B×Ltap
%   rcs_frame - 当前帧参考信号（用于均衡与 ISLR）
%   local_pn  - 当前帧本地 PN
%   opts      - 参数结构体
%
% 输出：
%   out.h_hat
%   out.req
%   out.best_islr
%   out.islr_hist
%   out.res_norm_hist
%   out.support

    y = y(:);
    rcs_frame = rcs_frame(:);

    [~, Ltap] = size(Theta);

    % 默认参数
    if ~isfield(opts, 'Kmax')
        opts.Kmax = min(Ltap, 200);
    end
    if ~isfield(opts, 'eps_islr')
        opts.eps_islr = -30;
    end
    if ~isfield(opts, 'main_half_width')
        opts.main_half_width = 2;
    end
    if ~isfield(opts, 'lambda_eq')
        opts.lambda_eq = 1e-3;
    end
    if ~isfield(opts, 'min_delta_islr')
        opts.min_delta_islr = 0.2;
    end
    if ~isfield(opts, 'patience')
        opts.patience = 5;
    end
    if ~isfield(opts, 'min_iter')
        opts.min_iter = 1;
    end

    % 初始化
    w = y;                          % 残差
    support = [];                   % 支持集
    h_hat = zeros(Ltap, 1);

    islr_hist = nan(opts.Kmax, 1);
    res_norm_hist = nan(opts.Kmax, 1);

    best_islr = inf;
    best_h = h_hat;
    best_req = rcs_frame;
    best_support = support;
    best_iter = 0;
    stop_reason = "Kmax";
    no_improve_count = 0;

    for k = 1:opts.Kmax

        % Step 1: 选最匹配原子
        proj = abs(Theta' * w);
        if ~isempty(support)
            proj(support) = -inf;
        end
        [~, zk] = max(proj);

        % Step 2: 更新支持集
        support = [support, zk];

        % Step 3: LS 求系数
        Omega = Theta(:, support);
        theta_k = (Omega' * Omega) \ (Omega' * y);

        % Step 4: 还原完整 h_hat
        h_hat = zeros(Ltap, 1);
        h_hat(support) = theta_k;

        % Step 5: 更新残差
        w = y - Omega * theta_k;

        % Step 6: 用当前信道均衡参考信号
        req = Equalize_By_Channel(rcs_frame, h_hat, opts.lambda_eq);

        % Step 7: 计算 ISLR
        [islr_db, ~, ~, ~] = Calc_ISLR_DTMB(req, opts.islr_opts);

        islr_hist(k) = islr_db;
        res_norm_hist(k) = norm(w);

        % 保存最好结果
        if islr_db < best_islr - opts.min_delta_islr
            best_islr = islr_db;
            best_h = h_hat;
            best_req = req;
            best_support = support;
            best_iter = k;
            no_improve_count = 0;
        elseif islr_db < best_islr
            best_islr = islr_db;
            best_h = h_hat;
            best_req = req;
            best_support = support;
            best_iter = k;
            no_improve_count = no_improve_count + 1;
        else
            no_improve_count = no_improve_count + 1;
        end

        % Step 8: 自适应停止
        if islr_db <= opts.eps_islr
            stop_reason = "epsilon";
            break;
        end

        if k >= opts.min_iter && no_improve_count >= opts.patience
            stop_reason = "no_improve";
            break;
        end
    end

    out = struct();
    out.h_hat = best_h;
    out.req = best_req;
    out.best_islr = best_islr;
    out.islr_hist = islr_hist(1:k);
    out.res_norm_hist = res_norm_hist(1:k);
    out.support = best_support;
    out.best_iter = best_iter;
    out.stop_iter = k;
    out.stop_reason = stop_reason;
end
