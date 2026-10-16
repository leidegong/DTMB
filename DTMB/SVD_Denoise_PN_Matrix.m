function [PN_denoised, svd_info] = SVD_Denoise_PN_Matrix(PN_matrix, rank_keep)
% SVD_Denoise_PN_Matrix
% 功能：
%   对 PN 矩阵做 SVD 分解，并只保留前 rank_keep 个奇异值成分，
%   重构去噪后的 PN 矩阵
%
% 输入：
%   PN_matrix  - 输入 PN 矩阵，每一行为一帧 PN
%   rank_keep  - 保留的有效秩
%
% 输出：
%   PN_denoised - SVD 去噪后的 PN 矩阵
%   svd_info    - 保存奇异值、能量占比等信息

    if nargin < 2
        rank_keep = 1;
    end

    if isempty(PN_matrix)
        error('输入 PN_matrix 为空，无法进行 SVD 去噪。');
    end

    % 转 double，避免某些情况下类型不一致
    PN_matrix = double(PN_matrix);

    % SVD 分解
    [U, S, V] = svd(PN_matrix, 'econ');

    singular_vals = diag(S);
    total_rank = length(singular_vals);

    if rank_keep > total_rank
        warning('rank_keep 大于矩阵可用秩，已自动截断为 %d', total_rank);
        rank_keep = total_rank;
    end

    % 只保留前 rank_keep 个奇异值
    S_denoised = zeros(size(S));
    S_denoised(1:rank_keep, 1:rank_keep) = S(1:rank_keep, 1:rank_keep);

    % 重构
    PN_denoised = U * S_denoised * V';

    % 能量统计
    energy_all = sum(singular_vals.^2);
    energy_keep = sum(singular_vals(1:rank_keep).^2);
    energy_ratio = energy_keep / (energy_all + eps);

    svd_info = struct();
    svd_info.U = U;
    svd_info.S = S;
    svd_info.V = V;
    svd_info.singular_vals = singular_vals;
    svd_info.rank_keep = rank_keep;
    svd_info.energy_ratio = energy_ratio;

    disp('    -> ✅ PN 矩阵 SVD 去噪完成');
    disp(['    -> 保留秩 rank = ', num2str(rank_keep)]);
    disp(['    -> 保留能量占比 = ', num2str(energy_ratio * 100, '%.2f'), '%']);
end