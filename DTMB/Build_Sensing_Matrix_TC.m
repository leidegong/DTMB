function Theta = Build_Sensing_Matrix_TC(local_pn, Ltap)
% Build_Sensing_Matrix_TC
% 功能：
%   根据本地 PN 序列构造 TC-AOMP 使用的感知矩阵 Theta
%
% 输入：
%   local_pn - 本地 PN 序列（列向量或行向量），长度 Lpn
%   Ltap     - 假设信道长度
%
% 输出：
%   Theta    - B × Ltap 感知矩阵，B = Lpn - Ltap + 1

    local_pn = local_pn(:);
    Lpn = length(local_pn);
    B = Lpn - Ltap + 1;

    if B <= 0
        error('Ltap 必须满足 Ltap <= Lpn');
    end

    Theta = zeros(B, Ltap);

    % 按论文式(2-17)构造
    % 每一列是本地 PN 的一个线性移位片段
    for col = 1:Ltap
        Theta(:, col) = local_pn(Ltap-col+1 : Lpn-col+1);
    end
end