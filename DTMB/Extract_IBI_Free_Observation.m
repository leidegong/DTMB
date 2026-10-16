function y = Extract_IBI_Free_Observation(pn_one_frame, Ltap)
% Extract_IBI_Free_Observation
% 功能：
%   从一帧 PN 序列中提取 IBI-free 区域作为观测向量
%
% 输入：
%   pn_one_frame - 1×Lpn 或 Lpn×1 的单帧 PN 序列
%   Ltap         - 假设信道长度
%
% 输出：
%   y            - B×1 观测向量，B = Lpn - Ltap + 1
%
% 说明：
%   论文中 IBI-free 区域写作：
%   y = [pn(Ltap-1), ..., pn(Lpn-1)]^T
%   这是 0 基索引。
%   MATLAB 中对应：
%   y = pn(Ltap : Lpn)

    pn_one_frame = pn_one_frame(:);
    Lpn = length(pn_one_frame);

    if Ltap > Lpn
        error('Ltap 不能大于 PN 长度。');
    end

    y = pn_one_frame(Ltap:Lpn);
end