function text = translate_stop_reason(reason)
    switch char(reason)
        case 'disabled'
            text = '未启用闭环';
        case 'max_iter'
            text = '达到最大迭代次数';
        case 'cfo_tol_before_compensation'
            text = '补偿前残留 CFO 已满足门限';
        case 'cfo_tol_after_compensation'
            text = '补偿后残留 CFO 满足门限';
        case 'small_cfo_improvement'
            text = '残留 CFO 改善量过小';
        case 'rejected_worse_cfo'
            text = '补偿后残留 CFO 变差，已拒绝该轮';
        otherwise
            text = char(reason);
    end
end


