function codeword = DTMB_BCH762_Encode(data)
% DTMB_BCH762_Encode 系统编码缩短BCH(762,752)，每列一个数据块。

    data = logical(data);
    if isvector(data)
        data = data(:);
    end
    if size(data, 1) ~= 752
        error('BCH 编码输入必须每列为 752 个信息比特。');
    end
    generator = logical([1 0 0 0 0 0 0 1 0 0 1]).';
    work = [data;false(10,size(data,2))];
    for idx = 1:752
        feedback = work(idx,:);
        if any(feedback)
            work(idx:idx+10,:) = xor(work(idx:idx+10,:), ...
                generator & feedback);
        end
    end
    codeword = [data;work(end-9:end,:)];
end
