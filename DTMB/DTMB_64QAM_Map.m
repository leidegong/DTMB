function symbols = DTMB_64QAM_Map(bits)
% DTMB_64QAM_Map 按GB 20600比特顺序[b0...b5]执行归一化64QAM映射。

    bits = double(bits(:) ~= 0);
    if mod(numel(bits), 6) ~= 0
        error('64QAM 映射比特数必须是 6 的整数倍。');
    end
    groups = reshape(bits, 6, []).';
    % Direct Gray-label lookup.  The previous ismember(...,'rows') calls
    % performed two table joins for every symbol block.
    level_lut = [-7;-5;-1;-3;7;5;1;3];
    idx_i = 1+groups(:,1)+2*groups(:,2)+4*groups(:,3);
    idx_q = 1+groups(:,4)+2*groups(:,5)+4*groups(:,6);
    symbols = complex(level_lut(idx_i),level_lut(idx_q))/sqrt(42);
end
