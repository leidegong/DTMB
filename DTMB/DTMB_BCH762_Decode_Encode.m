function [corrected_codewords, payload_bits, stats] = DTMB_BCH762_Decode_Encode(codewords)
% DTMB_BCH762_Decode_Encode 缩短BCH(762,752)单错纠正并系统重编码。
% 母码生成多项式为 g(x)=1+x^3+x^10。

    codewords = logical(codewords);
    if isvector(codewords)
        codewords = codewords(:);
    end
    if size(codewords, 1) ~= 762
        error('BCH 输入矩阵必须每列为 762 比特码字。');
    end

    persistent syndrome_position_lut
    if isempty(syndrome_position_lut)
        unit_errors = false(762,762);
        unit_errors(1:763:end) = true;
        keys = remainder_key(poly_remainder(unit_errors));
        syndrome_position_lut = zeros(1024,1,'int16');
        for bit_idx = 1:762
            lut_idx = double(keys(bit_idx))+1;
            if syndrome_position_lut(lut_idx) == 0
                syndrome_position_lut(lut_idx) = int16(bit_idx);
            else
                syndrome_position_lut(lut_idx) = int16(-1);
            end
        end
    end

    num_blocks = size(codewords, 2);
    received = codewords;
    syndrome_before = poly_remainder(received);
    has_syndrome = any(syndrome_before,1);
    keys = remainder_key(syndrome_before);
    positions = reshape(syndrome_position_lut(double(keys)+1),1,[]);
    corrected_position = zeros(num_blocks, 1);
    correctable = has_syndrome & positions > 0;
    uncorrectable = has_syndrome & positions <= 0;
    corrected_position(correctable) = double(positions(correctable));
    corrected_position(uncorrectable) = -1;
    columns = find(correctable);
    if ~isempty(columns)
        linear = sub2ind(size(received), ...
            double(positions(columns)),columns);
        received(linear) = ~received(linear);
    end
    payload_bits = received(1:752,:);
    corrected_codewords = DTMB_BCH762_Encode(payload_bits);

    stats = struct();
    stats.num_blocks = num_blocks;
    stats.syndrome_fail_before = sum(any(syndrome_before, 1));
    stats.single_error_corrected = sum(corrected_position > 0);
    stats.uncorrectable = sum(corrected_position < 0);
    stats.syndrome_fail_after = 0;
    stats.syndrome_fail_before_mask = has_syndrome.';
    stats.single_error_corrected_mask = corrected_position > 0;
    stats.uncorrectable_mask = corrected_position < 0;
    stats.corrected_position = corrected_position;
end

function remainder = poly_remainder(bits)
    generator = logical([1 0 0 0 0 0 0 1 0 0 1]).';
    work = logical(bits);
    if isvector(work), work = work(:); end
    for idx = 1:(size(work,1)-10)
        feedback = work(idx,:);
        if any(feedback)
            work(idx:idx+10,:) = xor(work(idx:idx+10,:), ...
                generator & feedback);
        end
    end
    remainder = work(end-9:end,:);
end

function key = remainder_key(remainder)
    key = zeros(1,size(remainder,2),'uint16');
    for idx = 1:10
        key = bitshift(key,1)+uint16(remainder(idx,:));
    end
end
