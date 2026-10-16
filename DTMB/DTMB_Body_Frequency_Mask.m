function mask = DTMB_Body_Frequency_Mask(include_tps, nfft)
% DTMB_Body_Frequency_Mask 生成帧体频域有效数据子载波掩码。
% 默认排除 36 个 TPS 子载波，用于 shared QAM scale、星座图和 EVM 评价。

    if nargin < 1 || isempty(include_tps)
        include_tps = false;
    end
    if nargin < 2 || isempty(nfft)
        nfft = 3780;
    end

    mask_unshifted = true(nfft, 1);
    if ~include_tps
        tps_z0 = [0,140,279,419,420,560,699,839,840,980, ...
                  1119,1259,1260,1400,1539,1679,1680,1820, ...
                  1959,2099,2100,2240,2379,2519,2520,2660, ...
                  2799,2939,2940,3080,3219,3359,3360,3500, ...
                  3639,3779];
        tps_y0 = map_z_positions_to_y_positions(tps_z0);
        tps_y0 = tps_y0(isfinite(tps_y0) & tps_y0 >= 0 & tps_y0 < nfft);
        mask_unshifted(tps_y0 + 1) = false;
    end

    mask = fftshift(mask_unshifted);
end

function y_pos0 = map_z_positions_to_y_positions(z_pos0)
    y_pos0 = nan(size(z_pos0));
    for i = 0:2
        for j = 0:2
            for k = 0:2
                for l = 0:1
                    for m = 0:1
                        for n = 0:4
                            for o = 0:6
                                dst = o * 540 + n * 108 + m * 54 + l * 27 + k * 9 + j * 3 + i;
                                src = i * 1260 + j * 420 + k * 140 + l * 70 + m * 35 + n * 7 + o;
                                y_pos0(z_pos0 == src) = dst;
                            end
                        end
                    end
                end
            end
        end
    end
end
