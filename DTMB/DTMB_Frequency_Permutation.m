function [z_to_y, tps_z_index] = DTMB_Frequency_Permutation()
% DTMB_Frequency_Permutation GB 20600-2006 3780载波频域交织置换。

    persistent permutation_cache tps_cache
    if isempty(permutation_cache)
        permutation_cache = zeros(3780, 1);
        for i = 0:2
            for j = 0:2
                for k = 0:2
                    for l = 0:1
                        for m = 0:1
                            for n = 0:4
                                for o = 0:6
                                    y0 = o*540 + n*108 + m*54 + l*27 + k*9 + j*3 + i;
                                    z0 = i*1260 + j*420 + k*140 + l*70 + m*35 + n*7 + o;
                                    permutation_cache(z0 + 1) = y0 + 1;
                                end
                            end
                        end
                    end
                end
            end
        end
        tps_cache = [0,140,279,419,420,560,699,839,840,980,1119,1259, ...
            1260,1400,1539,1679,1680,1820,1959,2099,2100,2240,2379, ...
            2519,2520,2660,2799,2939,2940,3080,3219,3359,3360,3500, ...
            3639,3779].' + 1;
    end
    z_to_y = permutation_cache;
    tps_z_index = tps_cache;
end
