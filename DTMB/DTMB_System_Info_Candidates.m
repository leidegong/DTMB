function candidates = DTMB_System_Info_Candidates()
% DTMB_System_Info_Candidates 64QAM的六种码率/时域交织系统信息字。

    words = [ ...
        "00100010011011010111010001110000"; ...
        "11011101100100101000101110001111"; ...
        "01000100000010110001001000010110"; ...
        "10111011111101001110110111101001"; ...
        "00010001010111100100011101000011"; ...
        "11101110101000011011100010111100"];
    rate_index = [1;1;2;2;3;3];
    rate_value = [0.4;0.4;0.6;0.6;0.8;0.8];
    interleaver_mode = [1;2;1;2;1;2];
    M = [240;720;240;720;240;720];

    candidates = repmat(struct(), 6, 1);
    for idx = 1:6
        bits32 = double(char(words(idx)) - '0').';
        bits36 = [ones(4, 1); bits32];
        signs = 2*bits36 - 1;
        candidates(idx).index = idx;
        candidates(idx).rate_index = rate_index(idx);
        candidates(idx).rate = rate_value(idx);
        candidates(idx).interleaver_mode = interleaver_mode(idx);
        candidates(idx).M = M(idx);
        candidates(idx).bits = bits36;
        candidates(idx).symbols = signs .* (1 + 1j) / sqrt(2);
        candidates(idx).word = char(words(idx));
    end
end
