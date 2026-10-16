function PN = PN_gen(frame_num)
%PN_GEN Return the PN945 sequence for one DTMB frame phase.
%   Cache the 200-sequence ROM.  The previous implementation evaluated the
%   200-by-945 literal script on every call, which dominated short-frame
%   demodulation/remodulation latency.
    persistent bipolarPnRom
    if isempty(bipolarPnRom)
        romFile = fullfile(fileparts(mfilename('fullpath')), ...
            'PN945_ROM_Bipolar.mat');
        if isfile(romFile)
            cached = load(romFile,'bipolar_pn_rom');
            bipolarPnRom = cached.bipolar_pn_rom;
        else
            PN_Rom_Create_945;
            bipolarPnRom = single(1 - 2*PN_rom);
        end
    end
    frameIndex = mod(frame_num,200) + 1;
    PN = bipolarPnRom(frameIndex,:);
end
