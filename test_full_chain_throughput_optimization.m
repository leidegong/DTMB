function report = test_full_chain_throughput_optimization()
%TEST_FULL_CHAIN_THROUGHPUT_OPTIMIZATION Exactness regression tests.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'DTMB'));
rng(20260903);

symbols = complex(randn(512,1),randn(512,1));
noise = 0.01+0.1*rand(512,1);
llrFast = DTMB_64QAM_Soft_Demap(symbols,noise);
llrReference = exhaustive_soft_demap(symbols,noise);
softMaxError = max(abs(llrFast-llrReference));
assert(softMaxError < 1e-10,'Optimized 64QAM LLRs changed.');

data = rand(752,200)>0.5;
validCodewords = DTMB_BCH762_Encode(data);
received = validCodewords;
positions = mod((1:size(data,2))*37,762)+1;
linear = sub2ind(size(received),positions,1:size(received,2));
received(linear) = ~received(linear);
[corrected,payload,stats] = DTMB_BCH762_Decode_Encode(received);
assert(isequal(corrected,validCodewords),'BCH corrected codewords changed.');
assert(isequal(payload,data),'BCH payload recovery changed.');
assert(stats.single_error_corrected == size(data,2));
assert(stats.uncorrectable == 0);

optimizedFile = fullfile(root,'DTMB_full_chain_optimization', ...
    'DTMB_Full_Throughput_Optimized.mat');
baselineFile = fullfile(root,'DTMB_groupwise_equalization', ...
    'DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat');
fullComparison = isfile(optimizedFile) && isfile(baselineFile);
fecMaskExact = NaN;
fecReconstructionExact = NaN;
fecValidFrames = NaN;
if fullComparison
    O = load(optimizedFile,'result');
    B = load(baselineFile,'result');
    n = min(O.result.num_frames,B.result.num_frames);
    maskO = O.result.fec_valid_frame_mask(1:n);
    maskB = B.result.fec_valid_frame_mask(1:n);
    fecMaskExact = isequal(maskO,maskB);
    common = maskO & maskB;
    fecReconstructionExact = isequal( ...
        O.result.canonical_reconstructed_frames(common,:), ...
        B.result.canonical_reconstructed_frames(common,:));
    fecValidFrames = nnz(maskO);
    assert(fecMaskExact && fecReconstructionExact);
end

report = struct('passed',true,'soft_llr_max_abs_error',softMaxError, ...
    'bch_blocks_tested',size(data,2), ...
    'full_comparison_available',fullComparison, ...
    'fec_valid_mask_exact',fecMaskExact, ...
    'fec_reconstruction_exact',fecReconstructionExact, ...
    'fec_valid_frames',fecValidFrames);
disp(report);
end

function llr = exhaustive_soft_demap(symbols,noise)
labels = dec2bin(0:63,6)-'0';
serialBits = labels(:,end:-1:1);
constellation = DTMB_64QAM_Map(reshape(serialBits.',[],1));
matrix = zeros(numel(symbols),6);
for bit = 1:6
    points0 = constellation(serialBits(:,bit)==0).';
    points1 = constellation(serialBits(:,bit)==1).';
    d0 = min(abs(symbols-points0).^2,[],2);
    d1 = min(abs(symbols-points1).^2,[],2);
    matrix(:,bit) = (d1-d0)./noise;
end
llr = reshape(matrix.',[],1);
end
