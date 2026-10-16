function comparison = compare_DTMB_groupwise_tracking(no_group_packet, ...
        group_packet, no_group_effect_packet, group_effect_packet, output_dir)
% compare_DTMB_groupwise_tracking 对比是否启用32帧分组时变信道跟踪。

    if nargin < 1 || isempty(no_group_packet)
        no_group_packet = fullfile('DTMB_no_groupwise_test', ...
            'DTMB_No_Groupwise_Demod_Remod_Reconstruction.mat');
    end
    if nargin < 2 || isempty(group_packet)
        group_packet = fullfile('DTMB_groupwise_equalization', ...
            'DTMB_Groupwise_Final_Demod_Remod_Reconstruction.mat');
    end
    if nargin < 3 || isempty(no_group_effect_packet)
        no_group_effect_packet = fullfile('DTMB_no_groupwise_test', 'effect', ...
            'DTMB_Effect_Test_Report.mat');
    end
    if nargin < 4 || isempty(group_effect_packet)
        group_effect_packet = fullfile('DTMB_final_effect_test', ...
            'DTMB_Effect_Test_Report.mat');
    end
    if nargin < 5 || isempty(output_dir)
        output_dir = 'DTMB_no_groupwise_test';
    end

    A = load(no_group_packet, 'result');
    B = load(group_packet, 'result');
    EA = load(no_group_effect_packet, 'report');
    EB = load(group_effect_packet, 'report');
    no_group = A.result;
    group = B.result;
    no_effect = EA.report;
    group_effect = EB.report;

    metric = ["LDPC通过块"; "LDPC通过率"; "BCH不可纠字"; "FEC有效帧"; ...
        "输入判决EVM中位数"; "重构EVM中位数"; "LDPC平均迭代"; ...
        "译码改比特率"; "重构相关系数中位数"; "重构ISLR中位数"; ...
        "重构PSLR中位数"];
    unit = ["blocks"; "%"; "words"; "frames"; "linear"; "linear"; ...
        "iterations"; "linear"; "linear"; "dB"; "dB"];
    no_group_value = [nnz(no_group.fec_stats.parity_ok); ...
        100*mean(no_group.fec_stats.parity_ok); ...
        sum(no_group.fec_stats.bch_uncorrectable); ...
        nnz(no_group.fec_valid_frame_mask); median(no_group.decision_evm_rms); ...
        no_effect.evm_reconstructed_median; mean(no_group.fec_stats.iterations); ...
        no_group.fec_stats.changed_bit_rate; no_effect.correlation_median; ...
        no_effect.islr_reconstructed_median_db; no_effect.pslr_reconstructed_median_db];
    group_value = [nnz(group.fec_stats.parity_ok); ...
        100*mean(group.fec_stats.parity_ok); ...
        sum(group.fec_stats.bch_uncorrectable); ...
        nnz(group.fec_valid_frame_mask); median(group.decision_evm_rms); ...
        group_effect.evm_reconstructed_median; mean(group.fec_stats.iterations); ...
        group.fec_stats.changed_bit_rate; group_effect.correlation_median; ...
        group_effect.islr_reconstructed_median_db; ...
        group_effect.pslr_reconstructed_median_db];
    difference_group_minus_no_group = group_value - no_group_value;
    comparison = table(metric, no_group_value, group_value, ...
        difference_group_minus_no_group, unit);
    writetable(comparison, fullfile(output_dir, 'groupwise_ablation_comparison.csv'));
    failed_blocks_no_group = find(~no_group.fec_stats.parity_ok);
    failed_blocks_group = find(~group.fec_stats.parity_ok);
    save(fullfile(output_dir, 'DTMB_Groupwise_Ablation_Report.mat'), ...
        'comparison', 'failed_blocks_no_group', 'failed_blocks_group', '-v7.3');

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1300 760]);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile;
    bar([no_group_value(1), group_value(1)]); hold on;
    yline(no_group.fec_stats.num_blocks, '--k', '总码块'); hold off; grid on;
    xticklabels({'无分组','32帧分组'}); ylabel('LDPC通过块数');
    title('LDPC校验结果');
    nexttile;
    bar([no_group_value(4), group_value(4)]); grid on;
    xticklabels({'无分组','32帧分组'}); ylabel('帧数');
    title('严格FEC回填帧');
    nexttile;
    bar([no_group_value(5), group_value(5); ...
         no_group_value(6), group_value(6)]); grid on;
    xticklabels({'均衡输入','重构输出'}); ylabel('EVM');
    legend({'无分组','32帧分组'}, 'Location', 'best');
    title('最近点EVM中位数');
    nexttile;
    yyaxis left;
    bar([1 2], [no_group_value(7), group_value(7)], 0.35); ylabel('平均迭代次数');
    yyaxis right;
    plot([1 2], 100*[no_group_value(8), group_value(8)], '-o', ...
        'LineWidth', 1.4); ylabel('译码改比特率 (%)');
    xticks([1 2]); xticklabels({'无分组','32帧分组'}); grid on;
    title('译码负担');
    exportgraphics(fig, fullfile(output_dir, 'groupwise_ablation_comparison.png'), ...
        'Resolution', 180);
    close(fig);

    disp(comparison);
    fprintf('无分组失败码块：'); fprintf('%d ', failed_blocks_no_group); fprintf('\n');
    fprintf('32帧分组失败码块：'); fprintf('%d ', failed_blocks_group); fprintf('\n');
end
