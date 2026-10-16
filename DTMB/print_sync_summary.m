function print_sync_summary(detected_mode, coarse_frame_starts, exact_frame_start, phase_index)
    fprintf('初始同步结果汇总：\n');
    fprintf('  帧头模式             ：%s\n', detected_mode);
    fprintf('  粗同步帧数           ：%d\n', length(coarse_frame_starts));
    fprintf('  第一帧粗起点         ：%.2f 采样点\n', coarse_frame_starts(1));
    fprintf('  第一帧精确起点       ：%.2f 采样点\n', exact_frame_start);
    fprintf('  第一帧 PN 相位号     ：%d\n', phase_index);
end


