function plot_required_sync_diagnostics(cfo_loop_info, PN_before_sync, PN_after_sync, ...
                                        PN_aligned, PN_denoised, svd_info)
    % 同步阶段必要诊断图集合：
    % 1. 闭环残留 CFO 收敛过程
    % 2. 帧间相参 FFT 残留 CFO 谱峰
    % 3. 时频同步前后帧头星座对比
    % 4. PN 星座同步效果汇总
    % 5. PN 矩阵 SVD 奇异值
    plot_cfo_loop_summary(cfo_loop_info);
    plot_cfo_fft_spectrum(cfo_loop_info.final_residual_cfo);
    plot_frame_head_constellation_sync_compare(PN_before_sync, PN_after_sync);
    plot_pn_constellation_summary(PN_after_sync, PN_aligned, PN_denoised, svd_info);
    plot_svd_singular_summary(svd_info);
end

function plot_cfo_loop_summary(cfo_loop_info)
    if ~isfield(cfo_loop_info, 'history') || isempty(cfo_loop_info.history)
        return;
    end

    hist = cfo_loop_info.history;
    iter = [hist.iteration];
    before = [hist.before_cfo_hz];
    after = [hist.after_cfo_hz];
    accepted = [hist.accepted];

    % 对应输出图：闭环残留 CFO 收敛过程
    % 作用：观察每轮残留 CFO 补偿前后是否逐步收敛到 0 Hz 附近。
    figure('Name', '闭环残留 CFO 收敛过程', 'Color', 'w');
    subplot(2,1,1);
    plot(iter, before, 'o-', 'LineWidth', 1.2); hold on;
    plot(iter, after, 's-', 'LineWidth', 1.2);
    yline(0, ':k');
    grid on;
    xlabel('闭环迭代次数');
    ylabel('残留 CFO / Hz');
    legend('补偿前', '补偿后', 'Location', 'best');
    title(sprintf('闭环 CFO 收敛过程，停止原因：%s', translate_stop_reason(cfo_loop_info.stop_reason)));

    subplot(2,1,2);
    bar(iter, [hist.applied_cfo_hz]);
    grid on;
    xlabel('闭环迭代次数');
    ylabel('本轮施加 CFO / Hz');
    title(sprintf('接受迭代：%d / %d，累计施加 CFO = %.6f Hz', ...
        sum(accepted), numel(accepted), cfo_loop_info.total_applied_cfo_hz));
end

function plot_cfo_fft_spectrum(cfo)
    if ~isstruct(cfo) || ~isfield(cfo, 'fft_info') || isempty(fieldnames(cfo.fft_info))
        return;
    end

    fft_info = cfo.fft_info;
    if ~isfield(fft_info, 'freq_axis') || isempty(fft_info.freq_axis)
        return;
    end

    freq_axis = fft_info.freq_axis(:);
    power_db = 10*log10(fft_info.power(:) / (max(fft_info.power(:)) + eps) + eps);
    view_half_width_hz = max(20, 5 * abs(fft_info.fft_cfo_hz));
    idx_show = abs(freq_axis - fft_info.fft_cfo_hz) <= view_half_width_hz;
    if sum(idx_show) < 10
        idx_show = true(size(freq_axis));
    end

    % 对应输出图：帧间相参 FFT 残留 CFO 谱峰
    % 作用：展示多帧 PN 相干积累后的残留 CFO 频谱峰值和细化估计结果。
    figure('Name', '帧间相参 FFT 残留 CFO 谱峰', 'Color', 'w');
    plot(freq_axis(idx_show), power_db(idx_show), 'LineWidth', 1.2); hold on;
    xline(fft_info.fft_cfo_hz, '--r', sprintf('FFT 峰值 %.4f Hz', fft_info.fft_cfo_hz));
    xline(cfo.res_cfo_hz, '--k', sprintf('细化结果 %.4f Hz', cfo.res_cfo_hz));
    grid on;
    xlabel('频率 / Hz');
    ylabel('归一化功率 / dB');
    title(sprintf('帧间相参 FFT：原始分辨率 %.4f Hz，零填充 bin 间隔 %.4f Hz', ...
        fft_info.native_resolution_hz, fft_info.fft_bin_spacing_hz));
end

function plot_frame_head_constellation_sync_compare(PN_before_sync, PN_after_sync)
    max_points = 40000;
    before_pts = pick_complex_points(PN_before_sync, max_points);
    after_pts = pick_complex_points(PN_after_sync, max_points);

    % 对应输出图：时频同步前后帧头星座对比
    % 作用：比较残留 CFO 闭环补偿前后 PN 帧头星座是否更稳定。
    figure('Name', '时频同步前后帧头星座对比', 'Color', 'w');
    subplot(1,2,1);
    scatter(real(before_pts), imag(before_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title('闭环精同步前帧头 PN 星座');

    subplot(1,2,2);
    scatter(real(after_pts), imag(after_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title('闭环精同步后帧头 PN 星座');
end

function plot_pn_constellation_summary(PN_raw, PN_aligned, PN_denoised, svd_info)
    max_points = 40000;
    raw_pts = pick_complex_points(PN_raw, max_points);
    aligned_pts = pick_complex_points(PN_aligned, max_points);
    denoised_pts = pick_complex_points(PN_denoised, max_points);

    % 对应输出图：PN 星座同步效果汇总
    % 作用：依次比较闭环同步后、逐帧公共相位补偿后、SVD 去噪后的 PN 星座。
    figure('Name', 'PN 星座同步效果汇总', 'Color', 'w');
    subplot(1,3,1);
    scatter(real(raw_pts), imag(raw_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title('闭环同步后');

    subplot(1,3,2);
    scatter(real(aligned_pts), imag(aligned_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    title('逐帧公共相位补偿后');

    subplot(1,3,3);
    scatter(real(denoised_pts), imag(denoised_pts), 3, '.');
    axis equal; grid on;
    xlabel('同相分量 I'); ylabel('正交分量 Q');
    if isfield(svd_info, 'enabled') && ~svd_info.enabled
        title('未启用 SVD 去噪');
    else
        title('SVD 去噪后');
    end
end

function pts = pick_complex_points(X, max_points)
    pts = X(:);
    if numel(pts) > max_points
        idx = round(linspace(1, numel(pts), max_points));
        pts = pts(idx);
    end
end

function plot_svd_singular_summary(svd_info)
    if isfield(svd_info, 'enabled') && ~svd_info.enabled
        return;
    end
    if ~isfield(svd_info, 'singular_vals') || isempty(svd_info.singular_vals)
        return;
    end

    s = svd_info.singular_vals(:);
    % 对应输出图：PN 矩阵 SVD 奇异值
    % 作用：查看 PN 矩阵奇异值衰减和当前 SVD 保留秩对应的累计能量。
    figure('Name', 'PN 矩阵 SVD 奇异值', 'Color', 'w');
    subplot(1,2,1);
    semilogy(s / (max(s) + eps), 'o-', 'LineWidth', 1.1);
    grid on;
    xlabel('奇异值序号');
    ylabel('归一化奇异值');
    title(sprintf('SVD 奇异值谱，保留秩 = %d', svd_info.rank_keep));

    subplot(1,2,2);
    energy = cumsum(s.^2) / (sum(s.^2) + eps);
    plot(energy, 'LineWidth', 1.2);
    yline(svd_info.energy_ratio, '--');
    grid on;
    xlabel('秩');
    ylabel('累计能量占比');
    title(sprintf('保留能量 = %.2f%%', 100 * svd_info.energy_ratio));
end

