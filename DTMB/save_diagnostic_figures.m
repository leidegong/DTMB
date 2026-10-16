function save_diagnostic_figures(cfg)
    if ~isfield(cfg, 'save_diagnostic_figures') || ~cfg.save_diagnostic_figures
        return;
    end

    if ~exist(cfg.figure_output_dir, 'dir')
        mkdir(cfg.figure_output_dir);
    end

    figs = findall(0, 'Type', 'figure');
    for k = 1:numel(figs)
        fig = figs(k);
        if ~isvalid(fig)
            continue;
        end
        fig_name = get(fig, 'Name');
        if isempty(fig_name)
            fig_name = sprintf('figure_%02d', k);
        end
        file_name = matlab.lang.makeValidName(fig_name);
        out_png = fullfile(cfg.figure_output_dir, [file_name, '.png']);
        try
            exportgraphics(fig, out_png, 'Resolution', 160);
        catch
            saveas(fig, out_png);
        end
    end

    fprintf('诊断图已保存到目录：%s\n', cfg.figure_output_dir);
end






