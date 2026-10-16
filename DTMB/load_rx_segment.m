function rx_raw = load_rx_segment(cfg)
    try
        info = whos('-file', cfg.filename, cfg.data_var_name);
    catch
        info = [];
    end

    total_samples_needed = cfg.frame_symbol_len * cfg.sps * cfg.target_frames;

    if ~isempty(info)
        raw_len = info.size(1);
        end_sample = min(cfg.start_sample + total_samples_needed - 1, raw_len);

        if cfg.start_sample > raw_len
            error('起始采样点 start_sample (%d) 超过原始数据长度 (%d)。', cfg.start_sample, raw_len);
        end

        m = matfile(cfg.filename);
        rx_raw = m.(cfg.data_var_name)(cfg.start_sample:end_sample, 1);
    else
        dataset_path = ['/' cfg.data_var_name];
        try
            h5_info = h5info(cfg.filename, dataset_path);
        catch ME
            error('无法在 %s 中找到 MAT/HDF5 变量 "%s"：%s', ...
                cfg.filename, cfg.data_var_name, ME.message);
        end

        h5_size = h5_info.Dataspace.Size;
        raw_len = h5_size(1);
        end_sample = min(cfg.start_sample + total_samples_needed - 1, raw_len);

        if cfg.start_sample > raw_len
            error('起始采样点 start_sample (%d) 超过原始数据长度 (%d)。', cfg.start_sample, raw_len);
        end

        raw_h5 = h5read(cfg.filename, dataset_path, [cfg.start_sample 1], [end_sample - cfg.start_sample + 1 1]);
        if isstruct(raw_h5) && isfield(raw_h5, 'real') && isfield(raw_h5, 'imag')
            rx_raw = complex(raw_h5.real, raw_h5.imag);
        else
            rx_raw = raw_h5;
        end
    end

    rx_raw = rx_raw(:);

    fprintf('读取采样点范围：%d : %d，共 %d 点\n', ...
        cfg.start_sample, end_sample, length(rx_raw));
end


