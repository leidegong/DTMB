function rx_out = apply_cfo_compensation(rx_in, Fs, cfo_hz)
    rx_in = rx_in(:);
    n = (0:length(rx_in)-1).';
    rx_out = rx_in .* exp(-1j * 2*pi*cfo_hz/Fs * n);
end


