function [PN_align, Frame_align, cpe_info] = align_pn_frame_common_phase(PN_matrix, Frame_matrix, phase_index)
    phase_index = phase_index(:);
    M = min([size(PN_matrix, 1), size(Frame_matrix, 1), numel(phase_index)]);

    PN_align = PN_matrix(1:M, :);
    Frame_align = Frame_matrix(1:M, :);
    phase_index = phase_index(1:M);

    phi = zeros(M, 1);
    amp = zeros(M, 1);

    for k = 1:M
        local_pn = PN_gen(phase_index(k)).';
        local_pn = local_pn(:);
        z = PN_align(k, :) * conj(local_pn);
        phi(k) = angle(z);
        amp(k) = abs(z);
        phase_corr = exp(-1j * phi(k));
        PN_align(k, :) = PN_align(k, :) * phase_corr;
        Frame_align(k, :) = Frame_align(k, :) * phase_corr;
    end

    cpe_info = struct();
    cpe_info.phi = phi;
    cpe_info.amp = amp;
    cpe_info.phase_index = phase_index;
end


