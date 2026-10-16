function s = rmfield_if_present(s, names)
    for k = 1:numel(names)
        if isfield(s, names{k})
            s = rmfield(s, names{k});
        end
    end
end


