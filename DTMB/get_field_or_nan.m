function val = get_field_or_nan(s, field_name)
    if isstruct(s) && isfield(s, field_name)
        val = s.(field_name);
    else
        val = nan;
    end
end


