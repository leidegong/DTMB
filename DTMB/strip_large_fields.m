function compact = strip_large_fields(pass)
    compact = pass;
    compact = rmfield_if_present(compact, {'PN', 'Frame'});
end


