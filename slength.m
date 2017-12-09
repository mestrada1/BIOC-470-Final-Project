function l = slength(S)
l = NaN;
if isstruct(S)
    l = 0;
    if ~isempty(S) && ~isempty(fieldnames(S))
        f = fields(S);
        nf = length(f);
        len = nan(nf, 1);
        for i = 1:nf
            f1 = getfield(S, f{i});
            if ischar(f1)
                len(i) = 1;
            else 
                len(i) = size(f1,1);
            end
        end
        ulen = unique(len);
        if length(ulen)==1
            l = ulen;
        else
            fprintf('Deprecated use of slength for structure with fields of nonuniform length\n');
            l = len(1);
        end
    end
end
end   