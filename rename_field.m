function S = rename_field(S, oldname, newname)

if iscell(oldname) && iscell(newname)
    if ~iscell(newname) || length(oldname)~=length(newname)
        error('lists must be same length');
    end
elseif ~iscell(oldname) && ~iscell(newname)
    oldname = {oldname};
    newname = {newname};
else 
    error('improper parameters');
end

flds = fieldnames(S);

for i = 1:length(oldname) 
    f = getfield(S, oldname{i});
    S = setfield(S, newname{i}, f);
    if ~strcmp(oldname{i}, newname{i})
        S = rmfield(S, oldname{i});
    end
    idx = find(strcmp(flds,oldname{i}));
    if length(idx)~=1
        error('unexpected behavior');
    end
    flds{idx} = newname{i};
end
S = order_fields_first(S, unique_keepords(flds)); % See order_fields_first code in repository
end