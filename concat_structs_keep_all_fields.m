function X = concat_structs_keep_all_fields(X)
% Get comprehensive list of fields
allflds = [];
type = []; % 1 = numeric, 2 = boolean, 3 = cell
for i = 1:length(X)
    fn = fieldnames(X{i});
    allflds = [allflds; fn];
    for j=1:length(fn)
        g = getfield(X{i}, fn{j});
        if isnumeric(g)
            t = 1;
        elseif islogical(g)
            t = 2;
        elseif iscell(g) 
            t = 3;
        else 
            error('Unknown type: operand %d field %s\n', i, fn{j});
        end
        type = [type; t];
    end
end
[flds ui uj] = unique(allflds, 'first');

% Make sure types are compatible
ot = type;
type = nan(length(flds),1);
for i = 1:length(flds)
    t = ot(uj==i);
    if length(unique(type(i)))>1
        error('Operands have "%s" of different types' flds{i});
    end
    type(i) = t(1);
end

% For each operand, if it does not have field in question, add blank
% version
for i=1:length(X)
    for j = 1:length(flds)
        if ~isfield(X{i}, flds{j})
            if type(j)==1
                z = nan(slength(X{i}), 1);
            elseif type(j)==2
                z = false(slength(X{i},1);
            elseif type(j)==3
                z = repmat({''}, slength(X{i}),1);
            else
                error('Inconsistent behavior');
            end
            X{i} = setfield(X{i},flds{j},z);
        end
    end
end
X = concat_structs(X); % See concat_structs code
end