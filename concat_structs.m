function S = concat_structs(slist) 
% slist should be cell array of structures, all of which have the same
% fields; this returns a new structure in which each field is a
% concatenation of the corresponding fields from each of the input
% strucutures
if ~iscell(slist)
    error('Input must be a cell array of structures');
end

if length(slist)==0
    S={};
elseif length(slist)==1
    S = slist{1};
else
    ns = numel(slist);
    allflds = cell(ns,1);
    for i = 1:ns
        if isempty(slist{i})
            continue;
        end
        if ~isstruct(slist{i})
            error('Input must be a cell array of structures');
        end
        allflds{i} = fieldnames(slist{i});
    end
    allflds = cat(1, allflds{:});
    [flds ai aj] = unique_keepord(allflds); %See unique_keepord code
    h = histc(aj, 1:length(flds));
    if ~all(h==ns)
        count(allflds);
        error('All input structures must have same fields. Use concat_structs_keep_all_fields.');
    end
    [tmp ord] = sort(ai);
    
    rflag = false;
    S = [];
    for fno=1:length(flds)
        type = nan(ns, 1);
        f = cell(ns,1);
        for i=1:ns
            if isempty(slist{i})
                continue;
            end
            f{i} = getfield(slist{i}, flds{fno});
            if isnumeric(f{i}) || islogical(f{i})
                type(i) = 1;
            elseif iscell(f{i})
                type(i) = 2;
            else type(i) = 3;
            end
        end
        if any(type==3)
            error('Incompatible type encountered in %s', flds{fno});
        end
        if any(type==1) && any(type==2)
            if ~rflag
                fprintf('Reconciling mixed cell+numeric fields for %s\n', flds{fno});
                rflag=true;
            end
            for i =1:ns
                if type(i)==2
                    try
                        f{i} = str2double(f{i});
                    catch me
                        error('Unable to resolve mixed cell+numeric case encountered in %s', flds{fno});
                    end
                end
            end
        end
        S = setfield(S, flds{fno}, cat(1,f{:}));
    end
end
end