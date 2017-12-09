function Y = nansub(X, idx, filler)
if numel(X)==2 && size(X,2)>1
    fprintf('Converting first argument to column vector\n');
    X = X';
end

if iscellstr(X) && size(X, 1) ==1 && size(X,2) >1
    X = X';
end 

if islogical(X)
    type = 0;
elseif isnumeric(X)
    type = 1;
elseif iscell(X)
    type = 2;
else
    error('Array type not supported');
end

if ~exist('filler', 'var')
    if type==0
        filler = false;
    elseif type==1
        filler = nan;
    elseif type==2
        filler = {''};
    else 
        error('Inconsistent behavior with "type"');
    end
end

if type==0
    if ~islogical(filler)
        error('Inappropriate filler for logical array');
    end
elseif type==1
    if ~isnumeric(filler)
        error('Inappropriate filler for numeric array');
    end
elseif type==2
    if ischar(filler)
        filler = {filler};
    end 
    if ~ischar(filler)
        error('Inappropriate filler for cell array');
    end
else
    error('Inconsistent behavior with "type"');
end

sz = size(X)
sz(1) = length(idx);
Y = repmat(filler, sz);
idx2 = find(~isnan(idx) & idx>=1 & idx<=length(X));
Y(idx2, :, :, :, :) = X(idx(idx2), :, :, :, :);
end