function [s order] = reorder_struct(s, order)
if nargin~=2
    error('reorder_struct(s, order)');
end
if islogical(order)
    order = find(order);
end
if ischar(order)
    if strcmpi(order, 'end')
        order = slength(s);
    else
        error('Invalid index parameter');
    end
end
if size(order,2) > 1 
    order = order';
end
nanflag = any(isnan(order));
fields = fieldnames(s);
nf = length(fields);

for i = 1:nf
    f = getfield(s, fields{i});
    if nanflag
        f = nansub(f, order); % See code for nansub
    else
        f = f(order, :, :, :, :, :, :, :, :, :);
    end
    s = setfield(s, fields{i}, f);
end
end