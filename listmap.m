function m = listmap(a, b)

if ischar(a)
    a={a};
end
if ischar(b)
    b={b};
end
m = nan(length(a), 1);
[a ai aj] = unique(a);
[c ia ib] = intersect(a,b);
for i = 1:length(c)
    if ~mod(i,1e5)
        fprintf('%d/%d ', i, length(c));
    end
    m(aj==ia(i)) = ib(i);
end
if length(c) >= 1e5
    fprintf('\n');
end
end