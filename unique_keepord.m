function [u ui uj] = unique_keepord(x, varargin)
if exist('varargin', 'var') && length(varargin)>=1 && ischar(varargin{1}) && (strcmpi(varargin{1}, 'first')||strcmpi(varargin{1}, 'last'))
    error('Do not specify "first" or "last" with this function. Default is set to "first")');
end

[u1 ui1 uj1] = unique(x, 'first', varargin{:});

[ui ord] = sort(ui1);
u = x(ui1(ord));
[tmp ord2] = sort(ord);
uj = ord2(uj1);
end