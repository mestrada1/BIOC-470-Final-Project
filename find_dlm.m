function d = find_dlm(fname,dlm)
 if ~exist('dlm', 'var') || isempty(dlm)
          dlm = [ char(9) ',|'];
      end
      f = fopen(fname, 'r');
      l = fgetl(f);
      for i = 1:length(dlm)
          h(i) = length(find(l==dlm(i)));
      end
      [hm, hi] = max(h);
      d = dlm(hi);
end