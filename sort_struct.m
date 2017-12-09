function [s2,ord]=sort_struct(s1,keyfield,order)
  if length(keyfield)==0
      return; 
  end
  if ~iscell(keyfield)
    keyfield = {keyfield};
  end
  if ~exist('order','var')
    order = repmat(1,length(keyfield),1);
  end
  if length(order) ~= length(keyfield)
    error('order and keyfield must have same number of elements');
  end
  
  orig_len = slength(s1);
  s2=s1;
  ord=(1:orig_len)';
  fields = fieldnames(s1);
  nf = length(fields);
  
  for k=length(keyfield):-1:1
    f = getfield(s2,keyfield{k});
    if length(f)<orig_len
        error('Attempted to sort on truncated field "%s"',keyfield{k}); 
    end
    if order(k)==1
      if isnumeric(f)
        [tmp ordi] = sortrows(f);
      else
        [tmp ordi] = sort(f);
      end
    elseif order(k)==-1
      if isnumeric(f)
        [tmp ordi] = sortrows(f,-1);
      else
        [tmp ordi] = sort(flip(f));
      end
    else
      error('Unknown order %d',order(k));
    end
    for i=1:nf
      f = getfield(s2,fields{i});
      f = f(ordi,:,:,:,:,:,:,:,:,:);
      s2 = setfield(s2,fields{i},f);
    end
    ord = ord(ordi);
  end
end