function S2 = keep_fields(S, flds)
if ischar(flds)
    flds = {flds};
end
S2=[];
for i=1:length(flds)
  if isempty(S)
    f = [];
  else
    f = getfield(S,flds{i});
  end
  S2=setfield(S2,flds{i},f);
end
end