function tab=read_table(fname,format,dlm,headerlines,varargin)
  fpos=0;
  if headerlines~=0
    if ischar(fname)
      f=fopen(fname,'r');
    else
      f=fname;
      fpos=ftell(f);
    end
    if length(dlm)~=1 
        tab.dlm = find_dlm(fname,dlm); % See find_dlm code in repository
    else
      tab.dlm=dlm;
    end
    if headerlines>0
      tab.headers=read_dlm(f,dlm,headerlines); % See read_dlm code in repository
    elseif headerlines==-1 
      headerlines=1;
      tab.headers=read_dlm(f,dlm,headerlines);
      tab.headers{1}=['EMPTY' tab.headers{1,:}];
    end
  else
    if ischar(fname)
      f=fopen(fname,'r');
    else
      f=fname;
      fpos=ftell(f);
    end
    tab.headers={};
    tab.dlm=dlm;
  end
  if isempty(format)
    if isempty(tab.headers)
      error('must have either format or headerlines');
    else
      format=[repmat('%s',1,length(tab.headers{end})) '\n'];   % (to allow for multiple header lines)
    end
  elseif iscell(format)
    if isempty(tab.headers)
      error('must have either format or headerlines');
    else
      if length(format)==1
        format=[repmat(format{1},1,length(tab.headers{end})) '\n'];
      else
        format=[format{1} repmat(format{3},1,length(tab.headers{end})-format{2}) '\n'];
      end
    end
  end

  if strcmp(format((end-1):end),'\n')
    format=format(1:(end-2));
  end
  
  verbose(['Reading file using format:' format],10); % See verbose code in repository
  fseek(f,fpos,'bof');
  tab.dat=textscan(f,format,'headerLines',headerlines,'delimiter',tab.dlm,varargin{:});
  fclose(f);
end