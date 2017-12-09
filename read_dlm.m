function [fl, fid] = read_dlm(fname,dlm,nlines)
% This function reads all/some of a delimited file into cell arrays of
% arrays; fl is a cell array of cell arrays of strings while fl{n} is a
% cell array stings containing the nth line of the file. FID is the fid fro
% the opened file, fname.
if nargin==1
    dlm=9;
  end
  
  if ischar(fname)
    fid=fopen(fname,'r');
  else
    fid=fname;
  end
  
  ln=1;
  if ~exist('nlines','var')
    nlines=Inf;
    do_close=1;
  else
    do_close=0;
  end
  
  had_output=0;
  fl={};
  while(ln<=nlines)
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    pos1=find(ismember(tline,dlm));
    if ~isempty(pos1)
        pos1=[ 0 pos1 length(tline)+1];
        for i=1:(length(pos1)-1)
        fl{i}=tline((pos1(i)+1):(pos1(i+1)-1));
        end
    else
        fl{1}=tline;
    end
    ln=ln+1;
    if mod(ln,1000)==0
      verbose(['...' num2str(ln)],30); % See verbose code in repository
      had_output=1;
    end
  end
  if do_close
    fclose(fid);
    fid=-1;
  end
  ln=ln-1;
  if had_output
    verbose([newline],30);
  end
end