function verbose(str,level,varargin)

  global VERBOSE_LEVEL
  global VERBOSE_FILE

  if nargin==1
    level=1;
  end
  
  if isempty(varargin)
    
    str = char(regexprep(cellstr(str),'%','%%'));
    if ~isempty(VERBOSE_LEVEL) && (level<=VERBOSE_LEVEL)
      fprintf(1,[str repmat('\n',size(str,1),1)]');  %escape the % to prevent line from commenting
      if ~isempty(VERBOSE_FILE)
        if ~exist(VERBOSE_FILE,'file')
          fid = fopen(VERBOSE_FILE,'w');
        else
          fid = fopen(VERBOSE_FILE,'a');
        end
        
        fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
        fprintf(fid,'\n');
        fclose(fid);
      end
    end
    
  else
    
    if ~isempty(VERBOSE_LEVEL) && (level<=VERBOSE_LEVEL)
      fprintf(1,str,varargin{:});  %escape the % to prevent line from commenting
      fprintf(1,'\n')
      if ~isempty(VERBOSE_FILE)
        fid = fopen(VERBOSE_FILE,'a');
        fprintf(fid,str,varargin{:});  %escape the % to prevent line from commenting
        fprintf(fid,'\n');
        fclose(fid);
      end
    end
  end
end