function S = load_struct(varargin)
% load_struct(filename, format, separator_character, header_lines,
% <other_parameters>, P)

%% Evaluate input parameters
args = varargin;
if length(args) < 1
    error('filename required');
end
filename = args{1};
if ~ischar(filename)
    error('first parameter should be filename as a character string');
end


% See if last parameter is a structure (P)
if isstruct(args{end})
    P = args{end};
    args = args(1:end-1);
else
    P = [];
end
P = impose_default_value(P, 'lowercase_fieldnames', false);
P = impose_default_value(P, 'ignore_poundsign_lines', true);
% See code for impose_default_value in repository

if length(args) >= 2
end
if length(args) >= 3
    if isempty(args{3})
        error('third parameter should be a separator character');
    end
    if isnumeric(args{3})
        error('header_lines is now fourth parameter...update call to load_struct');
    end
end
if length(args) >= 4
    if islogical(args{4})
        error('lowercase_fieldnames has been moved to P structure...update call to load_struct');
    end
    if ~isnumeric(args{4})
        error('fourth parameter should be number of header lines');
    end
end

%% Take care of comment lines
default_header_lines = 1;
% N.B. if table has comment lines at the beginning, increment header_lines
% to skip them
n_comment_lines = 0;
if P.ignore_poundsign_lines
    fid = fopen(args{1});
    while (true)
        x = fgetl(fid);
        if (isempty(x)) || (x(1)=='#') %Skipping empty lines or lines starting with '#'
            n_comment_lines = n_comment_lines + 1;
            continue
        elseif strncmp(x, 'Oncotator v', 11)
            fprintf('Un-poundsigned Oncotator header detected and skipped.\n');
            n_comment_lines = n_comment_lines + 1;
            continue
        else 
            break
        end
    end
    fclose(fid);
end

default_header_lines = 1 + n_comment_lines;
% Default args
if length(args) == 1
    args = [args {''}];
end
if length(args) == 2
    args = [args {char(9)}];
end
if length(args) == 3
    args = [args {default_header_lines}];
end
% Default whitespace
has_already = false;
for i = 4:length(args)
    if ischar(args{i}) && strcmpi(args{i}, 'whitespace')
        has_already=true;
        break;
    end
end
if ~has_already
    args = [args {'whitespace'} {'\b\r'}];
end

%% Load table
try
    table = read_table(args{:}); % See read_table code in repository
    nf = length(table.dat);
catch me
    in = fopen(args{1});
    t = fread(in, 'uint8=>char')';
    fclose(in);
    while ~isempty(t)
        if t(end) == 10
            t(end) = []; % remove trailing blank lines
        else 
            break;
        end
    end
    if isempty(t)
        q = {};
    else 
        q = split(t, char(10));
    end
    if isempty(q)
        fprintf('\n%s is a blank file\n', args{1});
        table = [];
        table.dlm = args{3};
        table.headers = {{}};
        table.dat = {{}};
        nf = 0;
    else
        disp(me)
        disp(me.message);
        error('Error loading struct file');
    end
end
if isempty(table.headers)
    table.headers{1} = cell(nf, 1);
    for f = 1:nf
        table.headers{1}{f} = sprintf('col%d', f);
    end
end
% Process header line

fields = table.headers{end};
if length(fields)~=nf
    fprintf('Header line has %d column headers instead of the expected %d:\n', length(fields), nf);
    fields{:}
    error('Unable to parse table header line.');
end

% Remove any characters except A-Z, a-z, 0-9, underscore from column headings
if P.lowercase_fieldnames
    fields = lower(fields);
end
fields = regexprep(fields, '\W', '');
fields_orig = fields;
fields = genvarname(fields_orig);

for f = 1:nf
    if strcmp(fields_orig(f), 'end')
        fields{f} = 'end';
        break
    end
    if strcmp(fields_orig(f), 'End')
        if P.lowercase_fieldnames
            fields{f} = 'end';
        else 
            fields{f} = 'End';
        end
        break;
    end
end

S = struct();
for f = 1:nf
    S.(fields(f)) = table.dat{f};
end
end