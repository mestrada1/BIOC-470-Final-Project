function ensure_writeable(fname) 
try 
    [path name ext] = fileparts(fname);
    if ~isempty(path) 
        if ~iscell(path)
            if ~exist(path, 'dir')
                mkdir(path)
            end
        else
            for i = 1:numel(path)
                if ~exist(path{i}, 'dir')
                    mkdir(path{i});
                end
            end
        end
    end
    testfile = [fname '.' rand() '.test'];
    out = fopen(testfile, 'wt');
    fwrite(out, 'test', 'char');
    fclose(out);
    delete(testfile);
catch me 
    error('%s is not writeable', fname);
end