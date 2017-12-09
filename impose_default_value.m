function S= impose_default_value(S, field, value, acceptable_values)
if ~isfield(S, field) || isempty(getfield(S, field))
    if ischar(value) && strcmp(value, '*required*')
        error('%s is a required field of P', field);
    else
        try
            S = setfield(S, field, value);
        catch me 
            fprintf('Error setting field "%s"\n', field);
            disp(me); disp(me.message);
        end
    end
end
if exist('acceptable_values', 'var')
    av = acceptable_values;
    v = getfield(S, field);
    if ischar(v) && isnumeric(av)
        v = str2double(v);
        S = setfield(S, field, v);
    end
    if ~ischar(v) && ischar(av)
        error('%s is assigned value of different type than accepted_values', field);
    end
    try 
       ism = ismember(v, av);
    catch me 
        error('%s is assigned value of different type than accepted_values', field);
    end
    if ~ism
        fprintf('Acceptable values for %s:\n', field);
        disp(av);
        fprintf('Attempted to set to:\n');
        disp(v);
        error('Invalid setting of %s', field);
    end 
end 
end