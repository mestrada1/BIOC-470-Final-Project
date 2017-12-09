function Mut_Detect_Significance_Preprocess(mut_file, coverage_file, covariate_file, output_file)
% Preprocessing code; Here I attempt to ensure that files are properly
% loaded and filter out any faulty arguments within file structures
params1 = 'mut_file, coverage_file, covariate_file, output_file';
if nargin < 4
    error(['usage: Mut_Detect_Preprocess(' params1 ')']); 
end

% Ensure that output_file exists and that it is writeable
[outpath] = fileparts(output_file);
if ~isempty(outpath) && ~exist(outpath, 'dir')
    mkdir(outpath);
end
ensure_writeable([output_file '.test.txt']); % See ensure_writeable code 

% Load mut file and ensure it has gene+patient field info
fprintf('Loading mut_file...\n');
M = load_struct(mut_file); % See load_struct code 

% Checking gene field information
if isfield(M, 'gene') && isfield(M, 'Hugo_Symbol')
    fprintf('Both "gene" and "Hugo_Symbol" present in mut_file');
elseif isfield(M, 'gene')
elseif isfield(M, 'Hugo_Symbol')
    M.gene = M.Hugo_Symbol;
else 
    error('mut_file is missing "gene" or "Hugo_Symbol" column');
end

% Checking patient field information
if isfield(M, 'patient') && isfield(M, 'Tumor_Sample_Barcode')
    fprintf('Both "patient" and "Tumor_Sample_Barcode" are present in mut_file');
elseif isfield(M, 'patient')
elseif isfield(M, 'Tumor_Sample_Barcode')
    M.patient = M.Tumor_Sample_Barcode;
else 
    error('mut_file is missing "patient" or "Tumor_Sample_Barcode" column');
end
if length(unique(M.patient)) < 2
    error('Need more than 1 patient for Mut_Detect_Significance Algorithm');
end

% Load coverage_file and ensure it has gene, effect, and categ information fields
fprintf('Loading coverage_file...\n');
C = load_struct_specify_string_cols(coverage_file, 1:3); % See load_struct_specify_string_cols code in repository
if ~isfield(C, 'gene')
    error('coverage_file is missing "gene" column');
end
if ~isfield(C, 'effect') && isfield(C, 'zone')
    C.effect = C.zone;
end
if ~isfield(C, 'effect')
    error('coverage_file is missing "effect" column');
end
C.effect = regexprep(C.effect, '^flank.*', 'noncoding');
if any(~ismember(unique(C.effect), {'noncoding', 'silent', 'nonsilent'}))
    error('"effect" must be one of noncoding/silent/nonsilent in coverage_file');
end
if ~isfield(C, 'categ')
    error('coverage_file is missing "categ" column');
end
ff = fieldnames(C); coverage_patient_names = ff(4:end);

% Ensure covariates file exists
if ~exist(covariate_file, 'file')
    error('covariates_file not found'); 
end

%% Process effect column in mut_file; remove 'is_coding' and 'is_silent' fields

fprintf('Processing mutation "effect"...\n');
if isfield(M, 'is_coding') || isfield(M, 'is_silent')
    f = {'is_coding', 'is_silent'};
    for i = 1:length(f)
        if isfield(M, f{i})
            M = rmfield(M, f{i}); 
        end
    end
end

if isfield(M, 'effect')
    fprintf('Will use pre-existing "effect" column.\n');
    M.effect = regexprep(M.effect, '^flank.*', 'noncoding');
    if any(~ismember(unique(M.effect), {'noncoding', 'silent', 'nonsilent', 'null'}))
        error('in mut_file, "effect" must be one of noncoding/silent/nonsilent/null.');
    end
else % Calculating effect column from Variant Classification/Type column 
    if ~isfield(M, 'Variant_Classification') && isfield(M, 'type')
        M.Variant_Classification = M.type;
    end
    if ~isfield(M, 'Variant_Classification')
        error('mut_file is missing Variant_Classification');
    end
end
%% Process categ column in mut_file

fprintf('Processing mutation "categ"...\n');

% See if coverage file is on the full 192 mutation set; 192 mutation set is
% expanded from 96 substitution types. E.g. A -> T mutation in N1_A_N2 (where
% N1 and N2 are flanking nucleotides) would be classified as two mutations
% if you consider this mutation on the transcribed AND non-transcribed
% strands

ucc = unique(C.categ);
names192 = generate192names(); % See code for generate192names in repository
coverage_is_on_full192 = (length(ucc)==192 && length(intersect(ucc,names192))==192);
categs_already_present = false;
% In case categ is already in M structure, ensure categories in C are the
% same as in M
if isfield(M, 'categ')
    mcc = unique(M.categ);
    ucc = unique(C.categ);
    if length(ucc)~=length(mcc) || any(~ismember(ucc, mcc)) || any(~ismember(mcc, ucc))
        fprintf('"categ" of mut_file does not match coverage_file. Ignoring it.\n');
        M = rmfield(M, 'categ');
    else
        categs_already_present = true;
    end 
end
% Potential methods to implement to decide what to do about categories
% Method 1: Use existing categories
% Method 2: Have only one category for non-nulls

if categ_flag==0 % Requested use of categories already present (categ_flag is an input for this preprocessing code)
    if ~categs_already_present
        error('when setting categ_flag==0, "categ" column must be already present in mut_file.');
    end
    method = 1;
elseif categ_flag==1
    method = 2; % Requested only one category for non-nulls
elseif isnan(categ_flag)
    % No method specified: decide what to do
    if categs_already_present
        method = 1; 
    else
        method = 2; % Will put all non-null in one category called "missense"
    end
end
if method==1
    fprintf('Will use the categories already present.\n');
    K=[];
    K.name = unique(M.categ);
elseif method==2
    fprintf('Will use two categories: missense and null+indel.\n');
    K1 = [];
    K1.left = {'ACGT'};
    K1.from = {'AC'};
    K1.change = {'in'};
    K1.right = {'ACGT'};
    K1.autoname = {'missense'};
    K1.name = {'missense'};
    K1.type = {'point'};
    K2 = [];
    K2.left = {'ACGT'};
    K2.from = {'AC'};
    K2.change = {'in'};
    K2.right = {'ACGT'};
    K2.autoname = {'null+indel'};
    K2.name = {'null+indel'};
    K2.type = {'non-point'};
    K = concat_structs_keep_all_fields({K1,K2}); % See concat_structs_keep_all_fields code in repository
    % Assign categories
    M.categ = nansub(K.name, 1+strcmp('null', M.effect)); % See nansub code in repository
    % Collapse coverage
    fprintf('Collapsing coverage...\n');
    [tmp tmp C.categ_idx] = unique(C.categ);
    C = sort_struct(C, {'gene', 'effect', 'categ_idx'});
    ug = unique(C.gene);
    ng = length(ug);
    ue = unique(C.effect);
    ne = length(ue);
    nk = slength(K);
    idx = find(C.categ_idx<=nk);
    C2 = reorder_struct(C,idx);
    C2.categ = nansub(K.name, C2.categ_idx);
    C2 = keep_fields(C2, {'gene', 'effect', 'categ'}); % See keep_fields code in repository
    np = length(coverage_patient_names);
    for p=1:np
        oldcov = reshape(C.(coverage_patient_names{p}), [192 ne ng]);
        newcov = repmat(sum(oldcov,1),2,1);
        C2.(coverage_patient_names{p}) = newcov(:);
    end
    C=C2;
    clear C2;
  end
  
  % Save output files
  
  fprintf('Writing preprocessed files.\n');
  
  % (1) mutation file
  f = {'newbase', 'chr_idx', 'triplet', 'yname', 'context65','newbase_idx','context65_right','triplet_middle'};
  for i = 1:length(f)
      if isfield(M, f{i})
          M = rmfield(M, f{i});
      end
  end
  save_struct(M,[output_filestem '.mutations.txt']); % See save_struct code in repository; this function saves the fully processed mut_file, as well as the coverage and categories files
  
  % (2) coverage file
  save_struct(C,[output_filestem '.coverage.txt']);
  
  % (3) categories file
  save_struct(K,[output_filestem '.categs.txt']);
  
  fprintf('Mut_Detect_Significance_Preprocess finished.\n');
  
end % end of Mut_Detect_Significance_Preprocess