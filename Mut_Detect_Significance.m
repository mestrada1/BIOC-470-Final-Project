function Mut_Detect_Significance(mut_file,coverage_file,covariate_file,output_file)

  if nargin~=4, error('usage: Mut_Detect_Significance(mut_file,coverage_file,covariate_file,output_file)'); end

  % First ensure output_file directory exists
  [outpath] = fileparts(output_file);
  if ~isempty(outpath) && ~exist(outpath,'dir')
      mkdir(outpath); 
  end
  
  % Load mut_file and make sure all required fields are present

  fprintf('Loading mut_file...\n');
  M = load_struct(mut_file); % See load_struct code in repository
  
  % Gene Field Check
  if isfield(M,'gene') && isfield(M,'Hugo_Symbol')
    fprintf('NOTE:  Both "gene" and "Hugo_Symbol" are present in mut_file.  Using "gene".\n');
  elseif isfield(M,'gene')
  elseif isfield(M,'Hugo_Symbol')
    M.gene = M.Hugo_Symbol;
  else
    error('mut_file lacks "gene" or "Hugo_Symbol" column.');
  end
  
  % Patient Field Check
  if isfield(M,'patient') && isfield(M,'Tumor_Sample_Barcode')
    fprintf('NOTE:  Both "patient" and "Tumor_Sample_Barcode" are present in mut_file.  Using "patient".\n');
  elseif isfield(M,'patient')
  elseif isfield(M,'Tumor_Sample_Barcode')
    M.patient = M.Tumor_Sample_Barcode;
  else
    error('mut_file lacks "patient" or "Tumor_Sample_Barcode" column.');
  end
  
  % Effect Field Check
  if isfield(M,'is_coding') || isfield(M,'is_silent')
    fprintf('Ignoring "is_coding" and "is_silent" fields. Looking for "effect" column.\n');
    f = {'is_coding', 'is_silent'};
    for i = 1:length(f)
        if isfield(M, f{i})
            M = rmfield(M, f{i}); 
        end
    end
  end
  if isfield(M,'effect')
    M.effect = regexprep(M.effect,'^flank.*','noncoding');
    if any(~ismember(unique(M.effect),{'noncoding','silent','nonsilent','null'}))
      error('in mut_file, "effect" must be one of noncoding/silent/nonsilent/null');
    end
  else
    error('mut_file lacks "effect" column.');
  end
  
  % Categ Field Check
  if ~isfield(M,'categ')
      error('mut_file lacks "categ" column.'); 
  end
  
  % Load coverage_file and covariate_file
  fprintf('Loading coverage_file...\n');
  C = load_struct_specify_string_cols(coverage_file,1:3); % See load_struct_specify_strings_cols code in repository 
  % gene, effect, categ are all strings
  G=[]; 
  G.gene = unique(C.gene); 
  ng=slength(G); % see slength code in repository
  fprintf('Loading covariate_file...\n');
  V = load_struct_specify_string_cols(covariate_file,1);  % gene is string
  f = fieldnames(V); 
  cvnames = f(2:end); 
  nv = length(cvnames);
  gidx = listmap(G.gene,V.gene); % See listmap code in repository
  for i=1:length(f)
    if strcmp(f{i},'gene')
        continue; 
    end
    G.(f{i}) = nansub(V.(f{i}),gidx); % See nansub code in repository
  end
  
  % Ensure that coverage_file has required fields
  if ~isfield(C,'gene')
      error('no "gene" column in coverage_file'); 
  end
  if ~isfield(C,'effect') && isfield(C,'zone')
      C = rename_field(C,'zone','effect'); % See rename_field code in repository
  end
  if ~isfield(C,'effect')
      error('no "effect" column in coverage_file'); 
  end
  C.effect = regexprep(C.effect,'^flank.*','noncoding');
  if any(~ismember(unique(C.effect),{'noncoding','silent','nonsilent'}))
    error('in coverage_file, "effect" must be one of noncoding/silent/nonsilent');
  end
  if ~isfield(C,'categ')
      error('no "categ" column in coverage_file'); 
  end
  f = fieldnames(C); 
  coverage_patient_names = f(4:end);
  
  % Remove any genes where there is no coverage for
  badgene = setdiff(M.gene,C.gene);
  if ~isempty(badgene)
    fprintf('NOTE:  %d/%d gene names could not be mapped to coverage information.  Excluding them.\n',length(badgene),length(unique(M.gene)));
    order = ismember(M.gene, badgene);
    if islogical(order)
        order = find(order);
    end
    M = reorder_struct(M, setdiff(1:slength(M), order));
  end

  % Ensure categories in C are same as in M
  bad = find(~ismember(M.categ,C.categ));
  if ~isempty(bad)
    fprintf('NOTE:  %d/%d mutations were outside the category set.  Excluding them.\n',length(bad),slength(M));
    if isfield(M,'ref_allele') && isfield(M,'newbase') && isfield(M,'type')
      is_probably_indel = strcmp('-',M.ref_allele(bad)) | strcmp('-',M.newbase(bad));
      is_probably_noncoding = grepmi('intron|utr|igr|flank',M.type(bad));
      nbad2 = bad(is_probably_indel&is_probably_noncoding);
      if nbad>0
        fprintf('           (%d of them are noncoding indels.)\n',length(bad2));
      end
    end
    if islogical(bad)
        bad = find(bad);
    end
    M = reorder_struct(M, setdiff(1:slength(M), bad));
  end
  if slength(M)==0
      error('No mutations left!\n');
  end

  % Map categories
  [K.name tmp C.categ_idx] = unique(C.categ);
  M.categ_idx = listmap(M.categ,K.name);
  ncat = slength(K);

  % Ensure that there is a null+indel category
  knum = str2double(K.name);
  if all(~isnan(knum))
    % all category names are numbers: assume null+indel is the last one
    [tmp null_categ] = max(knum);
  else
    null_categ = grepi('null|indel',K.name,1);
    if length(null_categ)==0
        error('Error: no null/indel category.\n'); 
    end
    if length(null_categ)>1
        error('Error: multiple null/indel categories.\n'); 
    end
  end

  % Ensure that C is sorted by the same gene order as in G
  C.gene_idx = listmap(C.gene,G.gene);
  C = sort_struct(C,'gene_idx'); % See sort_struct code in repository
  
  % Map genes
  M.gene_idx = listmap(M.gene,G.gene);
  
  % Regularize the sample name in the mutation table
  namebefore = M.patient;
  M.patient = regexprep(M.patient,'-Tumor$','');
  if any(~strcmp(namebefore,M.patient))
      fprintf('NOTE:  Trimming "-Tumor" from patient names.\n'); 
  end
  namebefore = M.patient;
  M.patient = regexprep(M.patient,'-','_');
  if any(~strcmp(namebefore,M.patient))
      fprintf('NOTE:  Converting "-" to "_" in patient names.\n'); 
  end

  pat=[]; 
  [pat.name tmp M.patient_idx] = unique(M.patient);
  pat.cov_idx = listmap(pat.name,coverage_patient_names);
  np = slength(pat);
  if np<2, error('Mut_Detect_Significance wont work with a single patient.\n'); 
  end

  % Generic coverage data given?
  generic_column_name = 'coverage';
  if length(coverage_patient_names)>1 || (length(coverage_patient_names)==1 && ~strcmpi(coverage_patient_names{1},generic_column_name))
    if any(strcmp(coverage_patient_names,generic_column_name))
        error('reserved name "%s" cannot appear in list of patient names',generic_column_name); 
    end
    if length(coverage_patient_names)~=length(unique(coverage_patient_names))
        error('patient names in coverage_file must be unique'); 
    end
    % Ensure all patients are accounted for in coverage_file
    if any(isnan(pat.cov_idx))
        error('some patients in mut_file are not accounted for in coverage_file'); 
    end
    generic_coverage_flag = false;
  else
    % using generic coverage
    pat.cov_idx(:) = 1;
    generic_coverage_flag = true;
  end
  
  % Construct n and N tables
  fprintf('Building n and N tables...\n');
  
  midx = strcmpi(M.effect,'silent');
  n_silent = hist3d(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  midx = strcmpi(M.effect,'nonsilent') | strcmpi(M.effect,'null');
  n_nonsilent = hist3d(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  midx = strcmpi(M.effect,'noncoding');
  n_noncoding = hist3d(M.gene_idx(midx),M.categ_idx(midx),M.patient_idx(midx),1,ng,1,ncat,1,np);
  
  N_silent = nan(ng,ncat,np);
  N_nonsilent = nan(ng,ncat,np);
  N_noncoding = nan(ng,ncat,np);
  
  for ci=1:ncat
    silent_idx = strcmpi(C.effect,'silent') & C.categ_idx==ci;
    nonsilent_idx = strcmpi(C.effect,'nonsilent') & C.categ_idx==ci;
    noncoding_idx = strcmpi(C.effect,'noncoding') & C.categ_idx==ci;
    for pi=1:np
      cpfld = coverage_patient_names{pat.cov_idx(pi)};
      N_silent(:,ci,pi) = C.(cpfld)(silent_idx);
      N_nonsilent(:,ci,pi) = C.(cpfld)(nonsilent_idx);
      N_noncoding(:,ci,pi) = C.(cpfld)(noncoding_idx);
    end
  end

  % Ensure all numbers are integers
  n_silent = round(n_silent);
  n_nonsilent = round(n_nonsilent);
  n_noncoding = round(n_noncoding);
  N_silent = round(N_silent);
  N_nonsilent = round(N_nonsilent);
  N_noncoding = round(N_noncoding);

  % Remove mutations with very low coverage
  n_silent(n_silent>N_silent) = 0;
  n_nonsilent(n_nonsilent>N_nonsilent) = 0;
  n_noncoding(n_noncoding>N_noncoding) = 0;
  
  % Checking totals
  tot_n_nonsilent = fullsum(n_nonsilent); % See fullsum code in repository
  tot_N_nonsilent = fullsum(N_nonsilent);
  tot_n_silent = fullsum(n_silent);
  tot_N_silent = fullsum(N_silent);
  tot_n_noncoding = fullsum(n_noncoding);
  tot_N_noncoding = fullsum(N_noncoding);
  tot_rate_nonsilent = tot_n_nonsilent/tot_N_nonsilent;
  tot_rate_silent = tot_n_silent/tot_N_silent;
  tot_rate_noncoding = tot_n_noncoding/tot_N_noncoding;
  tot_rate_coding = (tot_n_nonsilent+tot_n_silent)/(tot_N_nonsilent+tot_N_silent);
  
  min_tot_n_nonsilent = 50;
  min_tot_n_silent = 50;
  min_tot_n_noncoding = 50;
  min_rate_nonsilent = 1e-9;
  max_rate_nonsilent = 1e-3;
  min_rate_silent = 1e-9;
  max_rate_silent = 1e-3;
  min_rate_noncoding = 1e-9;
  max_rate_noncoding = 1e-3;
  max_abs_log2_difference_nonsilent_silent = 1.0;
  max_abs_log2_difference_noncoding_coding = 1.0;
  
  % Check if silent and nonsilent are ok, give warning otherwise
  if tot_n_nonsilent<min_tot_n_nonsilent || tot_n_silent<min_tot_n_silent
      error('not enough mutations to analyze'); 
  end
  if tot_rate_nonsilent<min_rate_nonsilent || tot_rate_nonsilent>max_rate_nonsilent
      error('nonsilent mutation rate out of range'); 
  end
  if tot_rate_silent<min_rate_silent || tot_rate_silent>max_rate_silent
      error('silent mutation rate out of range'); 
  end
  abs_log2_difference_nonsilent_silent = abs(log2(tot_rate_nonsilent/tot_rate_silent));
  if abs_log2_difference_nonsilent_silent>max_abs_log2_difference_nonsilent_silent
      fprintf('Warning: silent and nonsilent rates are too different.\n'); 
  end
  
  % Check if noncoding is ok, give warning otherwise and zero it all out
  ok = false;
  if tot_n_noncoding==0
    fprintf('NOTE:  no noncoding mutations.\n');
  else
    if tot_n_noncoding<min_tot_n_noncoding
      fprintf('Warning:  not enough noncoding mutations to analyze\n');
    else
      if tot_rate_noncoding<min_rate_noncoding || tot_rate_noncoding>max_rate_noncoding
        fprintf('Warning:  noncoding mutation rate out of range\n');
      else
        abs_log2_difference_noncoding_coding = abs(log2(tot_rate_noncoding/tot_rate_coding));
        if abs_log2_difference_noncoding_coding>max_abs_log2_difference_noncoding_coding
          fprintf('Warning:  coding and noncoding rates are too different\n');
        else
          ok = true;
        end
      end
    end
  end
  if ~ok
    fprintf('Zeroing out all noncoding mutations and coverage for rest of calculation.\n');
    n_noncoding(:) = 0;
    N_noncoding(:) = 0;
  end
  
  % Add total columns
  n_silent(:,end+1,:) = sum(n_silent,2);
  n_nonsilent(:,end+1,:) = sum(n_nonsilent,2);
  n_noncoding(:,end+1,:) = sum(n_noncoding,2);
  N_silent(:,end+1,:) = N_silent(:,null_categ,:); % Copy total coverage from null+indel coverage
  N_nonsilent(:,end+1,:) = N_nonsilent(:,null_categ,:);
  N_noncoding(:,end+1,:) = N_noncoding(:,null_categ,:);
  
  % Find total across patients and save in G
  G.N_nonsilent = sum(N_nonsilent(:,end,:),3);
  G.N_silent = sum(N_silent(:,end,:),3);
  G.N_noncoding = sum(N_noncoding(:,end,:),3);
  G.n_nonsilent = sum(n_nonsilent(:,end,:),3);
  G.n_silent = sum(n_silent(:,end,:),3);
  G.n_noncoding = sum(n_noncoding(:,end,:),3);
  
  % Process covariates
  
  fprintf('Processing covariates...\n');
  
  V = nan(ng,nv);
  for vi=1:nv
      V(:,vi) = G.(cvnames{vi}); 
  end
  
  % Convert covariate raw values to Z-scores
  Z = nan(ng,nv);
  for vi=1:nv
    missing = isnan(V(:,vi)) | isinf(V(:,vi));
    mn = mean(V(~missing,vi));
    sd = std(V(~missing,vi));
    Z(~missing,vi) = (V(~missing,vi)-mn)./sd;
  end
  
  % Find "bagels" (neighboring genes for BMR estimation) 
  
  fprintf('Finding bagels...  ');
  
  max_neighbors = 50;
  qual_min = 0.05;
  
  G.nnei = nan(ng,1); 
  G.x = nan(ng,1); 
  G.X = nan(ng,1);
  
  for gi=1:ng
      if ~mod(gi,1000)
          fprintf('%d/%d ',gi,ng); 
      end
    
    % Calculate distances from this gene
    df2 = bsxfun(@minus,Z,Z(gi,:)).^2;
    dist2 = nansum(df2,2)./sum(~isnan(df2),2);
    [tmp,ord] = sort(dist2); 
    ord = [gi;ord(ord~=gi)];
    
    % Expand bagel outward until quality falls below qual_min
    nfit=0; 
    Nfit=0;
    for ni=0:max_neighbors
      gidx = ord(ni+1);
      ngene = G.n_silent(gidx) + G.n_noncoding(gidx);
      Ngene = G.N_silent(gidx) + G.N_noncoding(gidx);
      if ni==0
          ngene0=ngene; 
          Ngene0=Ngene; 
      end
      nfit=nfit+ngene; 
      Nfit=Nfit+Ngene;
      
      % Compare the gene being added to the central gene
      hc = hyge2cdf(ngene,Ngene,ngene0,Ngene0); % See hyge2cdf code in repository
      qual_left = min(hc, 1-hc);
      qual = 2*qual_left;
      
      % Stop if this gene would drop quality below qual_min
      if ni>0 && qual<qual_min
          break; 
      end
      
      % Update gene's statistics
      G.nnei(gi) = ni;
      G.x(gi) = nfit; 
      G.X(gi) = Nfit;
      
    end % Next neighborhood size
  end, fprintf('\n'); % Next gene
  
  fprintf('Expanding to (x,X)_gcp...\n');
  
  n_gcp = n_nonsilent + n_silent + n_noncoding;
  N_gcp = N_nonsilent + N_silent + N_noncoding;
  
  n_cp = sum(n_gcp,1);
  N_cp = sum(N_gcp,1);
  
  n_c = sum(n_cp,3);
  N_c = sum(N_cp,3);
  mu_c = n_c./N_c;
  
  n_tot = n_c(end);
  N_tot = N_c(end);
  mu_tot = n_tot/N_tot;
  f_c = mu_c/mu_tot;
  f_Nc = N_c/N_tot;
  
  n_p = n_cp(:,end,:);
  N_p = N_cp(:,end,:);
  mu_p = n_p./N_p;
  f_p = mu_p/mu_tot;
  f_Np = N_p/mean(N_p);
  
  x_gcp = repmat(G.x,[1 ncat+1 np]); 
  X_gcp = repmat(G.X,[1 ncat+1 np]); % Last column = total
  x_gcp = bsxfun(@times,x_gcp,f_c.*f_Nc); 
  X_gcp = bsxfun(@times,X_gcp,f_Nc);
  x_gcp = bsxfun(@times,x_gcp,f_p.*f_Np); 
  X_gcp = bsxfun(@times,X_gcp,f_Np);
  
  % 2D Projection to calculate gene p-values
  
  fprintf('Calculating p-value using 2D Projection method...  ');
  
  null_score_boost = 3;
  min_effect_size = 1.25;
  convolution_numbins = 1000;
  
  G.p = nan(ng,1);
  
  for g=1:ng, 
      if ~mod(g,1000)
          fprintf('%d/%d ',g,ng); 
      end
      
    % For each sample, prioritize mutation categories according to how likely
    % it would be for this gene x sample to have a mutation there by chance.
    
    N = reshape(N_nonsilent(g,1:ncat,:),ncat,np)';
    n = reshape(n_nonsilent(g,1:ncat,:),ncat,np)';
    x = reshape(x_gcp(g,1:ncat,:),ncat,np)';
    X = reshape(X_gcp(g,1:ncat,:),ncat,np)';
    P0 = hyge2pdf(0,N,x,X);
    P1 = hyge2pdf(1,N,x,X);
    
    % Determine each patient's priority order of categories (according to P1)
    % left column of "priority" = least extreme category of mutation
    % right column of "priority" = most extreme category of mutation
    [tmp priority] = sort(P1,2,'descend');
    % sort the P arrays to put the columns in least->most priority order
    shft = (priority - repmat(1:ncat,np,1));
    map = reshape(1:(np*ncat),np,ncat);
    newmap = map + shft*np;
    P0 = P0(newmap);
    P1 = P1(newmap);
    P2 = 1-(P0+P1);  % P2 refers to 2 or more mutation in category probability
    P2(P2<0) = 0;
    
    % For each sample, compute probability that it would have been of each (2-dimensional) degree.
    % degree=(d1,d2), where d=0 (no mut) ..... ncat (most extreme mut)
    % d1 is the most extreme mutation (or no mutation) while d2 is the second most extreme mutation (or no mutation)
    % d1 can be 0-ncat; d2 can be 0-d1
    
    Pdeg = zeros(np,ncat+1,ncat+1);
    for d1=0:ncat
        for d2=0:d1
      % Must have 0 in any/all categories > d1
      p = prod(P0(:,d1+1:end),2);
      if (d1>0)  % And (if d1>0)
        if (d1==d2)
          % If d1==d2, must have 2+ in category d1
          p = p .* P2(:,d1);
        else
          % Else: must have exactly 1 in category d1
          % must be clear in any/all categories (d2+1) to (d1-1)
          % and (if d2>0) have (1 or 2+) in category d2
          p = p .* P1(:,d1);
          p = p .* prod(P0(:,d2+1:d1-1),2);
          if (d2>0)
            p = p .* (P1(:,d2)+P2(:,d2));
          end
        end
      end
      Pdeg(:,d1+1,d2+1) = p;
        end
    end

    %% Calculate score for a sample being of each possible degree
    Sdeg = zeros(np,ncat+1,ncat+1);
    for d1=1:ncat, for d2=0:d1
      if d1==d2
        p = P2(:,d1);
      else
        if d2>0
          p = P1(:,d1).*P1(:,d2);
        else
          p = P1(:,d1);
        end
      end
      Sdeg(:,d1+1,d2+1) = -log10(p);
    end,end

    % Null score boost (recommended by authors of Mut_Sig)
    priority2 = [zeros(np,1) priority];
    Sdeg(priority2==null_categ) = Sdeg(priority2==null_categ) + null_score_boost;
    
    % Determine actual (two-dimensional) degree and score for each sample and sum scores to get score_obs for gene
    
    degree = zeros(np,2);
    score_obs = 0;
    for p = 1:np
      i = 1;
      for d = ncat:-1:1
        c = priority(p,d);
        if i==1
          if n(p,c)>=2
            degree(p,:) = [d d];
            i = 3;
          elseif n(p,c)==1
            degree(p,i) = d;
            i=i+1;
          end
        elseif i==2
          if n(p,c)>=1
            degree(p,i) = d;
            i=i+1;
          end
        else 
          break
        end
      end
      score_sample = Sdeg(p,degree(p,1)+1,degree(p,2)+1);
      score_obs = score_obs + score_sample;
    end

    % Minimum effect size
    score_obs = score_obs / min_effect_size;
    
    % For zero score, don't do convolutions
    if score_obs<=0
        G.p(g)=1; 
        continue; 
    end
    
    % Compute p value for gene by convolutions
    
    numbins = convolution_numbins;
    binsize = score_obs / numbins;
    H = zeros(numbins,1);
    H(1) = 1;  % Initial condition: all probability is in first bin
    
    % Sequential convolution
    offset = min(numbins, round(Sdeg/binsize));
    ncols = (ncat+1)*(ncat+2)/2;
    newH = zeros(numbins,ncols);
    for p=1:np
      newH(:) = 0;
      col=1;
      for d1=0:ncat
          for d2=0:d1
            o = offset(p,d1+1,d2+1);
            newH(o+1:end,col) = Pdeg(p,d1+1,d2+1) .* H(1:end-o);
            col=col+1;
          end
      end
      H = sum(newH,2);
    end
    
    % Save p-value
    G.p(g) = max(0,1-sum(H));
    
  end
  fprintf('\n'); % Repeat for next gene
  
  % FDR
  G.q = calc_fdr_value(G.p); % See calc_fdr_value code in repository
  
  G = sort_struct(G,'p');
  save_struct(G,output_file);
  
  fprintf('Complete. Results are written to %s\n',output_file);
  
end