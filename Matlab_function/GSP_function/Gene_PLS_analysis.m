function [result_corr, result_VAR] = Gene_PLS_analysis(parcel,type,MSFD)

n_perm = 10000;
Tmap = MSFD;
n_ROI=size(Tmap,1);

% load expression data 
expression_gene = importdata(strcat('Data',filesep, parcel,'_gene_expression.xlsx'));
expre_data =  expression_gene.data;
expre_data = expre_data(1:n_ROI,2:end);
expre_name =  expression_gene.textdata(1,2:end)';

% load T map
Tmap_data =  Tmap(1:n_ROI,:);

%%% using the Moran spectral randomization to test the the statistical significance 
permID = load(strcat('Data',filesep, parcel,'_PermID.mat'));
permID = permID.permID;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLS analysis
X=expre_data;    % Predictors
Y=Tmap_data;     % Response variable
% z-score:
X=zscore(X);
Y=zscore(Y);

%perform full PLS and plot variance in Y explained
dim = 5;
[XL,YL,XS,YS,BETA,PCTVAR_real,MSE,stats]=plsregress(X,Y,dim);

% Loop through each PLS component and its corresponding variance explained
% Use fprintf to print each PLS component on a new line
for i = 1:dim
    fprintf('Variance Explained of PLS%d = %.4f\n', i, PCTVAR_real(2, i));
end

% %% Permutation testing based on spherical rotations, to account for spatial autocorrelation, of the T map (1000 times)
[rho,p]=corr(XS(:,1),Tmap_data,'type','Spearman');
p_spin = spin_perm_sphere_p(Tmap_data, XS(:,1), permID, 'spearman');
result_corr = [rho,p_spin];    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:
dim = 5;
out_Rsq = zeros(n_perm, dim);
Cumulative_Rsq = zeros(1, dim);
Cumulative_p = zeros(1, dim);
out_Expl = zeros(n_perm, dim);
Explain = zeros(1, dim);
Expl_p = zeros(1, dim);

parfor dim=1:dim
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    temp=cumsum(100*PCTVAR(2,1:dim));   % Cumulative interpretation rate
    Rsquared = temp(dim);    %real value
    
    Expl_vari = 100*PCTVAR(2,1:dim);
    Expl = Expl_vari(dim);
    %%A null model was produced below to prove that the reported transcriptome-imaging association was not a false positive
    Rsq = zeros(1, n_perm); 
    Expl_surr = zeros(1, n_perm); 
    for j=1:n_perm
         j
         Yp = Tmap_data(permID(:,j));
         [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);
         temp=cumsum(100*PCTVAR(2,1:dim));
         Rsq(j) = temp(dim);
         
         Expl_vari_surr = 100*PCTVAR(2,1:dim);
         Expl_surr(j) = Expl_vari_surr(dim);
    end
    dim
    out_Rsq(:,dim) = Rsq';
    Cumulative_Rsq(dim)=Rsquared;
    Cumulative_p(dim)=length(find(Rsq>=Rsquared))/n_perm;     % evaluates the P value of the true value compared to the null model
    
    out_Expl(:,dim) = Expl_surr';
    Explain(dim) = Expl;
    Expl_p(dim)=length(find(Expl_surr>=Expl))/n_perm;
end
result_VAR = [Explain;Expl_p];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output gene Rank results will be used for enrichment analysis
genes=expre_name; % this needs to be imported first
geneindex=1:size(expre_data,2);

% Bootstrap to get the gene list:
bootnum=n_perm;
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

pls_index1 = 1;
[R1,p1]=corr(XS(:,pls_index1),YS(:,pls_index1));

%align PLS components with desired direction for interpretability 
if R1(1,1)<0  
    stats.W(:,pls_index1)=-1*stats.W(:,pls_index1);
    XS(:,pls_index1)=-1*XS(:,pls_index1);
end

% choose the interest of component
[PLS1w,x1] = sort(stats.W(:,pls_index1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
%start bootstrap,
% The goal is to find the standard deviation
parfor i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; 
    Xr=X(myresample,:); 
    Yr=Y(myresample,:);
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data

    temp=stats.W(:,pls_index1);    
    newW=temp(x1);           
    if corr(PLS1w,newW)<0   
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');

%get bootstrap weights Z socre, the normalized weights of PLS
temp1=PLS1w./PLS1sw';

%order bootstrap weights (Z) and names of regions
[Z1 ,ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);

%%% Find P-value based on Z-score
P1_value_Pos = spm_z2p(Z1(Z1>=0));
P1_value_Neg = spm_z2p(abs(Z1(Z1<0)));
P1_value = [P1_value_Pos;P1_value_Neg];
correct_P1 = mafdr(P1_value,'BHFDR',1);  

P1_value_id = find(correct_P1<0.05);
P1_value_05 = Z1(P1_value_id);

%print out results
% later use first column of these csv files for pasting into GOrilla (for
% bootstrapped ordered list of genes) 
outdir = 'Results';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
filename = strcat(outdir, filesep,'Gene_PLS', ...
                  num2str(pls_index1), ...
                  '_geneWeights_', ...
                  type, '.csv');
fid1 = fopen(filename, 'w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f, %f\n', PLS1{i}, geneindex1(i), Z1(i),correct_P1(i));
end
fclose(fid1)

end

