clc
clear all

% add paths to necessary scripts
addpath(genpath('Matlab_function'))

nSurr = 19;
parcel = 'HCP_MMP';
type='SC_MPC_GD';

% connectome
Net = load(strcat('Data',filesep, parcel,'_Fusion_all.mat'));
net_param = Net.(type); 

% eigen decomposition of MSC
[U,LambdaL] = GSP_Laplacian(net_param);

% project fMRI signals into eigenmodes to find cutoff and split high and low frequencies
zX_RS = get_HCP_BOLD_fMRI(parcel);
[NN,Vlow, Vhigh]= Get_cut_off_freq(U,zX_RS);

% analyze graph signal¡ªfilted signal
[N_c,N_d,mean_c,mean_d] = GSanalysis(zX_RS,Vhigh,Vlow,U);

% generate Connectome-informed surrogates
[XrandS] = GSrandomozation_surrogates(zX_RS,U,nSurr);

% generate MSFD index
[MSFD_log,MSFD_thr] = Compute_MSFD(XrandS,U,Vlow,Vhigh,mean_d,mean_c,N_d,N_c);

% Relation to Macro and micro hierarchy organization
[hierachical] = Relationship_Hierarchy(MSFD_log,parcel);

% gene analyis
[Pls_gene, Pls_VAR_gene] = Gene_PLS_analysis(parcel,type,MSFD_log);