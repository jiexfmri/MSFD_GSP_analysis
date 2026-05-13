clc
clear all
% add paths to necessary scripts
addpath(genpath('Matlab_function'))

parcel = 'HCP_MMP';
N_ROI = 360;  % number of regions

% Normalisation: We peformed rank normalisation and rescaling to enforce the 
% same range of values for each modality, accounting for the variations in sparsity. 
% This helps to equalise the contribution of each modality, without enforcing an 
% arbitary threshold.      
ts_group = xlsread('Data/Connectome/HCP_MMP_group_TS.xlsx');
gd_group = xlsread('Data/Connectome/HCP_MMP_group_GD.xlsx');
mpc_group = xlsread('Data/Connectome/HCP_MMP_group_MPC.xlsx');

% tract strength (TS)  Normalisation
[~, idx]  = sort(ts_group(:), 'ascend');   % larger numbers are higher rank
ts_norm   = sort_back((1:length(ts_group(:)))', idx);
ts_norm   = reshape(ts_norm, [size(ts_group)]);
ts_norm   = (ts_norm - (numel(ts_group) - nnz(ts_group))) .* (ts_group>0);
unique(ts_norm)
    
% geodesic distance (GD)  Normalisation
this_gd     = 1./gd_group;  % invert distance matrix
[~, idx]    = sort(this_gd(:), 'ascend');
this_gd     = sort_back((1:length(this_gd(:)))', idx);
this_gd(this_gd==0) = nan;
gd_scale   = rescale(this_gd(:), 1, max(ts_norm(:)));
gd_norm    = reshape(gd_scale, [size(gd_group)]);
gd_norm(isnan(gd_norm)) = 0;
gd_norm = gd_norm-diag(diag(gd_norm));
    
% microstructure profile covariance (MPC)  Normalisation
[~, idx] = sort(mpc_group(:), 'ascend');  % larger numbers are higher rank
this_mpc = sort_back((1:length(mpc_group(:)))', idx);
this_mpc = rescale(this_mpc(:), 1, max(ts_norm(:)));
this_mpc(this_mpc==0) = nan;
mpc_scale  = rescale(this_mpc(:), 1, max(ts_norm(:)));
mpc_norm   = reshape(mpc_scale, [size(mpc_group)]);
mpc_norm(isnan(mpc_norm)) = 0;
mpc_norm = mpc_norm - diag(diag(mpc_norm));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fusion: Horizontal concatenation of matrices and production of a node-to-node affinity matrix using row-wise normalised angle similarity. 

% Store modalities in a structure for easy access
Connectome.SC = ts_group;
Connectome.GD = gd_norm;
Connectome.MPC = mpc_norm;

% Define combinations
combinations = {
    {'SC', 'MPC'}, ...
    {'SC', 'GD'}, ...
    {'MPC', 'GD'}, ...
    {'SC', 'MPC', 'GD'}
};

% Initialize results structure
results = struct();

% Loop through each combination
for i = 1:length(combinations)
    comb = combinations{i};
    mat_horz = [];
    
    % Concatenate matrices based on current combination
    for j = 1:length(comb)
        mat_horz = [mat_horz, Connectome.(comb{j})];
    end
    
    % Compute affinity matrix  
    affinity_matrix     = 1-squareform(pdist(mat_horz'.','cosine'));
    affinity_matrix(isnan(affinity_matrix)) = 0;
    norm_angle_matrix    = 1-acos(affinity_matrix)/pi;
    final_matrix = norm_angle_matrix -diag(diag(norm_angle_matrix));
    
    % Store result in the structure with dynamic field name
    field_name = [strjoin(comb, '_')];
    results.(field_name) = final_matrix;
end

results.SC = ts_group;
results.GD = gd_norm;
results.MPC = mpc_norm;

% Save all results
save('Data/Connectome/HCP_MMP_Fusion_all.mat', '-struct', 'results');