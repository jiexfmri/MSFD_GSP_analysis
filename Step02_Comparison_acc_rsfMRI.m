clc
clear all
% add paths to necessary scripts
addpath(genpath('Matlab_function'))

parcel = 'HCP_MMP';
n_ROI = 360;
n_modes = 200;
n_nulls = 1000;

% rest-fMRI
X_RS = load(strcat('Data',filesep, parcel,'_RSfMRI_example.mat'));
X_RS = X_RS.X_RS;

% connectome
Net = load(strcat('Data',filesep, parcel,'_Fusion_all.mat'));

networks = {'SC', 'MPC', 'GD','SC_MPC', 'SC_GD','MPC_GD','SC_MPC_GD'};
%networks = {'SC_MPC_GD'};

ACC_subj = struct(); AUC_subj = struct();
for i = 1:length(networks)
    net_name = networks{i}; 
    net_param = Net.(net_name);     
    [ACC_subj.(net_name), AUC_subj.(net_name)] = Recon_activity_FC_subj(net_param, X_RS, n_modes, n_nulls,parcel, net_name);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary AUC subject-level
task = {'activity_ratio','FC_Corr'};
for i = 1:length(task)
    acc = [(1:n_modes)'];
    auc = [(1:n_modes)']; 
    for j = 1:length(networks)
        acc = [acc, ACC_subj.(networks{j})(:, i+1)];
        auc = [auc, AUC_subj.(networks{j})(:, i+1)];
    end  
    outputFile = strcat('Results/recon_activity_FC/', parcel, '_Recon_rsfMRI_', task{i}, '_subj_Summary_seven.xlsx');
    xlswrite(outputFile, acc, strcat('A2:H', num2str(n_modes+1)));
    xlswrite(outputFile, auc, strcat('L2:S', num2str(n_modes+1)));
end
