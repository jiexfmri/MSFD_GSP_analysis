clc
clear all
% add paths to necessary scripts
addpath(genpath('Matlab_function'))

parcel = 'HCP_MMP';
n_ROI = 360;
n_modes = 200;
n_nulls = 1000;

% subject-level task maps
%  the seven key HCP task contrasts
X_task7 = load(strcat('Data',filesep, parcel,'_key_seven_task_zstat_example.mat'));
X_task7 = X_task7.zstat;

X_task47 = load(strcat('Data',filesep, parcel,'_alltask_zstat_example.mat'));
X_task47 = X_task47.zstat;

% connectome
Net = load(strcat('Data',filesep, parcel,'_Fusion_all.mat'));

% networks = {'SC', 'MPC', 'GD','SC_MPC', 'SC_GD','MPC_GD','SC_MPC_GD'};
networks = {'SC_MPC_GD'};

ACC_key_task7_subj = struct();
AUC_key_task7_subj = struct();
ACC_task47_subj = struct();
AUC_task47_subj = struct();
for i = 1:length(networks)
    net_name = networks{i}; 
    net_param = Net.(net_name);     
    [ACC_key_task7_subj.(net_name), AUC_key_task7_subj.(net_name)] = Recon_key_task7_activation_subj(net_param,X_task7,n_modes,n_nulls,parcel,net_name);
    [ACC_task47_subj.(net_name), AUC_task47_subj.(net_name)] = Recon_task47_activation_subj(net_param, X_task47, n_modes, n_nulls,parcel, net_name);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Summary AUC subject-level
acc_key = zeros(n_modes, length(networks)+1); 
auc_key = zeros(n_modes, length(networks)+1);
acc_key(:,1) = (1:n_modes)';
auc_key(:,1) = (1:n_modes)';

acc_all = zeros(n_modes, length(networks)+1); 
auc_all = zeros(n_modes, length(networks)+1);
acc_all(:,1) = (1:n_modes)';
auc_all(:,1) = (1:n_modes)';
for j = 1:length(networks)
    acc_key(:, j+1) = ACC_key_task7_subj.(networks{j})(:, 2);
    auc_key(:, j+1) = AUC_key_task7_subj.(networks{j})(:, 2);
    outputFile = strcat('Results/recon_key_task7/',parcel,'_Recon_mean_key_task7_subj_Summary.xlsx');
    xlswrite(outputFile, acc_key, strcat('A2:H', num2str(n_modes+1)));
    xlswrite(outputFile, auc_key, strcat('L2:S', num2str(n_modes+1)));
    
    acc_all(:, j+1) = ACC_task47_subj.(networks{j})(:, 2);
    auc_all(:, j+1) = AUC_task47_subj.(networks{j})(:, 2);
    outputFile = strcat('Results/recon_task47/',parcel,'_Recon_mean_task47_subj_Summary.xlsx');
    xlswrite(outputFile, acc_all, strcat('A2:H', num2str(n_modes+1)));
    xlswrite(outputFile, auc_all, strcat('L2:S', num2str(n_modes+1)));
end

