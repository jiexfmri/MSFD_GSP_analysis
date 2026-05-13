function [rewired_mean_task7,rewired_task7]= null_TS_rewired_MPC_GD_key_task7(XS_task,n_modes,n_nulls)

n_subjs = size(XS_task,3);
rewired_task7 = zeros(n_nulls,n_modes,7);
rewired_mean_task7 = zeros(n_nulls,n_modes);

% connectome
Net = load(strcat('Data/HCP_MMP_Fusion_all.mat'));
SC = Net.SC;
GD = Net.GD;
MPC = Net.MPC;

parfor null = 1:n_nulls
    [W_null_GD, ~] = null_model_und_sign(GD);  
    [W_null_MPC, ~] = null_model_und_sign(MPC);  
    % Concatenate matrices based on current combination
    mat_horz = [SC, W_null_GD, W_null_MPC];
     % Compute affinity matrix  
    affinity_matrix     = 1-squareform(pdist(mat_horz'.','cosine'));
    affinity_matrix(isnan(affinity_matrix)) = 0;
    norm_angle_matrix    = 1-acos(affinity_matrix)/pi;
    final_matrix = norm_angle_matrix -diag(diag(norm_angle_matrix));
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Recon key 7 task maps
    acc_task7 = zeros(n_modes,7,n_subjs);
    for s = 1:n_subjs
        X_task = XS_task(:,:,s);
        acc_task7(:,:,s) = Recon_task_activation(final_matrix,X_task,n_modes);
    end
    group_acc_task7 = mean(acc_task7,3);  
    rewired_task7(null,:,:) = group_acc_task7;
    rewired_mean_task7(null,:,:) = mean(group_acc_task7,2);
end

end
