function [rewired_mean_task47,rewired_task47]= null_TS_rewired_GD_task(XS_task,n_modes,n_nulls,task_num)

% XS_task = mean(XS_task,3);
n_subjs = size(XS_task,3);
rewired_task47 = zeros(n_nulls,n_modes,47);
rewired_mean_task47 = zeros(n_nulls,n_modes);

% connectome
Net = load(strcat('Data/HCP_MMP_Fusion_all.mat'));
SC = Net.SC;
GD = Net.GD;
MPC = Net.MPC;

parfor null = 1:n_nulls
    [W_null, ~] = null_model_und_sign(GD);   
    % Concatenate matrices based on current combination
    mat_horz = [SC, W_null];
    % Compute affinity matrix  
    affinity_matrix     = 1-squareform(pdist(mat_horz'.','cosine'));
    affinity_matrix(isnan(affinity_matrix)) = 0;
    norm_angle_matrix    = 1-acos(affinity_matrix)/pi;
    final_matrix = norm_angle_matrix -diag(diag(norm_angle_matrix)); 
    
    if task_num == 47
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Recon 47 task maps
        acc_task47 = zeros(n_modes,47,n_subjs);
        for s = 1:n_subjs
            X_task = XS_task(:,:,s);
            acc_task47(:,:,s) = Recon_task_activation(final_matrix,X_task,n_modes);
        end
        group_acc_task47 = mean(acc_task47,3);  
        rewired_task47(null,:,:) = group_acc_task47;
        rewired_mean_task47(null,:) = mean(group_acc_task47,2);
    end
end

end






