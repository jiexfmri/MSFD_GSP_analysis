function [spin_mean_task7,spin_task7] = null_spin_key_task7(W,XS_task,n_modes,n_nulls,parcel)

n_subjs = size(XS_task,3);
n_ROI = size(W,1);

% spin test
permID = load(strcat('Data/',parcel,'_PermID.mat'));
permID = permID.permID;

%% spin null model 
spin_task7 = zeros(n_nulls,n_modes,7);
spin_mean_task7 = zeros(n_nulls,n_modes);

parfor null = 1:n_nulls
    % Recon 7 task maps
    acc_task7 = zeros(n_modes,7,n_subjs);
    for s = 1:n_subjs
        X_task = XS_task(permID(:,null),:,s);   %spin null
               
        acc_task7(:,:,s) = Recon_task_activation(W,X_task,n_modes);
    end
    group_acc_task7 = mean(acc_task7,3);  
    spin_task7(null,:,:) = group_acc_task7;
    spin_mean_task7(null,:) = mean(group_acc_task7,2);
end


end






