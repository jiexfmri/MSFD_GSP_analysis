function [spin_mean_task47,spin_task47] = null_spin_task(W,XS_task,n_modes,n_nulls,parcel,task_num)

n_subjs = size(XS_task,3);
n_ROI = size(W,1);

% spin test
permID = load(strcat('Data/',parcel,'_PermID.mat'));
permID = permID.permID;

%% spin null model 
spin_task47 = zeros(n_nulls,n_modes,47);
spin_mean_task47 = zeros(n_nulls,n_modes);

parfor null = 1:n_nulls
        if task_num == 47
        %% Recon 47 task maps
        acc_task47 = zeros(n_modes,47,n_subjs);
        for s = 1:n_subjs
            X_task = XS_task(permID(:,null),:,s);   %spin null

            acc_task47(:,:,s) = Recon_task_activation(W,X_task,n_modes);
        end
        group_acc_task47 = mean(acc_task47,3);  
        spin_task47(null,:,:) = group_acc_task47;
        spin_mean_task47(null,:) = mean(group_acc_task47,2);
    end
end
end






