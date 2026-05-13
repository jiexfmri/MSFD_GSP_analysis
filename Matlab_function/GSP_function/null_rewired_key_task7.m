function [rewired_mean_task7,rewired_task7]= null_rewired_key_task7(W,XS_task,n_modes,n_nulls)

n_subjs = size(XS_task,3);
rewired_task7 = zeros(n_nulls,n_modes,7);
rewired_mean_task7 = zeros(n_nulls,n_modes);

parfor null = 1:n_nulls
    [W_null, ~] = null_model_und_sign(W);   
    %% Recon key 7 task maps
    acc_task7 = zeros(n_modes,7,n_subjs);
    for s = 1:n_subjs
        X_task = XS_task(:,:,s);
        acc_task7(:,:,s) = Recon_task_activation(W_null,X_task,n_modes);
    end
    group_acc_task7 = mean(acc_task7,3);  
    rewired_task7(null,:,:) = group_acc_task7;
    rewired_mean_task7(null,:) = mean(group_acc_task7,2);
end

end
