function [rewired_mean_task47,rewired_task47]= null_rewired_task(W,XS_task,n_modes,n_nulls,task_num)

n_subjs = size(XS_task,3);
rewired_task47 = zeros(n_nulls,n_modes,47);
rewired_mean_task47 = zeros(n_nulls,n_modes);

parfor null = 1:n_nulls
    [W_null, ~] = null_model_und_sign(W);   
    
    if task_num == 47
        %% Recon 47 task maps
        acc_task47 = zeros(n_modes,47,n_subjs);
        for s = 1:n_subjs
            X_task = XS_task(:,:,s);
            acc_task47(:,:,s) = Recon_task_activation(W_null,X_task,n_modes);
        end
        group_acc_task47 = mean(acc_task47,3);  
        rewired_task47(null,:,:) = group_acc_task47;
        rewired_mean_task47(null,:,:)=mean(group_acc_task47,2);
    end
end

end






