function [line_ACC_task47,line_AUC_task47]= Recon_task47_activation_subj(W,XS_task,n_modes,n_nulls,parcel,type)

file_name = strcat('Results/recon_task47/',type);
if ~exist(file_name,'dir')
    mkdir(file_name);
end

%%%% Recon 47 task
n_subjs = size(XS_task,3);
acc_task47 = zeros(n_modes,47,n_subjs);
parfor s = 1:n_subjs
    X_task = XS_task(:,:,s);
    acc_task47(:,:,s) = Recon_task_activation(W,X_task,n_modes);
end

n_subjs = size(XS_task,3);
acc_task7 = zeros(n_modes,7,n_subjs);
% average subj
save(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj','.mat'),'-v6','acc_task47');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% average subject
group_acc_task47 = mean(acc_task47,3);  
ACC_task47  = [(1:n_modes)',group_acc_task47];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_mean','.xlsx'),ACC_task47, strcat('A2:AV',num2str(n_modes+1)));

for mode= 1:n_modes
    for s = 1: 47
        task47_AUC(mode,s)=trapz(ACC_task47(1:mode,s+1))/n_modes;
    end
end
AUC_task47  = [(1:n_modes)',task47_AUC];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_mean','.xlsx'),AUC_task47,  strcat('AY2:CT',num2str(n_modes+1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average task47
acc_mean_task47 = squeeze(mean(acc_task47,2));  
mean_task47  = [(1:n_modes)',acc_mean_task47];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj','.xlsx'),mean_task47);

line_task47 = mean(acc_mean_task47,2);
line_ACC_task47  = [(1:n_modes)',line_task47];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_mean_all','.xlsx'),line_ACC_task47, strcat('A2:B',num2str(n_modes+1)));

for mode= 1:n_modes
     line_task47_AUC(mode,1)=trapz(line_task47(1:mode,1))/n_modes;
end
line_AUC_task47  = [(1:n_modes)',line_task47_AUC];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_mean_all','.xlsx'),line_AUC_task47, strcat('E2:F',num2str(n_modes+1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NULL mode

null_model = 1;
if null_model == 1

    disp('Running null model...')

    %% fusing the TS matrix with a random noise matrix of similar properties to MPC 
    [TS_rewired_mean_MPC_task47,TS_rewired_MPC_task47] = null_TS_rewired_MPC_task(XS_task,n_modes,n_nulls,47);
    % save data for R visual TASK 47
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:48) = reshape(TS_rewired_MPC_task47(:,i,:),n_nulls,47); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(TS_rewired_mean_MPC_task47(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_TS_rewired_MPC','.xlsx'),data,strcat('A2:AV',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj_null_TS_rewired_MPC','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    % Calculating p-value
    P_rewired47 = zeros(n_modes,47);
    P_rewired = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:47
            P_rewired47(i,j) = sum(TS_rewired_MPC_task47(:,i,j)>group_acc_task47(i,j))/n_nulls;  
        end
        P_rewired(i,1) = sum(TS_rewired_mean_MPC_task47(:,i)>acc_mean_task47(i,1))/n_nulls;
     end
    output_rewired47  = [(1:n_modes)',P_rewired47];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_TS_rewired_MPC_Pvalue','.xlsx'),output_rewired47, strcat('A2:AV',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47subj_null_TS_rewired_MPC_Pvalue','.xlsx'),output_rewired, strcat('A2:B',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fusing the TS matrix with a random noise matrix of similar properties to GD 
    [TS_rewired_mean_GD_task47,TS_rewired_GD_task47] = null_TS_rewired_GD_task(XS_task,n_modes,n_nulls,47);
    % save data for R visual TASK 47
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:48) = reshape(TS_rewired_GD_task47(:,i,:),n_nulls,47); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(TS_rewired_mean_GD_task47(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_TS_rewired_GD','.xlsx'),data,strcat('A2:AV',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj_null_TS_rewired_GD','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    % Calculating p-value
    P_rewired47 = zeros(n_modes,47);
    P_rewired = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:47
            P_rewired47(i,j) = sum(TS_rewired_GD_task47(:,i,j)>group_acc_task47(i,j))/n_nulls;  
        end
        P_rewired(i,1) = sum(TS_rewired_mean_GD_task47(:,i)>acc_mean_task47(i,1))/n_nulls;
     end
    output_rewired47  = [(1:n_modes)',P_rewired47];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_TS_rewired_GD_Pvalue','.xlsx'),output_rewired47, strcat('A2:AV',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47subj_null_TS_rewired_GD_Pvalue','.xlsx'),output_rewired, strcat('A2:B',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fusing the TS matrix with a random noise matrix of similar properties to MPC and GD
    [TS_rewired_mean_MPC_GD_task47,TS_rewired_MPC_GD_task47] = null_TS_rewired_MPC_GD_task(XS_task,n_modes,n_nulls,47);
    % save data for R visual TASK 47
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:48) = reshape(TS_rewired_MPC_GD_task47(:,i,:),n_nulls,47); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(TS_rewired_mean_MPC_GD_task47(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_TS_rewired_MPC_GD','.xlsx'),data,strcat('A2:AV',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj_null_TS_rewired_MPC_GD','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    % Calculating p-value
    P_rewired47 = zeros(n_modes,47);
    P_rewired = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:47
            P_rewired47(i,j) = sum(TS_rewired_MPC_GD_task47(:,i,j)>group_acc_task47(i,j))/n_nulls;  
        end
        P_rewired(i,1) = sum(TS_rewired_mean_MPC_GD_task47(:,i)>acc_mean_task47(i,1))/n_nulls;
     end
    output_rewired47  = [(1:n_modes)',P_rewired47];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_TS_rewired_MPC_GD_Pvalue','.xlsx'),output_rewired47, strcat('A2:AV',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47subj_null_TS_rewired_MPC_GD_Pvalue','.xlsx'),output_rewired, strcat('A2:B',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% rewired null model
    [rewired_mean_task47,rewired_task47] = null_rewired_task(W,XS_task,n_modes,n_nulls,47)
    % save data for R visual TASK 47
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:48) = reshape(rewired_task47(:,i,:),n_nulls,47); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(rewired_mean_task47(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_rewired','.xlsx'),data,strcat('A2:AV',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj_null_rewired','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spin null model
    [spin_mean_task47,spin_task47] = null_spin_task(W,XS_task,n_modes,n_nulls,parcel,47)
    % save data for R visual TASK 47
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:48) = reshape(spin_task47(:,i,:),n_nulls,47); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(spin_mean_task47(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_spin','.xlsx'),data,strcat('A2:AV',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj_null_spin','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculating p-value
    P_rewired47 = zeros(n_modes,47);P_spin47 = zeros(n_modes,47);
    P_rewired = zeros(n_modes,1);P_spin = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:47
            P_rewired47(i,j) = sum(rewired_task47(:,i,j)>group_acc_task47(i,j))/n_nulls;  
            P_spin47(i,j) = sum(spin_task47(:,i)>group_acc_task47(i,1))/n_nulls;
        end
        P_rewired(i,1) = sum(rewired_mean_task47(:,i)>acc_mean_task47(i,1))/n_nulls;  
        P_spin(i,1) = sum(spin_mean_task47(:,i)>acc_mean_task47(i,1))/n_nulls;
     end
    
    %%%%%%%%%%%%%%%%%%%%%%
    output_rewired47  = [(1:n_modes)',P_rewired47];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_rewired_Pvalue','.xlsx'),output_rewired47, strcat('A2:AV',num2str(n_modes+1)));
    output_spin47  = [(1:n_modes)',P_spin47];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task47_subj_null_spin_Pvalue','.xlsx'),output_spin47, strcat('A2:AV',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj_null_rewired_Pvalue','.xlsx'),output_rewired, strcat('A2:AV',num2str(n_modes+1)));
    output_spin  = [(1:n_modes)',P_spin];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task47_subj_null_spin_Pvalue','.xlsx'),output_spin, strcat('A2:B',num2str(n_modes+1)));
end

end






