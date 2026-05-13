function [line_ACC_task7,line_AUC_task7]= Recon_key_task7_activation_subj(W,XS_task,n_modes,n_nulls,parcel,type)

% save file path 
file_name = strcat('Results/recon_key_task7/',type);
if~exist(file_name,'dir')
    mkdir(file_name);
end

n_ROI = size(W,1);
n_subjs = size(XS_task,3);
RX_task7 = zeros(n_ROI,7,n_subjs);
acc_task7 = zeros(n_modes,7,n_subjs);
recon_task7 = zeros(n_ROI,n_modes,7,n_subjs);
parfor s = 1:n_subjs
    X_task = XS_task(:,:,s);
    % Working memory  % Gambling   % Motor % Language % Social  % Relational  % Emotion
    seven_task = X_task;
    %%%% Recon seven task
    [acc_task7(:,:,s),recon_task7(:,:,:,s)] = Recon_task_activation(W,seven_task,n_modes);
end
save(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task7_subj','.mat'),'-v6','acc_task7');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_acc_task7 = mean(acc_task7,3);  
ACC_task7  = [(1:n_modes)',group_acc_task7];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task7_subj_mean','.xlsx'),ACC_task7, strcat('A2:H',num2str(n_modes+1)));

for mode= 1:n_modes
    for s = 1: 7
        task7_AUC(mode,s)=trapz(ACC_task7(1:mode,s+1))/n_modes;
    end
end
AUC_task7  = [(1:n_modes)',task7_AUC];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task7_subj_mean','.xlsx'),AUC_task7, strcat('K2:R',num2str(n_modes+1)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% average task7
acc_mean_task7 = squeeze(mean(acc_task7,2));  
mean_task7  = [(1:n_modes)',acc_mean_task7];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_task7_subj','.xlsx'),mean_task7);

line_task7 = mean(acc_mean_task7,2);
line_ACC_task7  = [(1:n_modes)',line_task7];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task7_mean_all','.xlsx'),line_ACC_task7, strcat('A2:B',num2str(n_modes+1)));

for mode= 1:n_modes
     line_task7_AUC(mode,1)=trapz(line_task7(1:mode,1))/n_modes;
end
line_AUC_task7  = [(1:n_modes)',line_task7_AUC];
xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_task7_mean_all','.xlsx'),line_AUC_task7, strcat('E2:F',num2str(n_modes+1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
null_model = 0;
if null_model == 1

    disp('Running null model...')

    %% NULL mode
    %% fusing the TS matrix with a random noise matrix of similar properties to MPC 
    [TS_rewired_MPC_mean_task7,TS_rewired_MPC_task7] = null_TS_rewired_MPC_key_task7(XS_task,n_modes,n_nulls);
    % save data for R visual  TASK 7
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:8) = reshape(TS_rewired_MPC_task7(:,i,:),n_nulls,7); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(TS_rewired_MPC_mean_task7(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_TS_rewired_MPC','.xlsx'),data,strcat('A2:H',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7_subj_null_TS_rewired_MPC','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    % Calculating p-value
    P_rewired7 = zeros(n_modes,7);
    P_rewired = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:7
            P_rewired7(i,j) = sum(TS_rewired_MPC_task7(:,i,j)>group_acc_task7(i,j))/n_nulls;  
        end
        P_rewired(i,1) = sum(TS_rewired_MPC_mean_task7(:,i)>acc_mean_task7(i,1))/n_nulls; 
     end
    output_rewired7  = [(1:n_modes)',P_rewired7];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_TS_rewired_MPC_Pvalue','.xlsx'),output_rewired7, strcat('A2:H',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7subj_null_TS_rewired_MPC_Pvalue','.xlsx'),output_rewired, strcat('A2:B',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fusing the TS matrix with a random noise matrix of similar properties to GD 
    [TS_rewired_GD_mean_task7,TS_rewired_GD_task7] = null_TS_rewired_GD_key_task7(XS_task,n_modes,n_nulls);
    % save data for R visual  TASK 7
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:8) = reshape(TS_rewired_GD_task7(:,i,:),n_nulls,7); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(TS_rewired_GD_mean_task7(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_TS_rewired_GD','.xlsx'),data,strcat('A2:H',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7_subj_null_TS_rewired_GD','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    % Calculating p-value
    P_rewired7 = zeros(n_modes,7);
    P_rewired = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:7
            P_rewired7(i,j) = sum(TS_rewired_GD_task7(:,i,j)>group_acc_task7(i,j))/n_nulls;  
        end
        P_rewired(i,1) = sum(TS_rewired_GD_mean_task7(:,i)>acc_mean_task7(i,1))/n_nulls;
     end
    output_rewired7  = [(1:n_modes)',P_rewired7];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_TS_rewired_GD_Pvalue','.xlsx'),output_rewired7, strcat('A2:H',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7subj_null_TS_rewired_GD_Pvalue','.xlsx'),output_rewired, strcat('A2:B',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fusing the TS matrix with a random noise matrix of similar properties to MPC and GD 
    [TS_rewired_MPC_GD_mean_task7,TS_rewired_MPC_GD_task7] = null_TS_rewired_MPC_GD_key_task7(XS_task,n_modes,n_nulls);
    % save data for R visual  TASK 7
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:8) = reshape(TS_rewired_MPC_GD_task7(:,i,:),n_nulls,7); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(TS_rewired_MPC_GD_mean_task7(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_TS_rewired_MPC_GD','.xlsx'),data,strcat('A2:H',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7_subj_null_TS_rewired_MPC_GD','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    % Calculating p-value
    P_rewired7 = zeros(n_modes,7);
    P_rewired = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:7
            P_rewired7(i,j) = sum(TS_rewired_MPC_GD_task7(:,i,j)>group_acc_task7(i,j))/n_nulls;  
        end
        P_rewired(i,1) = sum(TS_rewired_MPC_GD_mean_task7(:,i)>acc_mean_task7(i,1))/n_nulls;
     end
    output_rewired7  = [(1:n_modes)',P_rewired7];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_TS_rewired_MPC_GD_Pvalue','.xlsx'),output_rewired7, strcat('A2:H',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7subj_null_TS_rewired_MPC_GD_Pvalue','.xlsx'),output_rewired, strcat('A2:B',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% rewired null model
    [rewired_mean_task7,rewired_task7] = null_rewired_key_task7(W,XS_task,n_modes,n_nulls);
    % save data for R visual  TASK 7
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:8) = reshape(rewired_task7(:,i,:),n_nulls,7); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(rewired_mean_task7(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_rewired','.xlsx'),data,strcat('A2:H',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7_subj_null_rewired','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spin null model
    [spin_mean_task7,spin_task7]  = null_spin_key_task7(W,XS_task,n_modes,n_nulls,parcel);
    % save data for R visual  TASK 7
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2:8) = reshape(spin_task7(:,i,:),n_nulls,7); 
    
        temp(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        temp(1+n_nulls*(i-1):n_nulls*i,2) = reshape(spin_mean_task7(:,i),n_nulls,1); 
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_spin','.xlsx'),data,strcat('A2:H',num2str(n_modes*n_nulls+1)))
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7_subj_null_spin','.xlsx'),temp,strcat('A2:B',num2str(n_modes*n_nulls+1)))
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculating p-value
    P_rewired7 = zeros(n_modes,7);P_spin7 = zeros(n_modes,7);
    P_rewired = zeros(n_modes,1);P_spin = zeros(n_modes,1);
    for i =1:n_modes
        for j = 1:7
            P_rewired7(i,j) = sum(rewired_task7(:,i,j)>group_acc_task7(i,j))/n_nulls;  
            P_spin7(i,j) = sum(spin_task7(:,i,j)>group_acc_task7(i,j))/n_nulls;
        end
        P_rewired(i,1) = sum(rewired_mean_task7(:,i)>acc_mean_task7(i,1))/n_nulls;  
        P_spin(i,1) = sum(spin_mean_task7(:,i)>acc_mean_task7(i,1))/n_nulls;
     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    output_rewired7  = [(1:n_modes)',P_rewired7];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_rewired_Pvalue','.xlsx'),output_rewired7, strcat('A2:H',num2str(n_modes+1)));
    output_spin7  = [(1:n_modes)',P_spin7];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_key_task7_subj_null_spin_Pvalue','.xlsx'),output_spin7, strcat('A2:H',num2str(n_modes+1)));
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7_subj_null_rewired_Pvalue','.xlsx'),output_rewired, strcat('A2:B',num2str(n_modes+1)));
    output_spin  = [(1:n_modes)',P_spin];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_corr_mean_key_task7_subj_null_spin_Pvalue','.xlsx'),output_spin, strcat('A2:B',num2str(n_modes+1)));
end

end

