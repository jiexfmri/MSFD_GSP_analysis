function [output_acc,output_AUC] = Recon_activity_FC_subj(W,X_RS,n_modes,n_nulls, parcel,type)

% save file path 
file_name = strcat('Results/recon_activity_FC/',type);
if~exist(file_name,'dir')
    mkdir(file_name);
end

%%% =================================================================
%   (1) Calculate connectome eigenmodes   
% ====================================================================
n_ROI = size(W,1);

U = Compute_Laplacian(W);

%%% =================================================================
%   (2) load fMRI time series
% ====================================================================
n_subjs=size(X_RS,3);
mean_data = mean(X_RS,2);
zX_RS = X_RS - mean_data;       % dmean mean centering
RX_RS = X_RS;                   % raw data for FC


%% human rs-fMRI data projected on the structural eigenmodes
%%% =================================================================
%  Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
% =========================================================================

%% reconstruction brain activity(BOLD-fMRI)
Norm_emp = zeros(n_ROI,n_subjs);
recon_activity_ratio = zeros(n_modes,n_subjs);

for s=1:n_subjs    
     for r=1:n_ROI      
          %total empirical energy
          Norm_emp(r,s)  = norm(zX_RS(r,:,s));  
     end
end

parfor mode= 1:n_modes
    for s=1:n_subjs    
        X_hat=U'*zX_RS(:,:,s);      
        basis=zeros(size(U));
        basis(:,1:mode)=U(:,1:mode);   
        X_all=basis*X_hat;      % reconstruct back full signal to check norm

        % norms of reconstructed BOLD-fMRI
        Norm_recon = zeros(n_ROI,1);
        for r=1:n_ROI
            Norm_recon(r,1)=norm(X_all(r,:));
        end
        % activity ratio
        recon_activity_ratio(mode,s) = mean(Norm_recon(:,1))/mean(Norm_emp(:,s));  
    end
end
activity_ratio = mean(recon_activity_ratio,2);

clear Norm_recon Norm_emp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconstruction FC
% find indices of all values above the diagonal
triu_ind = find(triu(ones(n_ROI,n_ROI),+1)==1); 

% FC of empirical signals
FCvec_emp = zeros(length(triu_ind), n_subjs);
recon_FC_R = zeros(n_modes, n_subjs);
for s = 1: n_subjs
    data_emp = RX_RS(:,:,s);
    FC_emp = squeeze(corr(data_emp'));
    FCvec_emp(:,s)=FC_emp(triu_ind);
end

parfor mode= 1:n_modes
    for s=1:n_subjs
        % Calculate reconstruction beta coefficients
        X_hat=U'*RX_RS(:,:,s); 
        % reconstructed BOLD-fMRI
        basis=zeros(size(U));
        basis(:,1:mode)=U(:,1:mode);   
        X_all=basis*X_hat;  
        
        FC_recon = squeeze(corr(X_all'));
        FCvec_recon = FC_recon(triu_ind);
        % calc r
        recon_FC_R(mode,s)=corr(FCvec_emp(:,s),FCvec_recon);   
    end    
end
FC_R = mean(recon_FC_R,2); 

% Reconstruction accuracy of each subject
data_out  = [(1:n_modes)',recon_activity_ratio];
xlswrite(strcat(file_name,filesep,parcel,'_',type,'_Recon_activity_ratio_subj','.xlsx'),data_out);

data_out  = [(1:n_modes)',recon_FC_R];
xlswrite(strcat(file_name,filesep,parcel,'_',type,'_Recon_FC_R_subj','.xlsx'),data_out);

% Mean reconstruction accuracy
output_acc  = [(1:n_modes)',activity_ratio,FC_R];
xlswrite(strcat(file_name,filesep,parcel,'_',type,'_Recon_rsfMRI_subj_mean','.xlsx'),output_acc, strcat('A2:C',num2str(n_modes+1)));

for mode= 1:n_modes
    activity_ratio_AUC(mode,1)=trapz(activity_ratio(1:mode,1))/n_modes;
    FC_R_AUC(mode,1)=trapz(FC_R(1:mode,1))/n_modes;
end
output_AUC  = [(1:n_modes)',activity_ratio_AUC,FC_R_AUC];
xlswrite(strcat(file_name,filesep,parcel,'_',type,'_Recon_rsfMRI_subj_mean','.xlsx'),output_AUC, strcat('F2:H',num2str(n_modes+1)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% performer null models

null_model = 0;
if null_model == 1
    disp('Running null model...')

    %% fusing the TS matrix with a random noise matrix of similar properties to MPC 
    [TS_rewired_MPC_activity_ratio, TS_rewired_MPC_FC_R]= null_TS_rewired_MPC_rsfMRI(X_RS,n_modes,n_nulls);
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2) = TS_rewired_MPC_activity_ratio(:,i); 
        data(1+n_nulls*(i-1):n_nulls*i,3) = TS_rewired_MPC_FC_R(:,i);  
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_TS_rewired_MPC','.xlsx'),data,strcat('A2:C',num2str(n_modes*n_nulls+1)))
    % Calculating p-value
    P_rewired = zeros(n_modes,3);
    for i =1:n_modes
        P_rewired(i,1) = sum(TS_rewired_MPC_activity_ratio(:,i)>activity_ratio(i,1))/n_nulls;
        P_rewired(i,2) = sum(TS_rewired_MPC_FC_R(:,i)> FC_R(i,1))/n_nulls;
    end
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_TS_rewired_MPC_Pvalue','.xlsx'),output_rewired, strcat('A2:C',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fusing the TS matrix with a random noise matrix of similar properties to GD 
    [TS_rewired_GD_activity_ratio, TS_rewired_GD_FC_R]= null_TS_rewired_GD_rsfMRI(X_RS,n_modes,n_nulls);
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2) = TS_rewired_GD_activity_ratio(:,i); 
        data(1+n_nulls*(i-1):n_nulls*i,3) = TS_rewired_GD_FC_R(:,i);  
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_TS_rewired_GD','.xlsx'),data,strcat('A2:C',num2str(n_modes*n_nulls+1)))
    % Calculating p-value
    P_rewired = zeros(n_modes,3);
    for i =1:n_modes
        P_rewired(i,1) = sum(TS_rewired_GD_activity_ratio(:,i)>activity_ratio(i,1))/n_nulls;
        P_rewired(i,2) = sum(TS_rewired_GD_FC_R(:,i)> FC_R(i,1))/n_nulls;
    end
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_TS_rewired_GD_Pvalue','.xlsx'),output_rewired, strcat('A2:C',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% fusing the TS matrix with a random noise matrix of similar properties to MPC and GD 
    [TS_rewired_MPC_GD_activity_ratio, TS_rewired_MPC_GD_FC_R]= null_TS_rewired_MPC_GD_rsfMRI(X_RS,n_modes,n_nulls);
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2) = TS_rewired_MPC_GD_activity_ratio(:,i); 
        data(1+n_nulls*(i-1):n_nulls*i,3) = TS_rewired_MPC_GD_FC_R(:,i);  
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_TS_rewired_MPC_GD','.xlsx'),data,strcat('A2:C',num2str(n_modes*n_nulls+1)))
    % Calculating p-value
    P_rewired = zeros(n_modes,3);
    for i =1:n_modes
        P_rewired(i,1) = sum(TS_rewired_MPC_GD_activity_ratio(:,i)>activity_ratio(i,1))/n_nulls;
        P_rewired(i,2) = sum(TS_rewired_MPC_GD_FC_R(:,i)> FC_R(i,1))/n_nulls;
    end
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_TS_rewired_MPC_GD_Pvalue','.xlsx'),output_rewired, strcat('A2:C',num2str(n_modes+1)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% rewired MSC null model
    [rewired_activity_ratio, rewired_FC_R]= null_rewired_rsfMRI(W,X_RS,n_modes,n_nulls);
    % save data for R visual
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2) = rewired_activity_ratio(:,i); 
        data(1+n_nulls*(i-1):n_nulls*i,3) = rewired_FC_R(:,i);  
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_rewired','.xlsx'),data,strcat('A2:C',num2str(n_modes*n_nulls+1)))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% spin null model
    [spin_activity_ratio, spin_FC_R]= null_spin_rsfMRI(W,X_RS,n_modes,n_nulls,parcel);
    % save data for R visual
    for i = 1: n_modes
        data(1+n_nulls*(i-1):n_nulls*i,1) = i;   
        data(1+n_nulls*(i-1):n_nulls*i,2) = spin_activity_ratio(:,i); 
        data(1+n_nulls*(i-1):n_nulls*i,3) = spin_FC_R(:,i);  
    end
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_spin','.xlsx'),data,strcat('A2:C',num2str(n_modes*n_nulls+1)))
    
    %% Calculating p-value
    P_rewired = zeros(n_modes,3);
    P_spin = zeros(n_modes,3);
    for i =1:n_modes
        P_rewired(i,1) = sum(rewired_activity_ratio(:,i)>activity_ratio(i,1))/n_nulls;
        P_rewired(i,2) = sum(rewired_FC_R(:,i)> FC_R(i,1))/n_nulls;
    
        P_spin(i,1) = sum(spin_activity_ratio(:,i)>activity_ratio(i,1))/n_nulls;
        P_spin(i,2) = sum(spin_FC_R(:,i)> FC_R(i,1))/n_nulls;
    end
    
    output_rewired  = [(1:n_modes)',P_rewired];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_rewired_Pvalue','.xlsx'),output_rewired, strcat('A2:C',num2str(n_modes+1)));
    
    output_spin  = [(1:n_modes)',P_spin];
    xlswrite(strcat(file_name,'/',parcel,'_',type,'_Recon_rsfMRI_subj_null_spin_Pvalue','.xlsx'),output_spin, strcat('A2:C',num2str(n_modes+1)));
end

end