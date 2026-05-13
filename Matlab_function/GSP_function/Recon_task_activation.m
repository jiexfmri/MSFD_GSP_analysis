function [output_acc,recon_task]= Recon_task_activation(W,X_RS,n_modes)

%%% =================================================================
%   (1) Calculate connectome eigenmodes   
% ====================================================================
n_ROI = size(W,1);

U = Compute_Laplacian(W);   %% Laplacian Decomposition

%%% =================================================================
%   (2) load fMRI time series
% ====================================================================
n_subjs=size(X_RS,2);
zX_RS = X_RS;

%% human tfMRI data projected on the structural eigenmodes
%%% =================================================================
%   (3) Calculate reconstruction accuracy using 1 to num_modes eigenmodes    
% =========================================================================
%% reconstruction brain activity(BOLD-fMRI)
recon_corr = zeros(n_modes,n_subjs);
recon_task = zeros(n_ROI,n_modes,n_subjs);
parfor s=1:n_subjs 
    for mode = 1:n_modes
   
        X_hat=U'*zX_RS(:,s);      
        basis=zeros(size(U));
        basis(:,1:mode)=U(:,1:mode);   
        X_all=basis*X_hat;      % reconstruct back full signal
        
        recon_task(:,mode,s) = X_all;
        % compute correlation 
        recon_corr(mode,s) = corr(X_all, zX_RS(:,s));
    end
end

output_acc  = recon_corr;
end