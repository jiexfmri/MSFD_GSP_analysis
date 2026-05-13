function [rewired_activity_ratio, rewired_FC_R]= null_rewired_rsfMRI(W,X_RS,n_modes,n_nulls);

% =========================================================================
%  Generate 1000 random undigraphs£º[W0, R] = null_model_und_sign(W)
%  further calculate reconstruction accuracy
% =========================================================================

%  load fMRI time series
n_ROI = size(W,1);
n_subjs=size(X_RS,3);
mean_data = mean(X_RS,2);
zX_RS = X_RS - mean_data;       % dmean mean centering  for activity
RX_RS = X_RS;                   % raw data for FC

rewired_activity_ratio = zeros(n_nulls,n_modes);
rewired_FC_R = zeros(n_nulls,n_modes);

Norm_emp = zeros(n_ROI,n_subjs);
for s=1:n_subjs    
     for r=1:n_ROI      
          %total empirical energy
          Norm_emp(r,s)  = norm(zX_RS(r,:,s));  
     end
end

% find indices of all values above the diagonal
triu_ind = find(triu(ones(n_ROI,n_ROI),+1)==1); 
% FC of empirical signals
FCvec_emp = zeros(length(triu_ind), n_subjs);
for s = 1: n_subjs
    data_emp = RX_RS(:,:,s);
    FC_emp = squeeze(corr(data_emp'));
    FCvec_emp(:,s)=FC_emp(triu_ind);
end

%% rewired null model 
parfor null = 1:n_nulls  
    [W_null, ~] = null_model_und_sign(W);   
    U = Compute_Laplacian(W_null);   %% Laplacian Decomposition

    %% reconstruction activity
    recon_activity_ratio = zeros(n_modes,n_subjs);
    for mode= 1:n_modes
         for s=1:n_subjs
             
            X_hat=U'*zX_RS(:,:,s);      
            basis=zeros(size(U));
            basis(:,1:mode)=U(:,1:mode);   
            X_all=basis*X_hat;      % reconstruct back full signal to check norm

            %% norms of reconstructed BOLD-fMRI
            Norm_recon = zeros(n_ROI,1);
            for r=1:n_ROI
                Norm_recon(r,1)=norm(X_all(r,:));
            end
            % Method 1  for every subject's activity energy
            % activity ratio
            recon_activity_ratio(mode,s) = mean(Norm_recon(:,1))/mean(Norm_emp(:,s));  
         end
    end
    activity_ratio = mean(recon_activity_ratio,2);
    rewired_activity_ratio(null,:) = activity_ratio';

    %% reconstruction FC
    recon_FC_R = zeros(n_modes, n_subjs);
    for mode= 1:n_modes
        for s=1:n_subjs
            X_hat=U'*RX_RS(:,:,s); % Calculate reconstruction beta coefficients
            basis=zeros(size(U)); 
            basis(:,1:mode)=U(:,1:mode);   
            X_all=basis*X_hat;  % reconstructed BOLD-fMRI

            FC_recon = squeeze(corr(X_all'));
            FCvec_recon = FC_recon(triu_ind);
            % calc r
            recon_FC_R(mode,s)=corr(FCvec_emp(:,s),FCvec_recon);   
        end    
    end
    FC_R = mean(recon_FC_R,2); 
    
    rewired_FC_R(null,:) = FC_R'; 
end

end