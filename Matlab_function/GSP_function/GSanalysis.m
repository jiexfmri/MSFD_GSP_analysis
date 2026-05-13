function [N_c,N_d,mean_c,mean_d] = GSanalysis(zX_RS,Vhigh,Vlow,U)

%GRAPH SIGNAL ANALYSIS 
data=zX_RS;             

n_ROI = size(data,1);
nsubjs_RS = size(data,3);

%% compute fMRI HF/LF portions   
clear X_hat X_c X_d X_m N_c N_d X_all 
for s=1:nsubjs_RS                
    X_hat(:,:,s)=U'*data(:,:,s);      %GFT
    X_c(:,:,s)=Vlow*X_hat(:,:,s);     % coupling filtered signal
    X_d(:,:,s)=Vhigh*X_hat(:,:,s);    % decoupling filtered signal
    X_all(:,:,s)=U*X_hat(:,:,s);      % reconstruct back full signal to check norm
    
    %% norms  of the weights
    for r=1:n_ROI
        N_c(r,s)=norm(X_c(r,:,s));    
        N_d(r,s)=norm(X_d(r,:,s));
    end
end
%% mean across subjects
mean_c=mean(N_c,2); %average coupling
mean_d=mean(N_d,2); %average decoupling

end
