function [XrandS] = GSrandomozation_surrogates(data,U,nSurr)


%% Create CC-informed graph signal surrogates by randomization of 
%%real harmonics Fourier coefficients
clear XrandS 

%% SPATIAL RANDOMIZATION
X_hat_rand=zeros(size(data,1),size(data,2),size(data,3),nSurr);   %nodes*times*subject*rand
for s=1:size(data,3)
    X=data(:,:,s);
    for n=1:nSurr
        %randomize sign of Fourier coefficients
        PHIdiag=round(rand(size(U,1),1));
        PHIdiag(PHIdiag==0)=-1;
        PHI=diag(PHIdiag);                  
        XrandS{s}(:,:,n)=U*PHI*(U'*X);         
    end
end



end

