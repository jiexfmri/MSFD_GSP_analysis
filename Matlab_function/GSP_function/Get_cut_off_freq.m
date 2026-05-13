function [NN,Vlow, Vhigh]=Get_cut_off_freq(U,zX_RS)

n_ROI = size(U,1);

%% Average energy spectral density of rs-fMRI data projected on the MSC eigenmode
clear X_hat_L 
for s=1:size(zX_RS,3)
    X_hat_L(:,:,s)=U'*zX_RS(:,:,s);
end
%%energy spectral density 
pow=abs(X_hat_L).^2;
PSD=squeeze(mean(pow,2)); %mean across time

%% compute cut-off frequency(low- and high frequency) 
%mean across subjects/epochs
mPSD=mean(PSD,2);
AUCTOT=trapz(mPSD(1:size(U,1))); %total area under the curve

i=0;
AUC=0;
while AUC<AUCTOT/2
    i=i+1;
    AUC=trapz(mPSD(1:i));
end

NN=i-1; %CUTOFF FREQUENCY : number of low frequency eigenvalues to consider in order to have the same energy as the high freq ones

%% split structural harmonics in high/low frequency
Vlow=zeros(size(U));
Vhigh=zeros(size(U));
Vlow(:,1:NN)=U(:,1:NN);  %low frequencies = coupled
Vhigh(:,NN+1:end)=U(:,NN+1:end);  %high frequencies= decoupled

end




