function [MSFD_log,MSFD_thr_log] = Compute_MSFD(XrandS,U,Vlow,Vhigh,mean_d,mean_c,N_d,N_c)

%% find coupled/decoupled components of surrogate signals and test significance of MSFD

%% 1) FILTER coupled/decoupled signal portions in surrogates
fdata=XrandS;        %MSC-informed graph signal surrogates 
n_ROI = size(U,1);

clear X_hat_surr X_c_surr X_d_surr X_m_surr N_c_surr N_d_surr
for s=1:size(fdata,2)   
    for i=1:size(fdata{1},3)  %1:19
        X_hat_surr{s}(:,:,i)=U'*fdata{s}(:,:,i);
        X_c_surr{s}(:,:,i)=Vlow*X_hat_surr{s}(:,:,i);     
        X_d_surr{s}(:,:,i)=Vhigh*X_hat_surr{s}(:,:,i);
        % norms  of the weights
        for r=1:size(fdata{1},1)
            N_c_surr(r,i,s)=norm(X_c_surr{s}(r,:,i));
            N_d_surr(r,i,s)=norm(X_d_surr{s}(r,:,i));
        end
    end
end


%% structural decoupling index structural decoupling index
%consider and test mean across subjects
MSFD_surr=N_d_surr./N_c_surr;                      % node*19*sunject
MSFD_surr_avgsurr=squeeze(mean(MSFD_surr,2));       
MSFD_surr_avgsurrsubjs=mean(MSFD_surr_avgsurr,2);    

mean_MSFD=mean_d./mean_c;          %empirical AVERAGE MSFD, group-level
MSFD=N_d./N_c;                     %emipirical individual MSFD, node*subject

%find threshold for max
%for every subject, max across surrogates
for s=1:size(MSFD_surr,3)
       max_MSFD_surr(:,s)=max(MSFD_surr(:,:,s)')';
end

%find threshold for min
%for every subject, in across surrogates
for s=1:size(MSFD_surr,3)
      min_MSFD_surr(:,s)=min(MSFD_surr(:,:,s)')';
end

%% select significant MSFD for each subject, across  surrogates
%individual thr, 
%for each subject, I threshold the ratio based on individual ratio's surrogate distribution 
for s=1:size(fdata,2) 
    significant_values_max(s)=size(find(MSFD(:,s)>max_MSFD_surr(:,s)),1);
    significant_values_min(s)=size(find(MSFD(:,s)<min_MSFD_surr(:,s)),1);
    MSFD_thr_max(:,s)=MSFD(:,s)>max_MSFD_surr(:,s);
    MSFD_thr_min(:,s)=MSFD(:,s)<min_MSFD_surr(:,s);
    detect_max=sum(MSFD_thr_max'); %amounts of detection per region
    detect_min=sum(MSFD_thr_min');
end

%%for every region, test across subjects 0.05, correcting for the number of tests (regions), 0.05/n_ROI
x=0:1:100;
y=binocdf_old(x,100,0.05,'upper');            
THRsubjects=x(min(find(y<0.05/n_ROI))); 
THRsubjects=floor(size(fdata,2)/100*THRsubjects)+1;

MSFD_sig_higher=detect_max>THRsubjects;
MSFD_sig_lower=detect_min>THRsubjects;

MSFD_sig_higher_positions=find(MSFD_sig_higher==1);
MSFD_sig_lower_positions=find(MSFD_sig_lower==1);

MSFD_sig_tot_positions=[find(MSFD_sig_higher==1),find(MSFD_sig_lower==1)];
MSFD_sig_tot_positions=sort(unique(MSFD_sig_tot_positions));

%%threshold empirical mean ratios
mean_MSFD_thr=ones(n_ROI,1);
mean_MSFD_thr(MSFD_sig_tot_positions)=mean_MSFD(MSFD_sig_tot_positions);


%% Fig surr MSFD pattern
saturate=1;
CC2=log2(MSFD_surr_avgsurrsubjs); 
%% adjust Cvalues for saturation (to eliminate outliers peaks)
if saturate
    thr=1;
    CC2new=CC2;
    CC2new(find(CC2>thr))=0;
    CC2new(find(CC2>thr))=max(CC2new);
    CC2new(find(CC2<-thr))=0;
    CC2new(find(CC2<-thr))=min(CC2new);
    CC2=CC2new;
end
MSFD_surr_log = CC2;

%% Fig mean_MSFD
saturate=1;
CC2=log2(mean_MSFD); 
if saturate
    thr=1;
    CC2new=CC2;
    CC2new(find(CC2>thr))=0;
    CC2new(find(CC2>thr))=max(CC2new);
    CC2new(find(CC2<-thr))=0;
    CC2new(find(CC2<-thr))=min(CC2new);
    CC2=CC2new;
end
MSFD_log = CC2;

%% significant decoupling and coupling regions
saturate=1;
CC2=log2(mean_MSFD_thr); 
if saturate
    thr=1;
    CC2new=CC2;
    CC2new(find(CC2>thr))=0;
    CC2new(find(CC2>thr))=max(CC2new);
    CC2new(find(CC2<-thr))=0;
    CC2new(find(CC2<-thr))=min(CC2new);
    CC2=CC2new;
end
MSFD_thr_log =  CC2;

end









