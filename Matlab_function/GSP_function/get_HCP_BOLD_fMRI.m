function zX_RS= get_HCP_BOLD_fMRI(parcel)

X_RS = load(strcat('Data',filesep, parcel,'_RSfMRI_example.mat'));
X_RS = X_RS.X_RS;


% Normalized fMRI timecourses
zX_RS=zscore(X_RS,0,2);

end
