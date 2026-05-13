function [rela_hierachical] = Relationship_Hierarchy(MSFD,parcel)

n_nulls = 10000;
n_ROI = size(MSFD,1);
data = importdata(strcat('Data',filesep, parcel,'_Hierarchy.xlsx'));
temp = data.data(:,2:end);
hier = temp;

load(strcat('Data',filesep, parcel,'_PermID.mat'))

rela_hierachical = zeros(size(hier,2),2);
for i=1:size(hier,2)
    [r(i),p(i)] = corr(MSFD,hier(:,i),'type','Spearman');  
    p_perm(i) = spin_perm_sphere_p(MSFD, hier(:,i), permID, 'spearman');
    rela_hierachical(i,:)=[r(i),p_perm(i)]; 
end

end



