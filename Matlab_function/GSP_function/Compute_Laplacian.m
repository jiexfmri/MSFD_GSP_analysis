function U = Compute_Laplacian(W)
                                     
n_ROI = size(W,2);  

sparse_data = W;

%%% We represent graph G by using the top 10% of connections per row 
% sparsity = 90;
% sparse_data(W < prctile(W,sparsity)) = 0; 

sparse_data(W < 0) = 0; 
A = sparse_data;
A = A-diag(diag(A));
A = (A+A')./2;

%% Symmetric Normalization of adjacency matrix
D=diag(sum(A,2)); %degree
Wsymm=D^(-1/2)*A*D^(-1/2);
Wnew=Wsymm;

%% compute normalized Laplacian
L=eye(n_ROI)-Wnew;
[U,LambdaL] = eig(L);
[LambdaL, IndL]=sort(diag(LambdaL));
U=U(:,IndL);
orth_U = U*U';  %Verify orthogonality 

end