function [Theta_new,Psi_new,log_det_Omega]=flip_flop_only_backtracking_per_column(S,T,Psi,Theta,beta)
 %If the regularisation parameter provided is only one dimensional, 
 %use the same regularization parameter for both rows and columns.
if(length(beta)<2)
        beta(2)=beta(1);
end
n = size(T, 2);     p = size(S, 2); %column number of T,S
[~, pd] = chol(T + beta(1)*eye(n)); 
    %%returns pd indicating whether (T+lambda(1)*eye(n)) is symmetric positive definite.
    %%pd=0-> (T+lambda(1)*eye(n)) is symmetric positive definite
    % ISSUE: depending on the relative size of n and p, either S or T will be low-rank.
if pd > 0
    error('T is not a proper covariance.')
end
[~, pd] = chol(S + beta(2)*eye(p));
if pd > 0
    error('S is not a proper covariance.')
end
%Initialise options for glmnet
options = glmnetSet; 
options.L1GPoptions.verbose = 0;
options.ltype='modified.Newton';
%eigen-decompose Psi and Theta
[U,L1]=eig(Psi);
%%produce L1: diagonal matrix L1 with eig(Psi) as diagonal
%%produce U: whose columns are corresponding eigenvectors
[V,L2]=eig(Theta);
%%produce L2: diagonal matrix L2 with eig(Theta) as diagonal
%%produce V: hwose columns are corresponding eigenvectors                                                  
% Regularization parameter.
options.lambda=beta(1);
%%\tilde{S_k}=S_k-(K-1)/K*tr(S_k)/d_k*I_{d_k} 
%%see teralasso.m line 74
S=S-0.5*mean(diag(S))*eye(size(S));
T=T-0.5*mean(diag(T))*eye(size(T));
Psi_new=sub_biglassoFastLoop_only_backtracking_per_column(Psi,U,L1,L2,T,options);
[~,L1new]=eig(Psi_new);
options.lambda=beta(2);
Theta_new = sub_biglassoFastLoop_no_eig_per_column(Theta,V,L2,L1new,S,options); 
% First optimise Psi with Theta fixed, then optimise Theta with Psi fixed.
[~,L2new]=eig(Theta_new);
L1=L1new;
L2=L2new;
lambda1=diag(L1);%n dim
lambda1(lambda1==0)=[];
lambda2=diag(L2);%p dim
lambda2(lambda2==0)=[];
log_det_Omega=0;
for i=1:n
        dum=log(abs(lambda1(i)+lambda2));
        log_det_Omega=log_det_Omega+sum(dum);
end   
end