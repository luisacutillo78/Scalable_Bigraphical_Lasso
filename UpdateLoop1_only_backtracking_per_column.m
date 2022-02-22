function [diffTheta,diffPsi,Psi,Theta,objectiveFunction]=UpdateLoop1_only_backtracking_per_column(S,T,beta,varargin)
if all(beta < 0)%% ARG lambda : the L1 shrinkage (regularisation) parameter.

error('regularization parameters must be non-negative.')
end
for k = 1:2:length(varargin)
    %%varargin is an input variable in a function definition statement 
    %%that enables the function to accept any number of input arguments.
    switch varargin{k}
        case 'maxIter'; maxIter = varargin{k+1};
        case 'thresh'; thresh = varargin{k+1};
    end
end

if ~exist('maxIter','var')
    maxIter = 100;           %%default maximum iterations                                               
    % Maximum number of glasso iterations.
end
if ~exist('thresh','var')
    thresh = 1e-4;           %%default threshold value                                               % Convergence threshold on W (estimated inverse-covariance).
end
converged = false; %%indicates whether it is converged
iter = 0;
objectiveFunction=[];
n=size(T,1);
p=size(S,2);
%Initialise Psi and Theta
Psi=T;
Theta=S;
while (~converged && (iter < maxIter))
    [Theta_new,Psi_new,log_det_Omega]=flip_flop_only_backtracking_per_column(S,T,Psi,Theta,beta);
    size_Theta=size(Theta,1);
    diffTheta=Theta_new-Theta - eye(size_Theta)*mean(diag((Theta_new - Theta)));  
    size_Psi=size(Psi,1);
    diffPsi=Psi_new-Psi - eye(size_Psi)*mean(diag((Psi_new - Psi)));
    iter = iter + 1;
    %%Determine the overall convergence via Frobenius norm of the
    %%difference between two updates
    nrm(iter)=norm(diffTheta,'fro')^2+norm(diffPsi,'fro')^2;
  if max(nrm(max(1,iter-3):iter)) < thresh
        converged = true;
  end
    %Store each update and objective function for tracking
    %PsiS{iter}{1}=Theta;
    %PsiS{iter}{2}=Psi;
    Theta=Theta_new;
    Psi=Psi_new;
    objectiveFunction(iter,1)=n*trace(Theta*S)+p*trace(Psi*T)-log_det_Omega+beta(2)*norm(Theta,1)+beta(1)*norm(Psi,1);
    objectiveFunction(iter,2)=p*trace(Psi*T)-0.5*log_det_Omega+beta(1)*norm(Psi,1);
    objectiveFunction(iter,3)=n*trace(Theta*S)-0.5*log_det_Omega+beta(2)*norm(Theta,1);
    fprintf('iter: %d diff: %1.16f\n', iter,nrm(iter))

end   
end