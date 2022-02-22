function [S,T,Y,Psi_true,Theta_true,Omega_true]=Simulate_data_theta3blocks(n,p)
%Here we generate two sparse precision matrices. One over covariates (cells, time points), and one over samples (genes, pixels).
s = RandStream('mcg16807','Seed', 1e+5);
RandStream.setGlobalStream(s)
%generate Theta_true with three blocks on diagonals
Theta_true=-1*blkdiag(ones(p/3,p/3),ones(p/3,p/3),ones(p/3,p/3))+randn(p,p)/10+2*eye(p);
%Convert Theta to the nearest positive definite matrix
Theta_true=nearestSPD(Theta_true);
Theta_true = Theta_true/mean(diag(Theta_true));

Psi_true = genspd(n, .1, 1, 2, true);  %generate a sparse positive definite matrix of size q and Nonzero/Zero ratio=0.1.
Psi_true = Psi_true/mean(diag(Psi_true)); 

Omega_true = kron(Psi_true, eye(p)) + kron(eye(n), Theta_true);                             % The Kronecker-sum (or Cartesian product) of Psi and Theta.
 % Sample (vectorised) random matrices from Gaussian with a structured Kronecker-sum precision.
figure(1), clf
subplot(131), imagesc(Psi_true), title('\Psi_0')
colorbar % Visualise precisions.    [W, W_dual, Psi_hat, Theta_hat,of] = biglasso(S, T, lambda(i), 'maxIter', 20, 'thresh', 1e-8, 'warmInit', warmInit);
subplot(132), imagesc(Theta_true), title('\Theta_0')
colorbar
subplot(133), imagesc(Omega_true), title('\Psi_0 \otimes \Theta_0')
colorbar
m=100;
Y = gaussSamp(pdinv(Omega_true), m);  
Y = Y - repmat(mean(Y), m, 1);   % Remove mean matrix.
Y(isnan(Y)) = 0;        % Remove NaNs.
Y = permute(reshape(Y', p, n, m), [2 1 3]); 

S=zeros(p);
T=zeros(n);% Compute sufficient statistics for iid matrices,
for i = 1:m  
         pd = makedist('Normal');
         PY=cdf(pd,Y(:,:,i));
         Y(:,:,i)=nbininv(PY,10,0.9);%Convert Y to count data 
         %  across rows (T).
         T = T+sin(pi/2*corr((Y(:,:,i))','Type','Kendall'));       
         %   across columns (S),
         S = S+sin(pi/2*corr(Y(:,:,i),'Type','Kendall'));%T1 + cov(Y(:,:,i)',1); 
end         
T = T / m;  S = S / m; 
T=nearestSPD(T);
S=nearestSPD(S);
figure(15)
subplot(121), imagesc(T), title('simT')
colorbar
subplot(122), imagesc(S), title('simS')
end