function [T,S,Y,Psi_true,Theta_true,Omega_true]=Simulate_Gaussian_data(n,p)
%Here we generate two sparse precision matrices. One over covariates (cells, time points), and one over samples (genes, pixels).
s = RandStream('mcg16807','Seed', 1e+5);
RandStream.setGlobalStream(s)
Theta_true = genspd(p, .1, 1, 2, true);
Theta_true = Theta_true/mean(diag(Theta_true));
Psi_true = genspd(n, .1, 1, 2, true);%generate a sparse positive definite matrix of size q and Nonzero/Zero ratio=0.1.
Psi_true = Psi_true/mean(diag(Psi_true)); 
Omega_true = kron(Psi_true, eye(p)) + kron(eye(n), Theta_true);                             % The Kronecker-sum (or Cartesian product) of Psi and Theta.
% Sample (vectorised) random matrices from Gaussian with a structured Kronecker-sum precision.
figure(1), clf
subplot(131), imagesc(Psi_true), title('\Psi_0')
colorbar % Visualise precisions. 
subplot(132), imagesc(Theta_true), title('\Theta_0')
colorbar
subplot(133), imagesc(Omega_true), title('\Psi_0 \otimes \Theta_0')
colorbar
m=100;
Y = gaussSamp(pdinv(Omega_true), m);  
Y = Y - repmat(mean(Y), m, 1); % Remove mean matrix.
Y(isnan(Y)) = 0;% Remove NaNs.
Y = permute(reshape(Y', p, n, m), [2 1 3]); 
S=zeros(p);
T=zeros(n);% Compute sufficient statistics for iid matrices,
for i = 1:m  
    S = S + cov(Y(:,:,i),1);%Covariance across columns (S),
    T = T + cov(Y(:,:,i)',1);%Covaraiance across rows (T)
end
T = T / m;  S = S / m; %Take the average over 100 synthetic samples
figure(15)
subplot(121), imagesc(T), title('simT')
colorbar
subplot(122), imagesc(S), title('simS')
end