function [S,T,Y,Psi_true,Theta_true,Omega_true]=Simulate_data(n,p)
%Here we generate two sparse precision matrices. One over covariates (cells, time points), and one over samples (genes, pixels).
s = RandStream('mcg16807','Seed', 1e+5);
RandStream.setGlobalStream(s)
Theta_true = genspd(p, .1, 1, 2, true);
Theta_true = Theta_true/mean(diag(Theta_true));
Psi_true = genspd(n, .1, 1, 2, true); %generate a sparse positive definite matrix of size q and Nonzero/Zero ratio=0.1.
Psi_true = Psi_true/mean(diag(Psi_true)); 
Omega_true = kron(Psi_true, eye(p)) + kron(eye(n), Theta_true);% The Kronecker-sum (or Cartesian product) of Psi and Theta.
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
Y = Y - repmat(mean(Y), m, 1); % Remove mean matrix.                                                              % Normalise matrices st each matrix has a unit norm.
Y(isnan(Y)) = 0; % Remove NaNs.
Y = permute(reshape(Y', p, n, m), [2 1 3]); 
S=zeros(p);
T=zeros(n);% Compute sufficient statistics for iid matrices,
ct=0;
cs=0;
for i = 1:m  
         pd = makedist('Normal');
         PY=cdf(pd,Y(:,:,i));
         %Convert it to count data
         Y(:,:,i)=nbininv(PY,10,0.9);
         dum=corr(Y(:,:,i)','Type','Kendall');
         vec=isnan(dum);
         %across rows (T),
            if(sum(sum(vec))==0)
                 T = T+sin(pi/2*dum);
                 ct=ct+1;
            end      
            %across columns (S),
            dum=corr(Y(:,:,i),'Type','Kendall');
            vec=isnan(dum);
            if(sum(sum(vec))==0)
               S = S+sin(pi/2*dum);   
               cs=cs+1;
            end
end
T = T / ct;  S = S / cs; 
figure(15)
subplot(121), imagesc(T), title('simT')
colorbar
subplot(122), imagesc(S), title('simS')
end