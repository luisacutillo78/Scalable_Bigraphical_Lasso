mypath='./';
plotdir='./plots/';
outputdir='./Output';
fullfile(mypath,'ndlutil/matlab');
addpath(fullfile(mypath,'ndlutil/matlab'));
addpath(fullfile(mypath,'/rca/matlab'));
addpath(fullfile(mypath,'glmnet_matlab'));
addpath(fullfile(mypath,'L1General'));
addpath(fullfile(mypath,'teralasso-master'));
beta=[0.007,0.01];
N=20:20:100;
for i=1:length(N)%1:length(lambda)
    p=N(i);
    n=N(i);
[S,T,Y,Psi_true,Theta_true,Omega_true]=Simulate_Gaussian_data(n,p);
tic
[W, W_dual, Psi_hat, Theta_hat,of] = biglasso(S, T, beta, 'maxIter', 100, 'thresh', 1e-4);
et_biglasso(i)=toc;
tic
[~,~,Psi_hatF, Theta_hatF,ofF] = UpdateLoop1_only_backtracking_per_column(S, T, beta, 'maxIter', 100, 'thresh', 1e-4);
et_scBiglasso(i)=toc;
SS{1}=T;
SS{2}=S;
ps=[size(T,1),size(S,2)];
tol = 1e-4;
maxiter = 100;
tic
[ PsiTheta ] = teralasso( SS,ps,'L1',0,tol,beta,maxiter);
et_TeraLasso(i)=toc;
Psi_hatT=PsiTheta{1};
Theta_hatT=PsiTheta{2};
filename=sprintf('Psi_biglasso_gaussSynth_p_%d_n_%d_lambda1_%d_lambda2_%d.txt',p,n,beta(1),beta(2));
writematrix(Psi_hat,filename);
filename=sprintf('Theta_biglasso_gaussSynth_p_%d_n_%d_lambda1_%d_lambda2_%d.txt',p,n,beta(1),beta(2));
writematrix(Theta_hat,filename);
filename=sprintf('Psi_scBiglasso_gaussSynth_p_%d_n_%d_lambda1_%d_lambda2_%d.txt',p,n,beta(1),beta(2));
writematrix(Psi_hatF,filename);
filename=sprintf('Theta_scBiglasso_gaussSynth_p_%d_n_%d_lambda1_%d_lambda2_%d.txt',p,n,beta(1),beta(2));
writematrix(Theta_hatF,filename);
filename=sprintf('Psi_Terralasso_gaussSynth_p_%d_n_%d_lambda1_%d_lambda2_%d.txt',p,n,beta(1),beta(2));
writematrix(Psi_hatT,filename);
filename=sprintf('Theta_Terralasso_gaussSynth_p_%d_n_%d_lambda1_%d_lambda2_%d.txt',p,n,beta(1),beta(2));
writematrix(Theta_hatT,filename);

[precision_Psi_hat(i),precision_Theta_hat(i),recall_Psi_hat(i),recall_Theta_hat(i),TPR_Psi_hat(i),TPR_Theta_hat(i),FPR_Psi_hat(i),FPR_Theta_hat(i),accuracy_Psi_hat(i),accuracy_Theta_hat(i)]=ROC_syntheticGaussian(Psi_true,Psi_hat,Theta_true,Theta_hat);
[precision_Psi_hatF(i),precision_Theta_hatF(i),recall_Psi_hatF(i),recall_Theta_hatF(i),TPR_Psi_hatF(i),TPR_Theta_hatF(i),FPR_Psi_hatF(i),FPR_Theta_hatF(i),accuracy_Psi_hatF(i),accuracy_Theta_hatF(i)]=ROC_syntheticGaussian(Psi_true,Psi_hatF,Theta_true,Theta_hatF);
[precision_Psi_hatT(i),precision_Theta_hatT(i),recall_Psi_hatT(i),recall_Theta_hatT(i),TPR_Psi_hatT(i),TPR_Theta_hatT(i),FPR_Psi_hatT(i),FPR_Theta_hatT(i),accuracy_Psi_hatT(i),accuracy_Theta_hatT(i)]=ROC_syntheticGaussian(Psi_true,Psi_hatT,Theta_true,Theta_hatT);
end
fileout=sprintf('et_biglasso.txt');
filename =fullfile(outputdir, fileout);
writematrix(et_biglasso,filename);
fileout=sprintf('et_scBiglasso.txt');
filename =fullfile(outputdir, fileout);
writematrix(et_scBiglasso,filename);
fileout=sprintf('et_TeraLasso.txt');
filename =fullfile(outputdir, fileout);
writematrix(et_TeraLasso,filename);
fileout=sprintf('accuracy_Psi_biglasso_synthetic_gaussian.txt');
filename =fullfile(outputdir, fileout);
writematrix(accuracy_Psi_hat,filename);
fileout=sprintf('accuracy_Psi_scBiglasso_synthetic_gaussian.txt');
filename =fullfile(outputdir, fileout);
writematrix(accuracy_Psi_hatF,filename);
fileout=sprintf('accuracy_Psi_TeraLasso_synthetic_gaussian.txt');
filename =fullfile(outputdir, fileout);
writematrix(accuracy_Psi_hatT,filename);
fileout=sprintf('accuracy_Theta_biglasso_synthetic_gaussian.txt');
filename =fullfile(outputdir, fileout);
writematrix(accuracy_Theta_hat,filename);
fileout=sprintf('accuracy_Theta_scBiglasso_synthetic_gaussian.txt');
filename =fullfile(outputdir, fileout);
writematrix(accuracy_Theta_hatF,filename);
fileout=sprintf('accuracy_Theta_TeraLasso_synthetic_gaussian.txt');
filename =fullfile(outputdir, fileout);
writematrix(accuracy_Theta_hatT,filename);
newcolors = [0.83 0.14 0.14
             1.00 0.54 0.00
             0.47 0.25 0.80
             0.25 0.80 0.54];
         
colororder(newcolors)

figure(9)
plot(accuracy_Psi_hat(:),et_biglasso(:),'-*','LineWidth',1.5)
hold on 
plot(accuracy_Psi_hatF(:),et_scBiglasso(:),'--o','LineWidth',1.5)
plot(accuracy_Psi_hatT,et_TeraLasso(:),'-.x','LineWidth',1.5);
xlabel('Accuracy')
ylabel('time','FontSize',11)
legend('biglasso','scBiglasso','Teralasso','Location','southwest','FontSize',11)
hold off

figure(10)
subplot(121)
plot(20:20:100,et_biglasso(:),'-*','LineWidth',1.5)
hold on
plot(20:20:100,et_scBiglasso(:),'--o','Linewidth',1.5)
hold off
xlabel('n=p')
ylabel('time','FontSize',11)
legend('biglasso','scBiglasso','Location','northwest','FontSize',11)
subplot(122)
plot(20:20:100,accuracy_Psi_hat(:),'-*','LineWidth',1.5)
hold on 
plot(20:20:100,accuracy_Psi_hatF(:),'--o','LineWidth',1.5)
hold off
xlabel('n=p')
ylabel('time','FontSize',11)
legend('biglasso','scBiglasso','Location','southeast','FontSize',11)

figure(11)
plot(20:20:100,et_biglasso(:),'-*','LineWidth',1.5)
hold on
plot(20:20:100,et_scBiglasso(:),'--o','LineWidth',1.5)
hold off
xlabel('n=p')
ylabel('time','FontSize',11)
legend('biglasso','scBiglasso','Location','northwest','FontSize',11)
saveas(gcf,'synthetic_gaussian_timing.pdf')