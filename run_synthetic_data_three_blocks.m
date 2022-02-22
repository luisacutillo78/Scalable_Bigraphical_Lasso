mypath='./';
plotdir='./plots/';
outputdir='./Output';
fullfile(mypath,'ndlutil/matlab');
addpath(fullfile(mypath,'ndlutil/matlab'));
addpath(fullfile(mypath,'/rca/matlab'));
addpath(fullfile(mypath,'glmnet_matlab'));
addpath(fullfile(mypath,'L1General'));

if ~isfolder(plotdir)
    mkdir(plotdir)
end

if ~isfolder(outputdir)
    mkdir(outputdir)
end
n=50;
p=45;
[S,T,Y,Psi_true,Theta_true,Omega_true]=Simulate_data_theta3blocks(n,p);
beta1=0.002;
beta2=0.0008;
nl=length(beta1);
for i=1:nl
    for j=1:nl
beta=[beta1(i),beta2(j)];
[diffTheta,diffPsi,Psi,Theta,objectiveFunction]=UpdateLoop1_only_backtracking_per_column(S,T,beta);

figure(1), clf
subplot(131), imagesc(Psi_true), title('\Psi_0')
colorbar % Visualise precisions. 
subplot(132), imagesc(Theta_true), title('\Theta_0')
colorbar
subplot(133), imagesc(Omega_true), title('\Psi_0 \oplus \Theta_0')
colorbar

figure(2), clf
subplot(121), imagesc(Psi_true), title('\Psi_0')
colorbar % Visualise precisions.
subplot(122), imagesc(Psi), title('\Psi')
colorbar

figure(3), 
subplot_tight(1,2,1,[0.1,0.05])
imagesc(Theta_true), title('\Theta_{true}')
colorbar % Visualise precisions. 
subplot_tight(1,2,2,[0.1,0.05])
imagesc(Theta), title('\Theta')
colorbar
saveas(gcf,'synthetic_theta3blocks.pdf')

Omega=kron(Psi,eye(size(Theta,1)))+kron(eye(size(Psi,1)),Theta);

figure(4), clf
subplot(121), imagesc(Omega_true), title('\Omega_{true}')
colorbar % Visualise precisions. 
subplot(122), imagesc(Omega), title('\Omega')
colorbar

fig=figure(5); clf
plot(objectiveFunction(:,1))
title('Overall obj fun')
pngfile=sprintf('./plots/Obj_n_%d_p_%d_beta1%f_beta2%f_n.pdf',n,p,beta1(i),beta2(j));
print(fig, '-dpdf', pngfile);

figure(6), clf
plot(objectiveFunction(:,2))
title('Obj Fun for \Psi')

figure(7), clf
plot(objectiveFunction(:,3))
title('Obj Fun for \Theta')

    end
end
