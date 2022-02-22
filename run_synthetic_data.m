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
n=40;
p=50;
[S,T,Y,Psi_true,Theta_true,Omega_true]=Simulate_data(n,p);
beta1=0.005:0.001:0.016;
beta2=0.005:0.001:0.016;
nl=length(beta1);
for i=1:nl
    for j=1:nl
beta=[beta1(i),beta2(j)];
[diffTheta,diffPsi,Psi,Theta,objectiveFunction]=UpdateLoop1_only_backtracking_per_column(S,T,beta);

figure(1), clf
subplot(121), imagesc(Psi_true), title('\Psi_0')
colorbar % Visualise precisions.
subplot(122), imagesc(Psi), title('\Psi')
colorbar

figure(2), clf
subplot(121), imagesc(Theta_true), title('\Theta_{true}')
colorbar % Visualise precisions. 
subplot(122), imagesc(Theta), title('\Theta')
colorbar

fig=figure(3); clf
plot(objectiveFunction(:,1))
title('Overall obj fun')

figure(4), clf
plot(objectiveFunction(:,2))
title('Obj Fun for \Psi')

figure(5), clf
plot(objectiveFunction(:,3))
title('Obj Fun for \Theta')

%binary transformation: only negative values are counted as connection
Psi_trueb=(Psi_true<0);
Psi_b=(Psi<0);
TP_Psi=sum(sum((Psi_trueb==Psi_b).*Psi_trueb));%only offdiagonals are counted
FN_Psi=sum(sum((Psi_trueb~=Psi_b).*Psi_trueb));
FP_Psi=sum(sum((Psi_trueb~=Psi_b).*(1-Psi_trueb)));
TN_Psi=sum(sum((Psi_trueb==Psi_b).*(1-Psi_trueb)));
TPR_Psi=TP_Psi/(TP_Psi+FN_Psi);
FPR_Psi=FP_Psi/(FP_Psi+TN_Psi);
FNR_Psi=FN_Psi/(TP_Psi+FN_Psi);

Theta_trueb=(Theta_true<0);
Theta_b=(Theta<0);

TP_Theta=sum(sum((Theta_trueb==Theta_b).*Theta_trueb));
FN_Theta=sum(sum((Theta_trueb~=Theta_b).*Theta_trueb));
FP_Theta=sum(sum((Theta_trueb~=Theta_b).*(1-Theta_trueb)));
TN_Theta=sum(sum((Theta_trueb==Theta_b).*(1-Theta_trueb)));
TPR_Theta=TP_Theta/(TP_Theta+FN_Theta);
FPR_Theta=FP_Theta/(FP_Theta+TN_Theta);
FNR_Theta=FN_Theta/(TP_Theta+FN_Theta);

precision_Psi=TP_Psi/(TP_Psi+FP_Psi);
recall_Psi=TP_Psi/(TP_Psi+FN_Psi);
accuracy_Psi=(TP_Psi+TN_Psi)/(TP_Psi+TN_Psi+FP_Psi+FN_Psi);
precision_Theta=TP_Theta/(TP_Theta+FP_Theta);
recall_Theta=TP_Theta/(TP_Theta+FN_Theta);
accuracy_Theta=(TP_Theta+TN_Theta)/(TP_Theta+TN_Theta+FP_Theta+FN_Theta);

BETA1(i,j)=beta1(i);
BETA2(i,j)=beta2(j);
PRECISION_Psi(i,j)=precision_Psi;
PRECISION_Theta(i,j)=precision_Theta;
RECALL_Psi(i,j)=recall_Psi;
RECALL_Theta(i,j)=recall_Theta;
FPR_PSI(i,j)=FPR_Psi;
FPR_THETA(i,j)=FPR_Theta;
FNR_PSI(i,j)=FNR_Psi;
FNR_THETA(i,j)=FNR_Theta;
TPR_PSI(i,j)=TPR_Psi;
TPR_THETA(i,j)=TPR_Theta;
ACCURACY_Psi(i,j)=accuracy_Psi;
ACCURACY_Theta(i,j)=accuracy_Theta;
log_likelihood_Psi=-(objectiveFunction(size(objectiveFunction,1),2)-beta(1)*norm(Psi));
log_likelihood_Theta=-(objectiveFunction(size(objectiveFunction,1),3)-beta(2)*norm(Theta));
gamma=0.5;
[ebic_Psi]=Calculate_EBIC(log_likelihood_Psi,Psi_b,Theta_b,gamma);
[ebic_Theta]=Calculate_EBIC(log_likelihood_Theta,Theta_b,Psi_b,gamma);
[bic_Psi]=Calculate_BIC(log_likelihood_Psi,Psi_b,Theta_b);
[bic_Theta]=Calculate_BIC(log_likelihood_Theta,Theta_b,Psi_b);
%%%%%%%%%%%%%plottings
EBIC_Psi(i,j)=ebic_Psi;
EBIC_Theta(i,j)=ebic_Theta;
BIC_Psi(i,j)=bic_Psi;
BIC_Theta(i,j)=bic_Theta;
    end
end
fileout=sprintf('beta1_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(BETA1,filename);

fileout=sprintf('beta2_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(BETA2,filename);

fileout=sprintf('precision_Psi_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(PRECISION_Psi,filename);

fileout=sprintf('precision_Theta_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(PRECISION_Theta,filename);

fileout=sprintf('recall_Psi_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(RECALL_Psi,filename);

fileout=sprintf('recall_Theta_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(RECALL_Theta,filename);

fileout=sprintf('FPR_Psi_n_%d_p_%d.xlsx',n,p);
filename =fullfile(outputdir, fileout);
writematrix(FPR_PSI,filename);

fileout=sprintf('FPR_Theta_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(FPR_THETA,filename);

fileout=sprintf('TPR_Psi_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(TPR_PSI,filename);

fileout=sprintf('TPR_Theta_n_%d_p_%d.xlsx',n,p);
filename =fullfile(outputdir,fileout);
writematrix(TPR_THETA,filename);

fileout=sprintf('accuracy_Psi_n_%d_p_%d.xlsx',n,p);
filename = fullfile(outputdir,fileout);
writematrix(ACCURACY_Psi,filename);

fileout=sprintf('accuracy_Theta_n_%d_p_%d.xlsx',n,p);
filename =fullfile(outputdir, fileout);
writematrix(ACCURACY_Theta,filename);

figure(6), clf
subplot_tight(1,3,1,[0.1,0.03])
imagesc(Psi_true), title('\Psi_0')
colorbar % Visualise precisions. 
subplot_tight(1,3,2,[0.1,0.03])
imagesc(Theta_true), title('\Theta_0')
colorbar
subplot_tight(1,3,3,[0.1,0.03])
imagesc(Omega_true), title('\Psi_0 \oplus \Theta_0')
colorbar
saveas(gcf,'Psi0_Theta0_Omega0.pdf')

figure(7)
subplot_tight(2,2,1,[0.14,0.08])
hold on;
plot(RECALL_Psi(:,8),PRECISION_Psi(:,8),'Color',[0.00 0.28 0.73],'LineStyle','-','Linewidth',0.7,'Marker','*');
plot(RECALL_Psi(:,1),PRECISION_Psi(:,1),'Color','r','LineStyle','--','Linewidth',0.8,'Marker','o')%[0.83 0.13 0.18]
hold off;
title("(a)",'FontSize',11)
xlabel('Recall \Psi','FontSize',11)
ylabel('Precision \Psi','FontSize',11)
legend('\beta_2=0.012','\beta_2=0.005','Location','southwest','FontSize',11)

subplot_tight(2,2,2,[0.14,0.08])
hold on;
plot(RECALL_Theta(8,:),PRECISION_Theta(8,:),'Color',[0.00 0.28 0.73],'LineStyle','-','Linewidth',0.7,'Marker','*');
plot(RECALL_Theta(1,:),PRECISION_Theta(1,:),'Color','r','LineStyle','--','Linewidth',0.8,'Marker','o')%[0.83 0.13 0.18]
hold off;
title("(b)",'FontSize',11)
xlabel('Recall \Theta','FontSize',11)
ylabel('Precision \Theta','FontSize',11)
legend('\beta_1=0.012','\beta_1=0.005','Location','southwest','FontSize',11)

subplot_tight(2,2,3,[0.14,0.08])
hold on;
plot(BETA1(:,6),ACCURACY_Psi(:,6),'Color',[0.00 0.5 0],'LineStyle','-.','Linewidth',0.7,'Marker','x');
plot(BETA2(6,:),ACCURACY_Theta(6,:),'Color',[0.35 0 0.7],'LineStyle',':','Linewidth',0.8,'Marker','s');%[0.83 0.13 0.18]
ylim([0.8 1])
hold off;
title("(c)",'FontSize',11)
xlabel('\beta','FontSize',11)
ylabel('Accuracy','FontSize',11)
legend('Accuracy of \Psi-\beta_1','Accuracy of \Theta-\beta_2','Location','southeast','FontSize',10)

subplot_tight(2,2,4,[0.14,0.08])
hold on;
plot(FPR_PSI(:,6),TPR_PSI(:,6),'Color',[0.00 0.5 0],'LineStyle','-.','Linewidth',0.7,'Marker','x');
plot(FPR_THETA(6,:),TPR_THETA(6,:),'Color',[0.35 0 0.7],'LineStyle',':','Linewidth',0.8,'Marker','s');%[0.83 0.13 0.18]
hold off;
title("(d)",'FontSize',11)
xlabel('FPR','FontSize',11)
ylabel('TPR','FontSize',11)
legend('\Psi','\Theta','Location','southeast','FontSize',11)
saveas(gcf,'sythetic_curve_scB.pdf')

x=(5:1:16)*10.^-3;
x1=(5:1:16)*10.^-3;%(2:2:16)*10.^-3;
figure(8);
subplot_tight(2,2,1,[0.14,0.1])
plot(RECALL_Psi(:,6),PRECISION_Psi(:,6),'Color',[0.00 0.5 0],'LineStyle','-.','Linewidth',0.7,'Marker','x')
str=string(BETA1(3:6,6));
str='\beta_1'+str;
text(RECALL_Psi(3:6,6),PRECISION_Psi(3:6,6),str,'FontSize',7)
ylim([0.55,1.05])
title("(a)",'FontSize',11)
xlabel('Recall_{\Psi}','FontSize',10)
ylabel('Precision_{\Psi}','FontSize',11)

subplot_tight(2,2,2,[0.14,0.1])
plot(RECALL_Theta(6,:),PRECISION_Theta(6,:),'Color',[0.35 0 0.7],'LineStyle',':','Linewidth',0.8,'Marker','s');%[0.83 0.13 0.18]
str=string(BETA2(6,3:5));
str='\beta_2'+str;
text(RECALL_Theta(6,3:5),PRECISION_Theta(6,3:5)+0.03,str,'FontSize',7)
ylim([0.55,1.05])
xlabel('Recall_{\Theta}','FontSize',10)
ylabel('Precision_{\Theta}','FontSize',11)
title("(b)",'FontSize',11)

subplot_tight(2,2,3,[0.14,0.1])
plot(x,BIC_Psi(:,6),'Color',[0.00 0.5 0],'LineStyle','-.','Linewidth',0.8,'Marker','x');
title("(c)",'FontSize',11)
%xticks(x1);
%xtickangle(90)
xlabel('\beta_1','FontSize',10)
ylabel('BIC_{\Psi}','FontSize',11)

subplot_tight(2,2,4,[0.14,0.1])
plot(x,BIC_Theta(6,:),'Color',[0.35 0 0.7],'LineStyle',':','Linewidth',0.7,'Marker','s');
title("(d)",'FontSize',11)
%xticks(x1)
%xtickangle(90)
xlabel('\beta_2','FontSize',10)
ylabel('BIC_{\Theta}','FontSize',11)
saveas(gcf,'Precision_Recall_BIC.pdf')
