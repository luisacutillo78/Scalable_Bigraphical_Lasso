addpath(genpath('../glmnet_matlab/')) % Friedman's LASSO.
addpath(genpath('../L1General/'))     % Schmidt's L1GeneralProjection.


%read stages and genes names
S1=load('./CCdata/S1.txt');
S2=load('./CCdata/S2.txt');
S3=load('./CCdata/S3.txt');
%Genecc=load('../CCdata/')

fid = fopen('./CCdata/genesCCrow.txt');
All = textscan(fid, '%s');
fclose(fid);
All=All{:};

fid = fopen('./CCdata/Nmythosis.txt');
C = textscan(fid, '%s');
fclose(fid);
C=C{:};%this starts from position 1
Myth=C(1:end);

fid = fopen('./CCdata/NcellCortex.txt');
C = textscan(fid, '%s');
fclose(fid);
C=C{:};
Cor=C(2:end);

fid = fopen('./CCdata/Nchromatin.txt');
C = textscan(fid, '%s');
fclose(fid);
C=C{:};
Chro=C(2:end);

fid = fopen('./CCdata/Nendosome.txt');
C = textscan(fid, '%s');
fclose(fid);
C=C{:};
Endos=C(2:end);

fid = fopen('./CCdata/Nnucleoplasm.txt');
C = textscan(fid, '%s');
fclose(fid);
C=C{:};
Nucl=C(2:end);

fid = fopen('./CCdata/NreplicationFork.txt');
C = textscan(fid, '%s');
fclose(fid);
C=C{:};
Fork=C(2:end);

%%%%%%%%SELECTING THE MYTHOSYS GENES
[common,iA,is]=intersect(All,Myth);%intersect(All,subset);

%%%%%%Look at all the Phases and just consider our subset example
[ncells,ngenes]=size(S1);
scaleg=1;
scalec=1;
subindg=iA;
subindc=(1:scalec:ncells);


Ng=length(subindg);
M=length(subindc);
BS1=S1(:,iA);
BS2=S2(:,iA);
BS3=S3(:,iA);
B=[BS1;BS2;BS3];
C=standardizeCols(B);
R=standardizeCols(B');
S =2*sin(pi/6*corr(C,'Type','Spearman'));%Genes correlation
T = 2*sin(pi/6*corr(R,'Type','Spearman'));%cells correlation

T=nearestSPD(T);%cells correlation
S=nearestSPD(S);%Genes correlation
lambda=[0.014,0.001];

tol = 1e-2;
maxIter = 50;
[diffTheta,diffPsi,Psi,Theta,objectiveFunction]=scBiglasso(S,T,lambda,'maxIter',maxIter,'thresh',tol);

%%%%%%%%%%%%%plottings
figure(1)
Psi_b=abs(Psi)>0;
imagesc(Psi_b);
title('Psi binary not zeros val')
colorbar

figure(2)
Theta_b=abs(Theta)>0;
imagesc(Theta_b);
title('Theta binary not zero val')
colorbar

figure(3)
Psi_n=Psi<0;
imagesc(Psi_n);
title('Psi binary neg val')
colorbar

figure(4)
Theta_n=Theta<0;
imagesc(Theta_n);
title('Theta binary neg val')
colorbar

figure(5)
subplot_tight(1,2,1,[0.1,0.05])
imagesc(Psi_n), title('\Psi')
colorbar % Visualise precisions.
subplot_tight(1,2,2,[0.1,0.05])
imagesc(Theta_n), title('\Theta')
colorbar
saveas(gcf,'gene_Psi_Theta_neg.pdf')

writematrix(Psi,'./Output/Psi.csv') 
writematrix(Psi_n,'./Output/Psi_n.csv') 
writematrix(Psi_b,'./Output/Psi_b.csv') 
writematrix(Theta,'./Output/Theta.csv') 
writematrix(Theta_n,'./Output/Theta_n.csv') 
writematrix(Theta_b,'./Output/Theta_b.csv') 
