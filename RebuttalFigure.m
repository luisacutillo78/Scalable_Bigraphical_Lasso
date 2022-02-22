addpath(genpath('../teralasso'))  
addpath(genpath('./glmnet_matlab/'))                                        % Friedman's LASSO.
addpath(genpath('./L1General/')) 
addpath(genpath('./ndlutil/matlab'))
addpath(genpath('./rca/matlab'))
addpath(genpath('./teralasso-master'))
%Script to illustrate numerical convergence of TeraLasso algorithm and
%scBiglasso
clear variables

id1=run_compareTeraScB(0.01);
id2=run_compareTeraScB(0.02);
id3=run_compareTeraScB(0.03);

f1=figure(id1);
%ax1 = gca; % get handle to axes of figure
f2=figure(id2);
%ax2 = gca; % get handle to axes of figure
f3=figure(id3);
%ax3 = gca; % get handle to axes of figure
