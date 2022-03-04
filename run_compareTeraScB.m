function id=run_compareTeraScB(beta)
type = 'MCP';%Type of regularization
rng(33);
opt = 1;
K = 4;
big = 1;
errNrm = 'fro';

ps = 40*ones(1,K);
N = 1;
a = 60;

tol = 1e-6;
lk = 1;
psdum = ps(1)^opt*ones(1,2);

for jj=1:length(beta)
    
lambdaScB=[beta(jj) beta(jj)];
lambdatera=lambdaScB;
range = .1+[.1 .3];
mn = 0;
for k = 1:K
    k_init = floor(ps(k));
    Psi0{k} = generate_cov(ps(k),k_init,range);
    mn = mn + mean(diag(Psi0{k}));
end

%Find eigenvalues

eigz = 0;
for k = 1:K
    [U{k},D] = eig(Psi0{k});
    if min(diag(D)) <= 0
        error('Non PSD Sigma');
    end
    eigz = eigz + kron(ones(prod(ps(1:k-1)),1),kron(diag(D),ones(prod(ps(k+1:end)),1)));
end
for k = 1:K
    S{k} = 0;
end
%Generate data using eigdecomp
v = randn(prod(ps),N)./(sqrt(eigz)*ones(1,N));
for k = 1:K
    for j = 1:prod(ps(1:k-1))
        for i = 1:prod(ps(k+1:end))
            
            ix = prod(ps(k:end))*(j-1) + (i + (0:ps(k)-1)*prod(ps(k+1:end)));
            
            v(ix,:) = U{k}*v(ix,:);
        end
    end
    
end

for n = 1:N
    
    
    X = permute(reshape(v(:,n),ps(end:-1:1)),K:-1:1);
    
    
    %Produce sample covariances
    for k = 1:K
        XX = permute(X,[1:k-1,k+1:K, k]);
        mat = reshape(XX,prod(ps)/ps(k),ps(k));
        S{k} = S{k} + (mat'*mat)/(prod(ps)/ps(k));
    end
    
end
for k = 1:K
    S{k} = S{k}/N;
    
end


maxiter = 200;


K = 2;
ps = ps(1)^opt*ones(1,K);

lambda = lambdatera; 
%lambda=lambda/2;


if opt == 1
    v = v(1:prod(ps),:);
end

%% Get S
for k = 1:K
    S{k} = 0;
end
for n = 1:N
    
    
    X = permute(reshape(v(:,n),ps(end:-1:1)),K:-1:1);
    
    %Produce sample covariances
    for k = 1:K
        XX = permute(X,[1:k-1,k+1:K, k]);
        mat = reshape(XX,prod(ps)/ps(k),ps(k));
        S{k} = S{k} + (mat'*mat)/(prod(ps)/ps(k));
    end
    
end
for k = 1:K
    S{k} = S{k}/N;
    
end

%% Estimate
tic

[ Psi,Psis ] = teralasso( S,ps,type,a,tol ,lambda,maxiter);

time2 = toc;
stat_err2 = 0;

for k = 1:K
    if opt == 2
        mm = kron(Psi0{(k-1)*2 + 1},eye(size(Psi0{(k-1)*2 + 2})))+ kron(eye(size(Psi0{(k-1)*2 + 1})),Psi0{(k-1)*2 + 2});
    else
        mm = Psi0{k};
    end
    diff = (mm - Psi{k}) - eye(ps(k))*mean(diag((mm - Psi{k})));
    
    if errNrm ~= inf
        stat_err2 = stat_err2 + prod(ps)/ps(k)*norm(diff,errNrm)^2;
    else %infinity norm
        stat_err2 = max(stat_err2, max(max(abs(diff))));
    end
end
for k = 1:K
    Psis{1}{k} = eye(ps(k));
end
for count = 1:length(Psis)
    nrm(count) = 0;
    for k = 1:K
        diff = (Psis{count}{k} - Psi{k}) - eye(ps(k))*mean(diag((Psis{count}{k} - Psi{k})));
        
        if errNrm ~= inf
            nrm(count) = nrm(count) + prod(ps)/ps(k)*norm(diff,errNrm)^2;
        else %infinity norm
            nrm(count) = max(nrm(count), max(max(abs(diff))));
        end
    end
end
if opt == 1
    nrm = nrm*100;
end

%% Plot

ix = find(nrm./nrm(1)< 1e-4,1,'first');
if N == 1
    ix = find(nrm./nrm(1)< .2e-2,1,'first');
end
ix = find(nrm < stat_err2,1,'first')+2;
id=figure;hold off;semilogy(sqrt(nrm(1:ix-1)./nrm(1)),'*-','LineWidth',1.5,'Color',[1.00 0.54 0.00]);
%hold on;semilogy(1:length(nrm(1:ix-1)),sqrt(stat_err2/nrm(1))*ones(1,length(nrm(1:ix-1))),'g');
ixx = find(nrm < stat_err2,1,'first');
text(1, 1.2*sqrt(stat_err2/nrm(1)),[num2str(ixx/length(nrm)*time2,3) ' Tera sec']);
xlabel('Iterations');
ylabel('Norm. $\|\Omega_t - \Omega^*\|_F$','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda=lambdaScB;

tic
[diffTheta,diffPsi,Psi1,Theta,Psis]=scBiglasso_compare(S{2},S{1},lambda,'maxIter',maxiter);

time2 = toc;
Psi{1}=Theta;
Psi{2}=Psi1;
stat_err2_ScB = 0;
for k = 1:K
    diff = (Psi0{k} - Psi{k}) - eye(ps(k))*mean(diag((Psi0{k} - Psi{k})));
    if errNrm ~= inf
        stat_err2_ScB = stat_err2_ScB + prod(ps)/ps(k)*norm(diff,errNrm)^2;
    else %infinity norm
        stat_err2_ScB = max(stat_err2_ScB, max(max(abs(diff))));
    end
end
for k = 1:K
    Psis{1}{k} = eye(ps(k));
end
nrm=[];
for count = 1:length(Psis)
    nrm(count) = 0;
    for k = 1:K
       diff = (Psis{count}{k} - Psi{k}) - eye(ps(k))*mean(diag((Psis{count}{k} - Psi{k})));
        if errNrm ~= inf
            nrm(count) = nrm(count) + prod(ps)/ps(k)*norm(diff,errNrm)^2;
        else %infinity norm
            nrm(count) = max(nrm(count), max(max(abs(diff))));
        end
    end
end
if opt==1
    nrm=nrm*100;
end
%%%%%----------end here
ix = find(nrm./nrm(1)< 1e-4,1,'first');
if N == 1
    ix = find(nrm./nrm(1)< .2e-2,1,'first');
end
ix = find(nrm < stat_err2_ScB,1,'first')+2;
hold on;semilogy(sqrt(nrm(1:ix-1)./nrm(1)),'x-','LineWidth',1.5,'Color',[0.83 0.14 0.14]);
ixx = find(nrm < stat_err2_ScB,1,'first');

text(1, 1.2*sqrt(stat_err2_ScB/nrm(1)),[num2str(ixx/length(nrm)*time2,3) ' ScB sec']);
xlabel('Iterations','FontSize',12);
ylabel('Norm. $\|\Omega_t - \Omega^*\|_F$','interpreter','latex','FontSize',12)
legend({'Tera K = 2, d_k = 40', 'ScB K = 2, d_k = 40'},'Location','southeast','FontSize',12);

title(sprintf('$\\beta_1=\\beta_2$ = %1.2f',lambdaScB(1)),'Interpreter','latex','FontSize',12)
saveas(gcf,sprintf('compareTera_scB_Beta%1.4f.pdf',lambdaScB(1)))

end
end

