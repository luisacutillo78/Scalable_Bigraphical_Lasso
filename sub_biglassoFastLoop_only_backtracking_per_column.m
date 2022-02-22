function [Psi] = sub_biglassoFastLoop_only_backtracking_per_column(Psi,U,L1,L2,T,options)         % Main subroutine: Provides updates of the optimised parameter.
p = size(L2, 2);    
n = size(Psi, 2);
options.order = -1;         % -1: L-BFGS (limited-memory), 1: BFGS (full-memory), 2: Newton
options.verbose = 0;
lambda2=diag(L2);%%eigenvalues of Theta
%%Find n diagonal blocks pxp.  Depending on each partition, the W
%%seperate n diagonal pxp blocks 
%%W is (np)x(np), inverse of (Psi Kronecker-sum Theta)
%%in 4 blocks.[ W11,   t(W1non1) ]
%%            [ W1non1, Wnon1non1]
%%for each i update the i-th column/row of Psi
A=zeros(n-1,n-1);

for i = 1:n
   rIndx_11 = [1 : (i-1)     (i + 1) : n];%length n-1 %%row index not i
   rIndx_t12 = rIndx_11;% [1:(i-1)   (i+1):n];
   U11 = U( rIndx_11, : );                                 
   t12 = T( rIndx_t12, i );  
   %%p*T(not i, i)+A(not i,not i).*Psi(not i, i)=zeros(n-1,n-1)
   v=(Psi(i,i)+lambda2).^(-1);%% 1/(psi(1,1)+lambda(2,j))
   lambda1=diag(L1);%%eigenvalues of Psi
    
    recip_lambdas = 1./(lambda1+lambda2');
    temp_mat = recip_lambdas*v;
    for k=1:n-1 %%%because each index lenghth is n-1(for "not i")
        %%calculate the block trace of column k
        uk=U11(k,:);%%the dimension here follows sum i=1 to n in Proposition 3.2
        A(:,k)=(uk.*U11)*(temp_mat); %the equation in Proposition 3.2
    end
     A=(1/p)*A;%A*psi12+p*t12=0
     psi12=Psi(rIndx_t12,i);
     fit=glmnet(A,-t12,'gaussian',options);
     psi12new=sparse(fit.beta);
     %Obtain the direction vector 
     %moving from the old psi12 to the new update
     diff=psi12-psi12new;
     step=diff;
     grad=2*(A'*A)*psi12+2*A'*t12;
     zeta=(min(lambda1))^2;
     %Apply backtracking line search
     Q=norm(A*psi12+t12,'fro')^2+grad'*(-step)+(0.5/zeta)*norm(-step,'fro')^2;
     Fnew=norm(A*psi12new+t12,'fro')^2;%A*psi12+t12=0
     %Check condition like Tera Lasso. Conduct backtracking line search
     tk=1;
     while (Fnew>Q)
            diff=psi12-psi12new;
            tkk=(1+sqrt(1+4*tk^2))/2;
            step=tk/tkk*diff;
            psi12new=shrinkL1(psi12,step);
            Q=norm(A*psi12+t12,'fro')^2+grad'*(-step)+(0.5/zeta)*norm(-step,'fro')^2;
            Fnew=norm(A*psi12new+t12,'fro')^2;
     end
     Psi(rIndx_t12,i)=psi12new;
     Psi(i,rIndx_t12)=psi12new'; 
     [U,L1]=eig(Psi);
end
end
