function [bic]=Calculate_BIC(log_likelihood,Psi_b,Theta_b)
n=size(Theta_b,1);
p=size(Psi_b,1);
Psi_b=double(Psi_b);
%Theta_b=double(Theta_b);
num_edge=round(sum(sum(Psi_b-diag(diag(Psi_b))))/2);
%binomial_efficient=nchoosek((p*(p+1))/2,num_edge);
bic=-2*log_likelihood+log(n)*num_edge;%+gamma*log(binomial_efficient);
end