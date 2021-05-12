%Function used to estimate Gelman-Rubin criterion.
%theta:4-D Tensor with multiple chains of posterior parameter estimates
%First dimension: Number of features
%Second dimension: Parameter dimensionality
%Third dimension: Number of samples per chain
%Fourth dimension: Number of chains
function [R_c] = PSRF(theta)
M=size(theta,4);%Num of chains
T=size(theta,3);%Length of chains

theta_j=mean(theta,3);
sig_j_sq=var(theta,0,3);

B=T*var(theta_j,0,4);%Between chain variance
W=mean(sig_j_sq,4);%Within chain variance
V_hat=(1-1/T)*W+(M+1)*B/(M*T);%pooled variance

d_hat=2*V_hat./var(V_hat,0,1);%degrees of freedom

nrm=(d_hat+3)./(d_hat+1);%Constant that accounts for degrees of freedom
R_c=sqrt(nrm.*(V_hat./W));%potential scale reduction factor
end

