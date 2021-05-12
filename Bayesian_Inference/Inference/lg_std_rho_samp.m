%Function builds conditional log posterior and the corresponding derivative
%of the covariance given its markov blanket for each genomic region. 

%param: vector with standard deviation for each epigenetic component and
%correlation (all unctostrained)
%SS: Is a vector with sufficient statistics for the standard deviations and correlation.
%a,b: parameters of the standard deviation priors
%alpha,beta: parameters of the correlation prior
%N_obs: number of observations used to build log-posterior

function [lg,glg] = lg_std_rho_samp(param,SS,a,b,alpha,beta,N_obs)
sig1=exp(param(1));
sig2=exp(param(2));
rr=tanh(param(3)/2);

lg=-(a(1)+1)*param(1)-(a(2)+1)*param(2)-b(1)/sig1-b(2)/sig2;%log-std prior
lg=lg+(alpha-1)*log(rr+1)+(beta-1)*log(1-rr);%log cor prior
lg=lg-(0.5/(1-rr^2))*(SS(1)/sig1^2+SS(2)/sig2^2-2*rr*SS(3)/(sig1*sig2));%llk pt.1
lg=lg-N_obs*(param(1)+param(2)+0.5*log(1-rr^2));%llk pt. 2


glg=zeros(size(param,2),1);
glg(1)=-N_obs/sig1-(-SS(1)/sig1+rr*SS(3)/sig2)/(sig1^2*(1-rr^2))-(a(1)+1)/sig1+b(1)/sig1^2;%scaled log-posterior derivative wrt. first standard deviation
glg(2)=-N_obs/sig2-(-SS(2)/sig2+rr*SS(3)/sig1)/(sig2^2*(1-rr^2))-(a(2)+1)/sig2+b(2)/sig2^2;%scaled log-posterior derivative wrt. second standard deviation
glg(3)=(-rr*(SS(1)/sig1^2+SS(2)/sig2^2-2*rr*SS(3)/(sig1*sig2))/(1-rr^2)+SS(3)/(sig1*sig2)+N_obs*rr)/(1-rr^2)+(alpha-1)/(1+rr)-(beta-1)/(1-rr); %scaled log-posterior derivative wrt. correlation
glg=glg.*[sig1;sig2;1/(1+cosh(param(3)))];%scaling gradient
end

