%Function builds conditional log posterior and the corresponding derivative
%of latent inflation given given its markov blanket for each genomic region. 

%param: posterior latent inflation (unconstrained)
%y1: observed expession counts
%x: latent expression  
%alpha,beta: prior inflation parameters

%Likelihood used in cases of Methylation-Expression and Accessibility-Expression Data 
function [llik,glik] = llik_ps_smp_zero_inf(param,y1,x,alpha,beta,sc)
%Build log posterior for latent inflation
llik=(alpha-1)*param-(alpha+beta-2)*log(1+exp(param));%prior
llik=llik+sum(log(exp(param)+exp(-sc(y1==0).*exp(x(y1==0)))))-size(y1,1)*log(1+exp(param));%evidence

%Build log posterior gradient for latent inflation
glik=sum(1./(1+exp(-sc(y1==0).*exp(x(y1==0)))*exp(-param)));
glik=glik+(alpha-1)-(alpha+beta-2+size(y1,1))./(1+exp(-param));
end



