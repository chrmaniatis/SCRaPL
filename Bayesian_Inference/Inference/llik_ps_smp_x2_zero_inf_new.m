%Function builds conditional log posterior and the corresponding derivative
%of latent expression (assuming inflation) given given its markov blanket for each genomic region. 

%param: posterior latent expression
%y2: observed expression
%m: latent mean sample for each epigenetic component
%sig_lt: latent standard deviation sample for each epigenetic component (unconstrained)
%cor_lt: latent correlation sample (unconstrained)
%x1: latent methylation/accessibility sample
%infl: latent inflation sample (unconstrained)
%sc: cell specific normalization

%Likelihood used in cases of Methylation-Expression and Accessibility-Expression Data 
function [llik,glik] = llik_ps_smp_x2_zero_inf_new(param,y2,m,sig_lt,cor_lt,x1,infl,sc)
Psi=cov_gen(1,sig_lt,cor_lt);
infl=1/(1+exp(-infl));    
    
m_p=m(2)+Psi(1,2)./Psi(1,1).*(x1-m(1));
Psi_p=Psi(2,2)-Psi(1,2)^2/Psi(1,1);

x=param;
%Build log posterior for latent expression
llik=zeros(size(x,1),1);
llik(y2==0)=log(infl+(1-infl)*exp(-sc(y2==0).*exp(x(y2==0))));
llik(y2>0)=-sc(y2>0).*exp(x(y2>0));
llik=llik-0.5.*(x-m_p-Psi_p*y2).^2/(Psi_p);
llik=sum(llik);

%Build log posterior gradient for latent expression wrt. 
%each latent expression variable
glik=zeros(size(x,1),1);
glik(y2==0)=-(1-infl)./(1-infl+infl*exp(sc(y2==0).*exp(x(y2==0)))).*sc(y2==0).*exp(x(y2==0));
glik(y2>0)=-sc(y2>0).*exp(x(y2>0));
glik=glik-(x-m_p-Psi_p*y2)/(Psi_p);
end
