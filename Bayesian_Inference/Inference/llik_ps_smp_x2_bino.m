%Function builds conditional log posterior and the corresponding derivative
%of latent methylation/accessibility given given its markov blanket for each genomic region. 

%param: posterior latent methylation/accessibility
%y1: observed methylation/accesibility
%m: latent mean sample for each epigenetic component
%sig_lt: latent standard deviation sample for each epigenetic component (unconstrained)
%cor_lt: latent correlation sample (unconstrained)
%x1: latent component number 1

%Likelihood used in case of Methylation-Accessibility Data 
function [llik,glik] = llik_ps_smp_x2_bino(param,y1,m,sig_lt,cor_lt,x1)
Psi=cov_gen(1,sig_lt,cor_lt);

m_p=m(2)+Psi(1,2)./Psi(1,1).*(x1-m(1));
Psi_p=Psi(2,2)-Psi(1,2)^2/Psi(1,1);

x=param;

tt=normcdf(x);
%Build log posterior for latent methylation/accessibility
llik=y1(:,1).*log(tt)+(y1(:,2)-y1(:,1)).*log(1-tt);
llik=llik-0.5.*(x-m_p).^2/Psi_p;
llik=sum(llik);

tt1=normpdf(x);
%Build log posterior gradient for latent methylation/accessibility wrt.
%each latent methylation/expression variable
glik=(tt1./(1-tt)).*(y1(:,1)./tt-y1(:,2));
glik=glik-(x-m_p)/Psi_p;
end

