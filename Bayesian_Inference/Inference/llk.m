%Builds likelihood for each genomic region from sampled data in the models that deals with methylation/expression or accessibility/expression 
%data assuming drop-out in expression data modeling. Since SCRaPL is a model with 1 hierarchy, we define the likelihood using
% two parts. One that corresponds to observation model and the other one to
% the first hierarchy. The likelihood is later used for tasks like model comparison.
%y_true: raw data matrix 
%x_lat: sampled latent states 
%m_gen: sampled latent means for each genomic region
%sig_gen: sampled latent standard deviations for each genomic region (log-transformed)
%cor_gen: sampled latent correlations for each genomic region (inv-tanh transformed)
%zero_inf: sampled zero-inflation for each genomic region (logit-transformed)
%norm_fact_all: cell specific normalization constant
function [lg_x_1,lg_x_2] = llk(y_true,x_lat,m_gen,sig_gen,cor_gen,zero_inf,norm_fact_all)
zinf=1./(1+exp(-zero_inf));
cr=tanh(cor_gen/2);
sg=exp(sig_gen);
m=m_gen;

lg_x_1=build_llk_inf(y_true,x_lat,zinf(x_lat(:,3)),norm_fact_all);
lg_x_2=build_lik_nrm(x_lat,m,sg,cr);

lg_x_1=[accumarray(lg_x_1(:,3),lg_x_1(:,1),[],@sum),accumarray(lg_x_1(:,3),lg_x_1(:,2),[],@sum)];
lg_x_2=accumarray(lg_x_2(:,2),lg_x_2(:,1),[],@sum);
end

function [llk] = build_lik_nrm(x,m,s,crr)
ind=squeeze(x(:,end,1));
m=m(ind,:);
s=s(ind,:);
crr=crr(ind,:);

z=(x(:,1,:)-m(:,1)).^2./s(:,1).^2 + (x(:,2)-m(:,2)).^2./s(:,2).^2 - 2*crr.*(x(:,1)-m(:,1)).*(x(:,2)-m(:,2))./(s(:,1).*s(:,2));
z=z./(2*(1-crr.^2));
llk=[-log(2*pi)-sum(log(s),2)-log(1-crr.^2)/2 - z,x(:,end,:)];
end

function [llk] = build_llk_inf(y,x,infl,nrm)
llk=x;
llk(:,1:2,:)=0;

infl=infl(y(:,4));

aa=0.001;
rt=[aa+(1-aa)*normcdf(x(:,1,:)),nrm.norm_fact(y(:,5)).*exp(x(:,2,:))];

llk(y(:,3)==0,2,:)=log(infl(y(:,3)==0,:,:)+(1-infl(y(:,3)==0,:,:)).*exp(-rt(y(:,3)==0,2,:)));
llk(y(:,3)>0,2,:)=log(1-infl(y(:,3)>0,1,:))+y(y(:,3)>0,3).*log(rt(y(:,3)>0,2,:))-rt(y(:,3)>0,2,:) -gammaln(y(y(:,3)>0,3)+1);

llk(:,1,:)=y(:,1).*log(rt(:,1,:))+(y(:,2)-y(:,1)).*log(1-rt(:,1,:));

end
