%Function used to compute signifincant features by using posterior corelations 
%for each feature and deciding criterion presented in SCRaPL paper.
%Decisions for significant features are done using posterior correlation
%samples in latent space. Essentially for combinateion of the proposed
%parameters alpha and gamma we estimate significant features (ie. features
%that have strong correlation). 
function [ft_summary,feature_ind_all,null_prob,ft_int] = EFDR_ft_det(cor,thresholds,prob,max_efdr)
gamma_lat=2*atanh(min(thresholds):0.005:max(thresholds));
alpha=min(prob):0.001:max(prob);%
ft_summary=zeros(size(alpha,2),5,size(gamma_lat,2));
null_prob=[];
parfor uuu=1:size(gamma_lat,2)
x_all1=-gamma_lat(uuu):0.005:gamma_lat(uuu);
f1_p=zeros(size(cor,1),size(x_all1,2));
p=zeros(size(cor,1),1);

for ii=1:size(cor,1)
   %Posterior correlation probability approximation for null hypothesis from samples and
   %integration for each feature
   y1 = ksdensity(cor(ii,:),x_all1,'Function','pdf');
   f1_p(ii,:)=y1;
   p(ii)=trapz(x_all1,y1);
end

p1=1-p;
EFDR=zeros(size(alpha,2),1);
EFNR=zeros(size(alpha,2),1);
feature=zeros(size(alpha,2),1);
feature_ind=zeros(size(cor,1),size(alpha,2));
null_prob=[null_prob,p1];

for jj=1:size(alpha,2)
    %Estimating EFDR and EFNR for each alpha and gamma combination
    zz1=p1>alpha(jj);
    zz2=p1<=alpha(jj);
    EFDR(jj)=sum(p.*zz1,1)./sum(zz1,1);
    EFNR(jj)=sum((1-p).*zz2,1)./sum(zz2,1);
    feature(jj,:)=sum(p1>alpha(jj),1);
    feature_ind(:,jj)=(p1>alpha(jj));
end

s=[tanh(gamma_lat(uuu)/2)*ones(size(alpha,2),1),alpha',EFDR,EFNR,feature];
s(s(:,3)>max_efdr,:)=0;%discarding comniations of alpha and gamma that yield efdr greater than the maximum threshold.
ft_summary(:,:,uuu)=s;
feature_ind_all{uuu}=sparse(feature_ind);%Converting logical arrays with significnat features to sparse arrays to save storage space.
uuu/size(gamma_lat,2)%progress
end
%Returning pair of alpha and gamma that yield maximum EFDR<10% and maximum
%number of features
ind_pair=zeros(2,size(gamma_lat,2));
[ind_pair(1,:),ind_pair(2,:)]=max(squeeze(ft_summary(:,end,:)),[],1);
[~,thrs_ind]=max(ind_pair(1,:));
alpha_ind=ind_pair(2,thrs_ind);

ft_int=[alpha_ind,thrs_ind];
end

