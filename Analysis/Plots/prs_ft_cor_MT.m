%Function used to estimate pearson correlation and p-value for each feature
%in methylation/accessibility-expression data (ie binomail-count data)
function [corr_all,p_t] = prs_ft_cor_MT(gns_ind,y_true)

n=histc(y_true(:,4),unique(y_true(:,4)));% number of observations per feature
x1=y_true(:,1)./y_true(:,2);%methylation/accessibility level 
x2=y_true(:,3);

ff_m=[accumarray(gns_ind,x1,[],@mean),accumarray(gns_ind,x2,[],@mean)];%per feature methylation/accessibility and expression mean
ff_s=[accumarray(gns_ind,x1,[],@std),accumarray(gns_ind,x2,[],@std)];%per feature methylation/accessibility and expression standard deviation

y_true1=[x1,x2];
y_true1=(y_true1-ff_m(gns_ind,:))./ff_s(gns_ind,:);
y_true1=y_true1(:,1).*y_true1(:,2);
corr_all=accumarray(gns_ind,y_true1,[],@sum)./(histc(gns_ind,unique(gns_ind))-1);%estimate per feature pearson correlation

t=corr_all.*sqrt((n-2)./(1-corr_all.^2));
p_t=2*tcdf(-abs(t),n-2);%estimate p-value
end

