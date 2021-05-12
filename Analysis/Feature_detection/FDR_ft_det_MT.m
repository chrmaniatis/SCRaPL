%Function used to compute signifincant features with pearson correlation.
%This function is designed for
%methylation/accesibility -expression data. (ie binomial-count data)
function [ft_summary,ft_det,fdr_s, cor_prs] = FDR_ft_det_MT(y_dat,nrm,fdr)

nrm=nrm.norm_fact;
y_dat(:,3)=y_dat(:,3)./nrm(:,1);
    
[cor_prs,p] = prs_ft_cor_MT(y_dat(:,4),y_dat);
fdr_s = mafdr(p,'BHFDR',true',"Showplot",false);

ft_det=fdr_s<fdr;
ft_summary=sum(ft_det);

end

