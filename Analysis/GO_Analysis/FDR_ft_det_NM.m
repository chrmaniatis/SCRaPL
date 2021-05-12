%Function used to compute signifincant features with pearson correlation.
%This function is designed for
%methylation - accesibility data. (ie binomial-binomial data)
function [ft_det,ft_summary,fdr_s,cor_prs,p] = FDR_ft_det_NM(y_dat,fdr)
    
[cor_prs,p] = prs_ft_cor_MT(y_dat(:,4),y_dat);
fdr_s = mafdr(p,'BHFDR',true',"Showplot",false);

ft_det=fdr_s<fdr;
ft_summary=sum(ft_det);

end

