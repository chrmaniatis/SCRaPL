%Buillds covariance matrix for each genomic region from gene specific
%standard deviations and correlations. 
%Warning: Contrary to other palces input needs to be constrained (ie
%positive or in (-1,1)
%sig: standard deviations for each genomic region 
%rho: correlation between epigenetic marks for each genomic region.
function [Psi] = Build_Mtrx(sig,rho)
Psi=zeros(2,2,size(sig,1));
for ii=1:size(sig,1)
    Psi(1,1,ii)=sig(ii,1)^2;
    Psi(2,2,ii)=sig(ii,2)^2;
    Psi(1,2,ii)=rho(ii)*sig(ii,1)*sig(ii,2);
    Psi(2,1,ii)=rho(ii)*sig(ii,1)*sig(ii,2);
    
end

end

