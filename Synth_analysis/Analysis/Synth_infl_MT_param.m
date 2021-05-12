%% Convergence to real parameters 
%In this script we investigate chain's convergence to real parameters in
%synthetic data
clc
clear all
close all

[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Inflation/'];%Folder with data

TypeOfExp='Cnst/';%Type of experiment
MetaDataFolder=append(DataFolder,'Meta/',TypeOfExp);%Meta data folder
RawDataFolder=append(DataFolder,'Full_Data/',TypeOfExp);%Raw data folder
ResultsFolder=append('Results/Inflation/',TypeOfExp);%Results folder

Num_data=5;%Number of datasets
lb=[0.05,0.10,0.25,0.50,0.75];
for ii=1:Num_data
File_Meta=sprintf(append(MetaDataFolder,'Synth_MT_infl_%0.2f_Meta.mat'),lb(ii));%Meta filename
File_Data=sprintf(append(RawDataFolder,'Synth_MT_full_infl_%0.2f.mat'),lb(ii));%Raw data filename
File_Results=sprintf(append(ResultsFolder,'Synth_MT_infl_res_%d'),ii);%Save filenaem


load(File_Meta,"cor_gen_col","m_gen_col","sig_gen_col","Sima_gene","y_true")%Load meta files
load(File_Data,"m_true","H","x_lat")%Load true mean and covariance 

sig_true=[sqrt(H(1,1,1)),sqrt(H(2,2,1))];
crr_true=squeeze(H(1,2,:)./sqrt(H(1,1,:).*H(2,2,:)));%cor_true;

[crr_prs,p_t] = prs_ft_cor_MT(y_true(:,4),y_true);%Estimate per feature correlation

aa=3;
crr_samp=tanh(cor_gen_col/2);%Posterior correlation samples
crr_samp=crr_samp(:,1:aa:end);%Improving latent correlation autocorrelation

sig_samp=exp(sig_gen_col);%Posterior variance samples
sig_samp=sig_samp(:,:,1:aa:end);%Improving latent mean autocorrelation

m_samp=m_gen_col(:,:,1:aa:end);%Improving latent mean autocorrelation

mse_cor=mean(crr_samp-crr_true,2);%Estimate difference of posterior mean from true correlation for each feature
mse_m=mean(m_samp-m_true,3);%Estimate difference of posterior mean from true mean for each feature
mse_sig=mean(sig_samp-sig_true,3);%Estimate difference of posterior mean from true standard deviation for each feature
 
Num_lags=30;%Maximum number of lags for auto correlation
auto_corr=zeros(size(crr_samp,1),Num_lags+1);
 
for jj=1:size(crr_samp,1)
    auto_corr(jj,:)=autocorr(crr_samp(jj,:),'NumLags',Num_lags) ;
end

save(File_Results,"mse_cor","mse_m","mse_sig","auto_corr","crr_samp","crr_prs","sig_samp","m_samp")%save workspace
end
