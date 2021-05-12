%% Convergence to real parameters 
%In this script we investigate chain's convergence to real parameters in
%synthetic data
clc
clear all
close all

[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Real/'];%Folder with data

TypeOfExp='Inflation/';%Type of experiment
MetaDataFolder=append(DataFolder,TypeOfExp,'Meta/Cnst/');%Meta data folder
RawDataFolder=append(DataFolder,TypeOfExp,'Data/Cnst/');%Raw data folder
ResultsFolder=append('Results/Real/',TypeOfExp,'Cnst/');%Results folder

Num_data=5;%Number of datasets
lb=[0.05,0.10,0.25,0.50,0.75];
for ii=1:Num_data
File_Meta=sprintf(append(MetaDataFolder,'Synth_real_infl_%0.2f_cnst_Meta.mat'),lb(ii));%Meta filename
File_Data=sprintf(append(RawDataFolder,'Synth_real_infl_%0.2f_cnst.mat'),lb(ii));%Raw data filename
File_Results=sprintf(append(ResultsFolder,'Synth_real_infl_res_%d'),ii);%Save filenaem


load(File_Meta,"cor_gen_col","m_gen_col","sig_gen_col","Sima_gene","y_true")%Load meta files
load(File_Data,"cor_true","x_lt")%Load true mean and covariance 

m_true=[accumarray(x_lt(:,3),x_lt(:,1),[],@mean),accumarray(x_lt(:,3),x_lt(:,2),[],@mean)];
sig_true=[accumarray(x_lt(:,3),x_lt(:,1),[],@std),accumarray(x_lt(:,3),x_lt(:,2),[],@std)];
crr_true= cor_true';
[crr_prs,p_t] = prs_ft_cor_MT(y_true(:,4),y_true);%Estimate per feature correlation

aa=3;
crr_samp=tanh(cor_gen_col/2);%Posterior correlation samples
crr_samp=crr_samp(:,1:aa:end);%Improving latent correlation autocorrelation

sig_samp=exp(sig_gen_col);%Posterior variance samples
sig_samp=sig_samp(:,:,1:aa:end);%Improving latent mean autocorrelation

m_samp=m_gen_col(:,:,1:aa:end);%Improving latent mean autocorrelation

ign_ind=mean(crr_samp==0,2)>0.8;%Indices of ignored genes during inference; 

mse_cor=mean(crr_samp,2)-crr_true;%Estimate difference of posterior mean from true correlation for each feature
mse_m=mean(m_samp,3)-m_true;%Estimate difference of posterior mean from true mean for each feature
mse_sig=mean(sig_samp,3)-sig_true;%Estimate difference of posterior mean from true standard deviation for each feature
 
mse_cor(ign_ind,:)=[];
mse_m(ign_ind,:)=[];
mse_sig(ign_ind,:)=[];

Num_lags=30;%Maximum number of lags for auto correlation
auto_corr=zeros(size(crr_samp,1),Num_lags+1);
 
for jj=1:size(crr_samp,1)
    auto_corr(jj,:)=autocorr(crr_samp(jj,:),'NumLags',Num_lags) ;
end

save(File_Results,"mse_cor","mse_m","mse_sig","auto_corr","crr_samp","crr_prs","sig_samp","m_samp")%save workspace
end
