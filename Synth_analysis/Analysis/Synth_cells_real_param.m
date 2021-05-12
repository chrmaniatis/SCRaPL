%% Convergence to real parameters 
%In this script we investigate chain's convergence to real parameters in
%synthetic data
clc
clear all
close all

[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Real/'];%Folder with data

TypeOfExp='Downsample_cells/';%Type of experiment
MetaDataFolder=append(DataFolder,TypeOfExp,'Meta/Beta/');%Meta data folder
RawDataFolder=append(DataFolder,TypeOfExp,'Data/Beta/');%Raw data folder
ResultsFolder=append('Results/Real/',TypeOfExp,'Beta/');%Results folder

Num_data=6;%Number of datasets
lb=[5,10,25,50,100,200];
for ii=1:Num_data
File_Meta=sprintf(append(MetaDataFolder,'Synth_real_cells_%d_beta_Meta.mat'),lb(ii));%Meta filename
File_Data=sprintf(append(RawDataFolder,'Synth_real_cells_%d_beta.mat'),lb(ii));%Raw data filename
File_Results=sprintf(append(ResultsFolder,'Synth_real_cells_res_%d'),ii);%Save filenaem


load(File_Meta,"cor_gen_col","m_gen_col","sig_gen_col","Sima_gene","y_true")%Load meta files
load(File_Data,"cor_true","x_lt")%Load true mean and covariance 

m_true=[accumarray(x_lt(:,3),x_lt(:,1),[],@mean),accumarray(x_lt(:,3),x_lt(:,2),[],@mean)];
sig_true=[accumarray(x_lt(:,3),x_lt(:,1),[],@std),accumarray(x_lt(:,3),x_lt(:,2),[],@std)];
crr_true=cor_true';
[crr_prs,p_t] = prs_ft_cor_MT(y_true(:,4),y_true);%Estimate per feature correlation

aa=3;
crr_samp=tanh(cor_gen_col/2);%Posterior correlation samples
crr_samp=crr_samp(:,1:aa:end);%Improving latent correlation autocorrelation

sig_samp=exp(sig_gen_col);%Posterior variance samples
sig_samp=sig_samp(:,:,1:aa:end);%Improving latent mean autocorrelation

m_samp=m_gen_col(:,:,1:aa:end);%Improving latent mean autocorrelation

ign_ind=mean(crr_samp==0,2)>0.8;%Indices of ignored genes during inference; 

mse_cor=mean(crr_samp-crr_true,2);%Estimate difference of posterior mean from true correlation for each feature
mse_m=mean(m_samp-m_true,3);%Estimate difference of posterior mean from true mean for each feature
mse_sig=mean(sig_samp-sig_true,3);%Estimate difference of posterior mean from true standard deviation for each feature

if sum(ign_ind>0)
    uu=mse_cor;
    uu(ign_ind)=[];
    md=median(uu);
    mse_cor(ign_ind)=md;
    
    uu=mse_m;
    uu(ign_ind,:)=[];
    md=median(uu,1);
    mse_m(ign_ind,1)=md(1);
    mse_m(ign_ind,2)=md(2);
    
    uu=mse_sig;
    uu(ign_ind,:)=[];
    md=median(uu,1);
    mse_sig(ign_ind,1)=md(1);
    mse_sig(ign_ind,2)=md(2);
end

Num_lags=30;%Maximum number of lags for auto correlation
auto_corr=zeros(size(crr_samp,1),Num_lags+1);
 
for jj=1:size(crr_samp,1)
    auto_corr(jj,:)=autocorr(crr_samp(jj,:),'NumLags',Num_lags) ;
end

save(File_Results,"mse_cor","mse_m","mse_sig","auto_corr","crr_samp","crr_prs","sig_samp","m_samp")%save workspace
end
