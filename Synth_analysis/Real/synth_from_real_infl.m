%% Synthetic from real data
%This script is used to generate synthetic data using latent expression
%data generated from variational auto-encoder (sc-VI). The ideal behind
%this experiment is to investigate the robustness of the model when latent
%state Gaussianity assumptions are dropped. For more information please
%look at SCRaPL paper
clc
clear all
close all

DataFolder="scVI_latent_expr";%Folder where latent expression data are placed
SaveFolder="Inflation/Data/Cnst";%Save folder
Filename="lt_exp_60.csv";%Latent expression rates filename.
SaveFilename="Synth_real_infl_0.75.mat";%Synthetic data filename.

T_lat=readtable(fullfile(DataFolder,Filename));%Load data
T_lat_y=T_lat(:,randperm(size(T_lat,2),300));%Permute columns to remove pattern created by the projections
x2=table2array(T_lat_y);%Convert latent expression data to array to make processing easier
x2=x2+0.001+0.06*rand(size(x2));%

clear T_lat T_lat_y
clc

cor_true=-0.8+0.2*rand(1,300);%2*betarnd(15,15,1,300)-1;%Generate random feature specific correlations that do not follow any standard distribution
x1=zeros(size(x2));
expre=zeros(size(x2));
cor_lt_obs=zeros(300,1);
scp=zeros(300,1);

for ii=1:size(cor_lt_obs,1)
    rr=log(x2(:,ii));%Latent expression
    m_meth=1+cor_true(ii)*(2/std(rr))*(rr-mean(rr));%conditional methylation mean
    sig_meth=sqrt(1-cor_true(ii)^2)*std(rr);%conditional methylation standard deviation
    x1(:,ii)=m_meth+sig_meth*randn(size(m_meth,1),1);%generate random methylation data with correlation cor_true 
    crr=corrcoef(x1(:,ii),log(x2(:,ii)));
    cor_lt_obs(ii)=crr(1,2);%make sure that correlation between methylation and expression is cor_true
    tt=rand(size(expre,1),1);
    scp(ii)=0.75;%Inflation probability
    expre(:,ii)=poissrnd(x2(:,ii));%Generate random poisson data
    expre(expre(:,ii)<quantile(expre(:,ii),0.95) & tt<scp(ii),ii)=0;%Censor them such that expression distribution follows zero-inflate Poisson.

end

norm_fact=ones(size(x1,1),2);
norm_fact=array2table(norm_fact);%Create nomralization constants (essentially ones file with ones) to mimic real data

dt=zeros(size(x1,1),size(x1,2),4);
dt(:,:,1)=x1;
dt(:,:,2)=expre;
dt(:,:,3)=repmat((1:size(x1,2)),size(x1,1),1);
dt(:,:,4)=repmat((1:size(x1,1))',1,size(x1,2));

%Tensor with latent data
dt_lt=zeros(size(x1,1),size(x1,2),4);
dt_lt(:,:,1)=x1;
dt_lt(:,:,2)=log(x2);
dt_lt(:,:,3)=repmat((1:size(x1,2)),size(x1,1),1);
dt_lt(:,:,4)=repmat((1:size(x1,1))',1,size(x1,2));

dt1=reshape(dt,size(dt,1)*size(dt,2),size(dt,3));
x_lt=reshape(dt_lt,size(dt_lt,1)*size(dt_lt,2),size(dt_lt,3));%matrix with latent methylation and expression data. (we save that to make comparisons with predictions)

Cap_CpG = randi([100 450],size(dt1,1),1);%Generate methylation coverage
y_true=[binornd(Cap_CpG,normcdf(dt1(:,1))),Cap_CpG,dt1(:,end-2:end)];%Generate methylation data

x_lt=sortrows(x_lt,3);
y_true=sortrows(y_true,4);

[a,~,c]=unique(y_true(:,4));
out_genes=[a,histc(c,a)];%Number of features per gene

save(strjoin([SaveFolder,SaveFilename],'/'),"out_genes","y_true","x_lt","cor_true","cor_lt_obs","norm_fact")