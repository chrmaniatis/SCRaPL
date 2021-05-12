%Computation of potential scale reduction factor to determine covergence
clear all
close all
clc

[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Meta_files/'];%Folder with data


aaa=3;%Thinning parameter
%Define tensors where posterior samples from different chains will go 
theta_cor=zeros(9480,1,450/aaa,5);
theta_zinf=zeros(9480,1,450/aaa,5);
theta_m=zeros(9480,2,450/aaa,5);
theta_sig=zeros(9480,2,450/aaa,5);

for ii=1:size(theta_cor,4)
    Filename=sprintf("gastrulation_2500_2500_neg0_inf_Meta_%d.mat",ii-1);%Data filename
    load(append(DataFolder,Filename),"cor_gen_col","m_gen_col","sig_gen_col","zero_inf_col")%Load data
    
    theta_cor(:,:,:,ii)=tanh(cor_gen_col(:,1:aaa:end)/2);
    theta_zinf(:,:,:,ii)=1./(1+exp(-zero_inf_col(:,1:aaa:end)));
    theta_m(:,:,:,ii)=m_gen_col(:,:,1:aaa:end);
    theta_sig(:,:,:,ii)=exp(sig_gen_col(:,:,1:aaa:end));
end

R_cor=PSRF(theta_cor);%scale reduction for correlation
R_zinf=PSRF(theta_zinf);%scale reduction for zero inflation
R_m=PSRF(theta_m);%scale reduction for mean
R_sig=PSRF(theta_sig);%scale reduction for standard deviation

save("Convergence_Files/gastrulation_conv.mat","R_cor","R_zinf","R_m","R_sig")