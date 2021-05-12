%% Model Comparison Bayesian
%In this script we use Deviance Information Criterion to compare between
%Poisson and zero-inflated poisson, in order to demonstrate the need for zero inflation 
%for most genes. For comparison we also use Akaike Information Criterion.

%Inflated model meta-data inference
clc
clear all

[parentdir,~,~]=fileparts(pwd);
DataFolder=[parentdir,'/Reconstructed'];

File=strcat('gastrulation_2500_2500_neg0_inf_Meta_0.mat');
load(fullfile(DataFolder,File),"cor_gen_col","m_gen_col","sig_gen_col","zero_inf_col","y_true","x_mn","norm_fact_all","lg_x_1","lg_x_2","out_genes")%Load Data

%Posterior mean for each latent parameter.
cor=mean(cor_gen_col,2);
mn=mean(m_gen_col,3);
sig=mean(sig_gen_col,3);
zinf=mean(zero_inf_col,2);
x_mn(:,3)=y_true(:,4);

%estimate mean deviance
[llik_mn1,llik_mn2]=llk_inf(y_true,x_mn,mn,sig,cor,zinf,norm_fact_all);
llk_mn=-2*(sum(llik_mn1,2)+llik_mn2);
llk=-2*(squeeze(sum(lg_x_1,2))+lg_x_2);

%compute DIC and AIC to compare models with inflation 
p_d=mean(llk,2)-llk_mn;
DIC_sp=2*p_d+llk_mn;
AIC_inf=min(llk,[],2)+2*(2*out_genes(:,2)+6);

clearvars -except DIC_sp DIC_gl p_d llk llk_mn p_d p_v AIC_inf

%% Non-inflated model meta-data inference
[parentdir,~,~]=fileparts(pwd);
DataFolder=[parentdir,'/Reconstructed'];

File=strcat('gastrulation_2500_2500_neg0_ninf_Meta_0.mat');
load(fullfile(DataFolder,File),"cor_gen_col","m_gen_col","sig_gen_col","y_true","x_mn","norm_fact_all","lg_x_1","lg_x_2","out_genes")%Load Data

%Posterior mean for each latent parameter.
cor=mean(cor_gen_col,2);
mn=mean(m_gen_col,3);
sig=mean(sig_gen_col,3);
x_mn(:,3)=y_true(:,4);

%estimate mean deviance
[llik_mn1_ninf,llik_mn2_ninf]=llk_ninf(y_true,x_mn,mn,sig,cor,norm_fact_all);
llk_mn_noinf=-2*(sum(llik_mn1_ninf,2)+llik_mn2_ninf);
llk_ninf=-2*(squeeze(sum(lg_x_1,2))+lg_x_2);

p_d_ninf=mean(llk_ninf,2)-llk_mn_noinf;
DIC_sp_ninf=2*p_d_ninf+llk_mn;
AIC_ninf=min(llk_ninf,[],2)+2*(2*out_genes(:,2)+5);

clearvars -except DIC_sp DIC_gl DIC_sp_ninf DIC_gl_ninf p_d llk llk_mn llk_ninf llk_mn_noinfp_d p_v p_d_ninf p_v_ninf AIC_inf AIC_ninf out_genes y_true

close all

aa=1;
uu=200;
D_sp=DIC_sp-DIC_sp_ninf;
D_aic=AIC_inf-AIC_ninf;

%Filter a small number of genes with very high DIC to make the plot more understandable. In the reporting these genes are included. 
D_aic=D_aic(abs(D_sp)<10000);
D_sp=D_sp(abs(D_sp)<10000);

%Plot figures
figure(1)
histogram(D_sp,uu)
hold on
histogram(D_aic,uu)
xline(0,'--r','LineWidth',1)
legend('DIC','AIC')
xlabel("Model Difference",'FontSize',25)
ylabel("Number of genes",'FontSize',25)
%%
name="Figures/M_selection_n";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all
