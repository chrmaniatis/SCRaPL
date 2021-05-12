%% Connecting subsets
%If the option for splitting raw data was initially applied, use this
%script to join meta files.
clc
clear all

DataFolder_Meta="Results/";%Data folder.Pick relevant from Result folders of this file.
SaveFolder_Rec="Results/";%Saving folder. 

load(append(DataFolder_Meta,"gastrulation_prom_2500_2500_part1_Meta.mat"),"Col_smps","ref_gen")% Load number of collected samples.
clear var gene_labels y_true
Num_data=10;%NUmber of feature subsets.
N_genes=size(ref_gen,1);%Total number of features

m_gen_col_all=zeros(N_genes,2,Col_smps);
m_prior_all=zeros(N_genes,2);
sig_gen_col_all=zeros(N_genes,2,Col_smps);
cor_gen_col_all=zeros(N_genes,Col_smps);
zero_inf_all=zeros(N_genes,Col_smps);
x_lat_all=[];
y_true_all=[];
out_genes_all=[(1:N_genes)',zeros(N_genes,1)];
for ii=1:Num_data
     name=sprintf(append(DataFolder_Meta,"/gastrulation_prom_2500_2500_part%d_Meta.mat"),ii);
     load(name)%Load subset.
     mapping=ismember(ref_gen,gene_labels);%Find mapping.
     
     min_gn_ind=find(mapping,1,'first');
     max_gn_ind=find(mapping,1,'last');
     
     m_gen_col_all(mapping,:,:)=m_gen_col;
     out_genes_all(mapping,2)=out_genes(:,2);
     m_prior_all(mapping,:)=m_prior;
     sig_gen_col_all(mapping,:,:)=sig_gen_col;
     cor_gen_col_all(mapping,:)=cor_gen_col;
     zero_inf_all(mapping,:)=zero_inf_col;
    
     y_true(:,end-1)=1000*(ii-1)+y_true(:,end-1);
     x_lat_all=[x_lat_all;x_lat];
     y_true_all=[y_true_all;y_true];
    
end
[~,~,y_true_ref(:,end-1)]=unique(y_true_ref(:,end-1));
x_lat_all(:,end)=y_true_ref(:,end-1);

%rename variables
x_lat=x_lat_all;
m_gen_col=m_gen_col_all;
m_prior=m_prior_all;
sig_gen_col=sig_gen_col_all;
cor_gen_col=cor_gen_col_all;
out_genes=out_genes_all;
gene_labels=ref_gen;
y_true=y_true_ref;

save(append(SaveFolder_Rec,"gastrulation_2500_2500.mat"),"y_true","x_lat","m_gen_col","m_prior","sig_gen_col","cor_gen_col","out_genes","gene_labels","cell_labels","a_gene","alpha_gene","alpha_inf","b_gene","beta_gene","beta_inf","Col_smps","norm_fact","Psi_0")