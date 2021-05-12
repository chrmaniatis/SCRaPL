%% Feature detection
%This script is used detect features using SCRaPL's posterrior correlation
%against standard approach (feature detection using pearson correlation).
clc
clear all
close all

%parpool(32)%Initialises pool of workers. If parallel computing toolbox is not available please comment this line.
type_of_data='MT';%Type of data. Available options MT(binomial-count data), NM (binomial-binomial data)
[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Meta_files'];%Folder with data

for ii=1:6
Filename=sprintf("gastrulation_2500_2500_neg%d_inf_strp.mat",ii-1);%Meta-file name.
Full_Filename_Inference = fullfile(DataFolder,Filename);%Full file directory

Filename_Meta = sprintf("Detection_Meta_files/gastrulation_2500_2500_neg%d_inf_strp_Meta.mat",ii-1);%Saving name in traget
load(Full_Filename_Inference,"cor_gen_col","y_true")%Load inference meta-file

thrs=[0.1,0.6];%gamma values for grod search
prob=[0.51,0.95];% alpha values for gird search
max_efdr=0.1;%Maximum EFDR
max_fdr=max_efdr;%%Maximum FDR
[ft_summary_efdr,ft_det_efdr,EFDR_null_prob,alpha_gamma_ind] = EFDR_ft_det(cor_gen_col,thrs,prob,max_efdr);%Feature detection using grid search and efdr
if strcmp(type_of_data,'MT') == 1
    
    load(Full_Filename_Inference,"norm_fact_all")%Load also normalisation constant from inference meta-file
    [ft_summary_prs,ft_det_prs,FDR_null_prob] = FDR_ft_det_MT(y_true,norm_fact_all,max_fdr);%Feature dection using pearson correlation

    
elseif strcmp(type_of_data,'NM')==1
    
    [ft_summary_prs,ft_det_prs,FDR_null_prob] = FDR_ft_det_NM(y_true,max_fdr);%Feature dection using pearson correlation
    
else
    error('Please use a valid input for type of data. Accepted inputs are MT and NM.');
end

save(Filename_Meta,"ft_summary_efdr","ft_summary_prs","ft_det_efdr","ft_det_prs","EFDR_null_prob","FDR_null_prob")%save features detected with alternative approaches
end
