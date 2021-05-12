%% Go Analysis
%In this script we prepare files for Gene Ontology(GO) analysis. The output of the script are three file. The first one is the table of starting genes, which will serve
%as background in GO anlaysis. The other two files are features detected
%using Pearson correlation and the decision criterion as presented in
%SCRaPL paper
clc
clear all
close all

type_of_data='MT';%Type of data. Available options MT(binomial-count data), NM (binomial-binomial data)
Filename="NMT_MT_prom_2500_neg0_Meta";%Name of inference meta-file
[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Meta_files'];%Folder with data
Full_File = fullfile(DataFolder,append(Filename,".mat"));%Full file directory

SaveFolder="GO_files";%Directory where output files will be saved
File_out_SCRaPL_imp=append(Filename,"_scrapl_imp.txt");%Saving name of file with important features detected with our approach
File_out_PRS_imp=append(Filename,"_prs_imp.txt");%Saving name of file with important features detected usinf prearson
File_out_SCRaPL_all=append(Filename,"_all_txt");%Saving name of file with starting pool of features
load(Full_File,'gene_labels',"cor_gen_col","y_true")%Load inference meta-file

%Before computing computing/loading statistically significant features, we
%have to make sure that posterior correlation estimates do not correlate
%with number of CpGs in festures of interest for choosen window, otherwise results might be
%misleading.
M_musc=readtable("Musculus_meta.csv");%Load numebr of CpGs for all genomic regions
M_musc=M_musc(ismember(M_musc.ens_id,gene_labels.Gene),:);%Discard genoimc regions that do not appear in our dataset

cor=tanh(cor_gen_col/2);

cor_cpg_den=corrcoef(median(cor,2),M_musc.cpg_den);%Correlation between posterior correlation estimtes and CpG density (Number of Cpgs divided by genomic region length)
cor_cpg=corrcoef(median(cor,2),M_musc.cpg_sum);%Correlation between posterior correlation estimtes and number of CpGs

if abs(cor_cpg(1,2))>0.1 || abs(cor_cpg_den(1,2))>0.1
   sprintf("It seems that posterior correlation estimates are correlated with CpGs in each feature, hence results might be misleading") 
end

%Detecting features by fixing gamma and varying alpha such that EFDR is
%less than 10%/Pearson correlation. Alternatively these could be loaded if
%relevant computations have been done.
thrs=[0.255,0.255];%Choosen value for parameter gamma
prob=[0.7,0.9];%Values of alpha for grid-search
max_efdr=0.1;%Maximum allowed EFDR
max_fdr=0.1;%Maximum allowed FDR
[ft_summary,feature_det_scrapl] =EFDR_ft_det(cor_gen_col,thrs,prob,max_efdr);%Grid search for alpha to achieve EFDR less than 10 %

gam_ind=1;
ft_summary(:,end,:)=ft_summary(:,end,:)/1000;%Divide detected features by 1000 to make them to bring them to similar scale with other recordings
ft_summary=squeeze(ft_summary(:,:,gam_ind));
[~,a_ind]=max(ft_summary(:,end));

feature_det_scrapl=feature_det_scrapl{1,gam_ind};
feature_det_scrapl=feature_det_scrapl(:,a_ind);

if strcmp(type_of_data,'MT') == 1
    load(Full_File,"norm_fact_all")%Load also normalisation constant from inference meta-file
    [~,ft_det_prs] = FDR_ft_det_MT(y_true,norm_fact,max_fdr);%Feature dection using pearson correlation
elseif strcmp(type_of_data,'NM')==1
    [~,ft_det_prs] = FDR_ft_det_NM(y_true,max_fdr);%Feature dection using pearson correlation
else
    error('Please use a valid input for type of data. Accepted inputs are MT and NM.');
end


imp_scrapl_features=gene_labels(feature_det_scrapl>0.5,:);%Subset of features identified as statistically significant using scrapl.
imp_prs_features=gene_labels(ft_det_prs>0.5,:);%Subset of features identified as statistically significant using pearson correlation.

%Save tables with features
writetable(cell2table(unique(gene_labels.Gene)),fullfile(SaveFolder,File_out_SCRaPL_all))
writetable(cell2table(unique(imp_scrapl_features.Gene)),fullfile(SaveFolder,File_out_SCRaPL_imp))
writetable(cell2table(unique(imp_prs_features.Gene)),fullfile(SaveFolder,File_out_PRS_imp))
