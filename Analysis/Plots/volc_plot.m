%% Volcano Plots
%We use this script to draw bayesian and standard volcano plot given parameter gamma.
clc
clear all
close all


type_of_data='MT';%Type of data. Available options MT(binomial-count data), NM (binomial-binomial data)
Filename="NMT_MT_prom_2500_neg0_Meta";%Name of inference meta-file
[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Meta_files'];%Folder with data
Full_File = fullfile(DataFolder,append(Filename,".mat"));%Full file directory

SaveFolder="Volc";%Directory where output files will be saved
File_out_volc=append(Filename,"_volc.pdf");%Saving name of volcano plots

load(Full_File,"cor_gen_col","y_true")%Load files
max_fdr=0.1;%Maximum allowed FDR
if strcmp(type_of_data,'MT') == 1
    load(Full_File,"norm_fact")%Load normalisation constant
    [ft_summary_prs,ft_det_prs,FDR_null_prob,corr_prs,p_prs] = FDR_ft_det_MT(y_true,norm_fact,max_fdr);%Feature dection using pearson correlation
elseif strcmp(type_of_data,'NM')==1
    [ft_summary_prs,ft_det_prs,FDR_null_prob,corr_prs,p_prs] = FDR_ft_det_NM(y_true,max_fdr);%Feature dection using pearson correlation
    
else
    error('Please use a valid input for type of data. Accepted inputs are MT and NM.');
end

thrs=[0.205,0.205];%Choosen value for parameter gamma
prob=[0.7,0.9];%Values of alpha for grid-search
max_efdr=0.1;%Maximum allowed EFDR
[ft_summary,feature_det_scrapl,p_scrapl] =EFDR_ft_det(cor_gen_col,thrs,prob,max_efdr);%Grid search for alpha to achieve EFDR less than 10 %

%Choose alpha to push EFDR below 10%
[~,index]=max(ft_summary(:,end));
bay_indices=feature_det_scrapl{1,1}(:,index);

ttt=[corr_prs,p_prs];
ttt1=[median(tanh(cor_gen_col/2),2),1-p_scrapl];

%Distinguish features to categories of interest
rr=(bay_indices+ft_det_prs)/2;
EFDR_imp_FDR_imp=find(rr==1);
EFDR_imp_FDR_nimp=find(rr==0.5 & bay_indices==1);
EFDR_nimp_FDR_imp=find(rr==0.5 & ft_det_prs==1);
EFDR_rest=find(rr==0);

figure(1)
%Bayesian Volcano plot
subplot(1,2,1)
scatter(ttt1(EFDR_rest,1),-log10(ttt1(EFDR_rest,2)),150,'fill','c')
hold on
scatter(ttt1(EFDR_imp_FDR_imp,1),-log10(ttt1(EFDR_imp_FDR_imp,2)),150,'fill','o')
scatter(ttt1(EFDR_imp_FDR_nimp,1),-log10(ttt1(EFDR_imp_FDR_nimp,2)),150,'fill','g')
scatter(ttt1(EFDR_nimp_FDR_imp,1),-log10(ttt1(EFDR_nimp_FDR_imp,2)),150,'fill','k')
xlabel('Posterior media correlation')
ylabel('-log10(p_b_a_y)')
title('Bayesian Volc. Plot')
set(gca,'FontSize',20)
legend({'Not signif.','EFDR signif. FDR signif.','EFDR signif. FDR not signif.','EFDR not signif. FDR signif.'},'FontSize',26)

%Standard Volcano plot
subplot(1,2,2)
scatter(ttt(EFDR_rest,1),-log10(ttt(EFDR_rest,2)),150,'fill','c')
hold on
scatter(ttt(EFDR_imp_FDR_imp,1),-log10(ttt(EFDR_imp_FDR_imp,2)),150,'fill','r')
scatter(ttt(EFDR_imp_FDR_nimp,1),-log10(ttt(EFDR_imp_FDR_nimp,2)),150,'fill','g')
scatter(ttt(EFDR_nimp_FDR_imp,1),-log10(ttt(EFDR_nimp_FDR_imp,2)),150,'fill','k')
xlabel('Pearson correlation')
ylabel('-log10(p_p_r_s)')
title('Standard Volc. Plot')
set(gca,'FontSize',20)
legend({'Not signif.','EFDR signif. FDR signif.','EFDR signif. FDR not signif.','EFDR not signif. FDR signif.'},'FontSize',16)

name=fullfile(SaveFolder,File_out_volc);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all
