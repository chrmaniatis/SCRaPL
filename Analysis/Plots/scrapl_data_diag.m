%% Pseudo-bulk and other diagnostics
%We use this script to demonstrate that our model can produce results similar to the outputs of pseudo-bulk apporaches and plot other diagnostics 
clc
clear all
close all

type_of_data='MT';%Type of data. Available options MT(binomial-count data), NM (binomial-binomial data)
Filename="NMT_MT_prom_2500_neg0_Meta";%Name of inference meta-file
[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Meta_files'];%Folder with data
Full_File = fullfile(DataFolder,append(Filename,".mat"));%Full file directory

SaveFolder="SCRaPL_Diagnostics";%Directory where output files will be saved
File_out_diag1=append(Filename,"_diagn1.pdf");%Saving name of diagnostic plots
File_out_diag2=append(Filename,"_diagn2.pdf");%Saving name of diagnostic plots
load(Full_File,"cor_gen_col","m_gen_col","sig_gen_col","y_true")%Load Files

Epi_comp_name_1="methylation";%Epigenetic component 1 name
Epi_comp_name_2="expression";%Epigenetic component 2 name

mm=mean(m_gen_col,3);%Posterior latent mean for each feature
ss=mean(exp(sig_gen_col),3);%Posterior latent standard devation for each feature
crr=mean(tanh(cor_gen_col/2),2);%Posterior mean correlation for each feature

if strcmp(Epi_comp_name_2,"expression")==1
     load(Full_File,"norm_fact_all")
     zz=norm_fact_all.norm_fact;
     y_true(:,3)=y_true(:,3)./zz(y_true(:,5));
     mm1=[accumarray(y_true(:,4),y_true(:,1)./y_true(:,2),[],@mean),accumarray(y_true(:,4),y_true(:,3),[],@mean)];
     crt_tr=corrcoef(y_true(:,1)./y_true(:,2),y_true(:,3));%General correlation
else
     mm1=[accumarray(y_true(:,5),y_true(:,1)./y_true(:,2),[],@mean),accumarray(y_true(:,5),y_true(:,3)./y_true(:,4),[],@mean)];
     crt_tr=corrcoef(y_true(:,1)./y_true(:,2),y_true(:,3)./y_true(:,4));%General correlation
end

crt_lt=corrcoef(mm(:,1),mm(:,2));%Correlation between latent means
crt_psb=corrcoef(mm1(:,1),mm1(:,2));%Pseudo-bulk correlation

sz=35;
figure(1)
subplot(1,3,1)
scatter(mm(:,1),mm(:,2),sz,'filled')
xlabel(append("Latent ",Epi_comp_name_1," mean"))
ylabel(append("Latent ",Epi_comp_name_2," mean"))
set(gca,'FontSize',19)
subplot(1,3,2)
scatter(mm(:,1),ss(:,1),sz,'filled')
xlabel(append("Latent ",Epi_comp_name_1," mean"))
ylabel(append("Latent ",Epi_comp_name_1," standard deviation"))
set(gca,'FontSize',19)
subplot(1,3,3)
scatter(mm(:,2),ss(:,2),sz,'filled')
xlabel(append("Latent ",Epi_comp_name_2," mean"))
ylabel(append("Latent ",Epi_comp_name_2," standard deviation"))
set(gca,'FontSize',19)

name=fullfile(SaveFolder,File_out_diag1);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(2)
subplot(1,3,1)
scatter(ss(:,1),ss(:,2),sz,'filled')
xlabel(append("Latent ",Epi_comp_name_1," standard deviation"))
ylabel(append("Latent ",Epi_comp_name_2," standard deviation"))
set(gca,'FontSize',19)
subplot(1,3,2)
scatter(ss(:,1),crr,sz,'filled')
xlabel(append("Latent ",Epi_comp_name_1," standard deviation"))
ylabel("Latent posterior correlation")
set(gca,'FontSize',19)
subplot(1,3,3)
scatter(ss(:,2),crr,sz,'filled')
xlabel(append("Latent ",Epi_comp_name_2," standard deviation"))
ylabel("Latent posterior correlation")
set(gca,'FontSize',19)

name=fullfile(SaveFolder,File_out_diag2);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all
