%% Raw data diagnostics
%We use this script to plot raw data diagnostics 
clc
clear all
close all

type_of_data='MT';%Type of data. Available options MT(binomial-count data), NM (binomial-binomial data)
Filename="NMT_MT_prom_2500_neg0_Meta";%Name of inference meta-file
[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Meta_files'];%Folder with data
Full_File = fullfile(DataFolder,append(Filename,".mat"));%Full file directory

SaveFolder="Raw_Diagnostics";%Directory where output files will be saved
File_out_diag1=append(Filename,"_diagn1.pdf");%Saving name of diagnostic plots 1
File_out_diag2=append(Filename,"_diagn2.pdf");%Saving name of diagnostic plots 2
File_out_diag3=append(Filename,"_diagn3.pdf");%Saving name of diagnostic plots 3
File_out_diag4=append(Filename,"_diagn4.pdf");%Saving name of diagnostic plots 4
File_out_diag5=append(Filename,"_diagn5.pdf");%Saving name of diagnostic plots 5
File_out_diag6=append(Filename,"_diagn6.pdf");%Saving name of diagnostic plots 6

load(Full_File,"y_true")%Load Files

Epi_comp_name_1="meth.";%Epigenetic component 1 name
Epi_comp_name_2="expr.";%Epigenetic component 2 name


if strcmp(type_of_data,'MT') == 1
    
    load(Full_File,"norm_fact")%Load normalisation constant from inference meta-file
     zz=norm_fact.norm_fact;
     y_true(:,3)=y_true(:,3)./zz(y_true(:,5));%Normalise expression data
     
     mn_raw=[accumarray(y_true(:,4),y_true(:,1)./y_true(:,2),[],@mean),accumarray(y_true(:,4),y_true(:,3),[],@mean)];%Methylation/Accessiility -Expression means
     st_raw=[accumarray(y_true(:,4),y_true(:,1)./y_true(:,2),[],@std),accumarray(y_true(:,4),y_true(:,3),[],@std)];%Methylation/Accessiility -Expression std
     cor_raw = prs_ft_cor_MT(y_true(:,4),y_true);%Per feature pearson correlation
     sec_raw=[accumarray(y_true(:,4),y_true(:,2),[],@mean),accumarray(y_true(:,4),y_true(:,3)==0,[],@mean)];%Per feature methylation/accessibility coverage, %of zeros
    
elseif strcmp(type_of_data,'NM')==1
    
     mn_raw=[accumarray(y_true(:,5),y_true(:,1)./y_true(:,2),[],@mean),accumarray(y_true(:,4),y_true(:,3)./y_true(:,4),[],@mean)];%Methylation-Accessiility means
     st_raw=[accumarray(y_true(:,5),y_true(:,1)./y_true(:,2),[],@std),accumarray(y_true(:,4),y_true(:,3)./y_true(:,4),[],@std)];%Methylation-Accessiility standard deviations
     cor_raw = prs_ft_cor_NM(y_true(:,5),y_true);%Per feature pearson correlation
     sec_raw=[accumarray(y_true(:,5),y_true(:,2),[],@mean),accumarray(y_true(:,5),y_true(:,4),[],@mean)];%Per feature methylation,accessibility coverage
    
else
    error('Please use a valid input for type of data. Accepted inputs are MT and NM.');
end

%Raw data diagnostic plots
figure(1)
subplot(1,3,1)
scatter(log(mn_raw(:,1)),log(mn_raw(:,2)))
xlabel(append("Raw ",Epi_comp_name_1," mean (log)"))
ylabel(append("Raw ",Epi_comp_name_2," mean (log)"))
set(gca,'FontSize',20)
subplot(1,3,2)
scatter(log(mn_raw(:,1)),log(st_raw(:,1)))
xlabel(append("Raw ",Epi_comp_name_1," mean (log)"))
ylabel(append("Raw ",Epi_comp_name_1," std (log)"))
set(gca,'FontSize',20)
subplot(1,3,3)
scatter(log(mn_raw(:,2)),log(st_raw(:,2)))
xlabel(append("Raw ",Epi_comp_name_2," mean (log)"))
ylabel(append("Raw ",Epi_comp_name_2," std (log)"))
set(gca,'FontSize',20)

name=fullfile(SaveFolder,File_out_diag1);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(2)
subplot(1,3,1)
histogram(cor_raw,30)
xlabel("Prs. cor.")
ylabel("# of features")
set(gca,'FontSize',20)

if strcmp(type_of_data,'MT') == 1
    subplot(1,3,2)
    histogram(sec_raw(:,1),30) 
    xlabel("Methylation/accessibility coverage")
    ylabel("# of features")
    set(gca,'FontSize',20)
    subplot(1,3,3)
    histogram(sec_raw(:,2),20) 
    xlabel("% of zeros")
    ylabel("# of features")
    set(gca,'FontSize',20)
else
    subplot(1,3,2)
    histogram(sec_raw(:,1),30) 
    xlabel("Methylation coverage")
    ylabel("# of features")
    set(gca,'FontSize',20)
    subplot(1,3,3)
    histogram(sec_raw(:,2),30) 
    xlabel("Accessibility coverage")
    ylabel("# of features")
    set(gca,'FontSize',20)
end

name=fullfile(SaveFolder,File_out_diag2);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(3)
subplot(1,3,1)
histogram(histc(y_true(:,4),unique(y_true(:,4))),25)
xlabel("% of observations")
ylabel("# of features")
set(gca,'FontSize',20)

if strcmp(type_of_data,'MT') == 1
    subplot(1,3,2)
    scatter(log(mn_raw(:,1)),sec_raw(:,1))
    xlabel(append("Raw ",Epi_comp_name_1," mean (log)"))
    ylabel("Meth./Acc. coverage")
    set(gca,'FontSize',20)
    subplot(1,3,3)
    scatter(log(mn_raw(:,2)),sec_raw(:,2))
    xlabel("Raw expr. mean (log)")
    ylabel("% of zeros")
    set(gca,'FontSize',20)
else
    subplot(1,3,2)
    scatter(log(mn_raw(:,1)),sec_raw(:,1))
    xlabel(append("Raw mathylation mean (log)"))
    ylabel("Meth. coverage")
    set(gca,'FontSize',20)
    subplot(1,3,3)
    scatter(log(mn_raw(:,2)),sec_raw(:,2))
    xlabel(append("Raw accessibility mean (log)"))
    ylabel("Acc. coverage")
    set(gca,'FontSize',20)
end

name=fullfile(SaveFolder,File_out_diag3);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(4)
subplot(1,3,1)
scatter(log(st_raw(:,1)),log(st_raw(:,2)))
xlabel(append("Raw ",Epi_comp_name_1," std (log)"))
ylabel(append("Raw ",Epi_comp_name_2," std (log)"))
set(gca,'FontSize',20)
subplot(1,3,2)
scatter(log(st_raw(:,1)),cor_raw)
xlabel(append("Raw ",Epi_comp_name_1," std (log)"))
ylabel(append("Raw Prs. cor."))
set(gca,'FontSize',20)
subplot(1,3,3)
scatter(log(st_raw(:,2)),cor_raw)
xlabel(append("Raw ",Epi_comp_name_2," std (log)"))
ylabel(append("Raw Prs. cor."))
set(gca,'FontSize',20)

name=fullfile(SaveFolder,File_out_diag4);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(5)

subplot(1,3,1)
scatter(log(mn_raw(:,1)),cor_raw)
xlabel(append("Raw ",Epi_comp_name_1," mean (log)"))
ylabel(append("Raw Prs. cor."))
set(gca,'FontSize',20)
subplot(1,3,2)
scatter(log(mn_raw(:,2)),cor_raw)
xlabel(append("Raw ",Epi_comp_name_2," mean (log)"))
ylabel(append("Raw Prs. cor."))
set(gca,'FontSize',20)
subplot(1,3,3)
if strcmp(type_of_data,'MT') == 1
    scatter(sec_raw(:,1),cor_raw)
    xlabel(append("Meth./Acc. coverage "))
    xlabel(append("Raw Prs. cor"))
    
    set(gca,'FontSize',20)
else
    scatter(sec_raw(:,1),cor_raw)
    xlabel(append("Meth. coverage "))
    xlabel(append("Raw Prs. cor"))
    set(gca,'FontSize',20)
end

name=fullfile(SaveFolder,File_out_diag5);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(6)
if strcmp(type_of_data,'MT') == 1
    subplot(1,2,1)
    scatter(sec_raw(:,2),cor_raw)
    xlabel(append("% of zeros "))
    ylabel(append("Raw Prs. cor"))
    set(gca,'FontSize',20)
    subplot(1,2,2)
    scatter(sec_raw(:,1),sec_raw(:,2))
    xlabel(append("Meth./Acc. coverage"))
    ylabel(append("% of zeros "))
    set(gca,'FontSize',20)
else
    subplot(1,2,1)
    scatter(sec_raw(:,2),cor_raw)
    xlabel(append("Acc. coverage "))
    ylabel(append("Raw Prs. cor"))
    set(gca,'FontSize',20)
    subplot(1,2,2)
    scatter(sec_raw(:,1),sec_raw(:,2))
    xlabel(append("Meth. coverage"))
    ylabel(append("Acc. coverage"))
    set(gca,'FontSize',20)
end
name=fullfile(SaveFolder,File_out_diag6);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

