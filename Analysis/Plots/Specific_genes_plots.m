%% Specific gene plots
%In this script we plot posterior correlation distribution and compare it
%with normalised raw data.
clc
clear all
close all

SaveFolder="Specific_genes_plots";%Save folder name

DataFolder="~/Desktop/SCRaPL/Analysis/Meta_files";%Data folder
DataFile="MT_enh_MT_25000_neg0_Meta.mat";%Data file name
load(fullfile(DataFolder,DataFile),"y_true","cor_gen_col","norm_fact_all","gene_labels")%Load file

gene_ind=[8390;7868;189;61;8406;9100;8528;7539;6671;8261];%Indices of genes of interest.
cor_biog_pros=tanh(cor_gen_col/2);%Posterior correlation for each gene.
aa=1;%Parameter for improving autocorrelation.

for ii=1:size(gene_ind,1)
yy=y_true(y_true(:,4)==gene_ind(ii),:);
nrm=norm_fact_all(y_true(:,4)==gene_ind(ii),:);
yy(:,3)=yy(:,3)./nrm.norm_fact;
figure(ii)
subplot(1,2,1)
scatter(yy(:,1)./yy(:,2),log(1+yy(:,3)),sz,yy(:,2),'filled')
xlabel('Raw Methylation')
ylabel('Log expression')
hold on 
colorbar
set(gca,'FontSize',35)
sfh2 = subplot(1,2,2);
sfh2.Position = sfh2.Position - [-0.13 0 0.07 0];
violin(cor(gene_ind(ii),:)')
ylabel('Posterior correlation')
xlabel(sprintf('Feature %d',ii))
set(gca,'FontSize',33,'XTick',[])

name=append(SaveFolder,"/",sprintf("Feature_%d_simple.pdf",gene_ind(ii)));
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all
end
