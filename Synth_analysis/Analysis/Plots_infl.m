%% Sythetic section: cell experiments
%In this script, we present the cell experiments for synthetic data.
clc
clear all
close all

addpath('ViolinPlot/');%Add path of violin plot function

DataFolder='Results/';
TypeOfExp='Inflation/Beta/';%Type of experiment %Real/
FullDataFolder=append(DataFolder,TypeOfExp);%Full data folder
PlotFolder=append('Plots/',TypeOfExp);%Folder with plots

PlotName_cor="Infl_exp_cor.pdf";%Correlation plot name
PlotName_mn="Infl_exp_lt_mn.pdf";%Latent mean plot name
PlotName_sig="Infl_exp_lt_sig.pdf";%Latent sigma plot name

Num_data=5;%Number of datasets
lb=[0.05,0.10,0.25,0.50,0.75];
Num_genes=300;%Number of features
mse_cor_all=[];
mse_mn_all=[];
mse_sig_all=[];
for ii=1:Num_data
File_Data=sprintf(append(FullDataFolder,'Synth_MT_infl_res_%d.mat'),ii);%Data filename
load(File_Data,"mse_cor","mse_m","mse_sig")%Load data generating parameters

mse_cor_all=[mse_cor_all,mse_cor];
mse_mn_all=[mse_mn_all,mse_m];
mse_sig_all=[mse_sig_all,mse_sig];
end
x_neg={'5%','10%','25%','50%','75%'};%

aa=25;
% Draw plot
figure(1)
violin(mse_cor_all)
hold on
xticklabels(x_neg)
xticks([1,2,3,4,5])%
set(gca,'FontSize',30)
xlabel('Inflation')
ylabel('Difference from truth')
set(gca,'FontSize',aa)
title('Correlation','FontSize',30)

% Save plot
name=fullfile(PlotFolder,PlotName_cor);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(2)
subplot(1,2,1)
violin(mse_sig_all(:,[1,3,5,7,9]))
xticklabels(x_neg)
xticks([1,2,3,4,5])
xlabel('Inflation')
ylabel('Difference from truth')
set(gca,'FontSize',aa)
title('Methylation std.','FontSize',30)
subplot(1,2,2)
violin(mse_sig_all(:,[2,4,6,8,10]))
xticklabels(x_neg)
xticks([1,2,3,4,5])
xlabel('Inflation')
ylabel('Difference from truth')
set(gca,'FontSize',aa)
title('Expression std.','FontSize',30)

name=fullfile(PlotFolder,PlotName_sig);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(3)
subplot(1,2,1)
violin(mse_mn_all(:,[1,3,5,7,9]))
xticklabels(x_neg)
xticks([1,2,3,4,5])
xlabel('Inflation')
ylabel('Difference from truth')
set(gca,'FontSize',aa)
title('Methylation lt. mean','FontSize',30)
subplot(1,2,2)
violin(mse_mn_all(:,[2,4,6,8,10]))
xticklabels(x_neg)
xticks([1,2,3,4,5])
xlabel('Inflation')
ylabel('Difference from truth')
set(gca,'FontSize',aa)
title('Expression lt. mean','FontSize',30)

name=fullfile(PlotFolder,PlotName_mn);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all