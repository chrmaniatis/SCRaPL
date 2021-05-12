%% Comparison between simulted data generated with zero-inflated and inflated poisson (pst. pred. check)

clc
clear all
close all

[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=[parentdir,'/Meta_files/'];%Folder with data

%% Load Inflated model

name=sprintf('NMT_MT_enh_25000_neg0_Meta1_Inf.mat');%Inflated poisson meta data
load(append(DataFolder,name),"m_gen_col","norm_fact","sig_gen_col","cor_gen_col","zero_inf_col","out_genes","y_true")%Load file
%Compute per feature posterior latent mean/covariance to simulate data raw
%expression 
nrm=norm_fact.norm_fact;
aa=3;
cor_inf=mean(tanh(cor_gen_col(:,1:aa:end)/2),2);
m_inf=mean(m_gen_col(:,:,1:aa:end),3);
sig_inf=mean(exp(sig_gen_col(:,:,1:aa:end)),3);
z_inf=mean(1./(1+exp(-zero_inf_col(:,1:aa:end))),2);
N_genes=size(cor_inf,1);

x_lat=zeros(size(y_true,1),3);
x_lat(:,end)=y_true(:,4);
y_sim_inf=y_true;
y_sim_inf(:,[1,3])=0*y_sim_inf(:,[1,3]); %Keeps methylation coverage same as in raw data.

%Simulate synthetic data using zero-inflated Poisson
for ii=1:N_genes
   H=diag(sig_inf(ii,:))*[1,cor_inf(ii);cor_inf(ii),1]*diag(sig_inf(ii,:));
   
   
   tt_x=mvnrnd(m_inf(ii,:),H,out_genes(ii,2));
   x_lat(x_lat(:,3)==ii,1:2)=tt_x; 
   tt=[normcdf(tt_x(:,1)),exp(tt_x(:,2))];
    
   
   tt_y=[binornd(y_sim_inf(y_sim_inf(:,4)==ii,2),tt(:,1)),poissrnd(nrm(y_sim_inf(y_sim_inf(:,4)==ii,end)).*tt(:,2))];
   tt_y(rand(out_genes(ii,2),1)<z_inf(ii),2)=0;
   y_sim_inf(y_sim_inf(:,4)==ii,[1,3])=tt_y; 
    
   ii/N_genes
end

%Clear variables
clearvars H inf x_lat tt_y tt_x tt m_gen_col sig_gen_col cor_gen_col cor_inf m_inf sig_inf zero_inf

%% Load Non-Inflated model
 
name=sprintf('NMT_MT_enh_25000_neg0_Meta_Noinf.mat');%No inflation Poisson data
load(append(DataFolder,name),"m_gen_col","sig_gen_col","cor_gen_col")%Load data

%Compute per feature posterior latent mean/covariance to simulate data raw
%expression 
cor_noinf=mean(tanh(cor_gen_col(:,1:aa:end)/2),2);
m_noinf=mean(m_gen_col(:,:,1:aa:end),3);
sig_noinf=mean(exp(sig_gen_col(:,:,1:aa:end)),3);

x_lat=zeros(size(y_true,1),3);
x_lat(:,end)=y_true(:,4);
y_sim_noinf=y_true;
y_sim_noinf(:,[1,3])=0*y_sim_noinf(:,[1,3]);

%Simulate synthetic data using Poisson
for ii=1:N_genes
   H=diag(sig_noinf(ii,:))*[1,cor_noinf(ii);cor_noinf(ii),1]*diag(sig_noinf(ii,:));
   
   tt_x=mvnrnd(m_noinf(ii,:),H,out_genes(ii,2));
   x_lat(x_lat(:,3)==ii,1:2)=tt_x; 
   tt=[normcdf(tt_x(:,1)),exp(tt_x(:,2))];
    
   
   tt_y=[binornd(y_sim_noinf(y_sim_noinf(:,4)==ii,2),tt(:,1)),poissrnd(nrm(y_sim_noinf(y_sim_noinf(:,4)==ii,end)).*tt(:,2))];
   y_sim_noinf(y_sim_noinf(:,4)==ii,[1,3])=tt_y; 
    
   ii/N_genes
end

%Clear variables
clearvars H x_lat tt_y tt_x tt m_gen_col sig_gen_col cor_gen_col cor_noinf m_noinf sig_noinf
%% Comparison plots
close all
%In this section we replicate some of the plots in figure 2 of Splatter
%(Zappia et al.)
mn_true=[accumarray(y_true(:,4),log2(1+y_true(:,1)./y_true(:,2)),[],@mean),accumarray(y_true(:,4),log2(1+y_true(:,3)),[],@mean)];
var_true=[accumarray(y_true(:,4),log2(1+y_true(:,1)./y_true(:,2)),[],@var),accumarray(y_true(:,4),log2(1+y_true(:,3)),[],@var)];

mn_inf=[accumarray(y_sim_inf(:,4),log2(1+y_sim_inf(:,1)./y_sim_inf(:,2)),[],@mean),accumarray(y_sim_inf(:,4),log2(1+y_sim_inf(:,3)),[],@mean)];
var_inf=[accumarray(y_sim_inf(:,4),log2(1+y_sim_inf(:,1)./y_sim_inf(:,2)),[],@var),accumarray(y_sim_inf(:,4),log2(1+y_sim_inf(:,3)),[],@var)];

mn_noinf=[accumarray(y_sim_noinf(:,4),log2(1+y_sim_noinf(:,1)./y_sim_noinf(:,2)),[],@mean),accumarray(y_sim_noinf(:,4),log2(1+y_sim_noinf(:,3)),[],@mean)];
var_noinf=[accumarray(y_sim_noinf(:,4),log2(1+y_sim_noinf(:,1)./y_sim_noinf(:,2)),[],@var),accumarray(y_sim_noinf(:,4),log2(1+y_sim_noinf(:,3)),[],@var)];

zer_gen_true=accumarray(y_true(:,4),y_true(:,3)==0,[],@mean);
zer_gen_inf=accumarray(y_sim_inf(:,4),y_sim_inf(:,3)==0,[],@mean);
zer_gen_noinf=accumarray(y_sim_noinf(:,4),y_sim_noinf(:,3)==0,[],@mean);

zer_cell_true=accumarray(y_true(:,5),y_true(:,3)==0,[],@mean);
zer_cell_inf=accumarray(y_sim_inf(:,5),y_sim_inf(:,3)==0,[],@mean);
zer_cell_noinf=accumarray(y_sim_noinf(:,5),y_sim_noinf(:,3)==0,[],@mean);

g_true_gene=repmat("Real data",N_genes,1);
g_inf_gen=repmat("Inf. data",N_genes,1);
g_noinf_gene=repmat("No Inf. data",N_genes,1);

N_cells=size(norm_fact,1);
g_true_cell=repmat("Real data",N_cells,1);
g_inf_cell=repmat("Inf. data",N_cells,1);
g_noinf_cell=repmat("No Inf. data",N_cells,1);

g=[g_true_gene;g_inf_gen;g_noinf_gene];
g_cell=[g_true_cell;g_inf_cell;g_noinf_cell];

mn=[mn_true;mn_inf;mn_noinf];
var_all=[var_true;var_inf;var_noinf];
zer_gen=[zer_gen_true;zer_gen_inf;zer_gen_noinf];
zer_cell=[zer_cell_true;zer_cell_inf;zer_cell_noinf];

figure(1)
boxplot(mn(:,2),g)
ylabel('Mean log2(CPM+1)')
title('Distribution of mean expression')
set(gca,'FontSize',30)

name="figure1.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(2)
boxplot(mn(:,2)-mean(mn(1:N_genes,2)),g)
ylabel('Rank difference mean log2(CPM+1)')
title('Difference in mean expression')
set(gca,'FontSize',30)

name="figure2.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(3)
boxplot(var_all(:,2),g)
ylabel('Var log2(CPM+1)')
title('Distribution of variance')
set(gca,'FontSize',30)

name="figure3.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(4)
boxplot(var_all(:,2)-mean(var_all(1:N_genes,2)),g)
ylabel('Rank difference var log2(CPM+1)')
title('Difference in variance')
set(gca,'FontSize',30)

name="figure4.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(5)
scatter(mn_true(:,2),var_true(:,2),'filled')
hold on
scatter(mn_inf(:,2),var_inf(:,2),'filled')
scatter(mn_noinf(:,2),var_noinf(:,2),'filled')
xlabel('Mean log2(CPM+1)')
ylabel('Var log2(CPM+1)')
title('Mean-variance relationship')
set(gca,'FontSize',30)
legend({'True','Inf.','No Inf.'},'FontSize',20)

name="figure5.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(6)
boxplot(zer_gen,g)
ylabel('Percentage zeros per gene')
title('Distribution of zeros per gene')
set(gca,'FontSize',30)

name="figure6.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(7)
boxplot(zer_gen-mean(zer_gen(1:N_genes)),g)
ylabel('Rank difference percentage zeros')
title('Difference in zeros per gene')
set(gca,'FontSize',30)

name="figure7.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(8)
boxplot(zer_cell,g_cell)
ylabel('Percentage zeros per cell')
title('Distribution of zeros per cell')
set(gca,'FontSize',30)

name="figure8.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all

figure(9)
boxplot(zer_cell-mean(zer_gen(1:N_cells)),g_cell)
ylabel('Rank difference percentage zeros')
title('Difference in zeros per cell')
set(gca,'FontSize',30)

name="figure9.pdf";
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all