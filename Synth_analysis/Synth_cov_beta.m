%Synthetic data generating script for coverage experiments where correlation is sampled from beta distribution
clc
clear all
close all

N_cells=60;%Number of cells
N_genes=300;%Number of genes
m_true=[1,4];%Data generating mean

Sig=[2,0;0,3];%Data generating standard deviation
cor=ones(2,2,N_genes);
cor(1,2,:)=2*betarnd(15,15,N_genes,1)-1;%Data generating correlation
cor(2,1,:)=cor(1,2,:);
H=zeros(size(cor));
for ii=1:size(H,3)
    H(:,:,ii)=Sig*cor(:,:,ii)*Sig;
end

zero_inf=0.1*ones(N_genes,1);%Zero inflation

x_lat=[];
y_true=[];
for ii=1:N_genes
    Cov_CpG = randi([50 250],N_cells,1);%Generate coverage
    hh=mvnrnd(m_true,H(:,:,ii),N_cells);%Generate latent means
    x_lat=[x_lat;hh,ii*ones(N_cells,1),(1:N_cells)'];%Sampling latent states
    tt=poissrnd(exp(hh(:,2)));
    tt(rand(size(tt,1),1)<zero_inf(ii))=0;
    y_true=[y_true;binornd(Cov_CpG,normcdf(hh(:,1))),Cov_CpG,tt,ii*ones(N_cells,1),(1:N_cells)'];%Raw data
end
norm_fact=ones(N_cells,1);
norm_fact=table(norm_fact);%Normalization constant

filename=sprintf('Synth_MT_cov_%d.mat',round(mean(y_true(:,2))));
filename_full=sprintf('Synth_MT_full_cov_%d.mat',round(mean(y_true(:,2))));

save(filename,'y_true','norm_fact') 
save(filename_full)
clearvars -except y_true norm_fact
rr=accumarray(y_true(:,4),y_true(:,1)./y_true(:,2),[],@mean);
rr1=accumarray(y_true(:,4),y_true(:,3),[],@mean);
