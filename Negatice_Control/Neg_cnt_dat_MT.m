%% Negative control dataset script
%This script takes a methylation/accessibility expression datset and
%outputs several negative control ones.
clc
clear all

Num_rep=6;%Number of negative control data at the end of the script.

for uuu=1:Num_rep
File_name_sv=sprintf('Meta_data/MT_enh_MT_25000_neg%d_new.mat',uuu);%save filename
File_name_inp=sprintf('Data/MT_enh_MT_25000_neg%d_new.mat',uuu);%input filename
load(File_name_inp)

y_true1=y_true;
norm_fact_all=norm_fact(y_true(:,5),:);

for ii=1:size(unique(y_true(:,4)),1)
    
    tt=y_true1(y_true(:,4)==ii,:);

    perm_met=randperm(size(tt,1));
    perm_rna=randperm(size(tt,1)); 

    y_true1(y_true1(:,4)==ii,1:2)=tt(perm_met,1:2);
    y_true1(y_true1(:,4)==ii,3)=tt(perm_rna,3);
    norm_fact_all(y_true1(:,4)==ii,1)=norm_fact(tt(perm_rna,5),1);
    
end
y_true=y_true1;
save(File_name_sv,'y_true','out_genes','norm_fact','gene_labels','cell_labels')
uuu/Num_rep%progress
end
