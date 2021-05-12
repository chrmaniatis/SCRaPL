%% Negative control dataset script
%This script takes a methylation accessibility datset and
%outputs several negative control ones.
clc
clear all

Num_rep=6;%Number of negative control data at the end of the script.

for uuu=1:Num_rep
File_name_sv=sprintf('Meta_data/MT_enh_MT_25000_0.75_neg%d_new.mat',uuu);%save filename
File_name_inp=sprintf('Data/MT_enh_MT_25000_0.75_neg%d_new.mat',uuu);%input filename
load(File_name_inp)

y_true1=y_true;
for ii=1:size(unique(y_true(:,5)),1)
    
    tt=y_true1(y_true(:,5)==ii,:);

    perm_met=randperm(size(tt,1));
    perm_acc=randperm(size(tt,1));  

    y_true1(y_true(:,5)==ii,1:2)=tt(perm_met,1:2);
    y_true1(y_true(:,5)==ii,3:4)=tt(perm_acc,3:4);
    
end
y_true=y_true1;
save(File_name_sv,'y_true','out_genes','gene_labels','cell_labels')
uuu/Num_rep%progress
end
