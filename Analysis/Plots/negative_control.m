%% Negative control Plots
%In this script we draw barcharts to compare detected feature numbers in original and negative control data. 
clc
clear all

type_of_data='MT';%Type of data. Available options MT(binomial-count data), NM (binomial-binomial data)
[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DetectFolder=[parentdir,'/Feature_detection/Detection_Meta_files/'];%Folder with detected features on meta files

File_basic="MT_enh_MT_25000";%Filename

SaveFolder="Negative_control";%Directory where output files will be saved
File_out_diag1=append(File_basic,"_neg_cnt_1.pdf");%Saving name of diagnostic plots 1
File_out_diag2=append(File_basic,"_neg_cnt_2.pdf");%Saving name of diagnostic plots 2

num_data=6;
col_ft_num_bay=zeros(num_data,1);
col_ft_num_prs=zeros(num_data,1);

for ii=1:num_data

    File_Dete=strcat(sprintf(append(DetectFolder,File_basic,"_neg%d_Meta.mat"),ii-1));
    load(File_Dete,"ft_summary_efdr","ft_summary_prs")

     ft_summary_efdr=ft_summary_efdr(:,:,32);
     num_ft=max(ft_summary_efdr(:,end));

     col_ft_num_prs(ii)=ft_summary_prs;
     col_ft_num_bay(ii)=num_ft;
end

x_or=categorical({'orig.'});
x_neg=categorical({'neg1','neg2','neg3','neg4','neg5'});
x_or=reordercats(x_or,{'orig.'});
x_neg=reordercats(x_neg,{'neg1','neg2','neg3','neg4','neg5'});

y_bay=col_ft_num_bay;
y_prs=col_ft_num_prs;

figure(1)
bar(x_or,y_bay(1))
hold on
bar(x_neg,y_bay(2:end))
xlabel('Dataset')
ylabel('Number of detected features')
title('Detected features using SCRaPL')
set(gca,'FontSize',30)
legend({'Original','Negative control'},'FontSize',20)

name=fullfile(SaveFolder,File_out_diag1);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all


figure(2)
bar(x_or,y_prs(1))
hold on
bar(x_neg,y_prs(2:end))
xlabel('Dataset')
ylabel('Number of detected features')
title('Detected features using Pearson correlation')
set(gca,'FontSize',30)
legend({'Original','Negative control'},'FontSize',20)

name=fullfile(SaveFolder,File_out_diag2);
h=gcf;
set(h,'PaperOrientation','landscape');
set(h,'PaperUnits','normalized');
set(h,'PaperPosition', [0 0 1 1]);
print(gcf, '-dpdf', name);
close all
