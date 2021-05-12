%% Joint  QC binomial-poisson data
%In this script we jointly QC binomial-poisson data. For very large
%datasets, there is the option of splitting it intno smale.
clc
clear all

type_of_data="MT";%Type of data. Available options MT/NT(binomial-count data), NM (binomial-binomial data).
genomic_region="prom";%Genomic region of interest. Available options are enh (enhancer) and prom (promoters).
[parentdir,~,~]=fileparts(pwd);%Get parent directory.
DataFolder=fullfile(append(parentdir,'/Join/',type_of_data,'/',genomic_region,'/'));%Folder with data. 
SaveFolder=append('Data_QC/',type_of_data,'/',genomic_region,'/');%Saving folder.

norm_fact = readtable(append(DataFolder,"norm_2500.csv"));
dt = readtable(fullfile(append(DataFolder,"prom_2500_2500.csv")));
save_name=strcat(append(SaveFolder,"gastrulation_2500_2500.mat"));

%Consider different columns for features/cells depending the the data and
%genomic region of interest.
if (type_of_data=="NT" || type_of_data=="MT" ) && (genomic_region=="prom")
    [gene_labels,~,gene_ind]=unique(dt(:,4));
    [cell_labels,~,cell_ind]=unique(dt(:,5));
    msg=1;
elseif type_of_data=="NM" && genomic_region=="prom"
    [gene_labels,~,gene_ind]=unique(dt(:,5));
    [cell_labels,~,cell_ind]=unique(dt(:,6));  
    msg=2;
elseif ((type_of_data=="NT" || type_of_data=="MT" ) && (genomic_region=="enh"))
    [gene_labels,~,gene_ind]=unique(dt(:,4:5));
    [cell_labels,~,cell_ind]=unique(dt(:,6));
    msg=3;
elseif type_of_data=="NM" && genomic_region=="enh"
    [gene_labels,~,gene_ind]=unique(dt(:,5:6));
    [cell_labels,~,cell_ind]=unique(dt(:,7));
    msg=4;
else
    msg = 'Please chech type of data and region of interest';
    error(msg)
end

%Relabel feature/cell columns using numbers.
dt_fin=dt;
dt_fin.feature=gene_ind;
dt_fin.cell=cell_ind;
dt_fin=table2array(dt_fin);

[a,~,c] = unique(dt_fin(:,end-1));
%Joint QC.
if msg==1 || msg ==3

gn_std=[histc(c,a),accumarray(c,dt_fin(:,3)==0,[],@mean),accumarray(c,dt_fin(:,1)./dt_fin(:,2),[],@std),accumarray(c,dt_fin(:,3),[],@std),a];
keep_gene=gn_std(gn_std(:,1)>4 & gn_std(:,2)<0.8 & gn_std(:,3)>0.0 & gn_std(:,4)>0 ,end);
  
else
gn_std=[histc(c,a),accumarray(c,dt_fin(:,1)./dt_fin(:,2),[],@std),accumarray(c,dt_fin(:,3)./dt_fin(:,4),[],@std),a];
keep_gene=gn_std(gn_std(:,1)>4 & gn_std(:,2)>0.0 & gn_std(:,3)>0.0 & gn_std(:,4)>0 ,end);   
    
end

%Remove features that do not pass joint QC and relabel features such that
%there are no missing indices between the first and the last.
keep_ind=ismember(c,keep_gene);
dt_fin=dt_fin(keep_ind,:);

map=[keep_gene,(1:size(keep_gene,1))'];

[~, loc] = ismember(dt_fin(:,end-1), map(:,1));
dt_fin(:,end-1) = map(loc,2);
out_genes=[unique(dt_fin(:,end-1)),histc(dt_fin(:,end-1),unique(dt_fin(:,end-1)))];
gene_labels=gene_labels(keep_gene,:);
y_true=dt_fin;

%Save results.
if msg==1 && msg==3
save(save_name,"y_true","norm_fact","out_genes","cell_labels","gene_labels")
else
save(save_name,"y_true","out_genes","cell_labels","gene_labels")    
end

%% Splitting data
%The code of this section splits features in a predefind number of subsets.
%This helps in cases where there are several low/medium size computational
%units.

Num_break=10;
Num_genes=size(gene_labels,1);
gene_ind=round(1:(Num_genes-1)/Num_break:Num_genes);
y_true_full=y_true;
ref_gen=gene_labels;
for ii=1:Num_break
   ind=[gene_ind(ii),gene_ind(ii+1)];
   if ii == Num_break 
        y_true=y_true_full(y_true_full(:,end-1) >= gene_ind(ii) & y_true_full(:,end-1) <= gene_ind(ii+1),:);
        gene_labels=gene_labels(unique(y_true(:,end-1)),1);
        out_genes=out_genes(unique(y_true(:,end-1)),:);
        [~,~,y_true(:,end-1)]=unique(y_true(:,end-1));
        out_genes(:,1)=unique(y_true(:,end-1));
   else
        y_true=y_true_full(y_true_full(:,end-1) >= gene_ind(ii) & y_true_full(:,end-1)< gene_ind(ii+1),:);
        gene_labels=gene_labels(unique(y_true(:,end-1)),1);
        out_genes=out_genes(unique(y_true(:,end-1)),:);
        [~,~,y_true(:,end-1)]=unique(y_true(:,end-1));
        out_genes(:,1)=unique(y_true(:,end-1));
   end
   name_part=sprintf(append(SaveFolder,"NAMEofChoice_part%d"),ii);
   if msg==1 && msg==3
    save(name_part,"y_true","norm_fact","out_genes","cell_labels","gene_labels","ind","ref_gen")
   else
    save(name_part,"y_true","out_genes","cell_labels","gene_labels","ind","ref_gen")    
   end
    
end