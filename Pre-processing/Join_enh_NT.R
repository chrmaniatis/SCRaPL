#This script joins accessibility and expression data for promoter regions.
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

io <- list()
io$base_dir   <- "~/Desktop/SCRaPL/Pre-processing"
io$exp_dir    <- paste0(io$base_dir, "/RNA/")
io$acc_dir    <- paste0(io$base_dir, "/Meth_Acc/acc/processed/unfiltered/out/")
io$out_dir    <- paste0(io$base_dir, "/Join/NT/enh/")
io$exp_data <-"SingleCellExperimentNMT.rds"
io$meta <- "sample_metadata.csv"
io$acc_data <- "enh_25000.csv"
io$met_rna_out <- "enh_25000.csv"
io$norm_out <- "norm_enh_25000.csv"


rna_sce_exp<- readRDS( paste0(io$exp_dir, io$exp_data))
rna_data <- assay(rna_sce_exp)
norm_fact <- norm_fact<-as.data.frame(rna_sce_exp@int_colData@listData[["size_factor"]]) %>% setnames(c("norm_fact"))
norm_fact[,"cell"] <- colnames(rna_data)

rna_data <- rna_data[order(rownames(rna_data)),order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

meta <- fread( paste0(io$exp_dir, io$meta) , showProgress = FALSE) #Lab QC. If not available comment all lines involving meta
meta <- meta[meta$pass_rnaQC==T & meta$pass_accQC==T ,] #Lab QC. If not available comment all lines involving meta
meta <- meta[order(meta$id_rna),]#Lab QC. If not available comment all lines involving meta

meta <- meta[is.element(meta$id_rna,colnames(rna_data)),]#Lab QC. If not available comment all lines involving meta
rna_data <- rna_data[,is.element(colnames(rna_data),meta$id_rna)]#Lab QC. If not available comment all lines involving meta
norm_fact <- norm_fact[is.element(norm_fact$cell,meta$id_rna),]#Lab QC. If not available comment all lines involving meta

colnames(rna_data) <- meta$sample#Lab QC. If not available comment all lines involving meta
norm_fact$cell <- meta$sample#Lab QC. If not available comment all lines involving meta
rna_data <- rna_data[,order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

acc_data <- fread(paste0(io$acc_dir,io$acc_data))
acc_data$cell <- gsub("~/SCRaPL/Pre-processing/Meth_Acc/acc/binarised/", "\\1", acc_data$cell)#This renaming is dataset-specific
acc_data$cell <- gsub(".tsv.gz", "\\1", acc_data$cell)
acc_data <- acc_data[is.element(acc_data$cell,meta$sample),] %>% setnames(c("acc_gpcs","total_gpcs","feature","Gene","cell"))
acc_data <- acc_data[order(acc_data$Gene,acc_data$cell),]

rna_data <- as.data.frame(as.table(as.matrix(rna_data))) %>% setnames(c("Gene","cell","expr"))
rna_data <- rna_data[order(rna_data$feature,rna_data$cell),]

acc_rna <- merge(acc_data,rna_data,allow.cartesian=TRUE)
acc_rna <- acc_rna[,c(3,4,6,1,5,2)]
acc_rna <- acc_rna[acc_rna$met_cpgs>-1,]
acc_rna <- acc_rna[order(acc_rna$feature,acc_rna$cell),]

fwrite(norm_fact, file = paste0(io$out_dir,io$norm_out), sep = "\t")
fwrite(met_rna, file = paste0(io$out_dir, io$met_rna_out) , sep = "\t")