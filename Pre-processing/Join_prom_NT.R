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
io$out_dir    <- paste0(io$base_dir, "/Join/NT/prom/")
io$exp_data <-"SingleCellExperimentNMT_filt.rds"
io$meta <- "sample_metadata.csv"
io$acc_data <- "prom_2500_2500.csv"
io$met_rna_out <- "prom_2500_2500.csv"
io$norm_out <- "norm_2500.csv"


rna_sce_exp<- readRDS( paste0(io$exp_dir, io$exp_data))
rna_data <- assay(rna_sce_exp)
norm_fact <- norm_fact<-as.data.frame(rna_sce_exp@int_colData@listData[["size_factor"]]) %>% setnames(c("norm_fact"))
norm_fact[,"cell"] <- colnames(rna_data)

rna_data <- rna_data[order(rownames(rna_data)),order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

meta <- fread( paste0(io$exp_dir, io$meta) , showProgress = FALSE)
meta <- meta[meta$pass_rnaQC==TRUE & meta$pass_accQC==TRUE ,]
meta <- meta[order(meta$id_rna),]


meta <- meta[is.element(meta$id_rna,colnames(rna_data)),]
rna_data <- rna_data[,is.element(colnames(rna_data),meta$id_rna)]
norm_fact <- norm_fact[is.element(norm_fact$cell,meta$id_rna),]

colnames(rna_data) <- meta$sample
norm_fact$cell <- meta$sample
rna_data <- rna_data[,order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

acc_data <- fread(paste0(io$met_dir,io$acc_data))
acc_data$cell <- gsub("~/SCRaPL/Pre-processing/Meth_Acc/acc/binarised/", "\\1", acc_data$cell)#This renaming is dataset-specific
acc_data$cell <- gsub(".tsv.gz", "\\1", acc_data$cell)
acc_data <- acc_data[is.element(acc_data$cell,meta$sample),]
acc_data <- acc_data[order(acc_data$feature,acc_data$cell),]

rna_data <- as.data.frame(as.table(as.matrix(rna_data))) %>% setnames(c("feature","cell","expr"))
rna_data <- rna_data[order(rna_data$feature,rna_data$cell),]

acc_rna <- merge(acc_data,rna_data)
acc_rna <- acc_rna[,c(3,4,5,1,2)]
acc_rna <- acc_rna[acc_rna$met_cpgs>-1,]
acc_rna <- acc_rna[order(acc_rna$feature,acc_rna$cell),]

fwrite(norm_fact, file = paste0(io$out_dir,io$norm_out), sep = "\t")
fwrite(met_rna, file = paste0(io$out_dir, io$met_rna_out) , sep = "\t")