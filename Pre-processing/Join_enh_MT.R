#This script joins methylation and expression data for promoter regions.
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
io$met_dir    <- paste0(io$base_dir, "/Meth_Acc/met/processed/unfiltered/out/")
io$out_dir    <- paste0(io$base_dir, "/Join/MT/enh/")
io$exp_data <-"SingleCellExperimentNMT_filt.rds"
io$meta <- "sample_metadata.csv"
io$met_data <- "enh_25000.csv"
io$met_rna_out <- "enh_25000.csv"
io$norm_out <- "norm_enh_25000.csv"


rna_sce_exp<- readRDS( paste0(io$exp_dir, io$exp_data))
rna_data <- assay(rna_sce_exp)
norm_fact <- norm_fact<-as.data.frame(rna_sce_exp@int_colData@listData[["size_factor"]]) %>% setnames(c("norm_fact"))
norm_fact[,"cell"] <- colnames(rna_data)

rna_data <- rna_data[order(rownames(rna_data)),order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

meta <- fread( paste0(io$exp_dir, io$meta) , showProgress = FALSE) #Lab QC. If not available comment all lines involving meta
meta <- meta[meta$pass_rnaQC==TRUE & meta$pass_metQC==TRUE ,]#Lab QC. If not available comment all lines involving meta
meta <- meta[order(meta$id_rna),]#Lab QC. If not available comment all lines involving meta

meta <- meta[is.element(meta$id_rna,colnames(rna_data)),]#Lab QC. If not available comment all lines involving meta
rna_data <- rna_data[,is.element(colnames(rna_data),meta$id_rna)]#Lab QC. If not available comment all lines involving meta
norm_fact <- norm_fact[is.element(norm_fact$cell,meta$id_rna),]#Lab QC. If not available comment all lines involving meta

colnames(rna_data) <- meta$sample#Lab QC. If not available comment all lines involving meta
norm_fact$cell <- meta$sample#Lab QC. If not available comment all lines involving meta
rna_data <- rna_data[,order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

met_data <- fread(paste0(io$met_dir,io$met_data))
met_data$cell <- gsub("~/SCRaPL/Pre-processing/Meth_Acc/met/binarised/", "\\1", met_data$cell)#This renaming is dataset-specific
met_data$cell <- gsub(".tsv.gz", "\\1", met_data$cell)
met_data <- met_data[is.element(met_data$cell,meta$sample),]  %>% setnames(c("met_cpgs","total_cpgs","feature","Gene","cell"))
met_data <- met_data[order(met_data$Gene,met_data$cell),]

rna_data <- as.data.frame(as.table(as.matrix(rna_data))) %>% setnames(c("Gene","cell","expr"))
rna_data <- rna_data[order(rna_data$Gene,rna_data$cell),]

met_rna <- merge(met_data,rna_data,allow.cartesian=TRUE)
met_rna <- met_rna[,c(3,4,6,5,1,2)]
met_rna <- met_rna[met_rna$met_cpgs>-1,]
met_rna <- met_rna[order(met_rna$Gene,met_rna$cell),]

fwrite(norm_fact, file = paste0(io$out_dir,io$norm_out), sep = "\t")
fwrite(met_rna, file = paste0(io$out_dir, io$met_rna_out) , sep = "\t")