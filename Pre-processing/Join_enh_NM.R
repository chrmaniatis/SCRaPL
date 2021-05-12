#This script joins accessibility and methylation data for promoter regions.
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
io$met_dir    <- paste0(io$base_dir, "/Meth_Acc/met/processed/unfiltered/out/")
io$acc_dir    <- paste0(io$base_dir, "/Meth_Acc/acc/processed/unfiltered/out/")
io$out_dir    <- paste0(io$base_dir, "/Join/NM/enh/")
io$meta <- "sample_metadata.csv"
io$acc_data <- "enh_25000.csv"
io$met_rna_out <- "enh_25000.csv"

#Lab QC
meta <- fread( paste0(io$exp_dir, io$meta) , showProgress = FALSE)
meta <- meta[meta$pass_accQC==TRUE & meta$pass_metQC==TRUE ,]
meta <- meta[order(meta$id_rna),]
meta <- meta[is.element(meta$id_rna,colnames(rna_data)),]


acc_data <- fread(paste0(io$acc_dir,io$acc_data))
acc_data$cell <- gsub("~/SCRaPL/Pre-processing/Meth_Acc/acc/binarised/", "\\1", acc_data$cell)#This renaming is dataset-specific
acc_data$cell <- gsub(".tsv.gz", "\\1", acc_data$cell)
acc_data <- acc_data[is.element(acc_data$cell,meta$sample),]#Quality control step (if no lab QC for cells does not exist please comment this line)
acc_data <- acc_data[order(acc_data$feature,acc_data$cell),] %>% setnames(c("acc_gpcs","total_gpcs","feature","Gene","cell"))

met_data <- fread(paste0(io$met_dir,io$acc_data))
met_data$cell <- gsub("~/SCRaPL/Pre-processing/Meth_Acc/met/binarised/", "\\1", met_data$cell)#This renaming is dataset-specific
met_data$cell <- gsub(".tsv.gz", "\\1", met_data$cell)
met_data <- acc_data[is.element(met_data$cell,meta$sample),]#Quality control step (if no lab QC for cells does not exist please comment this line)
met_data <- met_data[order(met_data$feature,met_data$cell),] %>% setnames(c("met_cpgs","total_cpgs","feature","Gene","cell"))

met_acc <- merge(met_data,acc_data)
met_acc <- met_acc[,c(4,5,6,7,1,2,3)]
met_acc <- met_acc[met_acc$met_cpgs>-1 & met_acc$acc_gpcs>-1,]
met_acc <- met_acc[order(met_acc$Gene,met_acc$cell),]

fwrite(met_acc, file = paste0(io$out_dir, io$met_rna_out) , sep = "\t")