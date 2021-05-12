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
io$out_dir    <- paste0(io$base_dir, "/Join/MT/prom/")
io$exp_data <-"SingleCellExperimentNMT_filt.rds"
io$meta <- "sample_metadata.csv"
io$met_data <- "prom_2500_2500.csv"
io$met_rna_out <- "prom_2500_2500.csv"
io$norm_out <- "norm_2500.csv"


rna_sce_exp<- readRDS( paste0(io$exp_dir, io$exp_data))
rna_data <- assay(rna_sce_exp)
norm_fact <- norm_fact<-as.data.frame(rna_sce_exp@int_colData@listData[["size_factor"]]) %>% setnames(c("norm_fact"))
norm_fact[,"cell"] <- colnames(rna_data)

rna_data <- rna_data[order(rownames(rna_data)),order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

meta <- fread( paste0(io$exp_dir, io$meta) , showProgress = FALSE)
meta <- meta[meta$pass_rnaQC==TRUE & meta$pass_metQC==TRUE ,]
meta <- meta[order(meta$id_rna),]


meta <- meta[is.element(meta$id_rna,colnames(rna_data)),]
rna_data <- rna_data[,is.element(colnames(rna_data),meta$id_rna)]
norm_fact <- norm_fact[is.element(norm_fact$cell,meta$id_rna),]

colnames(rna_data) <- meta$sample
norm_fact$cell <- meta$sample
rna_data <- rna_data[,order(colnames(rna_data))]
norm_fact <- norm_fact[order(norm_fact$cell),]

met_data <- fread(paste0(io$met_dir,io$met_data))
met_data$cell <- gsub("~/SCRaPL/Pre-processing/Meth_Acc/met/binarised/", "\\1", met_data$cell)#This renaming is dataset-specific
met_data$cell <- gsub(".tsv.gz", "\\1", met_data$cell)
met_data <- met_data[is.element(met_data$cell,meta$sample),]
met_data <- met_data[order(met_data$feature,met_data$cell),]

rna_data <- as.data.frame(as.table(as.matrix(rna_data))) %>% setnames(c("feature","cell","expr"))
rna_data <- rna_data[order(rna_data$feature,rna_data$cell),]

met_rna <- merge(met_data,rna_data)
met_rna <- met_rna[,c(3,4,5,1,2)]
met_rna <- met_rna[met_rna$met_cpgs>-1,]
met_rna <- met_rna[order(met_rna$feature,met_rna$cell),]

fwrite(norm_fact, file = paste0(io$out_dir,io$norm_out), sep = "\t")
fwrite(met_rna, file = paste0(io$out_dir, io$met_rna_out) , sep = "\t")