suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(SeuratData))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))

#This script is used to load and preprocess PBMC data.
io <- list()
io$base_dir   <- "/SCRaPL/SCRaPL_TF/Real/Results_atac_human_col/" #This line must change depending on the folder integration is performed.
io$raw_acc <- "human_acc_tmp_1.csv"
io$raw_rna <- "human_rna_tmp_1.csv"
io$nrm_acc <- "nrm_human_acc_tmp.csv"
io$nrm_rna <- "nrm_human_rna_tmp.csv"
io$feature_sel <- "features_select_1.csv"


feature_sel <- fread(paste0(io$base_dir,io$feature_sel),sep = ",",header = TRUE)
feature_sel <- feature_sel[,-c(1)]

pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")
pbmc.atac <- pbmc.atac[,pbmc.atac@meta.data$nCount_ATAC>23953]
cell_kp <- colnames(pbmc.atac@assays$ATAC@data)

sce_acc <- as.SingleCellExperiment(pbmc.atac)
clusters_acc <- quickCluster(sce_acc, min.size=50)
sce_acc <- computeSumFactors(sce_acc, clusters=clusters_acc)
sce_acc <- logNormCounts(sce_acc)

acc_nrm <- as.data.frame(sce_acc@colData)
acc_nrm <- acc_nrm[,-c(1,2,3,5)]

pbmc.atac <- FindVariableFeatures(pbmc.atac, nfeatures = 30000)
pbmc.atac <- subset(pbmc.atac, features = VariableFeatures(pbmc.atac))

acc_raw <- as.data.frame(as.matrix(pbmc.atac@assays$ATAC@counts))
pks_all <- as.data.table(rownames(acc_raw)) %>% setnames(c("peaks_all"))
pks_kpt <- as.data.table(pks_all$peaks_all[feature_sel$peak_ind]) %>% setnames(c("peaks"))
acc_raw <- acc_raw[feature_sel$peak_ind,]

fwrite(x=acc_raw,file = paste0(io$base_dir,io$raw_acc))
fwrite(x=acc_nrm,file = paste0(io$base_dir,io$nrm_acc))

rm(pbmc.atac,sce_acc,clusters_acc)

pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
pbmc.rna <- subset(x= pbmc.rna,cells=cell_kp)

sce_rna <- as.SingleCellExperiment(pbmc.rna)
clusters_rna <- quickCluster(sce_rna, min.size=50)
sce_rna <- computeSumFactors(sce_rna, clusters=clusters_rna)
sce_rna <- logNormCounts(sce_rna)

rna_nrm <- as.data.frame(sce_rna@colData)
rna_nrm <- rna_nrm[,-c(1,2,3,5)]

pbmc.rna <- FindVariableFeatures(pbmc.rna, nfeatures = 10000)
pbmc.rna <- subset(pbmc.rna, features = VariableFeatures(pbmc.rna))

rna_raw <- as.data.frame(as.matrix(pbmc.rna@assays$RNA@counts))
gns_all <- as.data.table(rownames(rna_raw)) %>% setnames(c("genes_all"))
gns_kpt <- gns_all[feature_sel$gene_ind,] %>% setnames(c("genes"))
rna_raw <- rna_raw[feature_sel$gene_ind,] 



fwrite(x=rna_raw,file = paste0(io$base_dir,io$raw_rna))
fwrite(x=rna_nrm,file = paste0(io$base_dir,io$nrm_rna))

fwrite(x=gns_kpt,file = paste0(io$base_dir,"gns_tmp_1.csv"))
fwrite(x=pks_kpt,file = paste0(io$base_dir,"pks_tmp_1.csv"))

rm(pbmc.rna,sce_rna)