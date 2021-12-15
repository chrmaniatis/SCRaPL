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

io <- list()
io$base_dir   <- "/home/christos/Desktop/ATAC-seq-trapnel/Files_new/Human/"

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

#To estimate Pearson correlation for choosing feature pairs please uncomment following lines
#rna_proc <- log(1+rna_raw/rna_nrm$sizeFactor)
#acc_proc <- log(1+acc_raw/acc_nrm$sizeFactor)

#corr <- cor(t(acc_proc),t(rna_proc))
#fwrite(x=corr,file = paste0(io$base_dir,"feature_cor_tmp.csv"))