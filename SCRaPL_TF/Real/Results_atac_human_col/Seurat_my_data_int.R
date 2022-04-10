library(Seurat)
library(Signac)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(tidyr)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(patchwork)
#This script is used to Run Seurat's integration analysis on data sampled from SCRaPL.
io <- list()
io$dt_ext <- "_4"
io$base_dir   <- "/SCRaPL/SCRaPL_TF/Real/Results_atac_human_col/"
io$data_acc <- paste0("Human_acc_seur_tmp",io$dt_ext,".csv")

acc <- read.table( paste0(io$base_dir,io$data_acc),sep = ",",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
colnames(acc) <- sub(".1","-1",colnames(acc))

pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
pbmc.atac@assays$ATAC@key <- "atac_"
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")

pbmc.atac <- pbmc.atac[,pbmc.atac@meta.data$nCount_ATAC>23953]
cell_kp <- colnames(pbmc.atac@assays$ATAC@data)

pbmc.atac <- pbmc.atac[rownames(acc),]
acc_seur <- CreateSeuratObject(log(1+acc),assay = "ATAC", project = "ATAC")#
pbmc.atac@assays[["ATAC"]]@counts@x <- acc_seur@assays[["ATAC"]]@counts@x
pbmc.atac@assays[["ATAC"]]@data@x <- acc_seur@assays[["ATAC"]]@data@x

io$data_exp <- paste0("Human_exp_seur_tmp",io$dt_ext,".csv")
rna <- read.table( paste0(io$base_dir,io$data_exp),sep = ",",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
#colnames(rna) <- sub(".1","-1",colnames(rna))
rna_seur <- CreateSeuratObject(counts = log(1+rna), project = "RNA" )#

pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
pbmc.rna <- subset(x= pbmc.rna,cells=cell_kp)
pbmc.rna <- pbmc.rna[rownames(rna_seur),]
pbmc.rna@assays[["RNA"]]@counts@x <- rna_seur@assays[["RNA"]]@counts@x
pbmc.rna@assays[["RNA"]]@data@x <- rna_seur@assays[["RNA"]]@data@x

# We exclude the first dimension as this is typically correlated with sequencing depth

#my_rna <- NormalizeData(my_rna)########################
pbmc.rna <- FindVariableFeatures(pbmc.rna, nfeatures = 2000)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)


# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc.atac) <- annotations

#my_acc <- RunTFIDF(my_acc)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

p1 <- DimPlot(pbmc.rna, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
print(p1 + p2)

# quantify gene activity
gene.activities_n <- GeneActivity(object=pbmc.atac, features = VariableFeatures(pbmc.rna))

# add gene activities as a new assay
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities_n)

# normalize gene activities
DefaultAssay(pbmc.atac) <- "ACTIVITY"

#my_acc <- NormalizeData(my_acc)###########################
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))

transfer.anchors_my <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),
                                                                                reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions_my <- TransferData(anchorset = transfer.anchors_my, refdata = pbmc.atac$seurat_annotations,
                                                                 weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions_my)
pbmc.atac$annotation_correct <- pbmc.atac$predicted.id == pbmc.atac$seurat_annotations

p1 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(pbmc.atac, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
print(p1 | p2)

predictions_my <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
predictions_my <- predictions_my/rowSums(predictions_my)  # normalize for number of cells in each cell type
predictions_my <- as.data.frame(predictions_my)

p1 <- ggplot(predictions_my, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
   low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct_my <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect_my <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data_my <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data_my, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect_my, ")"), paste0("TRUE (n = ", correct_my, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                     labels = c(paste0("FALSE (n = ", incorrect_my, ")"), paste0("TRUE (n = ", correct_my, ")"))) + xlab("Prediction Score")

print(p1 + p2)

genes.use_my <- VariableFeatures(pbmc.rna)
refdata_my <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use_my, ]

imputation_my <- TransferData(anchorset = transfer.anchors_my, refdata = refdata_my, weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)
pbmc.atac[["RNA"]] <- imputation_my

coembed_my <- merge(x = pbmc.rna, y = pbmc.atac)

coembed_my <- ScaleData(coembed_my , features = genes.use_my, do.scale = FALSE)
coembed_my  <- RunPCA(coembed_my , features = genes.use_my, verbose = FALSE)
coembed_my  <- RunUMAP(coembed_my , dims = 1:30)

coembed_my <- AddMetaData(object = coembed_my,metadata = coembed_my@meta.data$orig.ident,col.name = "Modality")
coembed_my <- AddMetaData(object = coembed_my,metadata = coembed_my@meta.data$seurat_annotations,col.name = "Seurat_Annotations")

p1 <- DimPlot(coembed_my, group.by = c("Modality", "Seurat_Annotations"))
print(p1)

