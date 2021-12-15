library(Seurat)
library(Signac)
library(SeuratData)
library(ggplot2)
library(cowplot)
library(tidyr)
library(data.table)
library(EnsDb.Hsapiens.v86)
library(patchwork)
#This script performs integration analysis
io <- list()
io$base_dir   <- "/home/christos/Desktop/ATAC-seq-trapnel/Files_new/Human/" #This line must change depending on the folder integration is performed.
io$data_acc <- "Human_acc_seur_tmp_1.csv"
io$data_exp <- "Human_exp_seur_tmp_1.csv"

acc <- read.table( paste0(io$base_dir,io$data_acc),sep = ",",header = TRUE,row.names = 1,stringsAsFactors = FALSE)
rna <- read.table( paste0(io$base_dir,io$data_exp),sep = ",",header = TRUE,row.names = 1,stringsAsFactors = FALSE)

pks <- as.data.table(rownames(acc)) %>% setnames(c("peaks"))
gns <- as.data.table(rownames(rna)) %>% setnames(c("genes"))

rm(acc,rna)

pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")
pbmc.atac <- pbmc.atac[,pbmc.atac@meta.data$nCount_ATAC>23953]
pbmc.atac <- pbmc.atac[pks$peaks,]
cell_kp <- colnames(pbmc.atac@assays$ATAC@data)

pbmc.rna <- LoadData("pbmcMultiome", "pbmc.rna")
pbmc.rna <- subset(pbmc.rna, seurat_annotations != "filtered")
pbmc.rna <- subset(x= pbmc.rna,cells=cell_kp)
pbmc.rna <- pbmc.rna[gns$genes,]

pbmc.rna <- NormalizeData(pbmc.rna)
pbmc.rna <- FindVariableFeatures(pbmc.rna,nfeatures = 2000)
pbmc.rna <- ScaleData(pbmc.rna)
pbmc.rna <- RunPCA(pbmc.rna)
pbmc.rna <- RunUMAP(pbmc.rna, dims = 1:30)

# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(pbmc.atac) <- annotations
# We exclude the first dimension as this is typically correlated with sequencing depth

pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")
pbmc.atac <- RunSVD(pbmc.atac)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

p1 <- DimPlot(pbmc.rna, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(pbmc.atac, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")

print(p1+p2)
# quantify gene activity
gene.activities <- GeneActivity(object=pbmc.atac, features = VariableFeatures(pbmc.rna))

# add gene activities as a new assay
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(pbmc.atac) <- "ACTIVITY"

pbmc.atac <- NormalizeData(pbmc.atac)
pbmc.atac <- ScaleData(pbmc.atac, features = rownames(pbmc.atac))
transfer.anchors <- FindTransferAnchors(reference = pbmc.rna, query = pbmc.atac, features = VariableFeatures(object = pbmc.rna),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = pbmc.rna$seurat_annotations,
                                     weight.reduction = pbmc.atac[["lsi"]], dims = 2:30)

pbmc.atac <- AddMetaData(pbmc.atac, metadata = celltype.predictions)
pbmc.atac$annotation_correct <- pbmc.atac$predicted.id == pbmc.atac$seurat_annotations

p1 <- DimPlot(pbmc.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(pbmc.atac, group.by = "seurat_annotations", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")

print(p1|p2)

predictions <- table(pbmc.atac$seurat_annotations, pbmc.atac$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                            low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(pbmc.atac$seurat_annotations == pbmc.atac$predicted.id))
incorrect <- length(which(pbmc.atac$seurat_annotations != pbmc.atac$predicted.id))
data <- FetchData(pbmc.atac, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
  geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
                                                                                                                    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
print(p1 + p2)

genes.use <- VariableFeatures(pbmc.rna)
refdata <- GetAssayData(pbmc.rna, assay = "RNA", slot = "data")[genes.use, ]
# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = pbmc.atac[["lsi"]],dims = 2:30)

pbmc.atac[["RNA"]] <- imputation
coembed <- merge(x = pbmc.rna, y = pbmc.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets

coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)

coembed <- AddMetaData(object = coembed,metadata = coembed@meta.data$orig.ident,col.name = "Modality")
coembed <- AddMetaData(object = coembed,metadata = coembed@meta.data$seurat_annotations,col.name = "Seurat_Annotations")


p1 <- DimPlot(coembed, group.by =  c("Modality", "Seurat_Annotations"))
print(p1)

