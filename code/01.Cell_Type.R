setwd("/thinker/tb269store/user/liuzhenyu/Project/BLCA_Neuro/BLCA_Neuro/results/01.Cell_Type")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

load("01.Cell_Type.RData")

# 1. Load data & QC ----
data.dir <- "/thinker/tb269store/user/liuzhenyu/Project/BLCA_Neuro/Summary/1_Cellranger_result"
samples <- list.dirs(data.dir, recursive = FALSE, full.names = FALSE)
mycolor <- list()

## 1.1. Load data ----
# counts
rc <- list()
rc <- lapply(samples, function(x) {
    data <- Read10X(data.dir = file.path(data.dir, x, "filtered_feature_bc_matrix"))
    colnames(data) <- paste(x, colnames(data), sep = "_")
    data
})
rc <- do.call(cbind, rc)
saveRDS(rc, file = "../../data/rc.rds")
# rc <- readRDS(file = "../../data/rc.rds")
rm(samples)

# metadata
sample.info <- data.frame(
    row.names = colnames(rc),
    "Cell" = colnames(rc),
    "Sample" = gsub("_A.+?*$|_T.+?*$|_C.+?*$|_G.+?*$", "", colnames(rc))
) %>% mutate(
    "Condition" = gsub("_00.+?*$", "", Sample)
)

# singleR
singleR <- read.table("../../data/cell_annot_MouseRNAseqData_label.main.txt", header = TRUE, sep = "\t")
singleR$Cell <- gsub("-.+?*$", "-1", singleR$Cell)
singleR$Cell <- paste(singleR$orig.ident, singleR$Cell, sep = "_")
table(singleR$orig.ident)
table(singleR$orig.ident) / table(sample.info$Sample)

## 1.2. QC for cells ----
table(singleR$Cell %in% sample.info$Cell)
rc <- rc[, singleR$Cell]
sample.info <- sample.info[singleR$Cell, ] %>%
    mutate(
        "percent.HB" = singleR$percent.HB,
        "percent.mito" = singleR$percent.mito,
        "Cell_Type_SingleR" = singleR$raw_main_celltype
    )

blca.obj <- CreateSeuratObject(
    counts = rc,
    project = "BLCA_Neuro",
    meta.data = sample.info
)
rm(singleR, rc)

mycolor$Condition <- ggsci::pal_npg("nrc")(3) %>%
    setNames(unique(sample.info$Condition))

pdf("Scatter.QC.pdf", 20, 12)
ggplot(sample.info, aes(x = nCount_RNA, y = nFeature_RNA, color = Condition)) +
    geom_point(show.legend = FALSE) +
    scale_color_manual(values = mycolor$Condition) +
    theme_classic() +
    facet_wrap(~Sample, ncol = 5) +
    labs(title = "Scatter plot of QC metrics")
dev.off()

pdf("Hist.QC.pdf", 10, 4)
patchwork::wrap_plots(
    ggplot(sample.info, aes(x = log10(nCount_RNA))) +
        geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
        theme_classic() +
        labs(title = "Histogram of log10 UMI counts"),
    ggplot(sample.info, aes(x = nFeature_RNA)) +
        geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
        theme_classic() +
        labs(title = "Histogram of gene counts"),
    ncol = 2
)
dev.off()

## 1.3. Reductions & clustering -----
blca.obj <- blca.obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(blca.obj, ndims = 50)
dev.off()

blca.obj <- RunUMAP(blca.obj, dims = 1:30)

blca.obj$Condition <- factor(blca.obj$Condition, levels = c("CTRL", "Capsaicin", "LC_60HDA"))
pdf("UMAP.Condition.pdf", 7, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Condition") +
    scale_color_manual(values = mycolor$Condition) +
    theme_classic() +
    labs(title = "UMAP plot colored by Condition")
dev.off()

pdf("UMAP.Condition.split.pdf", 12, 4)
DimPlot(blca.obj, reduction = "umap", group.by = "Condition") +
    scale_color_manual(values = mycolor$Condition) +
    theme_classic() +
    labs(title = "UMAP plot colored by Condition")
dev.off()

blca.obj <- FindNeighbors(blca.obj, dims = 1:30)
blca.obj <- FindClusters(blca.obj, resolution = 0.8)

mycolor$Cluster <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(blca.obj$seurat_clusters))) %>%
    setNames(levels(blca.obj$seurat_clusters))

pdf("UMAP.Clusters.pdf", 7, 6)
DimPlot(blca.obj, reduction = "umap", label = TRUE) +
    scale_color_manual(values = mycolor$Cluster) +
    theme_classic() +
    labs(title = "UMAP plot colored by Clusters")
dev.off()

pdf("UMAP.Clusters.split.Condition.pdf", 16, 6)
DimPlot(blca.obj, reduction = "umap", label = TRUE, split.by = "Condition") +
    scale_color_manual(values = mycolor$Cluster) +
    theme_classic() +
    labs(title = "UMAP plot colored by Clusters")
dev.off()

pdf("UMAP.Clusters.split.Sample.pdf", 20, 15)
DimPlot(blca.obj, reduction = "umap", label = TRUE, split.by = "Sample", ncol = 5) +
    scale_color_manual(values = mycolor$Cluster) +
    theme_classic() +
    labs(title = "UMAP plot colored by Clusters")
dev.off()

## 1.4. QC metrics ----
pdf("VlnPlot.QC.Cluster.pdf", 12, 6)
VlnPlot(blca.obj,
    group.by = "seurat_clusters", features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.HB"),
    ncol = 2, pt.size = 0, cols = mycolor$Cluster
) +
    theme_classic()
dev.off()

pdf("UMAP.QC.Cluster.pdf", 8, 6)
FeaturePlot(blca.obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.HB"),
    ncol = 2, cols = c("gray", "red")
)
dev.off()

# 2. Initial cell type annotation ----
## 2.1. SingleR annotation ----
table(blca.obj$Cell_Type_SingleR)
pdf("UMAP.SingleR.pdf", 7, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

## 2.2. Expression of marker genes ----
marker.all <- c(
    "Epcam", "Krt5", "Krt8", "Krt18", # epithelial cell
    "Pecam1", "Cdh5", # endothelial cell,
    "Ptprc", "Cd3d", "Cd4", "Cd8a", "Nkg7", "Ncam1", "Cd79a", "Ms4a1", "Igkc", "Ighg1", # leukocyte
    "Cd14", "Fcgr3", "Lyz2", # monocyte
    "Cd68", "Cd163", "Mrc1", "Adgre1", # macrophage
    "Itgax", "Fcer1a", "Cst3", # DC
    "S100a8", "S100a9", # neutrophils
    "Kit", "Hdc", # mast cell
    "Col1a2", "Dcn", "Pdgfra", # fibroblast,
    "Tubb3", "Nefl", "Map2", # neuron
    "S100b", "Sox10", # Schwann cell
    "Fabp7", "P2rx7" # Satellite glial cells
)
setdiff(marker.all, rownames(blca.obj))
marker.all[which(duplicated(marker.all))]

p <- FeaturePlot(blca.obj, features = marker.all, ncol = 6, cols = c("gray", "red"))
pdf("UMAP.MarkerGenes.pdf", 24, 24)
plot(p)
dev.off()
rm(p)

pdf("DotPlot.MarkerGenes.pdf", 8, 8)
DotPlot(blca.obj, features = marker.all) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

## 2.3. Label cell types ----
blca.obj$Cell_Type_1 <- "none"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(21))] <- "Epithelial"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(0, 1, 2, 4, 6, 11, 14, 18, 19, 28))] <- "Cancer"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(26, 33))] <- "Enothelial"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(5, 9, 13, 15, 32))] <- "T"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(17, 25))] <- "NK"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(29))] <- "B"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(7, 10, 20, 23, 24))] <- "Macrophage"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(8, 22, 31))] <- "Monocyte"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(16, 27))] <- "DC"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(3))] <- "Granulocyte"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(12, 30))] <- "Fibroblast"
table(blca.obj$seurat_clusters, blca.obj$Cell_Type_1)

Idents(blca.obj) <- "Cell_Type_1"

blca.obj$Cell_Type_1 <- factor(blca.obj$Cell_Type_1, levels = c(
    "Epithelial", "Enothelial", "T", "NK", "B",
    "Macrophage", "Monocyte", "DC", "Granulocyte",
    "Fibroblast", "Cancer"
))
mycolor$Cell_Type_1 <- c(ggsci::pal_npg()(10), "#FDB462FF") %>%
    setNames(levels(blca.obj$Cell_Type_1))

pdf("UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_1", label = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

pdf("UMAP.Cell_Type_1.split.Condition.pdf", 16, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_1", label = TRUE, , split.by = "Condition") +
    scale_color_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

# bigger categories
blca.obj$Cell_Type_0 <- "none"
blca.obj$Cell_Type_0[which(blca.obj$Cell_Type_1 %in% c("Epithelial", "Cancer"))] <- "Epithelial"
blca.obj$Cell_Type_0[which(blca.obj$Cell_Type_1 %in% c("Fibroblast", "Enothelial"))] <- "Stromal"
blca.obj$Cell_Type_0[which(blca.obj$Cell_Type_1 %in% c(
    "T", "B", "NK", "Macrophage",
    "Monocyte", "DC", "Granulocyte"
))] <- "Immune"
table(blca.obj$Cell_Type_0)

blca.obj$Cell_Type_0 <- factor(blca.obj$Cell_Type_0, levels = c("Epithelial", "Stromal", "Immune"))
mycolor$Cell_Type_0 <- c(ggsci::pal_jco()(4)[c(1, 2, 4)]) %>%
    setNames(levels(blca.obj$Cell_Type_0))

pdf("UMAP.Cell_Type_0.pdf", 7, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_0", label = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_0")
dev.off()

pdf("UMAP.Cell_Type_0.split.Condition.pdf", 16, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_0", label = TRUE, , split.by = "Condition") +
    scale_color_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_0")
dev.off()

## 2.4. Cell types cross conditions ----
plot.data <- as.data.frame(table(blca.obj$Cell_Type_1, blca.obj$Condition))
colnames(plot.data) <- c("Cell_Type_1", "Condition", "Count")

pdf("Barplot.Cell_Type_1.Condition.pdf", 8, 5)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = Cell_Type_1, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot.data, aes(x = Condition, y = Count, fill = Cell_Type_1)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_1) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

plot.data <- as.data.frame(table(blca.obj$Cell_Type_0, blca.obj$Condition))
colnames(plot.data) <- c("Cell_Type_0", "Condition", "Count")

pdf("Barplot.Cell_Type_0.Condition.pdf", 8, 5)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = Cell_Type_0, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot.data, aes(x = Condition, y = Count, fill = Cell_Type_0)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_0) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()
rm(plot.data)

# 3. Detailed cell type ----
## 3.1. pan-Lymphocyte ----
dir.create("sub1_Lymphocyte")

sub1.obj <- subset(blca.obj, subset = Cell_Type_1 %in% c("T", "B", "NK")) # 49423 cells

sub1.obj <- sub1.obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub1_Lymphocyte/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(sub1.obj, ndims = 50)
dev.off()

sub1.obj <- FindNeighbors(sub1.obj, dims = 1:20)
sub1.obj <- FindClusters(sub1.obj, resolution = 0.8)

sub1.obj <- RunUMAP(sub1.obj, dims = 1:20)

pdf("sub1_Lymphocyte/UMAP.Clusters.pdf", 7, 6)
DimPlot(sub1.obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub1_Lymphocyte/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(sub1.obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub1_Lymphocyte/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(sub1.obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

marker.sub <- list()
marker.sub[["Lymphocyte"]] <- c(
    "Ptprc",
    "Cd3d", "Cd3g", "Cd4", "Foxp3",
    "Cd8a", "Gzma", "Gzmb", "Havcr2", "Pdcd1",
    "Nkg7", "Ncam1", "Klrb1c",
    "Cd79a", "Ms4a1", "Igkc", "Sdc1",
    "Cxcr2", "S100a9", "Trem1",
    "Kit", "Hdc", "Ms4a2",
    "Siglech", "Flt3"
)
length(marker.sub[["Lymphocyte"]])
setdiff(marker.sub[["Lymphocyte"]], rownames(sub1.obj))

pdf("sub1_Lymphocyte/UMAP.MarkerGenes.pdf", 20, 18)
FeaturePlot(sub1.obj,
    features = marker.sub[["Lymphocyte"]],
    ncol = 5, cols = c("gray", "red"), raster = TRUE
)
dev.off()

pdf("sub1_Lymphocyte/DotPlot.MarkerGenes.pdf", 8, 6)
DotPlot(sub1.obj, features = marker.sub[["Lymphocyte"]]) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

# undefined cluster markers
marker.undefined.T_11 <- FindMarkers(sub1.obj,
    ident.1 = 11,
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker.undefined.T_11 <- marker.undefined.T_11[order(marker.undefined.T_11$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker.undefined.T_11), 20) # Neutrophil markers: S100a9, S100a8, Cxcr2, Trem1

marker.undefined.T_2022 <- FindMarkers(sub1.obj,
    ident.1 = c(20, 22),
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker.undefined.T_2022 <- marker.undefined.T_2022[order(marker.undefined.T_2022$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker.undefined.T_2022), 20) # Mast markers: Hdc, Tpsb2, Ms4a2

marker.undefined.T_17 <- FindMarkers(sub1.obj,
    ident.1 = c(17),
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker.undefined.T_17 <- marker.undefined.T_17[order(marker.undefined.T_17$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker.undefined.T_17), 20) # pDC markers: Siglech, Flt3

sub1.obj$Cell_Type_2 <- "none"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(9))] <- "T_CD4"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(0:2, 4:6, 8:10, 12:14, 18, 19))] <- "T_CD8"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(11))] <- "T_CD8_Neutrophil_like"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(20, 22))] <- "Mast"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(17))] <- "pDC"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(3, 16))] <- "Treg"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(7))] <- "NK"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(15))] <- "B"
sub1.obj$Cell_Type_2[which(sub1.obj$seurat_clusters %in% c(21))] <- "Plasma"

pdf("sub1_Lymphocyte/UMAP.Cell_Type_2.pdf", 7, 6)
DimPlot(sub1.obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE,
    cols = mycolor$Cell_Type_1 %>% unname()
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()
saveRDS(sub1.obj, file = "sub1_Lymphocyte/sub1.obj.rds")
rm(sub1.obj, marker.undefined.T_11, marker.undefined.T_2022, marker.undefined.T_17)

## 3.2. pan-Myeloid ----
dir.create("sub2_Myeloid")
sub2.obj <- subset(blca.obj, subset = Cell_Type_1 %in% c("Macrophage", "Monocyte", "DC", "Granulocyte")) # 274268 cells

sub2.obj <- sub2.obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub2_Myeloid/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(sub2.obj, ndims = 50)
dev.off()

sub2.obj <- FindNeighbors(sub2.obj, dims = 1:20)
sub2.obj <- FindClusters(sub2.obj, resolution = 0.8)

sub2.obj <- RunUMAP(sub2.obj, dims = 1:20)

pdf("sub2_Myeloid/UMAP.Clusters.pdf", 7, 6)
DimPlot(sub2.obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Myeloid/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(sub2.obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Myeloid/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(sub2.obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

marker.sub[["Myeloid"]] <- c(
    "Ptprc",
    "Cd14", "Fcgr3", "Lyz2", # monocyte
    "Cd68", "Cd163", "Mrc1", "Adgre1", # macrophage
    "Itgax", "Fcer1a", "Cst3", # DC
    "Siglech", "Flt3", # pDC
    "Cxcr2", "S100a8", "S100a9", "Trem1", # neutrophils
    "Kit", "Hdc", "Ms4a2" # mast cell

)
length(marker.sub[["Myeloid"]])
setdiff(marker.sub[["Myeloid"]], rownames(sub2.obj))

pdf("sub2_Myeloid/UMAP.MarkerGenes.pdf", 20, 14)
FeaturePlot(sub2.obj,
    features = marker.sub[["Myeloid"]],
    ncol = 5, cols = c("gray", "red"), raster = TRUE
)
dev.off()

pdf("sub2_Myeloid/DotPlot.MarkerGenes.pdf", 8, 6)
DotPlot(sub2.obj, features = marker.sub[["Myeloid"]]) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

saveRDS(sub2.obj, file = "sub2_Myeloid/sub2.obj.rds")
# save data ----
save.image("01.Cell_Type.RData")