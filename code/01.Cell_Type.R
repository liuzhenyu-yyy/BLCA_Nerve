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
marker.genes <- list()

marker.genes[["All"]] <- c(
    "Epcam", "Krt5", "Krt8", "Krt18", # epithelial cell
    "Pecam1", "Cdh5", # endothelial cell,
    "Ptprc", "Cd3d", "Cd3g", "Cd4", "Cd8a", "Nkg7", "Klrb1c", "Cd79a", "Ms4a1", "Sdc1", # leukocyte
    "Cd14", "Cd68", "Lyz2", # monocyte
    "Cd163", "Mrc1", "Adgre1", # macrophage
    "Itgax", "Fcer1a", "Cst3", # DC
    "S100a8", "S100a9", # neutrophils
    "Kit", "Hdc", # mast cell
    "Col1a2", "Dcn", "Pdgfra", # fibroblast,
    "Tubb3", "Nefl", "Map2", # neuron
    "S100b", "Sox10", # Schwann cell
    "Fabp7", "P2rx7", # Satellite glial cells,
    "Mki67", "Top2a" # cycling
)
setdiff(marker.genes$All, rownames(blca.obj))
marker.genes$All[which(duplicated(marker.genes$All))]

p <- FeaturePlot(blca.obj, features = marker.genes$All, ncol = 6, cols = c("gray", "red"))
pdf("UMAP.MarkerGenes.pdf", 24, 24)
plot(p)
dev.off()
rm(p)

pdf("DotPlot.MarkerGenes.pdf", 9, 8)
DotPlot(blca.obj, features = marker.genes$All) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

## 2.3. Label cell types ----
blca.obj$Cell_Type_1 <- "none"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(21))] <- "Bladder_epithelial"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(0, 1, 2, 4, 6, 11, 14, 18, 19, 28))] <- "Cancer"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(26, 33))] <- "Endothelial"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(5, 9, 13, 15, 25, 32))] <- "T"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(17))] <- "NK"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(29))] <- "B"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(7, 8, 10, 16, 20, 22, 23, 24, 27, 31))] <- "Myeloid"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(3))] <- "Granulocyte"
blca.obj$Cell_Type_1[which(blca.obj$seurat_clusters %in% c(12, 30))] <- "Fibroblast"
table(blca.obj$seurat_clusters, blca.obj$Cell_Type_1)

Idents(blca.obj) <- "Cell_Type_1"

blca.obj$Cell_Type_1 <- factor(blca.obj$Cell_Type_1, levels = c(
    "Bladder_epithelial", "Cancer",
    "Endothelial", "Fibroblast",
    "T", "B", "NK",
    "Granulocyte", "Myeloid"
))
mycolor$Cell_Type_1 <- c(ggsci::pal_npg()(9)) %>%
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

gene.selected <- c(
    "Epcam", "Krt8", "Krt18", # epithelial
    "Pecam1", "Cdh5", # endothelial
    "Col1a2", "Dcn", "Pdgfra", # fibroblast
    "Ptprc", # leukocyte
    "Cd3d", "Cd3g", "Cd79a", "Klrb1c", # lymphocyte
    "S100a8", "S100a9", # neutrophils
    "Cd14", "Cd68" # myeloid
)
blca.obj$Cell_Type_1
pdf("DotPlot.Cell_Type_1.pdf", 7, 3)
DotPlot(blca.obj, features = gene.selected, group.by = "Cell_Type_1") +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# bigger categories
blca.obj$Cell_Type_0 <- "none"
blca.obj$Cell_Type_0[which(blca.obj$Cell_Type_1 %in% c("Bladder_epithelial", "Cancer"))] <- "Epithelial"
blca.obj$Cell_Type_0[which(blca.obj$Cell_Type_1 %in% c("Fibroblast", "Endothelial"))] <- "Stromal"
blca.obj$Cell_Type_0[which(blca.obj$Cell_Type_1 %in% c("T", "B", "NK", "Myeloid", "Granulocyte"))] <- "Immune"
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

# fraction compare
plot.data <- as.data.frame(table(blca.obj$Cell_Type_0, blca.obj$Sample))
colnames(plot.data) <- c("Cell_Type_0", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(blca.obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_60HDA"))

pdf("Boxplot.Fraction.Cell_Type_0.pdf", 5, 3)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_0), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_60HDA"), c("CTRL", "LC_60HDA"))
    ) +
    facet_wrap(~Cell_Type_0, ncol = 3, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()

plot.data <- as.data.frame(table(blca.obj$Cell_Type_1, blca.obj$Sample))
colnames(plot.data) <- c("Cell_Type_1", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(blca.obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_60HDA"))

pdf("Boxplot.Fraction.Cell_Type_1.pdf", 8, 5)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_1), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_60HDA"), c("CTRL", "LC_60HDA"))
    ) +
    facet_wrap(~Cell_Type_1, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot.data, temp)

# 3. Immune compartment ----
##  3.1. Detailed cell typing ----
dir.create("sub1_Immune")

immune.obj <- subset(blca.obj, subset = Cell_Type_0 == "Immune") # 122993 cells

immune.obj <- immune.obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub1_Immune/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(immune.obj, ndims = 50)
dev.off()

immune.obj <- FindNeighbors(immune.obj, dims = 1:30)
immune.obj <- FindClusters(immune.obj, resolution = 1.0)

immune.obj <- RunUMAP(immune.obj, dims = 1:30)

pdf("sub1_Immune/UMAP.Clusters.pdf", 7, 6)
DimPlot(immune.obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub1_Immune/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(immune.obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

table(immune.obj$Cell_Type_SingleR)
pdf("sub1_Immune/UMAP.SingleR.pdf", 7, 6)
DimPlot(immune.obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

pdf("sub1_Immune/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(immune.obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

marker.genes[["Immune"]] <- c(
    "Ptprc",
    "Cd3d", "Cd3g", "Cd4", "Foxp3",
    "Cd8a", "Gzma", "Gzmb", "Havcr2", "Pdcd1",
    "Nkg7", "Ncam1", "Klrb1c",
    "Cd79a", "Ms4a1", "Igkc", "Sdc1",
    "Cd14", "Fcgr3", "Lyz2", "Ly6c1", # monocyte
    "Cd68", "Cd163", "Mrc1", "Adgre1", # macrophage
    "Itgax", "Fcer1a", "Cst3", # DC
    "Siglech", "Flt3", # pDC
    "Itgam", "Ly6g", "Ly6c2", # MDSC
    "Cxcr2", "S100a8", "S100a9", "Trem1", # neutrophils
    "Kit", "Hdc", "Ms4a2", # mast cell
    "Mki67", "Top2a" # cycling
)
length(marker.genes[["Immune"]])
setdiff(marker.genes[["Immune"]], rownames(immune.obj))
which(duplicated(marker.genes[["Immune"]]))

p <- FeaturePlot(immune.obj,
    features = marker.genes[["Immune"]],
    ncol = 6, cols = c("gray", "red"), raster = TRUE
)
pdf("sub1_Immune/UMAP.MarkerGenes.pdf", 24, 24)
plot(p)
dev.off()
rm(p)

pdf("sub1_Immune/DotPlot.MarkerGenes.pdf", 9, 7)
DotPlot(immune.obj, features = marker.genes[["Immune"]]) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

# undefined cluster markers
marker.undefined.C9 <- FindMarkers(immune.obj,
    ident.1 = 9,
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker.undefined.C9 <- marker.undefined.C9[order(marker.undefined.C9$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker.undefined.C9), 50) # Neurogenic markers: Ngf, Ncam1, Sox5, Tenm4

marker.undefined.C20 <- FindMarkers(immune.obj,
    ident.1 = c(20),
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker.undefined.C20 <- marker.undefined.C20[order(marker.undefined.C20$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker.undefined.C20), 50) #  Macrophage markers: C1qa Ms4a7 Ccl8

marker.undefined.C22 <- FindMarkers(immune.obj,
    ident.1 = 22,
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker.undefined.C22 <- marker.undefined.C22[order(marker.undefined.C22$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker.undefined.C22), 100) #  markers:

immune.obj$Cell_Type_2 <- "none"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(12, 16))] <- "T_CD4"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(29))] <- "T_CD4_Cycling"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(10, 33))] <- "Treg"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(0, 1, 18, 21, 28))] <- "T_CD8"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(11))] <- "T_CD8_Cycling"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(22))] <- "T_CD8_Neutrophil_like"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(13))] <- "NK"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(27))] <- "B"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(3, 4, 6, 8, 14, 19, 20))] <- "Macrophage"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(7, 24))] <- "Monocyte"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(15, 23, 26, 32))] <- "cDC"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(30))] <- "pDC"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(2, 5, 25))] <- "Neurophil"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(9))] <- "Myeloid_Neuro_like"
immune.obj$Cell_Type_2[which(immune.obj$seurat_clusters %in% c(17, 31))] <- "Myeloid_Cycling"

table(immune.obj$seurat_clusters, immune.obj$Cell_Type_2)

immune.obj$Cell_Type_2 <- factor(immune.obj$Cell_Type_2, levels = c(
    "T_CD4", "T_CD4_Cycling", "Treg",
    "T_CD8", "T_CD8_Cycling", "T_CD8_Neutrophil_like",
    "NK",
    "B",
    "Macrophage", "Monocyte",
    "cDC", "pDC",
    "Neurophil", "Myeloid_Neuro_like", "Myeloid_Cycling"
))
mycolor$Cell_Type_2 <- colorRampPalette(ggsci::pal_npg()(9))(length(levels(immune.obj$Cell_Type_2))) %>%
    setNames(levels(immune.obj$Cell_Type_2))

pdf("sub1_Immune/UMAP.Cell_Type_2.pdf", 13, 6)
DimPlot(immune.obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE,
    cols = mycolor$Cell_Type_2
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2") + NoLegend() +
    DimPlot(immune.obj,
        reduction = "umap", group.by = "Cell_Type_2", label = FALSE, raster = TRUE,
        cols = mycolor$Cell_Type_2
    ) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

## 3.2. Marker gene expression ----
Idents(immune.obj) <- immune.obj$Cell_Type_2
levels(Idents(immune.obj))

gene.selected <- c(
    "Cd3d", "Cd4", "Foxp3", "Cd8a", "Klrb1c", "Cd79a",
    "Cd68", "Mrc1", "Lrc2", "Cst3", "Siglech", "S100a8", "Hdc",
    "Ngf", "Dclk1", "Ncam1", "Mki67"
)

pdf("sub1_Immune/DotPlot.Cell_Type_2.pdf", 7, 4)
DotPlot(immune.obj, features = gene.selected) +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

rm(list = grep("marker.undefined", ls(), value = TRUE))
rm(gene.selected)

## 3.3. Percentage changes ----
# barplot composition
plot.data <- as.data.frame(table(immune.obj$Cell_Type_2, immune.obj$Condition))
colnames(plot.data) <- c("Cell_Type_2", "Condition", "Count")

pdf("sub1_Immune/Barplot.Cell_Type_2.Condition.pdf", 9, 5)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = Cell_Type_2, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot.data, aes(x = Condition, y = Count, fill = Cell_Type_2)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot.data <- as.data.frame(table(immune.obj$Cell_Type_2, immune.obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(immune.obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_60HDA"))

pdf("sub1_Immune/Boxplot.Fraction.pdf", 8, 7)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_60HDA"), c("CTRL", "LC_60HDA"))
    ) +
    facet_wrap(~Cell_Type_2, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction within immune compartment (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot.data, temp)

## 3.4. Differential analysis ----
dir.create("sub1_Immune/DEG")

unique(sample.info$Condition)
DEG.immune.Capsaicin <- list()
DEG.immune.LC_60HDA <- list()

for (cell_type in levels(immune.obj$Cell_Type_2)) {
    message(Sys.time(), ": Processing: ", cell_type)
    immune.sub <- subset(immune.obj, subset = Cell_Type_2 == cell_type)

    temp <- FindMarkers(immune.sub,
        ident.1 = "Capsaicin", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0.25, logfc.threshold = 0.25,
        max.cells.per.ident = 1000
    )
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Cell_Type_2 <- cell_type
    DEG.immune.Capsaicin[[cell_type]] <- temp

    temp <- FindMarkers(immune.sub,
        ident.1 = "LC_60HDA", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0.25, logfc.threshold = 0.25,
        max.cells.per.ident = 1000
    )
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Cell_Type_2 <- cell_type
    DEG.immune.LC_60HDA[[cell_type]] <- temp
}
rm(immune.sub, cell_type, temp)

lapply(names(DEG.immune.Capsaicin), function(x) {
    write.csv(DEG.immune.Capsaicin[[x]],
        file = paste0("sub1_Immune/DEG/DEG_Capsaicin_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})
lapply(names(DEG.immune.LC_60HDA), function(x) {
    write.csv(DEG.immune.LC_60HDA[[x]],
        file = paste0("sub1_Immune/DEG/DEG_LC_60HDA_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})

plot.data <- rbind(
    do.call(rbind, DEG.immune.Capsaicin) %>% mutate(Comparison = "Capsaicin_vs_CTRL"),
    do.call(rbind, DEG.immune.LC_60HDA) %>% mutate(Comparison = "LC_60HDA_vs_CTRL")
)
plot.data$DEG <- "Stable"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC > 0.25)] <- "Up"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC < -0.25)] <- "Down"
table(plot.data$Cell_Type_2, plot.data$Comparison, plot.data$DEG)

plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(immune.obj$Cell_Type_2))
pdf("sub1_Immune/DEG/Volcano.Capsaicin_vs_CTRL.pdf", 9, 7)
ggplot(plot.data %>% filter(Comparison == "Capsaicin_vs_CTRL"), aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Cell_Type_2), alpha = 0.5, size = 0.5) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    facet_wrap(~Cell_Type_2, ncol = 5) +
    theme_classic() +
    labs(title = "Volcano plot of DEGs in each cell type") +
    xlab("Average log2 fold change") +
    ylab("-log10 adjusted p value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    theme(legend.position = "none")
dev.off()

pdf("sub1_Immune/DEG/Volcano.LC_60HDA_vs_CTRL.pdf", 9, 7)
ggplot(plot.data %>% filter(Comparison == "LC_60HDA_vs_CTRL"), aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Cell_Type_2), alpha = 0.5, size = 0.5) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    facet_wrap(~Cell_Type_2, ncol = 5) +
    theme_classic() +
    labs(title = "Volcano plot of DEGs in each cell type") +
    xlab("Average log2 fold change") +
    ylab("-log10 adjusted p value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    theme(legend.position = "none")
dev.off()
rm(plot.data)

saveRDS(immune.obj, file = "sub1_Immune/immune.obj.rds")
rm(immune.obj)

# 4. Epithelial compartment ----
##  4.1. Clustering of cell states ----
dir.create("sub2_Cancer")
cancer.obj <- subset(blca.obj, subset = Cell_Type_1 %in% c("Cancer")) # 136537 cells

cancer.obj <- cancer.obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub2_Cancer/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(cancer.obj, ndims = 50)
dev.off()

cancer.obj <- RunUMAP(cancer.obj, dims = 1:30)

cancer.obj <- FindNeighbors(cancer.obj, dims = 1:30)
cancer.obj <- FindClusters(cancer.obj, resolution = 0.3)

library(ggraph)
pdf("sub2_Cancer/ClusteringTree.pdf", 8, 6)
clustree::clustree(cancer.obj, prefix = "RNA_snn_res.")
dev.off()

cancer.obj$seurat_clusters <- cancer.obj$RNA_snn_res.0.3
pdf("sub2_Cancer/UMAP.Clusters.pdf", 7, 6)
DimPlot(cancer.obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Cancer/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(cancer.obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Cancer/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(cancer.obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

table(cancer.obj$Cell_Type_SingleR)
pdf("sub2_Cancer/UMAP.SingleR.pdf", 7, 6)
DimPlot(cancer.obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

## 4.2. Cell states annotation ----
DEG.cancer.clusters <- FindAllMarkers(cancer.obj,
    group.by = "seurat_clusters",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
DEG.cancer.clusters <- DEG.cancer.clusters[order(
    DEG.cancer.clusters$cluster,
    0 - DEG.cancer.clusters$avg_log2FC
), ]

split(DEG.cancer.clusters$gene, DEG.cancer.clusters$cluster) %>% lapply(head, 50)

FeaturePlot(cancer.obj,
    features = c("Ptprc", "Krt18", "Krt8"),
    ncol = 2, cols = c("gray", "red"), raster = TRUE
)

cancer.obj$Cell_Type_2 <- "none"
cancer.obj$Cell_Type_2[which(cancer.obj$seurat_clusters %in% c(3, 5, 9))] <- "C1: Proliferative (Cell Cycle)"
cancer.obj$Cell_Type_2[which(cancer.obj$seurat_clusters %in% c(1, 4))] <- "C2: Quiescent & Stress-Adapted"
cancer.obj$Cell_Type_2[which(cancer.obj$seurat_clusters %in% c(0))] <- "C3: Basal & Immunogenic"
cancer.obj$Cell_Type_2[which(cancer.obj$seurat_clusters %in% c(6))] <- "C4: Invasive (Neurogenic/EMT)"
cancer.obj$Cell_Type_2[which(cancer.obj$seurat_clusters %in% c(2, 7, 8))] <- "C5: Evasive (Immune Mimic)"
table(cancer.obj$seurat_clusters, cancer.obj$Cell_Type_2)

cancer.obj$Cell_Type_2 <- factor(cancer.obj$Cell_Type_2, levels = c(
    "C1: Proliferative (Cell Cycle)", "C2: Quiescent & Stress-Adapted", "C3: Basal & Immunogenic",
    "C4: Invasive (Neurogenic/EMT)", "C5: Evasive (Immune Mimic)"
))

mycolor$Cell_Type_2 <- mycolor$Cell_Type_2[1:15]
mycolor$Cell_Type_2 <- c(
    mycolor$Cell_Type_2,
    setNames(ggsci::pal_jco()(5), levels(cancer.obj$Cell_Type_2))
)

pdf("sub2_Cancer/UMAP.Cell_Type_2.pdf", 8, 6)
DimPlot(cancer.obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE,
    cols = mycolor$Cell_Type_2
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

## 4.3. State marker genes ----
DEG.cancer.states <- FindAllMarkers(cancer.obj,
    group.by = "Cell_Type_2",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)

DEG.cancer.states <- DEG.cancer.states[order(
    DEG.cancer.states$cluster,
    0 - DEG.cancer.states$avg_log2FC
), ]

split(DEG.cancer.states$gene, DEG.cancer.states$cluster) %>% lapply(head, 50)

marker.genes[["Cancer"]] <- c(
    # C1: Proliferative (Cell Cycle)
    "Cenpf", "Cenpe", "Ttk", "Aurkb", "Brca1", "Exo1", "Polq", "Cdc20", "Cdc25c", "Mxd3",
    # C2: Quiescent & Stress-Adapted
    "Slc2a1", "Hmox1", "Gsta1", "Cryab", "Gdf15", "Adm", "Ftl1", "Hspb1", "Pgk1", "Ddit3",
    # C3: Basal & Immunogenic
    "H2-Ab1", "H2-Eb1", "H2-Aa", "C3", "C1ra", "C2", "Ager", "Lrrc26", "Pltp", "Il18",
    # C4: Invasive (Neurogenic/EMT)
    "Dclk1", "Ntrk3", "Robo2", "Efna5", "Sox5", "Lama3", "Adamts6", "Pdgfc", "Wwox", "Pbx1",
    # C5: Evasive (Immune Mimic)
    "Itk", "Fyn", "S100a9", "Ccl4", "C1qc", "Apoe", "Ptprc", "Itga4", "Rgs1", "Il1b"
)

length(marker.genes[["Cancer"]])
setdiff(marker.genes[["Cancer"]], rownames(cancer.obj))

cell.selected <- split(cancer.obj$Cell, cancer.obj$Cell_Type_2) %>%
    lapply(function(x) {
        sample(x, 500)
    }) %>%
    unlist() %>%
    unique()
pdf("sub2_Cancer/Heatmap.State.marker.pdf", 8, 7)
DoHeatmap(cancer.obj %>% ScaleData(., features = marker.genes[["Cancer"]]),
    cells = cell.selected,
    features = marker.genes[["Cancer"]],
    group.by = "Cell_Type_2",
    group.colors = mycolor$Cell_Type_2,
    size = 3,
    raster = TRUE
) +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "RdYlBu") %>% rev()) +
    theme(axis.text.x = element_blank())
dev.off()
rm(cell.selected)

gene.selected <- c(
    "Cenpf", "Cdk1", "Brca1", "Exo1", "Cdc20", # C1: Proliferative (Cell Cycle)
    "Slc2a1", "Hmox1", "Gdf15", "Ddit3", "Ftl1", # C2: Quiescent & Stress-Adapted
    "H2-Ab1", "Ager", "C3", "C1ra", "Il18", # C3: Basal & Immunogenic
    "Dclk1", "Ntrk3", "Robo2", "Lama3", "Adamts6", # C4: Invasive (Neurogenic/EMT)
    "S100a9", "C1qc", "Ptprc", "Gzmb", "Ccl4" # C5: Evasive (Immune Mimic)
)

pdf("sub2_Cancer/DotPlot.State.marker.pdf", 9, 3)
DotPlot(cancer.obj, features = gene.selected, group.by = "Cell_Type_2") +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sub2_Cancer/Violin.State.marker.pdf", 14, 18)
VlnPlot(cancer.obj,
    features = gene.selected, group.by = "Cell_Type_2",
    cols = mycolor$Cell_Type_2, ncol = 5, pt.size = 0
)
dev.off()
rm(gene.selected)

## 4.4. Percentage changes ----
# barplot composition
plot.data <- as.data.frame(table(cancer.obj$Cell_Type_2, cancer.obj$Condition))
colnames(plot.data) <- c("Cell_Type_2", "Condition", "Count")

pdf("sub2_Cancer/Barplot.Cell_Type_2.Condition.pdf", 9, 5)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = Cell_Type_2, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot.data, aes(x = Condition, y = Count, fill = Cell_Type_2)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot.data <- as.data.frame(table(cancer.obj$Cell_Type_2, cancer.obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(cancer.obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_60HDA"))

pdf("sub2_Cancer/Boxplot.Fraction.pdf", 8, 2.8)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.value", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_60HDA"), c("CTRL", "LC_60HDA"))
    ) +
    facet_wrap(~Cell_Type_2, ncol = 5, scales = "free") +
    labs(title = "Cell State composition across Samples") +
    ylab("Fraction of cancer cells (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot.data, temp)

## 4.5. DEG cross conditions ----
DEG.cancer.conditions <- list()

table(cancer.obj$Condition)

temp <- FindMarkers(cancer.obj,
    ident.1 = "Capsaicin", ident.2 = "CTRL",
    group.by = "Condition",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 10000
)

temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
temp$gene <- rownames(temp)
temp$Condition <- "Capsaicin"
DEG.cancer.conditions[["Capsaicin"]] <- temp
DEG.cancer.conditions[["Capsaicin"]]["Tubb3", ]
head(DEG.cancer.conditions[["Capsaicin"]])

temp <- FindMarkers(cancer.obj,
    ident.1 = "LC_60HDA", ident.2 = "CTRL",
    group.by = "Condition",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 10000
)

temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
temp$gene <- rownames(temp)
temp$Condition <- "LC_60HDA"
DEG.cancer.conditions[["LC_60HDA"]] <- temp

lapply(names(DEG.cancer.conditions), function(x) {
    write.csv(DEG.cancer.conditions[[x]],
        file = paste0("sub2_Cancer/DEG_", x, "_vs_CTRL_in_.csv"),
        row.names = FALSE
    )
})
rm(temp)

DEG.cancer.conditions[["LC_60HDA"]][c("Tubb3", "Ngf"), ]

# GO enrichment
library(clusterProfiler)
library(org.Mm.eg.db)

GO.res <- enrichGO(
    gene = DEG.cancer.conditions[["Capsaicin"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)

pdf("sub2_Cancer/GO.Capsaicin_vs_CTRL.Up.pdf", 7, 6)
dotplot(GO.res, showCategory = 20, label_format = 60) +
    theme_classic() +
    scale_color_viridis_c() +
    labs(title = "Up-regulated genes in Capsaicin vs CTRL")
dev.off()

GO.res <- enrichGO(
    gene = DEG.cancer.conditions[["Capsaicin"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)

pdf("sub2_Cancer/GO.Capsaicin_vs_CTRL.Down.pdf", 7, 6)
dotplot(GO.res, showCategory = 20, label_format = 60) +
    theme_classic() +
    labs(title = "Down-regulated genes in Capsaicin vs CTRL")
dev.off()

GO.res <- enrichGO(
    gene = DEG.cancer.conditions[["LC_60HDA"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)

pdf("sub2_Cancer/GO.LC_60HDA_vs_CTRL.Up.pdf", 7, 6)
dotplot(GO.res, showCategory = 20, label_format = 60) +
    theme_classic() +
    scale_color_viridis_c() +
    labs(title = "Up-regulated genes in LC_60HDA vs CTRL")
dev.off()

GO.res <- enrichGO(
    gene = DEG.cancer.conditions[["LC_60HDA"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)

pdf("sub2_Cancer/GO.LC_60HDA_vs_CTRL.Down.pdf", 7, 6)
dotplot(GO.res, showCategory = 20, label_format = 60) +
    theme_classic() +
    labs(title = "Down-regulated genes in LC_60HDA vs CTRL")
dev.off()

gc()
saveRDS(cancer.obj, file = "sub2_Cancer/cancer.obj.rds")
rm(cancer.obj)

# 5. Stromal compartment ----
##  5.1. Detailed cell typing ----
dir.create("sub3_Stromal")

stromal.obj <- subset(blca.obj, subset = Cell_Type_0 == "Stromal") # 12146 cells

stromal.obj <- stromal.obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub3_Stromal/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(stromal.obj, ndims = 50)
dev.off()

stromal.obj <- FindNeighbors(stromal.obj, dims = 1:30)
stromal.obj <- FindClusters(stromal.obj, resolution = 1.0)

stromal.obj <- RunUMAP(stromal.obj, dims = 1:30)

pdf("sub3_Stromal/UMAP.Clusters.pdf", 7, 6)
DimPlot(stromal.obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub3_Stromal/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(stromal.obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

table(stromal.obj$Cell_Type_SingleR)
pdf("sub3_Stromal/UMAP.SingleR.pdf", 7, 6)
DimPlot(stromal.obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

pdf("sub3_Stromal/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(stromal.obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

marker.genes[["Stromal"]] <- c(
    "Pecam1", "Cdh5", "Vwf", # endothelial
    "Prox1", "Pdpn", "Lyve1", # lymphatic endothelial
    "Cd34", # vascular
    "Col1a1", "Col1a2", "Dcn", # fibroblast
    "Lum", "Pdgfra",  # iCAF
    "Pdgfrb", "Rgs5", "Myh11", # myCAF
    "Acta2", "Tagln", "Myl9", # smooth muscle
    "Tubb3", "Nefl", "Map2", # neuron
    "S100b", "Sox10", # Schwann cell
    "Fabp7", "P2rx7", # Satellite glial cells,
    "Mki67", "Top2a" # cycling
)
length(marker.genes[["Stromal"]])
setdiff(marker.genes[["Stromal"]], rownames(stromal.obj))
which(duplicated(marker.genes[["Stromal"]]))

p <- FeaturePlot(stromal.obj,
    features = marker.genes[["Stromal"]],
    ncol = 6, cols = c("gray", "red"), raster = TRUE
)
pdf("sub3_Stromal/UMAP.MarkerGenes.pdf", 24, 18)
plot(p)
dev.off()
rm(p)

pdf("sub3_Stromal/DotPlot.MarkerGenes.pdf", 9, 7)
DotPlot(stromal.obj, features = marker.genes[["Stromal"]]) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

marker.undefined.C24 <- FindMarkers(stromal.obj,
    ident.1 = 24,
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker.undefined.C24 <- marker.undefined.C24[order(marker.undefined.C24$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker.undefined.C24), 100) #  Adipocyte  markers: Adipoq Pck1

stromal.obj$Cell_Type_2 <- "none"
stromal.obj$Cell_Type_2[which(stromal.obj$seurat_clusters %in% c(4, 16, 20, 23))] <- "Vescular_Endothelial"
stromal.obj$Cell_Type_2[which(stromal.obj$seurat_clusters %in% c(13, 25))] <- "Lymphatic_Endothelial"
stromal.obj$Cell_Type_2[which(stromal.obj$seurat_clusters %in% c(0:3, 5:11, 14, 15, 17:19))] <- "iCAF"
stromal.obj$Cell_Type_2[which(stromal.obj$seurat_clusters %in% c(12))] <- "myCAF"
stromal.obj$Cell_Type_2[which(stromal.obj$seurat_clusters %in% c(21))] <- "Smooth_Muscle"
stromal.obj$Cell_Type_2[which(stromal.obj$seurat_clusters %in% c(22))] <- "Schwann"
stromal.obj$Cell_Type_2[which(stromal.obj$seurat_clusters %in% c(24))] <- "Adipocyte"

table(stromal.obj$seurat_clusters, stromal.obj$Cell_Type_2)

stromal.obj$Cell_Type_2 <- factor(stromal.obj$Cell_Type_2, levels = c(
    "Vescular_Endothelial", "Lymphatic_Endothelial",
    "iCAF", "myCAF", "Smooth_Muscle", "Schwann", "Adipocyte"
))

mycolor$Cell_Type_2 <- mycolor$Cell_Type_2[1:20]
mycolor$Cell_Type_2 <- c(
    mycolor$Cell_Type_2,
    setNames(ggsci::pal_lancet()(7), levels(stromal.obj$Cell_Type_2))
)

pdf("sub3_Stromal/UMAP.Cell_Type_2.pdf", 13, 6)
DimPlot(stromal.obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE,
    cols = mycolor$Cell_Type_2
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2") + NoLegend() +
    DimPlot(stromal.obj,
        reduction = "umap", group.by = "Cell_Type_2", label = FALSE, raster = TRUE,
        cols = mycolor$Cell_Type_2
    ) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

## 5.2. Marker gene expression ----
Idents(stromal.obj) <- stromal.obj$Cell_Type_2
levels(Idents(stromal.obj))

gene.selected <- c(
    "Cdh5", "Pecam1", "Prox1", "Lyve1", "Col1a2", "Lum",
    "Pdgfra", "Rgs5", "Acta2", "Myh11", "Sox10", "S100b",
    "Adipoq", "Pck1"
)

pdf("sub3_Stromal/DotPlot.Cell_Type_2.pdf", 6, 2.6)
DotPlot(stromal.obj, features = gene.selected) +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

rm(list = grep("marker.undefined", ls(), value = TRUE))
rm(gene.selected, p)

## 5.3. Percentage changes ----
# barplot composition
plot.data <- as.data.frame(table(stromal.obj$Cell_Type_2, stromal.obj$Condition))
colnames(plot.data) <- c("Cell_Type_2", "Condition", "Count")

pdf("sub3_Stromal/Barplot.Cell_Type_2.Condition.pdf", 9, 5)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = Cell_Type_2, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot.data, aes(x = Condition, y = Count, fill = Cell_Type_2)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot.data <- as.data.frame(table(stromal.obj$Cell_Type_2, stromal.obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(stromal.obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_60HDA"))

pdf("sub3_Stromal/Boxplot.Fraction.pdf", 6.5, 4.8)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_60HDA"), c("CTRL", "LC_60HDA"))
    ) +
    facet_wrap(~Cell_Type_2, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction within immune compartment (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot.data, temp)

## 5.4. Differential analysis ----
dir.create("sub3_Stromal/DEG")

unique(sample.info$Condition)
DEG.stromal.Capsaicin <- list()
DEG.stromal.LC_60HDA <- list()

for (cell_type in levels(stromal.obj$Cell_Type_2)) {
    message(Sys.time(), ": Processing: ", cell_type)
    stromal.sub <- subset(stromal.obj, subset = Cell_Type_2 == cell_type)

    temp <- FindMarkers(stromal.sub,
        ident.1 = "Capsaicin", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0.25, logfc.threshold = 0.25,
        max.cells.per.ident = 1000
    )
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Cell_Type_2 <- cell_type
    DEG.stromal.Capsaicin[[cell_type]] <- temp

    temp <- FindMarkers(stromal.sub,
        ident.1 = "LC_60HDA", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0.25, logfc.threshold = 0.25,
        max.cells.per.ident = 1000
    )
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Cell_Type_2 <- cell_type
    DEG.stromal.LC_60HDA[[cell_type]] <- temp
}
rm(stromal.sub, cell_type, temp)

lapply(names(DEG.stromal.Capsaicin), function(x) {
    write.csv(DEG.stromal.Capsaicin[[x]],
        file = paste0("sub3_Stromal/DEG/DEG_Capsaicin_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})
lapply(names(DEG.stromal.LC_60HDA), function(x) {
    write.csv(DEG.stromal.LC_60HDA[[x]],
        file = paste0("sub3_Stromal/DEG/DEG_LC_60HDA_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})

plot.data <- rbind(
    do.call(rbind, DEG.stromal.Capsaicin) %>% mutate(Comparison = "Capsaicin_vs_CTRL"),
    do.call(rbind, DEG.stromal.LC_60HDA) %>% mutate(Comparison = "LC_60HDA_vs_CTRL")
)
plot.data$DEG <- "Stable"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC > 0.25)] <- "Up"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC < -0.25)] <- "Down"
table(plot.data$Cell_Type_2, plot.data$Comparison, plot.data$DEG)

plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(stromal.obj$Cell_Type_2))
pdf("sub3_Stromal/DEG/Volcano.Capsaicin_vs_CTRL.pdf", 7.5, 4.5)
ggplot(plot.data %>% filter(Comparison == "Capsaicin_vs_CTRL"), aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Cell_Type_2), alpha = 0.5, size = 0.5) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    facet_wrap(~Cell_Type_2, ncol = 4) +
    theme_classic() +
    labs(title = "Volcano plot of DEGs in each cell type") +
    xlab("Average log2 fold change") +
    ylab("-log10 adjusted p value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    theme(legend.position = "none")
dev.off()

pdf("sub3_Stromal/DEG/Volcano.LC_60HDA_vs_CTRL.pdf", 7.5, 4.5)
ggplot(plot.data %>% filter(Comparison == "LC_60HDA_vs_CTRL"), aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Cell_Type_2), alpha = 0.5, size = 0.5) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    facet_wrap(~Cell_Type_2, ncol = 4) +
    theme_classic() +
    labs(title = "Volcano plot of DEGs in each cell type") +
    xlab("Average log2 fold change") +
    ylab("-log10 adjusted p value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    theme(legend.position = "none")
dev.off()
rm(plot.data)

saveRDS(stromal.obj, file = "sub3_Stromal/stromal.obj.rds")
rm(stromal.obj)

# save data ----
save.image("01.Cell_Type.RData")
