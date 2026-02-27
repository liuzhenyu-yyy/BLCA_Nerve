# local
setwd("/Users/zhenyuliu/Research/Project/BLCA_Nerve/results/01.Cell_Type")
# server
setwd("/thinker/tb269store/user/liuzhenyu/Project/BLCA_Nerve/results/01.Cell_Type")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(parallel)

load("01.Cell_Type.RData")

# 1. Load data & QC ----
data_dir <- "/Users/zhenyuliu/Research/Project/BLCA_Nerve/data/CellRanger_Out"
samples <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
mycolor <- list()

## 1.1. Load data ----
# counts
rc <- list()
rc <- lapply(samples, function(x) {
    data <- Read10X(data_dir = file.path(data_dir, x, "filtered_feature_bc_matrix"))
    colnames(data) <- paste(x, colnames(data), sep = "_")
    data
})
lapply(rc, dim)

common_genes <- Reduce(intersect, lapply(rc, rownames))
rc <- lapply(rc, function(x) {
    x[common_genes, ]
})

rc <- do.call(cbind, rc)
saveRDS(rc, file = "../../data/rc.rds")
rc <- readRDS(file = "../../data/rc.rds")
rm(samples)

# metadata
rc[1:5, 1:5]
sample_info <- data.frame(
    row.names = colnames(rc),
    "Cell" = colnames(rc),
    "Sample" = gsub("_A.+?*$|_T.+?*$|_C.+?*$|_G.+?*$", "", colnames(rc))
) %>% mutate(
    "Condition" = gsub("_00.+?*$|_1.+?*$|_2.+?*$|_3.+?*$|1$|2$|3$", "", Sample)
)
sample_info$Sample <- gsub("_00", "_", sample_info$Sample)
sample_info$Sample <- gsub("60HDA", "6OHDA", sample_info$Sample)
sample_info$Sample <- gsub("RIme", "Rime_", sample_info$Sample)

sample_info$Condition <- gsub("60HDA", "6OHDA", sample_info$Condition)
sample_info$Condition <- gsub("RIme", "Rime", sample_info$Condition)

sample_info$Condition <- factor(sample_info$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA", "anti_PD1", "Rime"))
sample_info <- sample_info[order(sample_info$Condition), ]
sample_info$Sample <- factor(sample_info$Sample, levels = unique(sample_info$Sample))
table(sample_info$Condition)
table(sample_info$Sample)

rc <- rc[, sample_info$Cell]

# singleR_anno annotation
singler_anno_1 <- read.table("../../data/scRNAseq_1/cell_annot_MouseRNAseqData_label.main.csv",
    header = TRUE, sep = ","
)
colnames(singler_anno_1)
singler_anno_2 <- read.table("../../data/scRNAseq_2/cell_annot_MouseRNAseqData_label.main.csv",
    header = TRUE, sep = ","
)
colnames(singler_anno_2)

singler_anno <- rbind(singler_anno_1, singler_anno_2[, colnames(singler_anno_1)])
rm(singler_anno_1, singler_anno_2)

singler_anno$Cell <- gsub("-.+?*$", "-1", singler_anno$Cell)
singler_anno$Cell <- paste(singler_anno$orig.ident, singler_anno$Cell, sep = "_")
rownames(singler_anno) <- singler_anno$Cell

table(singler_anno$Cell %in% sample_info$Cell)
table(sample_info$Cell %in% singler_anno$Cell)
singler_anno[1:5, ]

## 1.2. QC for cells ----
sample_info$Cell_Type_SingleR <- "none"
sample_info[singler_anno$Cell, "Cell_Type_SingleR"] <- singler_anno$raw_main_celltype
sample_info[1:5, ]

blca_obj <- CreateSeuratObject(
    counts = rc,
    project = "BLCA_Neuro",
    meta.data = sample_info
)
rm(singler_anno, rc)

blca_obj$Percent_Mito <- PercentageFeatureSet(blca_obj, pattern = "^mt-")
quantile(blca_obj$Percent_Mito)

sample_info <- blca_obj@meta.data

mycolor$Condition <- ggsci::pal_npg("nrc")(5)[c(3, 2, 1, 4, 5)] %>%
    setNames(levels(sample_info$Condition))

pdf("Hist.QC.Metrics.pdf", 6, 4)
hist(blca_obj$nFeature_RNA, breaks = 100)
hist(log10(blca_obj$nCount_RNA), breaks = 100)
hist(blca_obj$Percent_Mito, breaks = 100)
dev.off()

cutoffs <- c(
    "nFeature_RNA" = 300,
    "nCount_RNA" = 1000,
    "Percent_Mito" = 10
)

pdf("Scatter.QC.pdf", 20, 12)
ggplot(sample_info, aes(x = nCount_RNA, y = nFeature_RNA, color = Condition)) +
    geom_point(show.legend = FALSE) +
    scale_color_manual(values = mycolor$Condition) +
    theme_classic() +
    geom_vline(xintercept = cutoffs["nCount_RNA"], linetype = "dashed", color = "red") +
    geom_hline(yintercept = cutoffs["nFeature_RNA"], linetype = "dashed", color = "red") +
    facet_wrap(~Condition, ncol = 3) +
    labs(title = "Scatter plot of QC metrics")
dev.off()

pdf("Hist.QC.Metrics.pdf", 15, 4)
patchwork::wrap_plots(
    ggplot(sample_info, aes(x = log10(nCount_RNA))) +
        geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
        geom_vline(xintercept = log10(cutoffs["nCount_RNA"]), linetype = "dashed", color = "red") +
        theme_classic() +
        labs(title = "Histogram of log10 UMI counts"),
    ggplot(sample_info, aes(x = nFeature_RNA)) +
        geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
        geom_vline(xintercept = cutoffs["nFeature_RNA"], linetype = "dashed", color = "red") +
        theme_classic() +
        labs(title = "Histogram of gene counts"),
    ggplot(sample_info, aes(x = Percent_Mito)) +
        geom_histogram(bins = 100, position = "identity", alpha = 0.5) +
        geom_vline(xintercept = cutoffs["Percent_Mito"], linetype = "dashed", color = "red") +
        theme_classic() +
        labs(title = "Histogram of mitochondrial percentage"),
    ncol = 3
)
dev.off()

blca_obj <- subset(blca_obj, subset = nFeature_RNA > cutoffs["nFeature_RNA"] &
    nCount_RNA > cutoffs["nCount_RNA"] &
    Percent_Mito < cutoffs["Percent_Mito"])
sample_info_all <- sample_info

## 1.3. Reductions & clustering -----
blca_obj <- blca_obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(blca_obj, ndims = 50)
dev.off()

blca_obj <- RunUMAP(blca_obj, dims = 1:30)

pdf("UMAP.Condition.pdf", 5, 4)
DimPlot(blca_obj, reduction = "umap", group.by = "Condition") +
    scale_color_manual(values = mycolor$Condition) +
    theme_classic() +
    labs(title = "UMAP plot colored by Condition")
dev.off()

blca_obj <- FindNeighbors(blca_obj, dims = 1:30)
blca_obj <- FindClusters(blca_obj, resolution = 0.8)

mycolor$Cluster <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(blca_obj$seurat_clusters))) %>%
    setNames(levels(blca_obj$seurat_clusters))

Idents(blca_obj) <- blca_obj$seurat_clusters
pdf("UMAP.Clusters.pdf", 5.5, 4)
DimPlot(blca_obj, reduction = "umap", label = TRUE, label.size = 3) +
    scale_color_manual(values = mycolor$Cluster) +
    theme_classic() +
    guides(color = guide_legend(ncol = 3)) +
    labs(title = "UMAP plot colored by Clusters")
dev.off()

pdf("UMAP.Clusters.split.Condition.pdf", 12, 8)
DimPlot(blca_obj, reduction = "umap", label = TRUE, split.by = "Condition", ncol = 3) +
    scale_color_manual(values = mycolor$Cluster) +
    theme_classic() +
    labs(title = "UMAP plot colored by Clusters")
dev.off()

## 1.4. QC metrics ----
blca_obj$Percent_Mito
pdf("VlnPlot.QC.Cluster.pdf", 8, 8)
VlnPlot(blca_obj,
    group.by = "seurat_clusters", features = c("nFeature_RNA", "nCount_RNA", "Percent_Mito"),
    ncol = 1, pt.size = 0, cols = mycolor$Cluster
) +
    theme_classic()
dev.off()

pdf("UMAP.QC.Cluster.pdf", 8, 6)
FeaturePlot(blca_obj,
    features = c("nFeature_RNA", "nCount_RNA", "Percent_Mito"),
    ncol = 2, cols = c("gray", "red")
)
dev.off()

# 2. Initial cell type annotation ----
## 2.1. SingleR annotation ----
table(blca_obj$Cell_Type_SingleR)
pdf("UMAP.SingleR.pdf", 7, 6)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

## 2.2. Expression of marker genes ----
marker_genes <- list()

marker_genes[["All"]] <- c(
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
setdiff(marker_genes$All, rownames(blca_obj))
marker_genes$All[which(duplicated(marker_genes$All))]

p <- FeaturePlot(blca_obj, features = marker_genes$All, ncol = 6, cols = c("gray", "red"))
pdf("UMAP.MarkerGenes.pdf", 24, 24)
plot(p)
dev.off()
rm(p)

pdf("DotPlot.MarkerGenes.pdf", 9, 8)
DotPlot(blca_obj, features = marker_genes$All) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

## 2.3. Label cell types ----
blca_obj$Cell_Type_1 <- "none"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(23))] <- "Bladder_epithelial"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(0:2, 7:9, 11, 15, 21, 29))] <- "Cancer"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(25, 36))] <- "Endothelial"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(4, 10, 16, 19, 30, 35))] <- "T"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(17))] <- "NK"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(31))] <- "B"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(5, 6, 13, 18, 22, 27, 28))] <- "Myeloid"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(3))] <- "Neutrophil"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(12, 24))] <- "Fibroblast"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(29, 30, 32, 33, 34, 20, 21, 14, 26))] <- "Doublets"
table(blca_obj$seurat_clusters, blca_obj$Cell_Type_1)

blca_obj$Cell_Type_1 <- factor(blca_obj$Cell_Type_1, levels = c(
    "Bladder_epithelial", "Cancer",
    "Endothelial", "Fibroblast",
    "T", "B", "NK",
    "Neutrophil", "Myeloid", "Doublets"
))

Idents(blca_obj) <- "Cell_Type_1"

mycolor$Cell_Type_1 <- c(ggsci::pal_npg()(length(levels(blca_obj$Cell_Type_1)))) %>%
    setNames(levels(blca_obj$Cell_Type_1))
mycolor$Cell_Type_1["Doublets"] <- "gray"

pdf("UMAP.Cell_Type_1.pdf", 5.5, 4)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_1", label = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

pdf("UMAP.Cell_Type_1.split.Condition.pdf", 16, 6)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_1", label = TRUE, split.by = "Condition") +
    scale_color_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

save.image("01.Cell_Type.RData")
saveRDS(mycolor, file = "mycolor.rds")

gene_selected <- c(
    "Epcam", "Krt8", "Krt18", # epithelial
    "Pecam1", "Cdh5", # endothelial
    "Col1a2", "Dcn", "Pdgfra", # fibroblast
    "Ptprc", # leukocyte
    "Cd3d", "Cd3g", "Cd79a", "Klrb1c", # lymphocyte
    "S100a8", "S100a9", # neutrophils
    "Cd14", "Cd68" # myeloid
)
blca_obj$Cell_Type_1
pdf("DotPlot.Cell_Type_1.pdf", 7, 3)
DotPlot(blca_obj, features = gene_selected, group.by = "Cell_Type_1") +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
rm(gene_selected)

# bigger categories
blca_obj$Cell_Type_0 <- "none"
blca_obj$Cell_Type_0[which(blca_obj$Cell_Type_1 %in% c("Bladder_epithelial", "Cancer"))] <- "Epithelial"
blca_obj$Cell_Type_0[which(blca_obj$Cell_Type_1 %in% c("Fibroblast", "Endothelial"))] <- "Stromal"
blca_obj$Cell_Type_0[which(blca_obj$Cell_Type_1 %in% c("T", "B", "NK", "Myeloid", "Granulocyte"))] <- "Immune"
table(blca_obj$Cell_Type_0)

blca_obj$Cell_Type_0 <- factor(blca_obj$Cell_Type_0, levels = c("Epithelial", "Stromal", "Immune"))
mycolor$Cell_Type_0 <- c(ggsci::pal_jco()(4)[c(1, 2, 4)]) %>%
    setNames(levels(blca_obj$Cell_Type_0))

pdf("UMAP.Cell_Type_0.pdf", 7, 6)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_0", label = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_0")
dev.off()

pdf("UMAP.Cell_Type_0.split.Condition.pdf", 16, 6)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_0", label = TRUE, , split.by = "Condition") +
    scale_color_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_0")
dev.off()

## 2.4. Cell types cross conditions ----
plot.data <- as.data.frame(table(blca_obj$Cell_Type_1, blca_obj$Condition))
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

plot.data <- as.data.frame(table(blca_obj$Cell_Type_0, blca_obj$Condition))
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
plot.data <- as.data.frame(table(blca_obj$Cell_Type_0, blca_obj$Sample))
colnames(plot.data) <- c("Cell_Type_0", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(blca_obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))

pdf("Boxplot.Fraction.Cell_Type_0.pdf", 5, 3)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_0), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_6OHDA"), c("CTRL", "LC_6OHDA"))
    ) +
    facet_wrap(~Cell_Type_0, ncol = 3, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()

plot.data <- as.data.frame(table(blca_obj$Cell_Type_1, blca_obj$Sample))
colnames(plot.data) <- c("Cell_Type_1", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(blca_obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))

pdf("Boxplot.Fraction.Cell_Type_1.pdf", 8, 5)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_1), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_6OHDA"), c("CTRL", "LC_6OHDA"))
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

immune_obj <- subset(blca_obj, subset = Cell_Type_0 == "Immune") # 122993 cells

immune_obj <- immune_obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub1_Immune/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(immune_obj, ndims = 50)
dev.off()

immune_obj <- FindNeighbors(immune_obj, dims = 1:30)
immune_obj <- FindClusters(immune_obj, resolution = 1.0)

immune_obj <- RunUMAP(immune_obj, dims = 1:30)

pdf("sub1_Immune/UMAP.Clusters.pdf", 7, 6)
DimPlot(immune_obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub1_Immune/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(immune_obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

table(immune_obj$Cell_Type_SingleR)
pdf("sub1_Immune/UMAP.SingleR.pdf", 7, 6)
DimPlot(immune_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

pdf("sub1_Immune/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(immune_obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

marker_genes[["Immune"]] <- c(
    "Ptprc",
    "Cd3d", "Cd3g", "Cd4", "Foxp3",
    "Cd8a", "Gzma", "Gzmb", "Havcr2", "Pdcd1",
    "Nkg7", "Ncam1", "Klrb1c",
    "Cd79a", "Ms4a1", "Igkc", "Sdc1",
    "Cd14", "Fcgr3", "Lyz2", "Ly6c1", # monocyte
    "Cd68", "Cd163", "Mrc1", "Adgre1", # macrophage
    "Itgax", "Fcer1a", "Cst3", # DC
    "Siglech", "Flt3", # pDC
    "Sox5", "Ngf", "Il34", # MDSC
    "Cxcr2", "S100a8", "S100a9", "Trem1", # neutrophils
    "Kit", "Hdc", "Ms4a2", # mast cell
    "Mki67", "Top2a" # cycling
)
length(marker_genes[["Immune"]])
setdiff(marker_genes[["Immune"]], rownames(immune_obj))
which(duplicated(marker_genes[["Immune"]]))

p <- FeaturePlot(immune_obj,
    features = marker_genes[["Immune"]],
    ncol = 6, cols = c("gray", "red"), raster = TRUE
)
pdf("sub1_Immune/UMAP.MarkerGenes.pdf", 24, 24)
plot(p)
dev.off()
rm(p)

pdf("sub1_Immune/DotPlot.MarkerGenes.pdf", 9, 7)
DotPlot(immune_obj,
    features = marker_genes[["Immune"]],
    group.by = "seurat_clusters"
) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

# undefined cluster markers
marker_undefined_c9 <- FindMarkers(immune_obj,
    ident.1 = 9,
    group.by = "seurat_clusters",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker_undefined_c9 <- marker_undefined_c9[order(marker_undefined_c9$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker_undefined_c9), 50) # Neurogenic markers: Ngf, Ncam1, Sox5, Tenm4 MDSC?

marker_undefined_c20 <- FindMarkers(immune_obj,
    ident.1 = c(20),
    group.by = "seurat_clusters",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker_undefined_c20 <- marker_undefined_c20[order(marker_undefined_c20$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker_undefined_c20), 50) #  Macrophage markers: C1qa Ms4a7 Ccl8

marker_undefined_c22 <- FindMarkers(immune_obj,
    ident.1 = 22,
    group.by = "seurat_clusters",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker_undefined_c22 <- marker_undefined_c22[order(marker_undefined_c22$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker_undefined_c22), 200) #  markers: CD8T + Myeloid, suppressed

immune_obj$Cell_Type_2 <- "none"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(12, 16))] <- "T_CD4"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(29))] <- "T_CD4_Cycling"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(10, 33))] <- "Treg"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(0, 1, 18, 21, 28))] <- "T_CD8"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(11))] <- "T_CD8_Cycling"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(22))] <- "T_CD8_Suppressed"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(13))] <- "NK"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(27))] <- "B"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(3, 4, 6, 8, 14, 19, 20))] <- "Macrophage"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(7, 24))] <- "Monocyte"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(15, 23, 26, 32))] <- "cDC"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(30))] <- "pDC"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(2, 5, 25))] <- "Neutrophil"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(9))] <- "MDSC_Neurogenic"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(17, 31))] <- "Myeloid_Cycling"

table(immune_obj$seurat_clusters, immune_obj$Cell_Type_2)

immune_obj$Cell_Type_2 <- factor(immune_obj$Cell_Type_2, levels = c(
    "T_CD4", "T_CD4_Cycling", "Treg",
    "T_CD8", "T_CD8_Cycling", "T_CD8_Suppressed",
    "NK",
    "B",
    "Macrophage", "Monocyte",
    "cDC", "pDC",
    "Neutrophil", "MDSC_Neurogenic", "Myeloid_Cycling"
))
mycolor$Cell_Type_2 <- colorRampPalette(ggsci::pal_npg()(9))(length(levels(immune_obj$Cell_Type_2))) %>%
    setNames(levels(immune_obj$Cell_Type_2))

pdf("sub1_Immune/UMAP.Cell_Type_2.pdf", 13, 6)
DimPlot(immune_obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE,
    cols = mycolor$Cell_Type_2
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2") + NoLegend() +
    DimPlot(immune_obj,
        reduction = "umap", group.by = "Cell_Type_2", label = FALSE, raster = TRUE,
        cols = mycolor$Cell_Type_2
    ) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

pdf("sub1_Immune/UMAP.Cell_Type_2.split.pdf", 12, 4)
DimPlot(immune_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition",
    label = FALSE, raster = TRUE,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()
dev.off()

p <- DimPlot(immune_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition",
    label = FALSE, raster = TRUE, pt.size = NULL,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()

pdf("sub1_Immune/UMAP.Cell_Type_2.split.density.pdf", 12, 4)
plot(p + stat_density_2d(aes_string(x = "umap_1", y = "umap_2", fill = "after_stat(level)"),
    linewidth = 0, geom = "density_2d_filled",
    n = 150, h = c(1.2, 1.2)
) +
    scale_fill_viridis_c(option = "D") +
    guides(fill = guide_colorbar(title = "Density")) +
    theme(legend.position = "right"))
dev.off()
rm(p)

## 3.2. Marker gene expression ----
Idents(immune_obj) <- immune_obj$Cell_Type_2
levels(Idents(immune_obj))

gene_selected <- c(
    "Cd3d", "Cd4", "Foxp3", "Cd8a", "Prf1", "Arg2", "Klrb1c", "Cd79a",
    "Cd68", "Mrc1", "Lrc2", "Cst3", "Siglech", "S100a8", "Hdc",
    "Sox5", "Dclk1", "Ncam1", "Mki67"
)

pdf("sub1_Immune/DotPlot.Cell_Type_2.pdf", 7, 4)
DotPlot(immune_obj, features = gene_selected) +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

rm(list = grep("marker.undefined", ls(), value = TRUE))
rm(gene_selected)

## 3.3. Percentage changes ----
# barplot composition
plot.data <- as.data.frame(table(immune_obj$Cell_Type_2, immune_obj$Condition))
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
plot.data <- as.data.frame(table(immune_obj$Cell_Type_2, immune_obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(immune_obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))

pdf("sub1_Immune/Boxplot.Fraction.pdf", 8, 7)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_6OHDA"), c("CTRL", "LC_6OHDA"))
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

saveRDS(immune_obj, file = "sub1_Immune/immune_obj.rds")
rm(immune_obj)

# 4. Epithelial compartment ----
##  4.1. Clustering of cell states ----
dir.create("sub2_Cancer")

cancer_obj <- subset(blca_obj, subset = Cell_Type_1 %in% c("Cancer")) # 136537 cells

cancer_obj <- cancer_obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub2_Cancer/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(cancer_obj, ndims = 50)
dev.off()

cancer_obj <- RunUMAP(cancer_obj, dims = 1:30)

cancer_obj <- FindNeighbors(cancer_obj, dims = 1:30)
cancer_obj <- FindClusters(cancer_obj, resolution = 0.3)

library(ggraph)
pdf("sub2_Cancer/ClusteringTree.pdf", 8, 6)
clustree::clustree(cancer_obj, prefix = "RNA_snn_res.")
dev.off()

cancer_obj$seurat_clusters <- cancer_obj$RNA_snn_res.0.3
pdf("sub2_Cancer/UMAP.Clusters.pdf", 7, 6)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Cancer/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Cancer/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

table(cancer_obj$Cell_Type_SingleR)
pdf("sub2_Cancer/UMAP.SingleR.pdf", 7, 6)
DimPlot(cancer_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

## 4.2. Cell states annotation ----
deg_cancer_clusters <- FindAllMarkers(cancer_obj,
    group.by = "seurat_clusters",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
deg_cancer_clusters <- deg_cancer_clusters[order(
    deg_cancer_clusters$cluster,
    0 - deg_cancer_clusters$avg_log2FC
), ]

split(deg_cancer_clusters$gene, deg_cancer_clusters$cluster) %>% lapply(head, 50)

FeaturePlot(cancer_obj,
    features = c("Ptprc", "Krt18", "Krt8"),
    ncol = 2, cols = c("gray", "red"), raster = TRUE
)

cancer_obj$Cell_Type_2 <- "none"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(3, 5, 9))] <- "C1: Proliferative (Cell Cycle)"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(1, 4))] <- "C2: Quiescent & Stress-Adapted"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(0))] <- "C3: Basal & Immunogenic"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(6))] <- "C4: Invasive (Neurogenic/EMT)"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(2, 7, 8))] <- "C5: Evasive (Immune Mimic)"
table(cancer_obj$seurat_clusters, cancer_obj$Cell_Type_2)

cancer_obj$Cell_Type_2 <- factor(cancer_obj$Cell_Type_2, levels = c(
    "C1: Proliferative (Cell Cycle)", "C2: Quiescent & Stress-Adapted", "C3: Basal & Immunogenic",
    "C4: Invasive (Neurogenic/EMT)", "C5: Evasive (Immune Mimic)"
))

mycolor$Cell_Type_2 <- mycolor$Cell_Type_2[1:15]
mycolor$Cell_Type_2 <- c(
    mycolor$Cell_Type_2,
    setNames(ggsci::pal_jco()(5), levels(cancer_obj$Cell_Type_2))
)

pdf("sub2_Cancer/UMAP.Cell_Type_2.pdf", 8, 6)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE,
    cols = mycolor$Cell_Type_2
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

## 4.3. State marker genes ----
deg_cancer_states <- FindAllMarkers(cancer_obj,
    group.by = "Cell_Type_2",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)

deg_cancer_states <- deg_cancer_states[order(
    deg_cancer_states$cluster,
    0 - deg_cancer_states$avg_log2FC
), ]

split(deg_cancer_states$gene, deg_cancer_states$cluster) %>% lapply(head, 50)

marker_genes[["Cancer"]] <- c(
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

length(marker_genes[["Cancer"]])
setdiff(marker_genes[["Cancer"]], rownames(cancer_obj))

cell_selected <- split(cancer_obj$Cell, cancer_obj$Cell_Type_2) %>%
    lapply(function(x) {
        sample(x, 500)
    }) %>%
    unlist() %>%
    unique()
pdf("sub2_Cancer/Heatmap.State.marker.pdf", 8, 7)
DoHeatmap(cancer_obj %>% ScaleData(., features = marker_genes[["Cancer"]]),
    cells = cell_selected,
    features = marker_genes[["Cancer"]],
    group.by = "Cell_Type_2",
    group.colors = mycolor$Cell_Type_2,
    size = 3,
    raster = TRUE
) +
    scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "RdYlBu") %>% rev()) +
    theme(axis.text.x = element_blank())
dev.off()
rm(cell_selected)

gene_selected <- c(
    "Cenpf", "Cdk1", "Brca1", "Exo1", "Cdc20", # C1: Proliferative (Cell Cycle)
    "Slc2a1", "Hmox1", "Gdf15", "Ddit3", "Ftl1", # C2: Quiescent & Stress-Adapted
    "H2-Ab1", "Ager", "C3", "C1ra", "Il18", # C3: Basal & Immunogenic
    "Dclk1", "Ntrk3", "Robo2", "Lama3", "Adamts6", # C4: Invasive (Neurogenic/EMT)
    "S100a9", "C1qc", "Ptprc", "Gzmb", "Ccl4" # C5: Evasive (Immune Mimic)
)

pdf("sub2_Cancer/DotPlot.State.marker.pdf", 9, 3)
DotPlot(cancer_obj, features = gene_selected, group.by = "Cell_Type_2") +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("sub2_Cancer/Violin.State.marker.pdf", 14, 18)
VlnPlot(cancer_obj,
    features = gene_selected, group.by = "Cell_Type_2",
    cols = mycolor$Cell_Type_2, ncol = 5, pt.size = 0
)
dev.off()
rm(gene_selected)

## 4.4. Percentage changes ----
# barplot composition
plot.data <- as.data.frame(table(cancer_obj$Cell_Type_2, cancer_obj$Condition))
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
plot.data <- as.data.frame(table(cancer_obj$Cell_Type_2, cancer_obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(cancer_obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))

pdf("sub2_Cancer/Boxplot.Fraction.pdf", 8, 2.8)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.value", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_6OHDA"), c("CTRL", "LC_6OHDA"))
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
deg_cancer_conditions <- list()

table(cancer_obj$Condition)

temp <- FindMarkers(cancer_obj,
    ident.1 = "Capsaicin", ident.2 = "CTRL",
    group.by = "Condition",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 10000
)

temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
temp$gene <- rownames(temp)
temp$Condition <- "Capsaicin"
deg_cancer_conditions[["Capsaicin"]] <- temp
deg_cancer_conditions[["Capsaicin"]]["Tubb3", ]
head(deg_cancer_conditions[["Capsaicin"]])

temp <- FindMarkers(cancer_obj,
    ident.1 = "LC_6OHDA", ident.2 = "CTRL",
    group.by = "Condition",
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 10000
)

temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
temp$gene <- rownames(temp)
temp$Condition <- "LC_6OHDA"
deg_cancer_conditions[["LC_6OHDA"]] <- temp

lapply(names(deg_cancer_conditions), function(x) {
    write.csv(deg_cancer_conditions[[x]],
        file = paste0("sub2_Cancer/deg_", x, "_vs_CTRL_in_.csv"),
        row.names = FALSE
    )
})
rm(temp)

deg_cancer_conditions[["LC_6OHDA"]][c("Tubb3", "Ngf"), ]

# GO enrichment
library(clusterProfiler)
library(org.Mm.eg.db)

go_res <- enrichGO(
    gene = deg_cancer_conditions[["Capsaicin"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)
go_res <- simplify(go_res)

pdf("sub2_Cancer/GO_Capsaicin_vs_CTRL.Up.pdf", 7, 6)
dotplot(go_res, showCategory = 20, label_format = 60) +
    theme_classic() +
    scale_color_viridis_c() +
    labs(title = "Up-regulated genes in Capsaicin vs CTRL")
dev.off()
write.table(go_res@result, "sub2_Cancer/GO_Capsaicin_vs_CTRL.Up.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
)

go_res <- enrichGO(
    gene = deg_cancer_conditions[["Capsaicin"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)
go_res <- simplify(go_res)

pdf("sub2_Cancer/GO_Capsaicin_vs_CTRL.Down.pdf", 7, 6)
dotplot(go_res, showCategory = 20, label_format = 60) +
    theme_classic() +
    labs(title = "Down-regulated genes in Capsaicin vs CTRL")
dev.off()
write.table(go_res@result, "sub2_Cancer/GO_Capsaicin_vs_CTRL.Down.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
)

go_res <- enrichGO(
    gene = deg_cancer_conditions[["LC_6OHDA"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)
go_res <- simplify(go_res)

pdf("sub2_Cancer/GO_LC_6OHDA_vs_CTRL.Up.pdf", 7, 6)
dotplot(go_res, showCategory = 20, label_format = 60) +
    theme_classic() +
    scale_color_viridis_c() +
    labs(title = "Up-regulated genes in LC_6OHDA vs CTRL")
dev.off()
write.table(go_res@result, "sub2_Cancer/GO_LC_6OHDA_vs_CTRL.Up.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
)

go_res <- enrichGO(
    gene = deg_cancer_conditions[["LC_6OHDA"]] %>%
        filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
        pull(gene),
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
)
go_res <- simplify(go_res)

pdf("sub2_Cancer/GO_LC_6OHDA_vs_CTRL.Down.pdf", 7, 6)
dotplot(go_res, showCategory = 20, label_format = 60) +
    theme_classic() +
    labs(title = "Down-regulated genes in LC_6OHDA vs CTRL")
dev.off()
write.table(go_res@result, "sub2_Cancer/GO_LC_6OHDA_vs_CTRL.Down.txt",
    sep = "\t", row.names = FALSE, quote = FALSE
)
rm(go_res)

gc()
saveRDS(cancer_obj, file = "sub2_Cancer/cancer_obj.rds")
rm(cancer_obj)

# 5. Stromal compartment ----
##  5.1. Detailed cell typing ----
dir.create("sub3_Stromal")

stromal_obj <- subset(blca_obj, subset = Cell_Type_0 == "Stromal") # 12146 cells

stromal_obj <- stromal_obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
    ScaleData() %>%
    RunPCA(verbose = FALSE)

pdf("sub3_Stromal/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(stromal_obj, ndims = 50)
dev.off()

stromal_obj <- FindNeighbors(stromal_obj, dims = 1:30)
stromal_obj <- FindClusters(stromal_obj, resolution = 1.0)

stromal_obj <- RunUMAP(stromal_obj, dims = 1:30)

pdf("sub3_Stromal/UMAP.Clusters.pdf", 7, 6)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = mycolor$Cluster, raster = TRUE, pt.size = 1.5
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub3_Stromal/UMAP.Clusters.split.pdf", 12, 4)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE, pt.size = 1.5
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

table(stromal_obj$Cell_Type_SingleR)
pdf("sub3_Stromal/UMAP.SingleR.pdf", 7, 6)
DimPlot(stromal_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE, pt.size = 1.5) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

pdf("sub3_Stromal/UMAP.Cell_Type_1.pdf", 7, 6)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE, pt.size = 1.5
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

marker_genes[["Stromal"]] <- c(
    "Pecam1", "Cdh5", "Vwf", # endothelial
    "Prox1", "Pdpn", "Lyve1", # lymphatic endothelial
    "Cd34", # vascular
    "Col1a1", "Col1a2", "Dcn", # fibroblast
    "Lum", "Pdgfra", # iCAF
    "Itga11",
    "Pdgfrb", "Rgs5", "Myh11", # myCAF
    "H2-Ab1", "H2-Eb1", "Cd74", # antigen-presenting CAF
    "Acta2", "Tagln", "Myl9", # smooth muscle
    "Tubb3", "Nefl", "Map2", # neuron
    "S100b", "Sox10", # Schwann cell
    "Fabp7", "P2rx7", # Satellite glial cells,
    "Mki67", "Top2a" # cycling
)
length(marker_genes[["Stromal"]])
setdiff(marker_genes[["Stromal"]], rownames(stromal_obj))
which(duplicated(marker_genes[["Stromal"]]))

p <- FeaturePlot(stromal_obj,
    features = marker_genes[["Stromal"]], pt.size = 1.5,
    ncol = 6, cols = c("gray", "red"), raster = TRUE
)
pdf("sub3_Stromal/UMAP.MarkerGenes.pdf", 24, 18)
plot(p)
dev.off()
rm(p)

pdf("sub3_Stromal/DotPlot.MarkerGenes.pdf", 9, 7)
DotPlot(stromal_obj,
    features = marker_genes[["Stromal"]],
    group.by = "seurat_clusters"
) +
    scale_color_viridis_c() +
    theme_classic() + coord_flip()
dev.off()

marker_undefined_c24 <- FindMarkers(stromal_obj,
    ident.1 = 24,
    min.pct = 0.25, logfc.threshold = 0.25,
    max.cells.per.ident = 1000
)
marker_undefined_c24 <- marker_undefined_c24[order(marker_undefined_c24$avg_log2FC, decreasing = TRUE), ]
head(rownames(marker_undefined_c24), 100) #  Adipocyte  markers: Adipoq Pck1

stromal_obj$Cell_Type_2 <- "none"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(4, 16, 20, 23))] <- "Vascular_Endothelial"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(13, 25))] <- "Lymphatic_Endothelial"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(0:3, 5, 7, 8, 17))] <- "iCAF"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(12))] <- "myCAF"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(6, 9, 10, 11, 14, 15, 18, 19))] <- "apCAF"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(21))] <- "Smooth_Muscle"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(22))] <- "Schwann"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(24))] <- "Adipocyte"

table(stromal_obj$seurat_clusters, stromal_obj$Cell_Type_2)

stromal_obj$Cell_Type_2 <- factor(stromal_obj$Cell_Type_2, levels = c(
    "Vascular_Endothelial", "Lymphatic_Endothelial",
    "iCAF", "myCAF", "apCAF", "Smooth_Muscle", "Schwann", "Adipocyte"
))

mycolor$Cell_Type_2 <- mycolor$Cell_Type_2[1:20]
mycolor$Cell_Type_2 <- c(
    mycolor$Cell_Type_2,
    setNames(ggsci::pal_lancet()(8), levels(stromal_obj$Cell_Type_2))
)

pdf("sub3_Stromal/UMAP.Cell_Type_2.pdf", 13, 6)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE, pt.size = 1.5,
    cols = mycolor$Cell_Type_2
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2") + NoLegend() +
    DimPlot(stromal_obj,
        reduction = "umap", group.by = "Cell_Type_2", label = FALSE, raster = TRUE, pt.size = 1.5,
        cols = mycolor$Cell_Type_2
    ) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

pdf("sub3_Stromal/UMAP.Cell_Type_2.split.pdf", 9, 3)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition",
    label = FALSE, raster = TRUE, pt.size = 1.5,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()
dev.off()

p <- DimPlot(stromal_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition",
    label = FALSE, raster = TRUE, pt.size = NULL,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()

pdf("sub3_Stromal/UMAP.Cell_Type_2.split.density.pdf", 12, 4)
plot(p + stat_density_2d(aes_string(x = "umap_1", y = "umap_2", fill = "after_stat(level)"),
    linewidth = 0, geom = "density_2d_filled",
    n = 150, h = c(1.2, 1.2)
) +
    scale_fill_viridis_c(option = "D") +
    guides(fill = guide_colorbar(title = "Density")) +
    theme(legend.position = "right"))
dev.off()
rm(p)

## 5.2. Marker gene expression ----
Idents(stromal_obj) <- stromal_obj$Cell_Type_2
levels(Idents(stromal_obj))

gene_selected <- c(
    "Cdh5", "Pecam1", "Prox1", "Lyve1", "Col1a2", "Lum",
    "Pdgfra", "Rgs5", "H2-Ab1", "Cd74", "Acta2", "Myh11", "Sox10", "S100b",
    "Adipoq", "Pck1"
)

pdf("sub3_Stromal/DotPlot.Cell_Type_2.pdf", 6.5, 2.6)
DotPlot(stromal_obj, features = gene_selected) +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

rm(list = grep("marker.undefined", ls(), value = TRUE))
rm(gene_selected, p)

## 5.3. Percentage changes ----
# barplot composition
plot.data <- as.data.frame(table(stromal_obj$Cell_Type_2, stromal_obj$Condition))
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
plot.data <- as.data.frame(table(stromal_obj$Cell_Type_2, stromal_obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(stromal_obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))

pdf("sub3_Stromal/Boxplot.Fraction.pdf", 6.5, 4.8)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_6OHDA"), c("CTRL", "LC_6OHDA"))
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

saveRDS(stromal_obj, file = "sub3_Stromal/stromal_obj.rds")
rm(stromal_obj)

# save data ----
save.image("01.Cell_Type.RData")
saveRDS(mycolor, file = "mycolor.rds")
