# local
setwd("/Users/zhenyuliu/Research/Project/BLCA_Nerve/results/01.Cell_Type")
# server
setwd("/thinker/tb269store/user/liuzhenyu/Project/BLCA_Nerve/results/01.Cell_Type")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(DoubletFinder)
library(parallel)

load("01.Cell_Type.RData")

# 1. Load data & QC ----
data_dir <- "/thinker/tb269store/user/liuzhenyu/Project/BLCA_Nerve/data/CellRanger_Out"
samples <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)
mycolor <- list()

## 1.1. Load data ----
# counts
rc_list <- list()
rc_list <- lapply(samples, function(x) {
    data <- Read10X(data_dir = file.path(data_dir, x, "filtered_feature_bc_matrix"))
    colnames(data) <- paste(x, colnames(data), sep = "_")
    data
})
lapply(rc_list, dim)

common_genes <- Reduce(intersect, lapply(rc_list, rownames))
rc_list <- lapply(rc_list, function(x) {
    x[common_genes, ]
})
saveRDS(rc_list, file = "../../data/rc_list.rds")
rc_list <- readRDS(file = "../../data/rc_list.rds")

rc <- do.call(cbind, rc_list)
saveRDS(rc, file = "../../data/rc.rds")
rc <- readRDS(file = "../../data/rc.rds")
rm(rc_list)

# metadata
cells <- colnames(rc)
sample_info <- data.frame(
    row.names = cells,
    "Cell" = cells,
    "Sample" = gsub("_A.+?*$|_T.+?*$|_C.+?*$|_G.+?*$", "", cells)
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
rm(cells)

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

table(singler_anno$Cell %in% sample_info$Cell) # 376867 cells
table(sample_info$Cell %in% singler_anno$Cell)
singler_anno[1:5, ]

## 1.2. QC for cells ----
sample_info$Cell_Type_SingleR <- "none"
sample_info[singler_anno$Cell, "Cell_Type_SingleR"] <- singler_anno$raw_main_celltype
sample_info[1:5, ]
rm(singler_anno)

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

## 1.3. Filter doublets -----
obj_list <- SplitObject(blca_obj, split.by = "Sample")

doublet_res <- lapply(obj_list, function(seu) {
    seu <- seu %>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA() %>%
        RunUMAP(dims = 1:15)

    sweep_res <- paramSweep(seu, PCs = 1:15, sct = FALSE)
    sweep_stats <- summarizeSweep(sweep_res, GT = FALSE)
    bcmvn <- find.pK(sweep_stats)
    optimal_pk <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

    # Estimate doublet rate based on cell number (~0.8% per 1000 cells)
    nexp_poi <- round(0.008 * (ncol(seu) / 1000) * ncol(seu))

    # Run DoubletFinder
    seu <- doubletFinder(seu, PCs = 1:15, pN = 0.25, pK = optimal_pk, nExp = nexp_poi, reuse.pANN = NULL, sct = FALSE)

    df_cols <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)
    seu$Doublet_Classification <- seu[[df_cols]]

    data.frame(
        row.names = colnames(seu),
        Cell = colnames(seu),
        Doublet_Classification = seu$Doublet_Classification
    )
})
doublet_res <- do.call(rbind, doublet_res)
rownames(doublet_res) <- gsub("^.+?\\.", "", rownames(doublet_res))

table(rownames(sample_info_all) %in% rownames(doublet_res))
sample_info_all$Doublet_Classification <- "none"
sample_info_all[rownames(doublet_res), "Doublet_Classification"] <- doublet_res$Doublet_Classification
table(sample_info_all$Doublet_Classification)

table(colnames(blca_obj) %in% rownames(doublet_res))
blca_obj$Doublet_Classification <- "none"
blca_obj@meta.data[rownames(doublet_res), "Doublet_Classification"] <- doublet_res$Doublet_Classification
table(blca_obj$Doublet_Classification)

blca_obj <- subset(blca_obj, subset = Doublet_Classification == "Singlet")
sample_info <- blca_obj@meta.data
rm(doublet_res, obj_list)

## 1.4. Reductions & clustering ----
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

mycolor$Cluster <- colorRampPalette(ggsci::pal_npg()(10))(length(unique(blca_obj$seurat_clusters))) %>%
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
mycolor$Cell_Type_SingleR <- colorRampPalette(brewer.pal(10, "Set3"))(length(unique(blca_obj$Cell_Type_SingleR))) %>%
    setNames(levels(blca_obj$Cell_Type_SingleR))
pdf("UMAP.SingleR.pdf", 5, 4)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, label.size = 3) +
    scale_color_manual(values = mycolor$Cell_Type_SingleR) +
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
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(26))] <- "Epithelial"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(0, 1, 5, 6, 8, 9, 12, 13, 21, 25, 31))] <- "Cancer"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(30, 33))] <- "Endothelial"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(3, 10, 15, 16, 18))] <- "T"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(19))] <- "NK"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(28))] <- "B"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(2, 7, 11, 17, 20, 23, 24, 32))] <- "Macro/Mono"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(17, 27, 34))] <- "DC"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(4))] <- "Granulocyte"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(14))] <- "Fibroblast"
blca_obj$Cell_Type_1[which(blca_obj$seurat_clusters %in% c(22, 29))] <- "Doublet"
table(blca_obj$seurat_clusters, blca_obj$Cell_Type_1)

blca_obj$Cell_Type_1 <- factor(blca_obj$Cell_Type_1, levels = c(
    "Epithelial", "Cancer",
    "Endothelial", "Fibroblast",
    "T", "B", "NK",
    "Granulocyte", "Macro/Mono", "DC",
    "Doublet"
))
table(blca_obj$Cell_Type_1)
Idents(blca_obj) <- "Cell_Type_1"

mycolor$Cell_Type_1 <- ggsci::pal_npg()(10) %>%
    setNames(levels(blca_obj$Cell_Type_1)[1:10])
mycolor$Cell_Type_1["Doublet"] <- "gray"

pdf("UMAP.Cell_Type_1.pdf", 5.5, 4)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_1", label = TRUE, label.size = 3) +
    scale_color_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

pdf("UMAP.Cell_Type_1.split.Condition.pdf", 12, 8)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_1", label = TRUE, split.by = "Condition", ncol = 3) +
    scale_color_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

gene_selected <- c(
    "Epcam", "Krt8", "Krt18", # epithelial
    "Pecam1", "Cdh5", # endothelial
    "Col1a2", "Dcn", # fibroblast
    "Ptprc", # leukocyte
    "Cd3d", "Cd3g", "Cd79a", "Klrb1c", # lymphocyte
    "S100a8", "S100a9", # neutrophils
    "Cd14", "Cd68", "Cst3"# myeloid
)
pdf("DotPlot.Cell_Type_1.pdf", 7, 3.5)
DotPlot(blca_obj,
    features = gene_selected, group.by = "Cell_Type_1",
    idents = setdiff(levels(blca_obj$Cell_Type_1), "Doublet")
) +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
rm(gene_selected)

# bigger categories
blca_obj$Cell_Type_0 <- "none"
blca_obj$Cell_Type_0[which(blca_obj$Cell_Type_1 %in% c("Epithelial", "Cancer"))] <- "Epithelial"
blca_obj$Cell_Type_0[which(blca_obj$Cell_Type_1 %in% c("Fibroblast", "Endothelial"))] <- "Stromal"
blca_obj$Cell_Type_0[which(blca_obj$Cell_Type_1 %in% c("T", "B", "NK", "Granulocyte", "Macro/Mono", "DC"))] <- "Immune"
blca_obj$Cell_Type_0[which(blca_obj$Cell_Type_1 %in% c("Doublet"))] <- "Doublet"
table(blca_obj$Cell_Type_0)

blca_obj$Cell_Type_0 <- factor(blca_obj$Cell_Type_0, levels = c("Epithelial", "Stromal", "Immune", "Doublet"))
mycolor$Cell_Type_0 <- c(ggsci::pal_jco()(4)[c(1, 2, 4)]) %>%
    setNames(levels(blca_obj$Cell_Type_0)[1:3])
mycolor$Cell_Type_0["Doublet"] <- "gray"

pdf("UMAP.Cell_Type_0.pdf", 5.5, 4)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_0", label = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_0")
dev.off()

pdf("UMAP.Cell_Type_0.split.Condition.pdf", 12, 8)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_0", label = TRUE, split.by = "Condition", ncol = 3) +
    scale_color_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_0")
dev.off()

## 2.4. Cell types cross conditions ----
plot_data <- as.data.frame(table(blca_obj$Cell_Type_1, blca_obj$Condition)) %>%
    filter(Var1 != "Doublet")
colnames(plot_data) <- c("Cell_Type_1", "Condition", "Count")

pdf("Barplot.Cell_Type_1.Condition.pdf", 8, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_1, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = Cell_Type_1)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_1) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

plot_data <- as.data.frame(table(blca_obj$Cell_Type_0, blca_obj$Condition)) %>%
    filter(Var1 != "Doublet")
colnames(plot_data) <- c("Cell_Type_0", "Condition", "Count")

pdf("Barplot.Cell_Type_0.Condition.pdf", 8, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_0, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = Cell_Type_0)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_0) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot_data <- as.data.frame(table(blca_obj$Cell_Type_0, blca_obj$Sample)) %>%
    filter(Var1 != "Doublet")
colnames(plot_data) <- c("Cell_Type_0", "Sample", "Count")
plot_data$Condition <- gsub("_[1-5]+?*$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(blca_obj$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("Boxplot.Fraction.Cell_Type_0.pdf", 5, 3)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_0), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 3, ref.group = "CTRL",
    ) +
    facet_wrap(~Cell_Type_0, ncol = 3, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_0) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()

plot_data <- as.data.frame(table(blca_obj$Cell_Type_1, blca_obj$Sample)) %>%
    filter(Var1 != "Doublet")
colnames(plot_data) <- c("Cell_Type_1", "Sample", "Count")
plot_data$Condition <- gsub("_[1-5]+?*$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(blca_obj$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("Boxplot.Fraction.Cell_Type_1.pdf", 8, 5)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_1), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 3, ref.group = "CTRL",
    ) +
    facet_wrap(~Cell_Type_1, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot_data, temp)

# Roe
mycolor$Roe <- colorRampPalette(brewer.pal(9, "RdBu"))(100) %>% rev()

roe <- Startrac::calTissueDist(
    blca_obj@meta.data,
    byPatient = FALSE,
    colname.cluster = "Cell_Type_0",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[setdiff(rownames(roe), "Doublet"), ]
pdf("Heatmap.Roe.Cell_Type_0.Condition.pdf", 4, 2.5)
pheatmap::pheatmap(
    roe,
    color = mycolor$Roe,
    breaks = seq(0, 2, length.out = 101),
    scale = "none",
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
)
dev.off()

roe <- Startrac::calTissueDist(
    blca_obj@meta.data,
    byPatient = FALSE,
    colname.cluster = "Cell_Type_1",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[setdiff(rownames(roe), "Doublet"), ]
pdf("Heatmap.Roe.Cell_Type_1.Condition.pdf", 4, 5)
pheatmap::pheatmap(
    roe,
    color = mycolor$Roe,
    breaks = seq(0, 2, length.out = 101),
    scale = "none",
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
)
dev.off()
rm(roe)

dir.create("rds")
saveRDS(blca_obj, file = "rds/blca_obj.rds")
# rm(blca_obj)

# 3. Immune compartment ----
##  3.1. Detailed cell typing ----
dir.create("sub1_Immune")

immune_obj <- subset(blca_obj, subset = Cell_Type_0 == "Immune") # 165931 cells
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

# saveRDS(immune_obj, file = "rds/immune_obj.rds")
# immune_obj <- readRDS("rds/immune_obj.rds")

pdf("sub1_Immune/UMAP.Clusters.pdf", 5.5, 4)
DimPlot(immune_obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    guides(color = guide_legend(ncol = 3)) +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub1_Immune/UMAP.Clusters.split.pdf", 10, 6)
DimPlot(immune_obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition", ncol = 3,
    label = TRUE, cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

table(immune_obj$Cell_Type_SingleR)
pdf("sub1_Immune/UMAP.SingleR.pdf", 5, 4)
DimPlot(immune_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_SingleR) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

pdf("sub1_Immune/UMAP.Cell_Type_1.pdf", 5, 4)
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

immune_obj$Cell_Type_2 <- "none"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(14, 17, 33))] <- "T_CD4"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(29))] <- "T_CD4_Cycling"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(10))] <- "Treg"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(0, 1, 9, 22, 11, 16))] <- "T_CD8"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(9))] <- "T_CD8_Cycling"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(18))] <- "NK"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(26))] <- "B"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(3, 4, 6, 13, 23))] <- "Macrophage"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(5, 8, 12, 25, 31, 19))] <- "Monocyte"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(15, 27, 24))] <- "cDC"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(30))] <- "pDC"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(34))] <- "Mast"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(2, 7, 20, 28))] <- "Neutrophil"
immune_obj$Cell_Type_2[which(immune_obj$seurat_clusters %in% c(32, 21))] <- "Myeloid_Cycling"

table(immune_obj$seurat_clusters, immune_obj$Cell_Type_2)

immune_obj$Cell_Type_2 <- factor(immune_obj$Cell_Type_2, levels = c(
    "T_CD4", "T_CD4_Cycling", "Treg", "T_CD8", "T_CD8_Cycling", "NK", "B",
    "Macrophage", "Monocyte", "cDC", "pDC", "Mast", "Neutrophil", "Myeloid_Cycling"
))

mycolor$Cell_Type_2 <- colorRampPalette(ggsci::pal_npg()(10))(length(levels(immune_obj$Cell_Type_2))) %>%
    setNames(levels(immune_obj$Cell_Type_2))

pdf("sub1_Immune/UMAP.Cell_Type_2.pdf", 9, 4)
DimPlot(immune_obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE, label.size = 3,
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

pdf("sub1_Immune/UMAP.Cell_Type_2.split.pdf", 8, 5)
DimPlot(immune_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition", ncol = 3,
    label = FALSE, raster = TRUE, pt.size = 2,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()
dev.off()

p <- DimPlot(immune_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition", ncol = 3,
    label = FALSE, raster = TRUE, pt.size = NULL,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()

pdf("sub1_Immune/UMAP.Cell_Type_2.split.density.pdf", 8, 5)
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
gene_selected <- c(
    "Cd3d", "Cd4", "Foxp3", "Cd8a", "Prf1", "Klrb1c", "Cd79a",
    "Cd68", "Mrc1", "Cst3", "Siglech", "Hdc", "S100a8", "Mki67"
)

pdf("sub1_Immune/DotPlot.Cell_Type_2.pdf", 6, 3.5)
DotPlot(immune_obj,  features = gene_selected, group.by = "Cell_Type_2") +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
rm(gene_selected)

## 3.3. Percentage changes ----
# barplot composition
plot_data <- as.data.frame(table(immune_obj$Cell_Type_2, immune_obj$Condition))
colnames(plot_data) <- c("Cell_Type_2", "Condition", "Count")

pdf("sub1_Immune/Barplot.Cell_Type_2.Condition.pdf", 9, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_2, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = Cell_Type_2)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot_data <- as.data.frame(table(immune_obj$Cell_Type_2, immune_obj$Sample))
colnames(plot_data) <- c("Cell_Type_2", "Sample", "Count")
plot_data$Condition <- gsub("_[0-5].+?*$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(immune_obj$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("sub1_Immune/Boxplot.Fraction.pdf", 8, 7)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 3, ref.group = "CTRL"
    ) +
    facet_wrap(~Cell_Type_2, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction within immune compartment (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot_data, temp)

# Roe
roe <- Startrac::calTissueDist(
    immune_obj@meta.data,
    byPatient = FALSE,
    colname.cluster = "Cell_Type_2",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[grep("Cycling|Doublet", rownames(roe), invert = TRUE), ]
pdf("sub1_Immune/Heatmap.Roe.Cell_Type_2.Condition.pdf", 4, 5.5)
pheatmap::pheatmap(
    roe,
    color = mycolor$Roe,
    breaks = seq(0, 2, length.out = 101),
    scale = "none",
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
)
dev.off()

saveRDS(immune_obj, file = "rds/immune_obj.rds")
rm(immune_obj, roe)

# 4. Epithelial compartment ----
##  4.1. Clustering of cell states ----
dir.create("sub2_Cancer")

cancer_obj <- subset(blca_obj, subset = Cell_Type_1 %in% c("Cancer")) # 161321 cells

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

# saveRDS(cancer_obj, file = "rds/cancer_obj.rds")
# cancer_obj <- readRDS("rds/cancer_obj.rds")

library(ggraph)
pdf("sub2_Cancer/ClusteringTree.pdf", 8, 6)
clustree::clustree(cancer_obj, prefix = "RNA_snn_res.")
dev.off()

cancer_obj$seurat_clusters <- cancer_obj$RNA_snn_res.0.3
pdf("sub2_Cancer/UMAP.Clusters.pdf", 5, 4)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,
    cols = mycolor$Cluster, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Cancer/UMAP.Clusters.split.pdf", 8, 6)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition",
    label = TRUE, cols = mycolor$Cluster, raster = TRUE, ncol = 3
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub2_Cancer/UMAP.Cell_Type_1.pdf", 5, 4)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "Cell_Type_1", label = TRUE,
    cols = mycolor$Cell_Type_1, raster = TRUE
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_1")
dev.off()

table(cancer_obj$Cell_Type_SingleR)
pdf("sub2_Cancer/UMAP.SingleR.pdf", 5, 4)
DimPlot(cancer_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_SingleR) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

## 4.2. Cell states annotation ----
table(cancer_obj$seurat_clusters)
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

pdf("sub2_Cancer/UMAP.Cluster.Markers.pdf", 8, 6)
FeaturePlot(cancer_obj,
    features = c("Ptprc", "Krt18", "Krt8"),
    ncol = 2, cols = c("gray", "red"), raster = TRUE
)
dev.off()

cancer_obj$Cell_Type_2 <- "none"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(4, 7))] <- "C1: Proliferative"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(0))] <- "C2: Hypoxic"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(1))] <- "C3: EMT"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(5))] <- "C4: Acute_Stress"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(2, 6))] <- "C5: Ribosomal_Adhesion"
cancer_obj$Cell_Type_2[which(cancer_obj$seurat_clusters %in% c(3, 8, 9))] <- "Doublet"
table(cancer_obj$seurat_clusters, cancer_obj$Cell_Type_2)

cancer_obj$Cell_Type_2 <- factor(cancer_obj$Cell_Type_2, levels = c(
    "C1: Proliferative", "C2: Hypoxic", "C3: EMT", "C4: Acute_Stress", "C5: Ribosomal_Adhesion", "Doublet"
))

mycolor$Cell_Type_2 <- c(
    mycolor$Cell_Type_2,
    setNames(ggsci::pal_jco()(6)[c(1:2, 4:6)], levels(cancer_obj$Cell_Type_2)[1:5])
)
mycolor$Cell_Type_2["Doublet"] <- "gray"

pdf("sub2_Cancer/UMAP.Cell_Type_2.pdf", 8, 6)
DimPlot(cancer_obj,
    reduction = "umap", group.by = "Cell_Type_2", label = TRUE, raster = TRUE,
    cols = mycolor$Cell_Type_2
) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

## 4.3. Percentage changes ----
# barplot composition
plot_data <- as.data.frame(table(cancer_obj$Cell_Type_2, cancer_obj$Condition)) %>%
    filter(Var1 != "Doublet")
colnames(plot_data) <- c("Cell_Type_2", "Condition", "Count")

pdf("sub2_Cancer/Barplot.Cell_Type_2.Condition.pdf", 9, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_2, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = Cell_Type_2)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot_data <- as.data.frame(table(cancer_obj$Cell_Type_2, cancer_obj$Sample)) %>%
    filter(Var1 != "Doublet")
colnames(plot_data) <- c("Cell_Type_2", "Sample", "Count")
plot_data$Condition <- gsub("_[1-5]+?*$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(cancer_obj$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("sub2_Cancer/Boxplot.Fraction.pdf", 8, 2.8)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 2.5, ref.group = "CTRL"
    ) +
    facet_wrap(~Cell_Type_2, ncol = 5, scales = "free") +
    labs(title = "Cell State composition across Samples") +
    ylab("Fraction of cancer cells (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot_data, temp)

# Roe
roe <- Startrac::calTissueDist(
    cancer_obj@meta.data,
    byPatient = FALSE,
    colname.cluster = "Cell_Type_2",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[grep("Cycling|Doublet", rownames(roe), invert = TRUE), ]
pdf("sub2_Cancer/Heatmap.Roe.Cell_Type_2.Condition.pdf", 5, 3)
pheatmap::pheatmap(
    roe,
    color = mycolor$Roe,
    breaks = seq(0, 2, length.out = 101),
    scale = "none",
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
)
dev.off()
rm(roe)

## 4.5. DEG cross conditions ----
table(cancer_obj$Condition)

deg_cancer_conditions <- lapply(names(mycolor$Condition)[-1], function(x) {
    temp <- FindMarkers(cancer_obj,
        ident.1 = x, ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0.25, logfc.threshold = 0.25,
        max.cells.per.ident = 10000
    )

    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Condition <- x
    temp
})
names(deg_cancer_conditions) <- names(mycolor$Condition)[-1]

lapply(names(deg_cancer_conditions), function(x) {
    write.csv(deg_cancer_conditions[[x]],
        file = paste0("sub2_Cancer/deg_", x, "_vs_CTRL_in_.csv"),
        row.names = FALSE, quote = FALSE
    )
})
rm(temp)

# volcano plots
plot_data <- do.call(rbind, deg_cancer_conditions)

plot_data$DEG <- "Stable"
plot_data$DEG[which(plot_data$p_val_adj < 0.05 & plot_data$avg_log2FC > 0.25)] <- "Up"
plot_data$DEG[which(plot_data$p_val_adj < 0.05 & plot_data$avg_log2FC < -0.25)] <- "Down"
table(plot_data$Condition, plot_data$DEG)

plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition)[-1])

label_data <- table(plot_data$Condition, plot_data$DEG) %>%
    as.data.frame()
colnames(label_data) <- c("Condition", "DEG", "Count")
label_data <- label_data %>%
    filter(DEG != "Stable") %>%
    mutate("x" = ifelse(DEG == "Up", 2, -2), label = paste0(DEG, ": ", Count))

pdf("sub2_Cancer/Volcano.Condition_vs_CTRL.pdf", 7, 3)
ggplot(plot_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = Condition), alpha = 0.5, size = 0.5) +
    geom_text(
        data = label_data,
        aes(x = x, y = 350, label = label), inherit.aes = FALSE, size = 3
    ) +
    scale_color_manual(values = mycolor$Condition) +
    facet_wrap(~Condition, ncol = 5) +
    theme_classic() +
    labs(title = "Volcano plot of DEGs in each Condition") +
    xlab("Average log2 fold change") +
    ylab("-log10 adjusted p value") +
    xlim(c(-3, 3)) + ylim(c(0, 380)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    theme(legend.position = "none")
dev.off()
rm(plot_data, label_data)

# GO enrichment
library(clusterProfiler)
library(org.Mm.eg.db)

go_deg_cancer_conditions_up <- lapply(deg_cancer_conditions, function(x) {
    temp <- x %>%
        filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
        pull(gene)
    go_res <- enrichGO(
        gene = temp, OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05
    )
    if (is.null(go_res)) {
        return(NULL)
    }
    simplify(go_res)
})

go_deg_cancer_conditions_down <- lapply(deg_cancer_conditions, function(x) {
    temp <- x %>%
        filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
        pull(gene)
    go_res <- enrichGO(
        gene = temp, OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05
    )
    if (is.null(go_res)) {
        return(NULL)
    }
    simplify(go_res)
})

lapply(names(go_deg_cancer_conditions_up), function(x) {
    if (is.null(go_deg_cancer_conditions_up[[x]])) {
        return(NULL)
    }
    write.csv(go_deg_cancer_conditions_up[[x]]@result,
        file = paste0("sub2_Cancer/GO_", x, "_vs_CTRL_up.csv"),
        row.names = TRUE
    )
    pdf(paste0("sub2_Cancer/GO_", x, "_vs_CTRL_up.pdf"), 7, 5)
    p <- dotplot(go_deg_cancer_conditions_up[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("Cancer - Upregulated in ", x, sep = ""))
    plot(p)
    dev.off()
})

lapply(names(go_deg_cancer_conditions_down), function(x) {
    if (is.null(go_deg_cancer_conditions_down[[x]])) {
        return(NULL)
    }
    write.csv(go_deg_cancer_conditions_down[[x]]@result,
        file = paste0("sub2_Cancer/GO_", x, "_vs_CTRL_down.csv"),
        row.names = TRUE
    )
    pdf(paste0("sub2_Cancer/GO_", x, "_vs_CTRL_down.pdf"), 7, 5)
    p <- dotplot(go_deg_cancer_conditions_down[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("Cancer - Downregulated in ", x, sep = ""))
    plot(p)
    dev.off()
})

saveRDS(cancer_obj, file = "rds/cancer_obj.rds")
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
stromal_obj <- readRDS("rds/stromal_obj.rds")

pdf("sub3_Stromal/ElbowPlot.PCA.pdf", 6, 4)
ElbowPlot(stromal_obj, ndims = 50)
dev.off()

stromal_obj <- FindNeighbors(stromal_obj, dims = 1:30)
stromal_obj <- FindClusters(stromal_obj, resolution = 1.0)

stromal_obj <- RunUMAP(stromal_obj, dims = 1:30)

pdf("sub3_Stromal/UMAP.Clusters.pdf", 5.5, 4)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,
    cols = mycolor$Cluster, raster = TRUE, pt.size = 1.5
) +
    theme_classic() +
    guides(color = guide_legend(ncol = 3)) +
    labs(title = "UMAP plot colored by clusters")
dev.off()

pdf("sub3_Stromal/UMAP.Clusters.split.pdf", 8, 6)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "seurat_clusters", split.by = "Condition", ncol = 3,
    label = TRUE, cols = mycolor$Cluster, raster = TRUE, pt.size = 1.5
) +
    theme_classic() +
    labs(title = "UMAP plot colored by clusters")
dev.off()

table(stromal_obj$Cell_Type_SingleR)
pdf("sub3_Stromal/UMAP.SingleR.pdf", 5, 4)
DimPlot(stromal_obj, reduction = "umap", group.by = "Cell_Type_SingleR", label = TRUE, raster = TRUE, pt.size = 1.5) +
    scale_color_manual(values = brewer.pal(9, "Set3")) +
    theme_classic() +
    labs(title = "UMAP plot colored by SingleR annotation")
dev.off()

pdf("sub3_Stromal/UMAP.Cell_Type_1.pdf", 5, 4)
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
    "Adipoq", "Pck1", # adipocyte
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
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(8, 23))] <- "Vascular_Endothelial"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(2, 11, 17, 19, 22))] <- "Lymphatic_Endothelial"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(0, 1, 3:7, 9, 12:16, 18))] <- "iCAF"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(10))] <- "myCAF"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(21))] <- "Smooth_Muscle"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(20))] <- "Schwann"
stromal_obj$Cell_Type_2[which(stromal_obj$seurat_clusters %in% c(24))] <- "Adipocyte"

table(stromal_obj$seurat_clusters, stromal_obj$Cell_Type_2)

stromal_obj$Cell_Type_2 <- factor(stromal_obj$Cell_Type_2, levels = c(
    "Vascular_Endothelial", "Lymphatic_Endothelial",
    "iCAF", "myCAF", "Smooth_Muscle", "Schwann", "Adipocyte"
))

mycolor$Cell_Type_2 <- c(
    mycolor$Cell_Type_2,
    setNames(ggsci::pal_lancet()(8), levels(stromal_obj$Cell_Type_2))
)

pdf("sub3_Stromal/UMAP.Cell_Type_2.pdf", 9, 4)
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

pdf("sub3_Stromal/UMAP.Cell_Type_2.split.pdf", 8, 5)
DimPlot(stromal_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition",
    label = FALSE, raster = TRUE, pt.size = 3, ncol = 3,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()
dev.off()

p <- DimPlot(stromal_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition",
    label = FALSE, raster = TRUE, pt.size = NULL, ncol = 3,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()

pdf("sub3_Stromal/UMAP.Cell_Type_2.split.density.pdf", 8, 5)
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
    "Cdh5", "Pecam1", "Prox1", "Lyve1", "Cd34", "Col1a2", "Lum",
    "Pdgfra", "Rgs5", "Acta2", "Myh11", "Sox10", "S100b",
    "Adipoq", "Pck1"
)

pdf("sub3_Stromal/DotPlot.Cell_Type_2.pdf", 6.5, 2.6)
DotPlot(stromal_obj, features = gene_selected) +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

rm(gene_selected, p)

## 5.3. Percentage changes ----
# barplot composition
plot_data <- as.data.frame(table(stromal_obj$Cell_Type_2, stromal_obj$Condition))
colnames(plot_data) <- c("Cell_Type_2", "Condition", "Count")

pdf("sub3_Stromal/Barplot.Cell_Type_2.Condition.pdf", 9, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_2, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = Cell_Type_2)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot_data <- as.data.frame(table(stromal_obj$Cell_Type_2, stromal_obj$Sample))
colnames(plot_data) <- c("Cell_Type_2", "Sample", "Count")
plot_data$Condition <- gsub("_[0-5].+?*$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(stromal_obj$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("sub3_Stromal/Boxplot.Fraction.pdf", 6.5, 4.8)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 2.5, ref.group = "CTRL"
    ) +
    facet_wrap(~Cell_Type_2, ncol = 4, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction within immune compartment (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot_data, temp)

saveRDS(stromal_obj, file = "rds/stromal_obj.rds")
rm(stromal_obj)

# save data ----
save.image("01.Cell_Type.RData")
saveRDS(mycolor, file = "rds/mycolor.rds")
