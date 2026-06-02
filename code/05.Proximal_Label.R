# local
setwd("/Users/zhenyuliu/Research/Project/BLCA_Nerve/results/05.Proximal_Label")
# server
setwd("/thinker/tb269store/user/liuzhenyu/Project/BLCA_Nerve/results/05.Proximal_Label")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(DoubletFinder)
library(parallel)
library(cli)

mycolor <- readRDS("../rds/mycolor.rds")

load("05.Proximal_Label.RData")

# 1. load data ----
data_dir <- "/Users/zhenyuliu/Research/Project/BLCA_Nerve/data/CellRanger_Proximal"
samples <- list.dirs(data_dir, recursive = FALSE, full.names = FALSE)

rc_list <- list()
rc_list <- lapply(samples, function(x) {
    data <- Read10X(data.dir = file.path(data_dir, x, "filtered_feature_bc_matrix"))
    colnames(data) <- paste(x, colnames(data), sep = "_")
    data
})
lapply(rc_list, dim)

common_genes <- Reduce(intersect, lapply(rc_list, rownames))
rc_list <- lapply(rc_list, function(x) {
    x[common_genes, ]
})
saveRDS(rc_list, file = "../../data/rc_list_Proximal.rds")
rc_list <- readRDS(file = "../../data/rc_list_Proximal.rds")

rc <- do.call(cbind, rc_list)
saveRDS(rc, file = "../../data/rc_Proximal.rds")
rc <- readRDS(file = "../../data/rc_Proximal.rds")
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
table(sample_info$Condition)

sample_info$Sample <- gsub("CAP_", "Capsaicin_", sample_info$Sample)
sample_info$Sample <- gsub("LC_6OHDA", "6OHDA", sample_info$Sample)
sample_info$Sample <- gsub("PBS_", "CTRL_", sample_info$Sample)

sample_info$Condition <- gsub("CAP", "Capsaicin", sample_info$Condition)
sample_info$Condition <- gsub("LC_6OHDA", "6OHDA", sample_info$Condition)
sample_info$Condition <- gsub("PBS", "CTRL", sample_info$Condition)

sample_info$Condition <- factor(sample_info$Condition,
    levels = c("CTRL", "Capsaicin", "6OHDA")
)

sample_info <- sample_info[order(sample_info$Condition), ]
sample_info$Sample <- factor(sample_info$Sample, levels = unique(sample_info$Sample))

table(sample_info$Condition)
table(sample_info$Sample)
rm(cells, samples, common_genes)

## 1.2. QC for cells ----
blca_proximal <- CreateSeuratObject(
    counts = rc,
    project = "BLCA_Neuro",
    meta.data = sample_info
)
blca_proximal$Percent_Mito <- PercentageFeatureSet(blca_proximal, pattern = "^mt-")

cutoffs <- c(
    "nFeature_RNA_low" = 1000,
    "nFeature_RNA_high" = 9000,
    "nCount_RNA" = 1000,
    "Percent_Mito" = 5
)

blca_proximal <- subset(blca_proximal, subset = nFeature_RNA > cutoffs["nFeature_RNA_low"] &
    nFeature_RNA < cutoffs["nFeature_RNA_high"] &
    nCount_RNA > cutoffs["nCount_RNA"] &
    Percent_Mito < cutoffs["Percent_Mito"])
sample_info <- blca_proximal@meta.data

## 1.3. Filter doublets -----
obj_list <- SplitObject(blca_proximal, split.by = "Sample")

doublet_res <- lapply(obj_list, function(seu) {
    message(
        col_yellow(format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        " Processing sample: ",
        col_green(unique(seu$Sample))
    )
    seu <- seu %>%
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA()

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

    gc()
    data.frame(
        row.names = colnames(seu),
        Cell = colnames(seu),
        Doublet_Classification = seu$Doublet_Classification
    )
})

doublet_res <- do.call(rbind, doublet_res)
rownames(doublet_res) <- gsub("^.+?\\.", "", rownames(doublet_res))

table(rownames(sample_info) %in% rownames(doublet_res))
sample_info$Doublet_Classification <- "none"
sample_info[rownames(doublet_res), "Doublet_Classification"] <- doublet_res$Doublet_Classification
table(sample_info$Doublet_Classification)

table(colnames(blca_proximal) %in% rownames(doublet_res))
blca_proximal$Doublet_Classification <- "none"
blca_proximal@meta.data[rownames(doublet_res), "Doublet_Classification"] <- doublet_res$Doublet_Classification
table(blca_proximal$Doublet_Classification)

blca_proximal <- subset(blca_proximal, subset = Doublet_Classification == "Singlet")
sample_info <- blca_proximal@meta.data
rm(doublet_res, obj_list)

# 2. Transfer label from previous rds ----
## 2.1. Run Seurat label transfer ----
blca_proximal <- NormalizeData(blca_proximal)
blca_obj <- readRDS("../rds/blca_obj.rds")

blca_anchors <- FindTransferAnchors(
    reference = blca_obj, query = blca_proximal, dims = 1:30,
    reference.reduction = "pca"
)

predictions <- TransferData(anchorset = blca_anchors, refdata = blca_obj$Cell_Type_2, dims = 1:30)
blca_proximal$Cell_Type_2 <- predictions$predicted.id
table(blca_proximal$Cell_Type_2)

predictions <- TransferData(anchorset = blca_anchors, refdata = blca_obj$Cell_Type_1, dims = 1:30)
blca_proximal$Cell_Type_1 <- predictions$predicted.id
table(blca_proximal$Cell_Type_1)

saveRDS(blca_anchors, file = "./blca_anchors.rds")
rm(predictions, blca_anchors)

## 2.2. Check label transfer results ----
blca_proximal <- readRDS("../rds/blca_proximal.rds")
gene_selected <- c(
    "Epcam", "Krt8", "Krt18", # epithelial
    "Pecam1", "Cdh5", # endothelial
    "Col1a2", "Dcn", # fibroblast
    "Ptprc", # leukocyte
    "Cd3d", "Cd3g", "Cd79a", "Klrb1c", # lymphocyte
    "S100a8", "S100a9", # neutrophils
    "Cd14", "Cd68", "Cst3" # myeloid
)
blca_proximal$Cell_Type_1 <- factor(blca_proximal$Cell_Type_1,
    levels = names(mycolor$Cell_Type_1)
)
Idents(blca_proximal) <- blca_proximal$Cell_Type_1
pdf("DotPlot.Cell_Type_1.pdf", 7, 3.5)
DotPlot(blca_proximal,
    features = gene_selected, group.by = "Cell_Type_1",
    idents = setdiff(levels(blca_proximal$Cell_Type_1), "Others")
) +
    scale_color_viridis_c() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
rm(gene_selected)

# 3. Fraction of cells ----
## 3.1. Labeled groups ----
plot_data <- as.data.frame(table(blca_proximal$Cell_Type_1, blca_proximal$Condition)) %>%
    filter(Var1 != "Others")
colnames(plot_data) <- c("Cell_Type_1", "Condition", "Count")

pdf("Barplot.Cell_Type_1.Condition.pdf", 8, 4)
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

plot_data <- as.data.frame(table(blca_proximal$Cell_Type_2, blca_proximal$Condition)) %>%
    filter(Var1 != "Others")
colnames(plot_data) <- c("Cell_Type_2", "Condition", "Count")

pdf("Barplot.Cell_Type_2.Condition.pdf", 10, 6)
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

plot_data <- as.data.frame(table(blca_proximal$Cell_Type_1, blca_proximal$Sample)) %>%
    filter(Var1 != "Others")
colnames(plot_data) <- c("Cell_Type_1", "Sample", "Count")

pdf("Barplot.Cell_Type_1.Sample.pdf", 8, 4)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_1, y = Count, fill = Sample)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_1 %>% unname()) +
        theme_classic() +
        labs(title = "Cell Type composition across Samples") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Sample, y = Count, fill = Cell_Type_1)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_1) +
        theme_classic() +
        labs(title = "Cell Type composition across Samples") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

plot_data <- as.data.frame(table(blca_proximal$Cell_Type_2, blca_proximal$Sample)) %>%
    filter(Var1 != "Others")
colnames(plot_data) <- c("Cell_Type_2", "Sample", "Count")

pdf("Barplot.Cell_Type_2.Sample.pdf", 10, 6)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_2, y = Count, fill = Sample)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2 %>% unname()) +
        theme_classic() +
        labs(title = "Cell Type composition across Samples") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Sample, y = Count, fill = Cell_Type_2)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Cell_Type_2) +
        theme_classic() +
        labs(title = "Cell Type composition across Samples") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

rm(plot_data)

# fraction compare
plot_data <- as.data.frame(table(blca_proximal$Cell_Type_1, blca_proximal$Sample)) %>%
    filter(Var1 != "Others")
colnames(plot_data) <- c("Cell_Type_1", "Sample", "Count")
plot_data$Condition <- gsub("_[1-5]+?*$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(blca_proximal$Sample)
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

pdf("Boxplot.Fraction.Cell_Type_1.Merge.pdf", 6, 3)
ggplot(
    plot_data %>% filter(Cell_Type_1 != "Cancer"),
    aes(x = Cell_Type_1, y = Fraction, fill = Condition)
) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    scale_fill_manual(values = mycolor$Condition) +
    ggpubr::stat_compare_means(
        aes(group = Condition),
        label = "p.signif",
        label.y.npc = 0.98,
        size = 3
    ) +
    theme_classic() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
rm(plot_data, temp)

# Roe
roe <- Startrac::calTissueDist(
    blca_proximal@meta.data,
    byPatient = FALSE,
    colname.cluster = "Cell_Type_1",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[setdiff(rownames(roe), c("Others")), ]
pdf("Heatmap.Roe.Cell_Type_1.Condition.pdf", 3, 4)
pheatmap::pheatmap(
    roe,
    color = mycolor$Roe,
    # breaks = seq(0, 2, length.out = 101),
    scale = "none",
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE
)
dev.off()
rm(roe)

## 3.2. compare label & unlabel ----
common.cols <- intersect(colnames(blca_proximal@meta.data), colnames(blca_obj@meta.data))
sample.info.merge <- rbind(
    blca_proximal@meta.data[, common.cols] %>% mutate(Dataset = "Proximal"),
    blca_obj@meta.data[, common.cols] %>% mutate(Dataset = "Unlabeled")
)
sample.info.merge$Dataset <- factor(sample.info.merge$Dataset,
    levels = c("Unlabeled", "Proximal")
)
rm(common.cols)

plot_data <- as.data.frame(table(
    sample.info.merge$Cell_Type_1,
    sample.info.merge$Condition, sample.info.merge$Dataset
)) %>%
    filter(Var1 != "Others")
colnames(plot_data) <- c("Cell_Type_1", "Condition", "Dataset", "Count")

pdf("CMP.Barplot.Cell_Type_1.Condition.pdf", 10, 4)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = Cell_Type_1, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        facet_wrap(~Dataset, ncol = 2) +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = Cell_Type_1)) +
        geom_bar(stat = "identity", position = "fill") +
        facet_wrap(~Dataset, ncol = 2) +
        scale_fill_manual(values = mycolor$Cell_Type_1) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
sample.info.merge$Sample <- paste(sample.info.merge$Dataset, sample.info.merge$Sample, sep = "_")

plot_data <- as.data.frame(table(
    sample.info.merge$Cell_Type_1,
    sample.info.merge$Sample, sample.info.merge$Dataset
)) %>%
    filter(Var1 != "Others")

colnames(plot_data) <- c("Cell_Type_1", "Sample", "Dataset", "Count")
plot_data$Condition <- gsub("_[1-5]+?*$|Unlabeled_|Proximal_", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(sample.info.merge$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("CMP.Boxplot.Fraction.Cell_Type_1.pdf", 10, 5)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Dataset), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 3, ref.group = "CTRL",
    ) +
    facet_wrap(~Cell_Type_1, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    # scale_fill_manual(values = mycolor$Cell_Type_1) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf("CMP.Boxplot.Fraction.Cell_Type_1.Merge.pdf", 6, 6)
ggplot(
    plot_data %>% filter(Cell_Type_1 != "Cancer"),
    aes(x = Cell_Type_1, y = Fraction, fill = Dataset)
) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    facet_wrap(~Condition, ncol = 1) +
    # scale_fill_manual(values = mycolor$Condition) +
    ggpubr::stat_compare_means(
        aes(group = Dataset),
        label = "p.signif",
        label.y.npc = 0.98,
        size = 3
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
rm(plot_data, temp)

# Roe
sample.info.merge$Dataset_Condition <- paste(sample.info.merge$Dataset, sample.info.merge$Condition, sep = "_")
roe <- Startrac::calTissueDist(
    sample.info.merge,
    byPatient = FALSE,
    colname.cluster = "Cell_Type_1",
    colname.patient = "Sample",
    colname.tissue = "Dataset_Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[setdiff(rownames(roe), c("Others")), ]
pdf("CMP.Heatmap.Roe.Cell_Type_1.Condition.pdf", 4.5, 5)
pheatmap::pheatmap(
    roe[, 6:1],
    color = mycolor$Roe,
    # breaks = seq(0, 2, length.out = 101),
    scale = "none",
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE
)
dev.off()
rm(roe, temp)

saveRDS(blca_proximal, file = "../rds/blca_proximal.rds")
rm(blca_proximal, blca_obj)

save.image("05.Proximal_Label.RData")
