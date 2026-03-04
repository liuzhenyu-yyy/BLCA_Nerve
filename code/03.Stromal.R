setwd("/Users/zhenyuliu/Research/Project/BLCA_Nerve/results/03.Stromal")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(parallel)

stromal_obj <- readRDS("../01.Cell_Type/rds/stromal_obj.rds")
mycolor <- readRDS("../01.Cell_Type/rds/mycolor.rds")

# load("03.Stromal.RData")

# 1. Differential analysis of each cell type ----
## 1.1. Identify DEGs ----
table(stromal_obj$Cell_Type_2, stromal_obj$Condition)
DEG_list_conditions <- mclapply(levels(stromal_obj$Condition)[-1], function(x) {
    DEG_list <- list()

    for (cell_type in levels(stromal_obj$Cell_Type_2)[1:4]) {
        message(Sys.time(), ": Differential analysis: ", x, " vs Control, ", cell_type)

        stromal_sub <- subset(stromal_obj, subset = Cell_Type_2 == cell_type)

        temp <- FindMarkers(stromal_sub,
            ident.1 = x, ident.2 = "CTRL",
            group.by = "Condition",
            min.pct = 0.25, logfc.threshold = 0.25
        )
        temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
        temp$gene <- rownames(temp)
        temp$Cell_Type_2 <- cell_type
        temp$Condition <- paste0(x, "_vs_CTRL")

        temp$DEG <- "Stable"
        temp$DEG[which(temp$p_val_adj < 0.05 & temp$avg_log2FC > 0.25)] <- "Up"
        temp$DEG[which(temp$p_val_adj < 0.05 & temp$avg_log2FC < -0.25)] <- "Down"
        DEG_list[[cell_type]] <- temp
    }

    DEG_list
}, mc.cores = 4)
names(DEG_list_conditions) <- levels(stromal_obj$Condition)[-1]
lapply(DEG_list_conditions, names)

lapply(names(DEG_list_conditions), function(x) {
    dir.create(paste0("DEG_", x))
    lapply(names(DEG_list_conditions[[x]]), function(cell_type) {
        write.csv(DEG_list_conditions[[x]][[cell_type]],
            file = paste0("DEG_", x, "/DEG_", x, "_in_", cell_type, ".csv"),
            row.names = FALSE
        )
    })
})
rm(temp, stromal_sub, cell_type, x)

saveRDS(DEG_list_conditions, file = "DEG_list_conditions.rds")
DEG_list_conditions <- readRDS("DEG_list_conditions.rds")

# volcano plot
lapply(names(DEG_list_conditions), function(x) {
    plot_data <- do.call(rbind, DEG_list_conditions[[x]]) %>% mutate(Comparison = paste0(x, "_vs_CTRL"))
    plot_data$Cell_Type_2 <- factor(plot_data$Cell_Type_2, levels = names(mycolor$Cell_Type_2)[21:24])

    x_max <- max(abs(plot_data$avg_log2FC), na.rm = TRUE)
    y_max <- max(-log10(plot_data$p_val_adj), na.rm = TRUE)

    label_data <- table(plot_data$Cell_Type_2, plot_data$Comparison, plot_data$DEG) %>%
        as.data.frame()
    colnames(label_data) <- c("Cell_Type_2", "Comparison", "DEG", "Count")
    label_data <- label_data %>%
        filter(DEG != "Stable") %>%
        mutate(
            "x" = ifelse(DEG == "Up", x_max * 0.4, -x_max * 0.4), "y" = y_max * 0.9,
            "label" = paste0(DEG, ": ", Count)
        )

    p <- ggplot(plot_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = Cell_Type_2), alpha = 0.5, size = 0.5) +
        geom_text(
            data = label_data,
            aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3
        ) +
        scale_color_manual(values = mycolor$Cell_Type_2) +
        facet_wrap(~Cell_Type_2, ncol = 5) +
        theme_classic() +
        labs(title = paste0("Volcano plot of DEGs in each cell type, ", x, "vs CTRL")) +
        xlab("Average log2 fold change") +
        ylab("-log10 adjusted p value") +
        xlim(-x_max * 0.7, x_max * 0.7) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
        theme(legend.position = "none")
    pdf(paste0("Volcano.", x, "_vs_CTRL.pdf"), 8, 3.5)
    plot(p)
    dev.off()
})

## 1.2. GO enrichment analysis ----
DEG_list_genes <- lapply(DEG_list_conditions, function(x) {
    lapply(x, function(y) {
        list(
            up = y %>%
                filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
                pull(gene),
            down = y %>%
                filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
                pull(gene)
        )
    })
})

GO_list <- mclapply(names(DEG_list_genes), function(x) {
    lapply(names(DEG_list_genes[[x]]), function(cell_type) {
        message(Sys.time(), ": GO enrichment: ", x, " vs Control, ", cell_type)

        gene_up <- DEG_list_genes[[x]][[cell_type]]$up
        gene_down <- DEG_list_genes[[x]][[cell_type]]$down

        GO_up <- enrichGO(
            gene = gene_up, OrgDb = org.Mm.eg.db,
            keyType = "SYMBOL", ont = "BP",
            pAdjustMethod = "BH", qvalueCutoff = 0.05
        )
        GO_down <- enrichGO(
            gene = gene_down, OrgDb = org.Mm.eg.db,
            keyType = "SYMBOL", ont = "BP",
            pAdjustMethod = "BH", qvalueCutoff = 0.05
        )

        list(
            up = if (is.null(GO_up)) NULL else simplify(GO_up),
            down = if (is.null(GO_down)) NULL else simplify(GO_down)
        )
    }) %>% setNames(names(DEG_list_genes[[x]]))
}, mc.cores = 4)
names(GO_list) <- names(DEG_list_genes)

# output GO results
lapply(names(GO_list), function(x) {
    dir.create(paste0("GO_", x))
    lapply(names(DEG_list_conditions[[x]]), function(cell_type) {
        if (!is.null(GO_list[[x]][[cell_type]]$up)) {
            write.csv(GO_list[[x]][[cell_type]]$up@result,
                file = paste0("GO_", x, "/GO_", x, "_vs_CTRL_", cell_type, "_UP.csv"),
                row.names = FALSE
            )
        }
        if (!is.null(GO_list[[x]][[cell_type]]$down)) {
            write.csv(GO_list[[x]][[cell_type]]$down@result,
                file = paste0("GO_", x, "/GO_", x, "_vs_CTRL_", cell_type, "_DOWN.csv"),
                row.names = FALSE
            )
        }
    })
})
lapply(names(GO_list), function(x) {
    lapply(names(DEG_list_conditions[[x]]), function(cell_type) {
        if (!is.null(GO_list[[x]][[cell_type]]$up)) {
            p <- dotplot(GO_list[[x]][[cell_type]]$up, showCategory = 20, label_format = 60) +
                theme_classic() +
                labs(title = paste(x, " vs CTRL in ", cell_type, " (UP)", sep = ""))
            pdf(paste0("GO_", x, "/Dot.GO_", x, "_vs_CTRL_", cell_type, "_UP.pdf"), 7, 5)
            plot(p)
            dev.off()
        }
        if (!is.null(GO_list[[x]][[cell_type]]$down)) {
            p <- dotplot(GO_list[[x]][[cell_type]]$down, showCategory = 20, label_format = 60) +
                theme_classic() +
                labs(title = paste(x, " vs CTRL in ", cell_type, " (DOWN)", sep = ""))
            pdf(paste0("GO_", x, "/Dot.GO_", x, "_vs_CTRL_", cell_type, "_DOWN.pdf"), 7, 5)
            plot(p)
            dev.off()
        }
    })
})

## 1.3. expression of CGRP receptor ----
pdf("UMAP_Ramp1.split.pdf", 16, 3.5)
FeaturePlot(
    stromal_obj,
    features = c("Ramp1"), cols = c("lightgrey", "red"),
    split.by = "Condition", order = TRUE, raster = TRUE, pt.size = 3
)
dev.off()

pdf("Violin_Ramp1.CAF.pdf", 8, 3.5)
VlnPlot(
    stromal_obj,
    features = c("Ramp1"),
    split.by = "Condition", cols = mycolor$Condition
)
dev.off()

stromal_obj$Cell_Type_2_Condition <- paste0(stromal_obj$Cell_Type_2, "_", stromal_obj$Condition)
stromal_obj$Cell_Type_2_Condition <- factor(
    stromal_obj$Cell_Type_2_Condition,
    levels = c(
        paste0(levels(stromal_obj$Cell_Type_2), "_", "CTRL"),
        paste0(levels(stromal_obj$Cell_Type_2), "_", "Capsaicin"),
        paste0(levels(stromal_obj$Cell_Type_2), "_", "LC_6OHDA"),
        paste0(levels(stromal_obj$Cell_Type_2), "_", "anti_PD1"),
        paste0(levels(stromal_obj$Cell_Type_2), "_", "Rime")
    )
)

pdf("DotPlot_Ramp1.pdf", 7, 7)
DotPlot(
    stromal_obj,
    features = c("Ramp1"),
    group.by = "Cell_Type_2_Condition"
) + scale_color_viridis_c()
dev.off()

# 2. iCAF functions ----
## 2.1. identify iCAF subsets ----
dir.create("iCAF")
iCAF_obj <- subset(stromal_obj, subset = Cell_Type_2 == "iCAF") # 7254 cells
rm(stromal_obj)

iCAF_obj <- iCAF_obj %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(dims = 1:20) %>%
    RunUMAP(dims = 1:20)

iCAF_obj <- FindClusters(iCAF_obj, resolution = 0.3)

pdf("iCAF/UMAP_iCAF.Subclusters.pdf", 5, 4)
DimPlot(iCAF_obj, raster = TRUE, pt.size = 2, group.by = "seurat_clusters",
 label = TRUE, cols = mycolor$Cluster) +
    labs(title = "iCAF Subclusters")
dev.off()

pdf("iCAF/UMAP_iCAF.Subclusters.split.pdf", 9, 6)
DimPlot(iCAF_obj,
    split.by = "Condition", raster = TRUE, ncol = 3,
    cols = mycolor$Cluster, pt.size = 3, label = FALSE
) +
    labs(title = "iCAF Subclusters")
dev.off()

pdf("iCAF/UMAP_iCAF.Ramp1.split.pdf", 16, 3.2)
FeaturePlot(
    iCAF_obj,
    features = c("Ramp1"), cols = c("lightgrey", "red"),
    split.by = "Condition", order = TRUE, raster = TRUE, pt.size = 3
)
dev.off()

## 2.2. functional state of iCAFs ----
marker_iCAF_subset <- FindAllMarkers(
    iCAF_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

marker_iCAF_subset %>%
    group_by(cluster) %>%
    slice_max(n = 50, order_by = avg_log2FC) %>%
    as.data.frame() %>%
    split(.$cluster) %>%
    lapply(function(x) x$gene)

icaf_state_names_v2 <- c(
    "0" = "c01_Matrix_Smoc2_Cxcl12",
    "1" = "c02_Niche_Car3_Hhip",
    "2" = "Others",
    "3" = "Others",
    "4" = "c03_Inflammatory_Ptx3_Cxcl13",
    "5" = "c04_IFN_Ifit3_Cxcl10",
    "6" = "Others",
    "7" = "Others",
    "8" = "Others"
)
iCAF_obj <- RenameIdents(iCAF_obj, icaf_state_names_v2)
iCAF_obj$iCAF_subtype <- Idents(iCAF_obj) %>%
    as.character() %>%
    factor(levels = c(
        "c01_Matrix_Smoc2_Cxcl12", "c02_Niche_Car3_Hhip",
        "c03_Inflammatory_Ptx3_Cxcl13", "c04_IFN_Ifit3_Cxcl10", "Others"
    ))

mycolor$iCAF_subtype <- setNames(
    c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "gray"),
    levels(iCAF_obj$iCAF_subtype)
)

pdf("iCAF/UMAP_iCAF.Subtype.pdf", 7, 4)
DimPlot(iCAF_obj,
    raster = TRUE, pt.size = 2,
    label = FALSE, cols = mycolor$iCAF_subtype
) +
    labs(title = "iCAF Subtypes")
dev.off()

pdf("iCAF/UMAP_iCAF.Subtype.split.pdf", 12, 3)
DimPlot(iCAF_obj,
    group.by = "iCAF_subtype", split.by = "Condition", raster = TRUE,
    pt.size = 3, cols = mycolor$iCAF_subtype
) +
    labs(title = "iCAF Subtypes")
dev.off()

## 2.3. compositional changes of iCAF subtypes ----
identical(colnames(iCAF_obj), colnames(iCAF_obj))

plot_data <- as.data.frame(table(iCAF_obj$iCAF_subtype, iCAF_obj$Condition))
colnames(plot_data) <- c("iCAF_subtype", "Condition", "Count")

pdf("iCAF/Barplot.iCAF_subtype.Condition.pdf", 8, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = iCAF_subtype, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = iCAF_subtype)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$iCAF_subtype) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot_data <- as.data.frame(table(iCAF_obj$iCAF_subtype, iCAF_obj$Sample))
colnames(plot_data) <- c("iCAF_subtype", "Sample", "Count")
plot_data$Condition <- gsub("_[0-9]$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(iCAF_obj$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("iCAF/Boxplot.Fraction.iCAF_subtype.pdf", 6, 6)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = iCAF_subtype), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 2.5, ref.group = "CTRL",
    ) +
    facet_wrap(~iCAF_subtype, ncol = 3, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$iCAF_subtype) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot_data, temp)

# roe
roe <- Startrac::calTissueDist(
    iCAF_obj@meta.data,
    byPatient = FALSE,
    colname.cluster = "iCAF_subtype",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[setdiff(rownames(roe), "Doublet"), ]
pdf("iCAF/Heatmap.Roe.iCAF_subtype.Condition.pdf", 4, 2.5)
pheatmap::pheatmap(
    roe,
    color = mycolor$Roe,
    breaks = seq(0, 2, length.out = 101),
    scale = "none",
    display_numbers = TRUE,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
)
dev.off()
rm(roe)

## save data
gc()
rm(stromal_obj)
save.image("03.Stromal.RData")
