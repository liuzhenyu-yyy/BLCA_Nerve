setwd("/Users/zhenyuliu/Research/Project/BLCA_Nerve/results/03.Stromal")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(parallel)

stromal_obj <- readRDS("../rds/stromal_obj.rds")
mycolor <- readRDS("../rds/mycolor.rds")

# load("03.Stromal.RData")

# 1. Differential analysis of each cell type ----
## 1.1. Identify DEGs ----
table(stromal_obj$Cell_Type_2, stromal_obj$Condition)
DEG_list_conditions <- mclapply(levels(stromal_obj$Condition)[-1], function(x) {
    DEG_list <- list()

    for (cell_type in levels(stromal_obj$Cell_Type_2)) {
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
    plot_data$Cell_Type_2 <- factor(plot_data$Cell_Type_2, levels = levels(stromal_obj$Cell_Type_2))

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
        facet_wrap(~Cell_Type_2, ncol = 4) +
        theme_classic() +
        labs(title = paste0("Volcano plot of DEGs in each cell type, ", x, " vs CTRL")) +
        xlab("Average log2 fold change") +
        ylab("-log10 adjusted p value") +
        xlim(-x_max * 0.7, x_max * 0.7) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
        theme(legend.position = "none")
    pdf(paste0("Volcano.", x, "_vs_CTRL.pdf"), 6.5, 4.5)
    plot(p)
    dev.off()
})

# only iCAF
plot.data <- do.call(rbind, DEG_list_conditions$Capsaicin) %>%
    filter(Cell_Type_2 == "iCAF") %>%
    mutate(Comparison = paste0("Capsaicin_vs_CTRL"))
plot.data$label <- ifelse(plot.data$gene %in% c("Pdgfd", "Runx2", "Gas6"), plot.data$gene, NA)

pdf(paste0("Volcano.iCAF.Capsaicin_vs_CTRL.pdf"), 6.5, 4.5)
ggplot(plot.data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(aes(color = is.na(plot.data$label) %>% as.numeric() %>% as.character()), alpha = 0.5, size = 0.5) +
    scale_color_manual(values = c("0" = "red", "1" = "gray")) +
    ggrepel::geom_text_repel(
        aes(label = label),
        na.rm = TRUE,
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = "grey50",
        size = 3
    ) +
    theme_classic() +
    xlab("Average log2 fold change") +
    ylab("-log10 adjusted p value") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
    theme(legend.position = "none")
dev.off()
rm(plot.data)

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

# plot selected cell types
dir.create("GO_Capsaicin_Selected")

# iCAF

term.selected <- c("GO:0030098", "GO:0051963", "GO:0045787", "GO:0071772", "GO:0048009", "GO:0050670")
plot.data <- GO_list$Capsaicin$iCAF$up@result[term.selected, ] %>%
    arrange(FoldEnrichment) %>%
    mutate(
        Description = stringr::str_wrap(Description, width = 42),
        Description = factor(Description, levels = Description),
        logP = -log10(p.adjust)
    )
pdf("GO_Capsaicin_Selected/Dot.GO_Capsaicin_in_iCAF_UP.pdf", 6, 2.5)
ggplot(plot.data, aes(x = FoldEnrichment, y = Description)) +
    geom_segment(aes(x = 0, xend = FoldEnrichment, yend = Description),
        color = "grey80", linewidth = 0.7
    ) +
    geom_point(aes(fill = logP),
        shape = 21, size = 4, color = "grey80", stroke = 0.4
    ) +
    scale_fill_gradient(
        low = "#bcd8c3", high = "#367a80",
        name = expression(-log[10]("adjusted P"))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
        title = "iCAF Upregulated genes",
        subtitle = "Capsaicin vs. CTRL",
        x = "Fold enrichment", y = NULL
    ) +
    theme_classic(base_size = 11)
dev.off()

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

plot_data <- DotPlot(
    stromal_obj,
    features = c("Ramp1"),
    group.by = "Cell_Type_2",
    split.by = "Condition",
    cols = mycolor$Condition
)$data
plot_data$Cell_Type_2 <- gsub(paste0(names(mycolor$Condition), collapse = "|"), "", plot_data$id)
plot_data$Cell_Type_2 <- gsub("_$", "", plot_data$Cell_Type_2)
plot_data$Condition <- gsub(paste0(names(mycolor$Cell_Type_2), collapse = "|"), "", plot_data$id)
plot_data$Condition <- gsub("^_", "", plot_data$Condition)
plot_data$Cell_Type_2 <- factor(plot_data$Cell_Type_2, levels = names(mycolor$Cell_Type_2))
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("DotPlot_Ramp1.pdf", 5, 3.5)
ggplot(plot_data, aes(x = Condition, y = Cell_Type_2)) +
    geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
    scale_fill_viridis_c() +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
    )
dev.off()

## 1.4. compare DEG & GO across conditions ----
names(DEG_list_genes)
dir.create("Compare_Conditions")
condition_selected <- c("Capsaicin", "anti_PD1", "Rime")

DEG_list_gene_cell_type <- lapply(names(DEG_list_genes$Capsaicin), function(cell_typex) {
    list(
        up = sapply(names(DEG_list_genes), function(condition) {
            DEG_list_genes[[condition]][[cell_typex]]$up
        }) %>%
            setNames(names(DEG_list_genes)) %>%
            .[condition_selected] %>%
            .[sapply(., length) > 0],
        down = sapply(names(DEG_list_genes), function(condition) {
            DEG_list_genes[[condition]][[cell_typex]]$down
        }) %>%
            setNames(names(DEG_list_genes)) %>%
            .[condition_selected] %>%
            .[sapply(., length) > 0]
    )
}) %>%
    setNames(names(DEG_list_genes$Capsaicin))
DEG_list_gene_cell_type$Vascular_Endothelial$up

# compare DEGs
lapply(names(DEG_list_gene_cell_type)[2:4], function(x) {
    temp <- DEG_list_gene_cell_type[[x]]$up
    if (length(temp) > 1) {
        pdf(paste0("Compare_Conditions/Venn_DEGs_Up_", x, ".pdf"), 5, 5)
        plot(Vennerable::Venn(temp),
            doWeights = TRUE,
            show = list(Faces = FALSE, DarkMatter = FALSE)
        )
        dev.off()
    }

    temp <- DEG_list_gene_cell_type[[x]]$down
    if (length(temp) > 1) {
        pdf(paste0("Compare_Conditions/Venn_DEGs_Down_", x, ".pdf"), 5, 5)
        plot(Vennerable::Venn(temp),
            doWeights = TRUE,
            show = list(Faces = FALSE, DarkMatter = FALSE)
        )
        dev.off()
    }
})

temp <- DEG_list_gene_cell_type$iCAF$down
temp$Capsaicin %>%
    setdiff(temp$anti_PD1) %>%
    setdiff(temp$Rime) %>%
    sort()
rm(temp)

# compare GO
names(DEG_list_gene_cell_type)
GO_list_compare <- mclapply(c("iCAF", "Vascular_Endothelial"), function(x) {
    list(
        up = compareCluster(DEG_list_gene_cell_type[[x]]$up,
            fun = "enrichGO",
            OrgDb = org.Mm.eg.db,
            keyType = "SYMBOL",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05
        ) %>% simplify(),
        down = compareCluster(DEG_list_gene_cell_type[[x]]$down,
            fun = "enrichGO",
            OrgDb = org.Mm.eg.db,
            keyType = "SYMBOL",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05
        ) %>% simplify()
    )
}, mc.cores = 2) %>%
    setNames(c("iCAF", "Vascular_Endothelial"))

lapply(names(GO_list_compare), function(x) {
    if (!is.null(GO_list_compare[[x]]$up)) {
        pdf(paste0("Compare_Conditions/GO_Compare_Up_", x, ".pdf"), 7, 6)
        p <- dotplot(GO_list_compare[[x]]$up, showCategory = 15, label_format = 80) +
            theme_classic() +
            labs(title = paste("GO of up-regulated genes in ", x, sep = ""))
        plot(p)
        dev.off()
    }
    if (!is.null(GO_list_compare[[x]]$down)) {
        pdf(paste0("Compare_Conditions/GO_Compare_Down_", x, ".pdf"), 7, 6)
        p <- dotplot(GO_list_compare[[x]]$down, showCategory = 15, label_format = 80) +
            theme_classic() +
            labs(title = paste("GO of down-regulated genes in ", x, sep = ""))
        plot(p)
        dev.off()
    }
})
rm(p, plot_data, label_data, x, cell_typex, temp)

# # 2. iCAF functions ----
# ## 2.1. identify iCAF subsets ----
# dir.create("iCAF")
# iCAF_obj <- subset(stromal_obj, subset = Cell_Type_2 == "iCAF") # 7254 cells
# rm(stromal_obj)

# iCAF_obj <- iCAF_obj %>%
#     NormalizeData() %>%
#     FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
#     ScaleData() %>%
#     RunPCA(npcs = 30) %>%
#     FindNeighbors(dims = 1:20) %>%
#     RunUMAP(dims = 1:20)

# iCAF_obj <- FindClusters(iCAF_obj, resolution = 0.3)

# pdf("iCAF/UMAP_iCAF.Subclusters.pdf", 5, 4)
# DimPlot(iCAF_obj, raster = TRUE, pt.size = 2, group.by = "seurat_clusters",
#  label = TRUE, cols = mycolor$Cluster) +
#     labs(title = "iCAF Subclusters")
# dev.off()

# pdf("iCAF/UMAP_iCAF.Subclusters.split.pdf", 9, 6)
# DimPlot(iCAF_obj,
#     split.by = "Condition", raster = TRUE, ncol = 3,
#     cols = mycolor$Cluster, pt.size = 3, label = FALSE
# ) +
#     labs(title = "iCAF Subclusters")
# dev.off()

# pdf("iCAF/UMAP_iCAF.Ramp1.split.pdf", 16, 3.2)
# FeaturePlot(
#     iCAF_obj,
#     features = c("Ramp1"), cols = c("lightgrey", "red"),
#     split.by = "Condition", order = TRUE, raster = TRUE, pt.size = 3
# )
# dev.off()

# ## 2.2. functional state of iCAFs ----
# marker_iCAF_subset <- FindAllMarkers(
#     iCAF_obj,
#     only.pos = TRUE,
#     min.pct = 0.25,
#     logfc.threshold = 0.25
# )

# marker_iCAF_subset %>%
#     group_by(cluster) %>%
#     slice_max(n = 100, order_by = avg_log2FC) %>%
#     as.data.frame() %>%
#     split(.$cluster) %>%
#     lapply(function(x) x$gene)

# Idents(iCAF_obj) <- iCAF_obj$seurat_clusters
# icaf_state_names_v2 <- c(
#     "0" = "c01_Matrix_Smoc2_Cxcl12",
#     "1" = "c02_Niche_Car3_Hhip",
#     "2" = "c09_Epithelial_Krt8_Lama3",
#     "3" = "c06_Myeloid_Adgre1_Mrc1",
#     "4" = "c03_Inflammatory_Ptx3_Cxcl13",
#     "5" = "c04_IFN_Ifit3_Cxcl10",
#     "6" = "c08_Lymphoid_Cd3g_Gzmb",
#     "7" = "c07_Vascular_Pecam1_Vwf",
#     "8" = "c05_Mesothelial_Upk3b_Wt1"
# )
# iCAF_obj <- RenameIdents(iCAF_obj, icaf_state_names_v2)
# iCAF_obj$iCAF_subtype <- Idents(iCAF_obj) %>%
#     as.character() %>%
#     factor(levels = c(
#         "c01_Matrix_Smoc2_Cxcl12", "c02_Niche_Car3_Hhip",
#         "c03_Inflammatory_Ptx3_Cxcl13", "c04_IFN_Ifit3_Cxcl10",
#         "c05_Mesothelial_Upk3b_Wt1", "c06_Myeloid_Adgre1_Mrc1",
#         "c07_Vascular_Pecam1_Vwf", "c08_Lymphoid_Cd3g_Gzmb", "c09_Epithelial_Krt8_Lama3"
#     ))

# mycolor$iCAF_subtype <- setNames(
#     ggsci::pal_npg()(9),
#     levels(iCAF_obj$iCAF_subtype)
# )

# pdf("iCAF/UMAP_iCAF.Subtype.pdf", 7, 4)
# DimPlot(iCAF_obj,
#     raster = TRUE, pt.size = 2, group.by = "iCAF_subtype",
#     label = FALSE, cols = mycolor$iCAF_subtype
# ) +
#     labs(title = "iCAF Subtypes")
# dev.off()

# pdf("iCAF/UMAP_iCAF.Subtype.split.pdf", 12, 3)
# DimPlot(iCAF_obj,
#     group.by = "iCAF_subtype", split.by = "Condition", raster = TRUE,
#     pt.size = 3, cols = mycolor$iCAF_subtype
# ) +
#     labs(title = "iCAF Subtypes")
# dev.off()

# ## 2.3. compositional changes of iCAF subtypes ----
# identical(colnames(iCAF_obj), colnames(iCAF_obj))

# plot_data <- as.data.frame(table(iCAF_obj$iCAF_subtype, iCAF_obj$Condition))
# colnames(plot_data) <- c("iCAF_subtype", "Condition", "Count")
# plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

# pdf("iCAF/Barplot.iCAF_subtype.Condition.pdf", 9, 5)
# patchwork::wrap_plots(
#     ggplot(plot_data, aes(x = iCAF_subtype, y = Count, fill = Condition)) +
#         geom_bar(stat = "identity", position = "fill") +
#         scale_fill_manual(values = mycolor$Condition) +
#         theme_classic() +
#         labs(title = "Cell Type composition across Conditions") +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#     ggplot(plot_data, aes(x = Condition, y = Count, fill = iCAF_subtype)) +
#         geom_bar(stat = "identity", position = "fill") +
#         scale_fill_manual(values = mycolor$iCAF_subtype) +
#         theme_classic() +
#         labs(title = "Cell Type composition across Conditions") +
#         theme(axis.text.x = element_text(angle = 45, hjust = 1)),
#     ncol = 2
# )
# dev.off()

# # fraction compare
# plot_data <- as.data.frame(table(iCAF_obj$iCAF_subtype, iCAF_obj$Sample))
# colnames(plot_data) <- c("iCAF_subtype", "Sample", "Count")
# plot_data$Condition <- gsub("_[0-9]$", "", plot_data$Sample)
# table(plot_data$Condition)
# table(plot_data$Sample)

# temp <- table(iCAF_obj$Sample)
# plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
# plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

# pdf("iCAF/Boxplot.Fraction.iCAF_subtype.pdf", 8, 5.5)
# ggplot(plot_data, aes(x = Condition, y = Fraction)) +
#     geom_boxplot(aes(fill = iCAF_subtype), outlier.shape = NA, alpha = 0.5) +
#     geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
#     ggpubr::stat_compare_means(
#         label = "p.signif", label.y.npc = 0.98, size = 2.5, ref.group = "CTRL",
#     ) +
#     facet_wrap(~iCAF_subtype, ncol = 5, scales = "free") +
#     labs(title = "Cell Type composition across Samples") +
#     ylab("Fraction (%)") +
#     scale_fill_manual(values = mycolor$iCAF_subtype) +
#     theme_classic() +
#     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#     NoLegend()
# dev.off()
# rm(plot_data, temp)

# # roe
# roe <- Startrac::calTissueDist(
#     iCAF_obj@meta.data,
#     byPatient = FALSE,
#     colname.cluster = "iCAF_subtype",
#     colname.patient = "Sample",
#     colname.tissue = "Condition",
#     method = "chisq",
#     min.rowSum = 0
# )
# roe <- roe[setdiff(rownames(roe), "Doublet"), ]
# pdf("iCAF/Heatmap.Roe.iCAF_subtype.Condition.pdf", 4, 3.5)
# pheatmap::pheatmap(
#     roe,
#     color = mycolor$Roe,
#     breaks = seq(0, 2, length.out = 101),
#     scale = "none",
#     display_numbers = TRUE,
#     cluster_rows = FALSE,
#     cluster_cols = FALSE,
# )
# dev.off()
# rm(roe)

## save data
gc()
rm(stromal_obj)
save.image("03.Stromal.RData")
