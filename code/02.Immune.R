setwd("/thinker/tb269store/user/liuzhenyu/Project/BLCA_Neuro/BLCA_Neuro/results/02.Immune")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

immune.obj <- readRDS("../01.Cell_Type/sub1_Immune/immune.obj.rds")
mycolor <- readRDS("../01.Cell_Type/mycolor.rds")

# 1. Differential analysis of each cell type ----
dir.create("DEG_Cell_Type")

## 1.1. Identify DEGs ----
DEG.immune.Capsaicin <- list()
DEG.immune.LC_6OHDA <- list()

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
        ident.1 = "LC_6OHDA", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0.25, logfc.threshold = 0.25,
        max.cells.per.ident = 1000
    )
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Cell_Type_2 <- cell_type
    DEG.immune.LC_6OHDA[[cell_type]] <- temp
}
rm(immune.sub, cell_type, temp)

dir.create("DEG_Cell_Type/tables")
lapply(names(DEG.immune.Capsaicin), function(x) {
    write.csv(DEG.immune.Capsaicin[[x]],
        file = paste0("DEG_Cell_Type/tables/DEG_Capsaicin_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})
lapply(names(DEG.immune.LC_6OHDA), function(x) {
    write.csv(DEG.immune.LC_6OHDA[[x]],
        file = paste0("DEG_Cell_Type/tables/DEG_LC_6OHDA_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})
save(DEG.immune.Capsaicin, DEG.immune.LC_6OHDA, file = "DEG_Cell_Type/DEG.immune.rda")

# volcano plot
plot.data <- rbind(
    do.call(rbind, DEG.immune.Capsaicin) %>% mutate(Comparison = "Capsaicin_vs_CTRL"),
    do.call(rbind, DEG.immune.LC_6OHDA) %>% mutate(Comparison = "LC_6OHDA_vs_CTRL")
)
plot.data$DEG <- "Stable"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC > 0.25)] <- "Up"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC < -0.25)] <- "Down"
table(plot.data$Cell_Type_2, plot.data$Comparison, plot.data$DEG)

plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(immune.obj$Cell_Type_2))
pdf("DEG_Cell_Type/Volcano.Capsaicin_vs_CTRL.pdf", 9, 7)
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

pdf("DEG_Cell_Type/Volcano.LC_6OHDA_vs_CTRL.pdf", 9, 7)
ggplot(plot.data %>% filter(Comparison == "LC_6OHDA_vs_CTRL"), aes(x = avg_log2FC, y = -log10(p_val_adj))) +
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

## 1.2. GO enrichment analysis ----
GO.DEG.Capsaicin.Up <- lapply(DEG.immune.Capsaicin, function(x) {
    temp <- x %>%
        filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
        pull(gene)
    GO.res <- enrichGO(
        gene = temp, OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05
    )
    if (is.null(GO.res)) {
        return(NULL)
    }
    simplify(GO.res)
})

GO.DEG.Capsaicin.Down <- lapply(DEG.immune.Capsaicin, function(x) {
    temp <- x %>%
        filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
        pull(gene)
    GO.res <- enrichGO(
        gene = temp, OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05
    )
    if (is.null(GO.res)) {
        return(NULL)
    }
    simplify(GO.res)
})

GO.DEG.LC_6OHDA.Up <- lapply(DEG.immune.LC_6OHDA, function(x) {
    temp <- x %>%
        filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
        pull(gene)
    GO.res <- enrichGO(
        gene = temp, OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05
    )
    if (is.null(GO.res)) {
        return(NULL)
    }
    simplify(GO.res)
})

GO.DEG.LC_6OHDA.Down <- lapply(DEG.immune.LC_6OHDA, function(x) {
    temp <- x %>%
        filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
        pull(gene)
    GO.res <- enrichGO(
        gene = temp, OrgDb = org.Mm.eg.db,
        keyType = "SYMBOL", ont = "BP",
        pAdjustMethod = "BH", qvalueCutoff = 0.05
    )
    if (is.null(GO.res)) {
        return(NULL)
    }
    simplify(GO.res)
})

# output GO results
dir.create("DEG_Cell_Type/GO_Capasaicin_Up")
lapply(names(GO.DEG.Capsaicin.Up), function(x) {
    if (is.null(GO.DEG.Capsaicin.Up[[x]])) {
        return(NULL)
    }
    write.csv(GO.DEG.Capsaicin.Up[[x]]@result,
        file = paste0("DEG_Cell_Type/GO_Capasaicin_Up/GO_Capsaicin_vs_CTRL_Up_in_", x, ".csv"),
        row.names = TRUE
    )
    pdf(paste0("DEG_Cell_Type/GO_Capasaicin_Up/GO_Capsaicin_vs_CTRL_Up_in_", x, ".pdf"), 7, 5)
    p <- dotplot(GO.DEG.Capsaicin.Up[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("Capsaicin vs CTRL - Upregulated in ", x, sep = ""))
    plot(p)
    dev.off()
})

dir.create("DEG_Cell_Type/GO_Capasaicin_Down")
lapply(names(GO.DEG.Capsaicin.Down), function(x) {
    if (is.null(GO.DEG.Capsaicin.Down[[x]])) {
        return(NULL)
    }
    write.csv(GO.DEG.Capsaicin.Down[[x]]@result,
        file = paste0("DEG_Cell_Type/GO_Capasaicin_Down/GO_Capsaicin_vs_CTRL_Down_in_", x, ".csv"),
        row.names = TRUE
    )
    pdf(paste0("DEG_Cell_Type/GO_Capasaicin_Down/GO_Capsaicin_vs_CTRL_Down_in_", x, ".pdf"), 7, 5)
    p <- dotplot(GO.DEG.Capsaicin.Down[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("Capsaicin vs CTRL - Downregulated in ", x, sep = ""))
    plot(p)
    dev.off()
})

dir.create("DEG_Cell_Type/GO_LC_6OHDA_Up")
lapply(names(GO.DEG.LC_6OHDA.Up), function(x) {
    if (is.null(GO.DEG.LC_6OHDA.Up[[x]])) {
        return(NULL)
    }
    write.csv(GO.DEG.LC_6OHDA.Up[[x]]@result,
        file = paste0("DEG_Cell_Type/GO_LC_6OHDA_Up/GO_LC_6OHDA_vs_CTRL_Up_in_", x, ".csv"),
        row.names = TRUE
    )
    pdf(paste0("DEG_Cell_Type/GO_LC_6OHDA_Up/GO_LC_6OHDA_vs_CTRL_Up_in_", x, ".pdf"), 7, 5)
    p <- dotplot(GO.DEG.LC_6OHDA.Up[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("LC_6OHDA vs CTRL - Upregulated in ", x, sep = ""))
    plot(p)
    dev.off()
})

dir.create("DEG_Cell_Type/GO_LC_6OHDA_Down")
lapply(names(GO.DEG.LC_6OHDA.Down), function(x) {
    if (is.null(GO.DEG.LC_6OHDA.Down[[x]])) {
        return(NULL)
    }
    write.csv(GO.DEG.LC_6OHDA.Down[[x]]@result,
        file = paste0("DEG_Cell_Type/GO_LC_6OHDA_Down/GO_LC_6OHDA_vs_CTRL_Down_in_", x, ".csv"),
        row.names = TRUE
    )
    pdf(paste0("DEG_Cell_Type/GO_LC_6OHDA_Down/GO_LC_6OHDA_vs_CTRL_Down_in_", x, ".pdf"), 7, 5)
    p <- dotplot(GO.DEG.LC_6OHDA.Down[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("LC_6OHDA vs CTRL - Downregulated in ", x, sep = ""))
    plot(p)
    dev.off()
})

## 1.3. GSEA on specific cell types ----
dir.create("DEG_Cell_Type/GSEA_Capsaicin")

# run GSEA
table(immune.obj$Cell_Type_2, immune.obj$Condition)
cell_type <- c("T_CD8")

gseGO.Capsaicin <- sapply(c("T_CD4", "T_CD8", "Treg", "Neutrophil"), function(cell_type) {
    message(Sys.time(), ": Processing: ", cell_type)
    immune.sub <- subset(immune.obj, subset = Cell_Type_2 == cell_type)

    message(Sys.time(), ": Identify DEGs: ", cell_type)
    temp <- FindMarkers(immune.sub,
        ident.1 = "Capsaicin", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0, logfc.threshold = 0,
        max.cells.per.ident = 1000
    )
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Cell_Type_2 <- cell_type

    message(Sys.time(), ": Running GSEA: ", cell_type)
    gene.list <- temp$avg_log2FC %>% setNames(temp$gene)

    gseGO.res <- gseGO(
        geneList = gene.list,
        OrgDb = org.Mm.eg.db,
        ont = "BP",
        keyType = "SYMBOL",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
    )
    gseGO.res
    # simplify(gseGO.res)
})

# output GSEA results
lapply(gseGO.Capsaicin, function(x) nrow(x@result))

lapply(names(gseGO.Capsaicin) %>% setdiff("T_CD8"), function(x) {
    if (is.null(gseGO.Capsaicin[[x]])) {
        return(NULL)
    }
    write.csv(gseGO.Capsaicin[[x]]@result,
        file = paste0("DEG_Cell_Type/GSEA_Capsaicin/GSEA_Capsaicin_in_", x, ".csv"),
        row.names = TRUE
    )
    pdf(paste0("DEG_Cell_Type/GSEA_Capsaicin/Dot.GSEA_Capsaicin_in_", x, ".pdf"), 7, 5)
    p <- dotplot(gseGO.Capsaicin[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("GSEA - Capsaicin vs CTRL in ", x, sep = ""))
    plot(p)
    dev.off()

})

# gseaplot2 - CD4 T
gseGO.Capsaicin$T_CD4@result %>%
    filter(grepl("T cell", Description)) %>%
    arrange(-NES) %>%
    dplyr::select(Description, NES, pvalue, p.adjust)
term.selected <- c("GO:0050798", "GO:2000406", "GO:0045582", "GO:0002294", "GO:0002369", "GO:0050870")

p.list <- lapply(term.selected, function(term) {
    gseaplot2(gseGO.Capsaicin$T_CD4,
        geneSetID = term,
        color = c("#e1423b"),
        title = paste(
            gseGO.Capsaicin$T_CD4@result[term, ]$Description, "\n",
            "NES =", round(gseGO.Capsaicin$T_CD4@result[term, ]$NES, 2),
            ", p.adj =", signif(gseGO.Capsaicin$T_CD4@result[term, ]$p.adjust, 2)),
    )
})
length(p.list)

pdf("DEG_Cell_Type/GSEA_Capsaicin/GSEA_Capsaicin_in_T_CD4_gseaplot.pdf", 5, 4)
for (i in seq_along(p.list)) {
    print(p.list[[i]])
}
dev.off()

# gseaplot2 - Neutrophil
gseGO.Capsaicin$Neutrophil@result %>%
    filter(grepl("T cell", Description)) %>%
    arrange(NES) %>%
    dplyr::select(Description, NES, pvalue, p.adjust)
term.selected <- c("GO:0050777", "GO:0002820", "GO:0002838")

p.list <- lapply(term.selected, function(term) {
    gseaplot2(gseGO.Capsaicin$Neutrophil,
        geneSetID = term,
        color = c("#458ec3"),
        title = paste(
            gseGO.Capsaicin$Neutrophil@result[term, ]$Description, "\n",
            "NES =", round(gseGO.Capsaicin$Neutrophil@result[term, ]$NES, 2),
            ", p.adj =", signif(gseGO.Capsaicin$Neutrophil@result[term, ]$p.adjust, 2)),
        )
})
length(p.list)

pdf("DEG_Cell_Type/GSEA_Capsaicin/GSEA_Capsaicin_in_Neutrophil_gseaplot.pdf", 5, 4)
for (i in seq_along(p.list)) {
    print(p.list[[i]])
}
dev.off()
rm(p.list, term.selected, i)

# 2. T cell states ----
## 2.1. Project T cell states ----
library(ProjecTILs)

table(immune.obj$Cell_Type_2)
T.obj <- subset(immune.obj, subset = Cell_Type_2 %in% c("T_CD4", "T_CD4_Cycling", "Treg", "T_CD8", "T_CD8_Cycling"))
rm(immune.obj)

ref <- load.reference.map()

T.obj.projected <- Run.ProjecTILs(T.obj,
    ref = ref, filter.cells = FALSE,
    split.by = "Sample",
    ncores = 20
)

levels(ref$functional.cluster)
mycolor$T_subtype <- c(
    "#edbe2a", "#A58AFF", "#53B400", "#F8766D",
    "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd"
) %>% setNames(levels(ref$functional.cluster))

p.list <- list()
sample_names <- (levels(T.obj$Condition))
for (sample_i in seq_along(sample_names)) {
    sample <- sample_names[sample_i]
    p.list[[sample_i]] <- plot.projection(ref,
        T.obj.projected[, T.obj.projected$Condition == sample],
        linesize = 0.3, pointsize = 0.6,
        raster = TRUE
    ) + ggtitle(sample) + theme(aspect.ratio = 1)
}
p.list$ref <- DimPlot(ref,
    group.by = "functional.cluster", label = TRUE,
    cols = mycolor$T_subtype, raster = TRUE
) +
    ggtitle("Reference") + theme(aspect.ratio = 1)

pdf("Projection_of_T_cells.pdf", 12, 10)
patchwork::wrap_plots(p.list, ncol = 2)
dev.off()
rm(p.list, sample, sample_i, sample_names)

table(T.obj.projected$functional.cluster)

p.list <- plot.states.radar(ref = ref, T.obj.projected, min.cells = 20, return.as.list = TRUE)
length(p.list)

pdf("Radar.T.Projection.pdf", 12, 10)
patchwork::wrap_plots(p.list[1:9], ncol = 3)
dev.off()

rm(p.list, Hs2Mm.convert.table, ref)

## 2.2. Composition of T cell states ----
T.obj.projected$Condition <- factor(T.obj.projected$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))
T.obj.projected$functional.cluster <- factor(T.obj.projected$functional.cluster,
    levels = names(mycolor$T_subtype)
)
plot.data <- as.data.frame(table(T.obj.projected$functional.cluster, T.obj.projected$Condition))
colnames(plot.data) <- c("functional.cluster", "Condition", "Count")

pdf("Barplot.T_subtype.Condition.pdf", 8, 5)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = functional.cluster, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot.data, aes(x = Condition, y = Count, fill = functional.cluster)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$T_subtype) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot.data <- as.data.frame(table(T.obj.projected$functional.cluster, T.obj.projected$Sample))
colnames(plot.data) <- c("T_subtype", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(T.obj.projected$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))

pdf("Boxplot.Fraction.T_subtype.pdf", 8, 4.7)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = T_subtype), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_6OHDA"), c("CTRL", "LC_6OHDA"))
    ) +
    facet_wrap(~T_subtype, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$T_subtype) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot.data, temp)

## 2.3. T state marker expression ----
identical(colnames(T.obj.projected), colnames(T.obj))
T.obj$T_state <- T.obj.projected$functional.cluster

gene.selected <- c(
    "Foxp3", "Il2ra", "Il10", "Tgfb1", # Treg
    "Gzmb", "Prf1", "Ifng", "Tnf", "Fasl", # Cytotoxicity
    "Pdcd1", "Lag3", "Havcr2", "Tigit", "Ctla4", "Tox" # exhaustion
)

plot.data <- DotPlot(T.obj,
    features = gene.selected,
    group.by = "Condition"
)$data
plot.data$group <- rep(rep(c("Treg(Supression)", "Cytotoxicity", "Exhaustion"), times = c(4, 5, 6)),
    times = 3
)

pdf("Dot.T_state.markers.Condition.pdf", 7, 6)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
        scale_fill_viridis_c() +
        facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
        theme_classic() +
        labs(x = "Gene", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle("Expression of T cell state markers"),
    ggplot(plot.data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
        scale_fill_viridis_c() +
        facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
        theme_classic() +
        labs(x = "Gene", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle("Scaled expression of T cell state markers"),
    ncol = 1
)
dev.off()

plot.data <- DotPlot(T.obj,
    features = gene.selected, scale = FALSE,
    group.by = "Condition", split.by = "T_state", cols = mycolor$Cluster %>% unname()
)$data
plot.data$group <- "none"
plot.data[plot.data$features.plot %in% c("Foxp3", "Il2ra", "Il10", "Tgfb1"), "group"] <- "Treg(Supression)"
plot.data[plot.data$features.plot %in% c("Gzmb", "Prf1", "Ifng", "Tnf", "Fasl"), "group"] <- "Cytotoxicity"
plot.data[plot.data$features.plot %in% c("Pdcd1", "Lag3", "Havcr2", "Tigit", "Tox", "Ctla4"), "group"] <- "Exhaustion"
table(plot.data$group)

plot.data$Condition <- gsub("_.+?$", "", plot.data$id)
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC"))
plot.data$T_state <- gsub("Capsaicin_|CTRL_|LC_6OHDA_", "", plot.data$id)
plot.data$T_state <- factor(plot.data$T_state, levels = names(mycolor$T_subtype))
plot.data <- plot.data[!is.na(plot.data$T_state), ]

pdf("Dot.T_state.markers.Condition.Split.pdf", 6, 7)
ggplot(plot.data, aes(x = features.plot, y = Condition)) +
    geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
    scale_fill_viridis_c() +
    facet_grid(cols = vars(group), rows = vars(T_state), scales = "free_x", space = "free_x") +
    theme_classic() +
    labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of T cell state markers across Conditions")
dev.off()

table(T.obj$T_state)
p.list <- lapply(levels(T.obj$T_state), function(state) {
    plot.data <- DotPlot(T.obj %>% subset(T_state == state),
        features = gene.selected,
        group.by = "Condition"
    )$data
    plot.data$group <- rep(rep(c("Treg(Supression)", "Cytotoxicity", "Exhaustion"), times = c(4, 5, 6)),
        times = 3
    )

    ggplot(plot.data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
        scale_fill_viridis_c() +
        facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
        theme_classic() +
        labs(x = "Gene", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
        ) +
        ggtitle(paste("Expression in", state))
})
pdf("Dot.T_state.markers.Condition.Split2.pdf", 10, 10)
patchwork::wrap_plots(p.list, ncol = 2)
dev.off()
rm(plot.data, state, p.list)

## 2.4. T state program ----
library(SignatuR)
data(SignatuR)

programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)
rm(SignatuR)

T.obj <- AddModuleScore(T.obj, features = programs)
colnames(T.obj@meta.data)[17:23] <- names(programs)

T.obj@meta.data %>%
    group_by(Condition) %>%
    summarise(across(all_of(names(programs)[5:7]), list(mean = mean)))

plot.data <- DotPlot(T.obj,
    features = names(programs) %>% setdiff(c("HeatShock", "IFN")),
    group.by = "Condition"
)$data

pdf("Dot.T_state.Signature.Condition.pdf", 6, 5)
patchwork::wrap_plots(
    ggplot(plot.data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
        scale_fill_viridis_c() +
        theme_classic() +
        labs(x = "Signature", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle("T cell state signatures"),
    ggplot(plot.data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
        scale_fill_viridis_c() +
        theme_classic() +
        labs(x = "Signature", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle("Scaled T cell state signatures"),
    ncol = 1
)
dev.off()

plot.data <- DotPlot(T.obj,
    features = names(programs) %>% setdiff(c("HeatShock", "IFN")), scale = TRUE,
    group.by = "Condition", split.by = "T_state", cols = mycolor$Cluster %>% unname()
)$data

plot.data$Condition <- gsub("_.+?$", "", plot.data$id)
plot.data$T_state <- gsub("Capsaicin_|CTRL_|LC_6OHDA_", "", plot.data$id)
plot.data$T_state <- factor(plot.data$T_state, levels = names(mycolor$T_subtype))
plot.data <- plot.data[!is.na(plot.data$T_state), ]

pdf("Dot.T_state.Signature.Condition.Split.pdf", 6, 7)
ggplot(plot.data, aes(x = features.plot, y = Condition)) +
    geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
    scale_fill_viridis_c() +
    facet_grid(rows = vars(T_state), scales = "free_x") +
    theme_classic() +
    labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of T cell state markers across Conditions")
dev.off()

table(T.obj$T_state)
p.list <- lapply(levels(T.obj$T_state), function(state) {
    plot.data <- DotPlot(T.obj %>% subset(T_state == state),
        features = names(programs) %>% setdiff(c("HeatShock", "IFN")),
        group.by = "Condition"
    )$data

    ggplot(plot.data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
        scale_fill_viridis_c() +
        theme_classic() +
        labs(x = "Gene", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
        ) +
        ggtitle(paste(state))
})
pdf("Dot.T_state.Signatures.Condition.Split2.pdf", 10, 6.5)
patchwork::wrap_plots(p.list, ncol = 3)
dev.off()
rm(plot.data, state, p.list)

Idents(T.obj)
pdf("Vln.Exhausion.pdf", 12, 7)
VlnPlot(T.obj,
    idents = names(mycolor$T_subtype),
    features = programs$Tcell.exhaustion,
    group.by = "T_state", split.by = "Condition",
    pt.size = 0, ncol = 3, cols = mycolor$Condition
)
dev.off()

# save data ----
gc()
save.image("02.Immune.RData")
load("02.Immune.RData")
