setwd("/thinker/tb269store/user/liuzhenyu/Project/BLCA_Neuro/BLCA_Neuro/results/03.Stromal")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

stromal.obj <- readRDS("../01.Cell_Type/sub3_Stromal/stromal.obj.rds")
mycolor <- readRDS("../01.Cell_Type/mycolor.rds")

# 1. Differential analysis of each cell type ----
dir.create("DEG_Cell_Type")

## 1.1. Identify DEGs ----
load("DEG_Cell_Type/DEG.stromal.rda")
DEG.stromal.Capsaicin <- list()
DEG.stromal.LC_6OHDA <- list()

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
        ident.1 = "LC_6OHDA", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0.25, logfc.threshold = 0.25,
        max.cells.per.ident = 1000
    )
    temp <- temp[order(temp$avg_log2FC, decreasing = TRUE), ]
    temp$gene <- rownames(temp)
    temp$Cell_Type_2 <- cell_type
    DEG.stromal.LC_6OHDA[[cell_type]] <- temp
}
rm(stromal.sub, cell_type, temp)

dir.create("DEG_Cell_Type/tables")
lapply(names(DEG.stromal.Capsaicin), function(x) {
    write.csv(DEG.stromal.Capsaicin[[x]],
        file = paste0("DEG_Cell_Type/tables/DEG_Capsaicin_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})
lapply(names(DEG.stromal.LC_6OHDA), function(x) {
    write.csv(DEG.stromal.LC_6OHDA[[x]],
        file = paste0("DEG_Cell_Type/tables/DEG_LC_6OHDA_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})

plot.data <- rbind(
    do.call(rbind, DEG.stromal.Capsaicin) %>% mutate(Comparison = "Capsaicin_vs_CTRL"),
    do.call(rbind, DEG.stromal.LC_6OHDA) %>% mutate(Comparison = "LC_6OHDA_vs_CTRL")
)
plot.data$DEG <- "Stable"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC > 0.25)] <- "Up"
plot.data$DEG[which(plot.data$p_val_adj < 0.05 & plot.data$avg_log2FC < -0.25)] <- "Down"
table(plot.data$Cell_Type_2, plot.data$Comparison, plot.data$DEG)

plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(stromal.obj$Cell_Type_2))
pdf("DEG_Cell_Type/Volcano.Capsaicin_vs_CTRL.pdf", 7.5, 4.5)
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

pdf("DEG_Cell_Type/Volcano.LC_6OHDA_vs_CTRL.pdf", 7.5, 4.5)
ggplot(plot.data %>% filter(Comparison == "LC_6OHDA_vs_CTRL"), aes(x = avg_log2FC, y = -log10(p_val_adj))) +
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

## 1.2. GO enrichment analysis ----
GO.DEG.Capsaicin.Up <- lapply(DEG.stromal.Capsaicin, function(x) {
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

GO.DEG.Capsaicin.Down <- lapply(DEG.stromal.Capsaicin, function(x) {
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

GO.DEG.LC_6OHDA.Up <- lapply(DEG.stromal.LC_6OHDA, function(x) {
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

GO.DEG.LC_6OHDA.Down <- lapply(DEG.stromal.LC_6OHDA, function(x) {
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
