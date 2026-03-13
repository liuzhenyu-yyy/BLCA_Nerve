setwd("/Users/zhenyuliu/Research/Project/BLCA_Nerve/results/02.Immune")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(parallel)

immune_obj <- readRDS("../01.Cell_Type/rds/immune_obj.rds")
mycolor <- readRDS("../01.Cell_Type/rds/mycolor.rds")

# load("02.Immune.RData")

# 1. Differential analysis of each cell type ----
## 1.1. Identify DEGs ----
table(immune_obj$Cell_Type_2, immune_obj$Condition)

DEG_list_conditions <- mclapply(levels(immune_obj$Condition)[-1], function(x) {
    DEG_list <- list()

    for (cell_type in levels(immune_obj$Cell_Type_2)) {
        message(Sys.time(), ": Differential analysis: ", x, " vs Control, ", cell_type)

        immune_sub <- subset(immune_obj, subset = Cell_Type_2 == cell_type)

        temp <- FindMarkers(immune_sub,
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
names(DEG_list_conditions) <- levels(immune_obj$Condition)[-1]
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

rm(temp, immune_sub, cell_type, x)

saveRDS(DEG_list_conditions, file = "DEG_list_conditions.rds")
DEG_list_conditions <- readRDS("DEG_list_conditions.rds")

# volcano plot
lapply(names(DEG_list_conditions), function(x) {
    plot_data <- do.call(rbind, DEG_list_conditions[[x]]) %>% mutate(Comparison = paste0(x, "_vs_CTRL"))
    plot_data$Cell_Type_2 <- factor(plot_data$Cell_Type_2, levels = names(mycolor$Cell_Type_2)[1:14])

    x_max <- max(abs(plot_data$avg_log2FC), na.rm = TRUE)

    label_data <- table(plot_data$Cell_Type_2, plot_data$Comparison, plot_data$DEG) %>%
        as.data.frame()
    colnames(label_data) <- c("Cell_Type_2", "Comparison", "DEG", "Count")
    label_data <- label_data %>%
        filter(DEG != "Stable") %>%
        mutate("x" = ifelse(DEG == "Up", x_max * 0.4, -x_max * 0.4), label = paste0(DEG, ": ", Count))

    p <- ggplot(plot_data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = Cell_Type_2), alpha = 0.5, size = 0.5) +
        geom_text(
            data = label_data,
            aes(x = x, y = 250, label = label), inherit.aes = FALSE, size = 3
        ) +
        scale_color_manual(values = mycolor$Cell_Type_2) +
        facet_wrap(~Cell_Type_2, ncol = 5) +
        theme_classic() +
        labs(title = paste0("Volcano plot of DEGs in each cell type, ", x, " vs Control")) +
        xlab("Average log2 fold change") +
        ylab("-log10 adjusted p value") +
        xlim(-x_max * 0.7, x_max * 0.7) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "#ec8a8a", linewidth = 0.5) +
        theme(legend.position = "none")
    pdf(paste0("Volcano.", x, "_vs_CTRL.pdf"), 9, 7)
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

## 1.3. GSEA on specific cell types ----
dir.create("GSEA_Capsaicin")

# run GSEA
table(immune_obj$Cell_Type_2, immune_obj$Condition)
cell_type <- c("T_CD8")

gseGO_Capsaicin <- sapply(c("T_CD4", "T_CD8", "Treg", "Neutrophil"), function(cell_type) {
    message(Sys.time(), ": Processing: ", cell_type)
    immune_sub <- subset(immune_obj, subset = Cell_Type_2 == cell_type)

    message(Sys.time(), ": Identify DEGs: ", cell_type)
    temp <- FindMarkers(immune_sub,
        ident.1 = "Capsaicin", ident.2 = "CTRL",
        group.by = "Condition",
        min.pct = 0, logfc.threshold = 0,
        max.cells.per.ident = 2000
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
lapply(gseGO_Capsaicin, function(x) nrow(x@result))

lapply(names(gseGO_Capsaicin) %>% setdiff("T_CD8"), function(x) {
    if (is.null(gseGO_Capsaicin[[x]])) {
        return(NULL)
    }
    write.csv(gseGO_Capsaicin[[x]]@result,
        file = paste0("GSEA_Capsaicin/GSEA_Capsaicin_in_", x, ".csv"),
        row.names = TRUE
    )
    pdf(paste0("GSEA_Capsaicin/Dot.GSEA_Capsaicin_in_", x, ".pdf"), 7, 5)
    p <- dotplot(gseGO_Capsaicin[[x]], showCategory = 20, label_format = 60) +
        theme_classic() +
        labs(title = paste("GSEA - Capsaicin vs CTRL in ", x, sep = ""))
    plot(p)
    dev.off()
})

# gseaplot2 - CD4 T
gseGO_Capsaicin$T_CD4@result %>%
    filter(grepl("T cell", Description)) %>%
    arrange(-NES) %>%
    dplyr::select(Description, NES, pvalue, p.adjust)
term_selected <- c("GO:0050798", "GO:2000406", "GO:0045582", "GO:0002294", "GO:0050870")

p_list <- lapply(term_selected, function(term) {
    gseaplot2(gseGO_Capsaicin$T_CD4,
        geneSetID = term,
        color = c("#e1423b"),
        title = paste(
            gseGO_Capsaicin$T_CD4@result[term, ]$Description, "\n",
            "NES =", round(gseGO_Capsaicin$T_CD4@result[term, ]$NES, 2),
            ", p.adj =", signif(gseGO_Capsaicin$T_CD4@result[term, ]$p.adjust, 2)
        ),
    )
})
length(p_list)

pdf("GSEA_Capsaicin/GSEA_Capsaicin_in_T_CD4_gseaplot.pdf", 5, 4)
for (i in seq_along(p_list)) {
    print(p_list[[i]])
}
dev.off()

# gseaplot2 - Neutrophil
gseGO_Capsaicin$Neutrophil@result %>%
    filter(grepl("T cell", Description)) %>%
    arrange(NES) %>%
    dplyr::select(Description, NES, pvalue, p.adjust)
term_selected <- c("GO:0050777", "GO:0002820", "GO:0002838")

p_list <- lapply(term_selected, function(term) {
    gseaplot2(gseGO_Capsaicin$Neutrophil,
        geneSetID = term,
        color = c("#458ec3"),
        title = paste(
            gseGO_Capsaicin$Neutrophil@result[term, ]$Description, "\n",
            "NES =", round(gseGO_Capsaicin$Neutrophil@result[term, ]$NES, 2),
            ", p.adj =", signif(gseGO_Capsaicin$Neutrophil@result[term, ]$p.adjust, 2)
        ),
    )
})
length(p_list)

pdf("GSEA_Capsaicin/GSEA_Capsaicin_in_Neutrophil_gseaplot.pdf", 5, 4)
for (i in seq_along(p_list)) {
    print(p_list[[i]])
}
dev.off()
rm(p_list, term_selected, i)

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
DEG_list_gene_cell_type$B$down

# compare DEGs
lapply(names(DEG_list_gene_cell_type), function(x) {
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

temp <- DEG_list_gene_cell_type$Treg$up
temp$Capsaicin %>%
    setdiff(temp$anti_PD1) %>%
    setdiff(temp$Rime) %>%
    sort()

# compare GO
temp <- compareCluster(DEG_list_gene_cell_type$T_CD8$up,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
) %>% simplify()
dotplot(temp) + theme_classic() + labs(title = "GO comparison of up-regulated genes in T_CD8")
rm(temp)

names(DEG_list_gene_cell_type)
GO_list_compare <- mclapply(c("T_CD4", "T_CD8", "Treg", "Neutrophil"), function(x) {
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
}, mc.cores = 4) %>%
    setNames(c("T_CD4", "T_CD8", "Treg", "Neutrophil"))

lapply(names(GO_list_compare), function(x) {
    if (!is.null(GO_list_compare[[x]]$up)) {
        pdf(paste0("Compare_Conditions/GO_Compare_Up_", x, ".pdf"), 7, 5)
        p <- dotplot(GO_list_compare[[x]]$up, showCategory = 10, label_format = 80) +
            theme_classic() +
            labs(title = paste("GO of up-regulated genes in ", x, sep = ""))
        plot(p)
        dev.off()
    }
    if (!is.null(GO_list_compare[[x]]$down)) {
        pdf(paste0("Compare_Conditions/GO_Compare_Down_", x, ".pdf"), 7, 5)
        p <- dotplot(GO_list_compare[[x]]$down, showCategory = 10, label_format = 80) +
            theme_classic() +
            labs(title = paste("GO of down-regulated genes in ", x, sep = ""))
        plot(p)
        dev.off()
    }
})

# 2. T cell states ----
## 2.1. Project T cell states ----
library(ProjecTILs)

table(immune_obj$Cell_Type_2)
T_obj <- subset(immune_obj,
    subset = Cell_Type_2 %in% c("T_CD4", "T_CD4_Cycling", "Treg", "T_CD8", "T_CD8_Cycling")
)
saveRDS(T_obj, file = "T_obj.rds")
rm(immune_obj)

ref <- load.reference.map()

T_obj.projected <- Run.ProjecTILs(T_obj,
    ref = ref, filter.cells = FALSE, k = 30, split.by = "Sample",
    ncores = 1
)
gc()
table(T_obj.projected$functional.cluster, T_obj.projected$Cell_Type_2)

levels(ref$functional.cluster)
mycolor$T_subtype <- c(
    "#edbe2a", "#A58AFF", "#53B400", "#F8766D",
    "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd"
) %>% setNames(levels(ref$functional.cluster))

p_list <- list()
sample_names <- (levels(T_obj$Condition))
for (sample_i in seq_along(sample_names)) {
    sample <- sample_names[sample_i]
    p_list[[sample_i]] <- plot.projection(ref,
        T_obj.projected[, T_obj.projected$Condition == sample],
        linesize = 0.3, pointsize = 0.6,
        raster = TRUE
    ) + ggtitle(sample) + theme(aspect.ratio = 1)
}
p_list$ref <- DimPlot(ref,
    group.by = "functional.cluster", label = TRUE,
    cols = mycolor$T_subtype, raster = TRUE
) +
    ggtitle("Reference") + theme(aspect.ratio = 1)

pdf("Projection_of_T_cells.pdf", 12, 10)
patchwork::wrap_plots(p_list, ncol = 2)
dev.off()
rm(p_list, sample, sample_i, sample_names)

table(T_obj.projected$functional.cluster)

p_list <- plot.states.radar(ref = ref, T_obj.projected, min.cells = 20, return.as.list = TRUE)
length(p_list)

pdf("Radar.T.Projection.pdf", 9, 8)
patchwork::wrap_plots(p_list[1:9], ncol = 3)
dev.off()
rm(p_list, Hs2Mm.convert.table, ref)

## 2.2. Composition of T cell states ----
T_obj.projected$Condition <- factor(T_obj.projected$Condition, levels = names(mycolor$Condition))
T_obj.projected$functional.cluster <- factor(T_obj.projected$functional.cluster,
    levels = names(mycolor$T_subtype)
)

identical(colnames(T_obj.projected), colnames(T_obj))
T_obj$T_state <- T_obj.projected$functional.cluster

plot_data <- as.data.frame(table(T_obj.projected$functional.cluster, T_obj.projected$Condition))
colnames(plot_data) <- c("functional.cluster", "Condition", "Count")

pdf("Barplot.T_subtype.Condition.pdf", 8, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = functional.cluster, y = Count, fill = Condition)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$Condition) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ggplot(plot_data, aes(x = Condition, y = Count, fill = functional.cluster)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(values = mycolor$T_subtype) +
        theme_classic() +
        labs(title = "Cell Type composition across Conditions") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)),
    ncol = 2
)
dev.off()

# fraction compare
plot_data <- as.data.frame(table(T_obj.projected$functional.cluster, T_obj.projected$Sample))
colnames(plot_data) <- c("T_subtype", "Sample", "Count")
plot_data$Condition <- gsub("_[0-9]$", "", plot_data$Sample)
table(plot_data$Condition)
table(plot_data$Sample)

temp <- table(T_obj.projected$Sample)
plot_data$Fraction <- plot_data$Count / temp[plot_data$Sample] * 100
plot_data$Condition <- factor(plot_data$Condition, levels = names(mycolor$Condition))

pdf("Boxplot.Fraction.T_subtype.pdf", 8, 4.7)
ggplot(plot_data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = T_subtype), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 2.5, ref.group = "CTRL",
    ) +
    facet_wrap(~T_subtype, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction (%)") +
    scale_fill_manual(values = mycolor$T_subtype) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot_data, temp, T_obj.projected)

# roe
roe <- Startrac::calTissueDist(
    T_obj@meta.data,
    byPatient = FALSE,
    colname.cluster = "T_state",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
roe <- roe[setdiff(rownames(roe), "Doublet"), ]
pdf("Heatmap.Roe.T_state.Condition.pdf", 4, 4.5)
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

## 2.3. T state marker expression ----
gene_selected <- c(
    "Foxp3", "Il2ra", "Il10", "Tgfb1", # Treg
    "Gzmb", "Prf1", "Ifng", "Tnf", "Fasl", # Cytotoxicity
    "Pdcd1", "Lag3", "Havcr2", "Tigit", "Ctla4", "Tox" # exhaustion
)

Idents(T_obj) <- T_obj$T_state
pdf("Vln.T_state.marker.pdf", 12, 9)
VlnPlot(T_obj,
    idents = names(mycolor$T_subtype),
    features = gene_selected,
    group.by = "Condition",
    pt.size = 0, ncol = 4, cols = mycolor$Condition
)
dev.off()

pdf("Vln.T_state.marker.split.pdf", 16, 12)
VlnPlot(T_obj,
    idents = names(mycolor$T_subtype),
    features = gene_selected,
    split.by = "Condition",
    pt.size = 0, ncol = 4, cols = mycolor$Condition
)
dev.off()

plot_data <- DotPlot(T_obj %>% subset(!is.na(T_state)),
    features = gene_selected,
    group.by = "Condition"
)$data
plot_data$group <- rep(rep(c("Treg(Supression)", "Cytotoxicity", "Exhaustion"), times = c(4, 5, 6)),
    times = 5
)

pdf("Dot.T_state.markers.Condition.pdf", 7, 5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
        scale_fill_viridis_c() +
        facet_grid(cols = vars(group), scales = "free_x", space = "free_x") +
        theme_classic() +
        labs(x = "Gene", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle("Expression of T cell state markers"),
    ggplot(plot_data, aes(x = features.plot, y = id)) +
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

table(T_obj$Cell_Type_2)
pdf("Dot.T_PDCD1.Condition.pdf", 4, 5)
DotPlot(T_obj %>% subset(!is.na(T_state)),
    features = "Pdcd1",
    group.by = "Condition"
) + scale_color_viridis_c() +
    theme_classic() + labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Pdcd1 in T cells")
dev.off()

pdf("Dot.T_CD8_PDCD1.Condition.pdf", 4, 5)
DotPlot(T_obj %>% subset(Cell_Type_2 == "T_CD8" & !is.na(T_state)),
    features = "Pdcd1",
    group.by = "Condition"
) + scale_color_viridis_c() +
    theme_classic() + labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Pdcd1 in CD8 T cells")
dev.off()

pdf("Dot.T_CD4_PDCD1.Condition.pdf", 4, 5)
DotPlot(T_obj %>% subset(Cell_Type_2 == "T_CD4" & !is.na(T_state)),
    features = "Pdcd1",
    group.by = "Condition"
) + scale_color_viridis_c() +
    theme_classic() + labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Pdcd1 in CD4 T cells")
dev.off()

pdf("Vln.T_CD8_PDCD1.Condition.pdf", 5, 3)
VlnPlot(T_obj %>% subset(Cell_Type_2 == "T_CD8" & !is.na(T_state)),
    features = "Pdcd1", col = mycolor$Condition,
    group.by = "Condition",  pt.size = 0,
) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Pdcd1 in CD8 T cells")
dev.off()

table(T_obj$T_state)
p_list <- lapply(levels(T_obj$T_state), function(state) {
    plot_data <- DotPlot(T_obj %>% subset(T_state == state),
        features = gene_selected,
        group.by = "Condition"
    )$data
    plot_data$group <- rep(rep(c("Treg(Supression)", "Cytotoxicity", "Exhaustion"), times = c(4, 5, 6)),
        times = nrow(plot_data) / 15
    )

    ggplot(plot_data, aes(x = features.plot, y = id)) +
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
pdf("Dot.T_state.markers.Condition.Split.pdf", 12, 12)
patchwork::wrap_plots(p_list, ncol = 2)
dev.off()
rm(plot_data, state, p_list)

## 2.4. T state program ----
library(SignatuR)
data(SignatuR)

programs <- GetSignature(SignatuR$Mm$Programs)
names(programs)
rm(SignatuR)

T_obj <- AddModuleScore(T_obj, features = programs)
colnames(T_obj@meta.data)[17:23] <- names(programs)

T_obj@meta.data %>%
    group_by(Condition) %>%
    summarise(across(all_of(names(programs)[5:7]), list(mean = mean)))

plot_data <- DotPlot(T_obj %>% subset(!is.na(T_state)),
    features = names(programs) %>% setdiff(c("HeatShock", "IFN")),
    group.by = "Condition"
)$data

pdf("Dot.T_state.Signature.Condition.pdf", 9, 2.5)
patchwork::wrap_plots(
    ggplot(plot_data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
        scale_fill_viridis_c() +
        theme_classic() +
        labs(x = "Signature", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle("T cell state signatures"),
    ggplot(plot_data, aes(x = features.plot, y = id)) +
        geom_point(aes(size = pct.exp, fill = avg.exp.scaled), shape = 21) +
        scale_fill_viridis_c() +
        theme_classic() +
        labs(x = "Signature", y = "Condition") +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        ggtitle("Scaled T cell state signatures"),
    ncol = 2
)
dev.off()

table(T_obj$T_state)
p_list <- lapply(levels(T_obj$T_state), function(state) {
    plot_data <- DotPlot(T_obj %>% subset(T_state == state),
        features = names(programs) %>% setdiff(c("HeatShock", "IFN")),
        group.by = "Condition"
    )$data

    ggplot(plot_data, aes(x = features.plot, y = id)) +
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
pdf("Dot.T_state.Signatures.Condition.Split.pdf", 8, 8)
patchwork::wrap_plots(p_list, ncol = 3)
dev.off()
rm(plot_data, state, p_list)

# save data ----
gc()
save.image("02.Immune.RData")
