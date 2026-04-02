setwd("/Users/zhenyuliu/Research/Project/BLCA_Nerve/results/04.Conmmunication")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(CellChat)
library(parallel)

# load("04.Communication.RData")

# 1. Integrate cell types ----
blca_obj <- readRDS("../01.Cell_Type/rds/blca_obj.rds")
mycolor <- readRDS("../01.Cell_Type/rds/mycolor.rds")

## 1.1. merge cell types from sub-objects ----
blca_obj$Cell_Type_2 <- "none"

immune_obj <- readRDS("../01.Cell_Type/rds/immune_obj.rds")
blca_obj@meta.data[immune_obj$Cell, "Cell_Type_2"] <- immune_obj$Cell_Type_2 %>% as.character()
rm(immune_obj)

cancer_obj <- readRDS("../01.Cell_Type/rds/cancer_obj.rds")
blca_obj@meta.data[cancer_obj$Cell, "Cell_Type_2"] <- cancer_obj$Cell_Type_1 %>% as.character()
rm(cancer_obj)

stromal_obj <- readRDS("../01.Cell_Type/rds/stromal_obj.rds")
blca_obj@meta.data[stromal_obj$Cell, "Cell_Type_2"] <- stromal_obj$Cell_Type_2 %>% as.character()
rm(stromal_obj)

blca_obj@meta.data %>%
    filter(Cell_Type_2 == "none") %>%
    pull(Cell_Type_1) %>%
    table()
blca_obj@meta.data[blca_obj$Cell_Type_1 == "Epithelial", "Cell_Type_2"] <- "Epithelial"
blca_obj@meta.data[blca_obj$Cell_Type_1 == "Doublet", "Cell_Type_2"] <- "Doublet"

table(blca_obj$Cell_Type_2, blca_obj$Cell_Type_1)

mycolor$Cell_Type_2["Epithelial"] <- mycolor$Cell_Type_1["Epithelial"]
mycolor$Cell_Type_2["Cancer"] <- mycolor$Cell_Type_1["Cancer"]

table(blca_obj$Cell_Type_2 %in% names(mycolor$Cell_Type_2))

temp <- intersect(names(mycolor$Cell_Type_2), unique(blca_obj$Cell_Type_2))
temp <- setdiff(temp, "Doublet")
temp <- c(temp, "Doublet")

blca_obj$Cell_Type_2 <- factor(blca_obj$Cell_Type_2, levels = temp)
rm(temp)

pdf("UMAP.Cell_Type_2.pdf", 7, 3.6)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_2", label = TRUE, label.size = 2.5) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

pdf("UMAP.Cell_Type_2.NoLabel.pdf", 7, 3.6)
DimPlot(blca_obj, reduction = "umap", group.by = "Cell_Type_2", label = FALSE) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

pdf("UMAP.Cell_Type_2.split.Condition.pdf", 10, 5.5)
DimPlot(blca_obj,
    reduction = "umap", group.by = "Cell_Type_2", label = FALSE,
    split.by = "Condition", ncol = 3
) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

p <- DimPlot(blca_obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition", ncol = 3,
    label = FALSE, raster = TRUE, pt.size = NULL,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()

pdf("UMAP.Cell_Type_2.split.density.pdf", 10, 5.5)
plot(p + stat_density_2d(aes_string(x = "umap_1", y = "umap_2", fill = "after_stat(level)"),
    linewidth = 0, geom = "density_2d_filled",
    n = 150, h = c(1.2, 1.2)
) +
    scale_fill_viridis_c(option = "D") +
    guides(fill = guide_colorbar(title = "Density")) +
    theme(legend.position = "right"))
dev.off()
rm(p)

## 1.2. cell types cross conditions ----
plot.data <- as.data.frame(table(blca_obj$Cell_Type_2, blca_obj$Condition))
colnames(plot.data) <- c("Cell_Type_2", "Condition", "Count")

pdf("Barplot.Cell_Type_2.Condition.pdf", 11, 5)
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
plot.data <- as.data.frame(table(blca_obj$Cell_Type_2, blca_obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_[0-9]$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(blca_obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = names(mycolor$Condition))

pdf("Boxplot.Fraction.Cell_Type_2.pdf", 12, 9)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 0.98, size = 2.5, ref.group = "CTRL") +
    facet_wrap(~Cell_Type_2, ncol = 6, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction of all cells (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot.data, temp)

# Ro/e
colnames(blca_obj@meta.data)
Roe <- Startrac::calTissueDist(
    blca_obj@meta.data,
    byPatient = FALSE,
    colname.cluster = "Cell_Type_2",
    colname.patient = "Sample",
    colname.tissue = "Condition",
    method = "chisq",
    min.rowSum = 0
)
anno.row <- data.frame(
    Cell_Type_2 = rownames(Roe),
    "Compartment" = rep(c("Immune", "Stromal", "Epithelial", "Doublet"), times = c(14, 7, 2, 1))
)
cell_type_selected <- table(blca_obj$Cell_Type_2) %>%
    .[. >= 100] %>%
    names()

pdf("Heatmap.Roe.Cell_Type_2.Condition.pdf", 4.6, 5.5)
pheatmap::pheatmap(
    Roe[cell_type_selected, ],
    scale = "none",
    color = mycolor$Roe,
    annotation_row = anno.row %>%
        filter(Cell_Type_2 %in% cell_type_selected) %>%
        tibble::column_to_rownames("Cell_Type_2"),
    annotation_colors = list(
        Compartment = c("Immune" = "#FC8D62", "Stromal" = "#66C2A5", "Epithelial" = "#8DA0CB", "Doublet" = "gray")
    ),
    display_numbers = TRUE,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    cutree_rows = 3,
    fontsize = 8,
    fontsize_number = 6.2,
    main = "Ro/e of Cell_Type_2"
)
dev.off()
rm(Roe, anno.row, cell_type_selected)

## 1.3. expression of selected genes ----
dir.create("Expression_Selected")
blca_obj <- readRDS("blca_obj.rds")
Idents(blca_obj) <- blca_obj$Cell_Type_2

# neurotrophic factor genes
gene.selected <- c(
    # Neurotrophin Family
    "Ngf", "Ntrk1", "Bdnf", "Ntrk2", "Ntf3", "Ntrk3", "Ngfr",
    # GDNF Family
    "Gdnf", "Gfra1", "Nrtn", "Gfra2", "Artn", "Gfra3", "Pspn", "Gfra4", "Ret"
)

pdf("Expression_Selected/Violin.Expression.Neuro_factor.pdf", 20, 12)
VlnPlot(blca_obj,
    features = gene.selected, cols = mycolor$Cell_Type_2,
    group.by = "Cell_Type_2",
    pt.size = 0, ncol = 4
)
dev.off()

pdf("Expression_Selected/Violin.Expression.Neuro_factor.split.pdf", 24, 24)
VlnPlot(blca_obj,
    features = gene.selected, split.by = "Condition",
    cols = mycolor$Condition,
    pt.size = 0, ncol = 3
)
dev.off()

pdf("Expression_Selected/Dot.Expression.Neuro_factor.pdf", 7, 5)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = FALSE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Neurotrophic Factors")
dev.off()

pdf("Expression_Selected/Dot.Expression.Neuro_factor.scaled.pdf", 7, 5)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = TRUE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Scaled Expression of Neurotrophic Factors")
dev.off()

blca_obj$Cell_Type_2_Condition <- paste0(blca_obj$Cell_Type_2, "-", blca_obj$Condition)

plot.data <- DotPlot(blca_obj,
    features = gene.selected, scale = FALSE,
    group.by = "Cell_Type_2_Condition"
)$data

plot.data$Condition <- gsub("^.+?-", "", plot.data$id)
plot.data$Condition <- factor(plot.data$Condition, levels = names(mycolor$Condition))
plot.data$Cell_Type_2 <- gsub("-.+?*$", "", plot.data$id)
plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(blca_obj$Cell_Type_2))
# plot.data <- plot.data[!is.na(plot.data$Cell_Type_2), ]
colnames(plot.data)

pdf("Expression_Selected/Dot.Expression.Neuro_factor.split.pdf", 14, 9)
ggplot(plot.data, aes(x = features.plot, y = Condition)) +
    geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
    scale_fill_viridis_c() +
    facet_wrap(~Cell_Type_2, scales = "free_x", ncol = 5) +
    theme_classic() +
    labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Neuronal Factors")
dev.off()
rm(plot.data, gene.selected)

# growth factors
gene.selected <- c(
    "Egf", "Hbegf", "Egfr", # Epidermal Growth Factor (EGF) family
    "Ngf", "Bdnf", "Ntrk1", # Nerve Growth Factor (NGF) / Neurotrophin family
    "Igf1", "Igf2", "Igf1r", "Igf2r", # Insulin-like Growth Factor (IGF) family
    "Fgf1", "Fgf2", "Fgf7", "Fgf13", "Fgfr1", "Fgfr2", # Fibroblast Growth Factor (FGF) family
    "Tgfa", "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2" # Transforming Growth Factor (TGF) family
)

grep("gf", rownames(blca_obj), value = TRUE) %>% sort()

pdf("Expression_Selected/Violin.Expression.Growth_factor.pdf", 20, 20)
VlnPlot(blca_obj,
    features = gene.selected, cols = mycolor$Cell_Type_2,
    group.by = "Cell_Type_2",
    pt.size = 0, ncol = 4
)
dev.off()

pdf("Expression_Selected/Violin.Expression.Growth_factor.split.pdf", 28, 35)
VlnPlot(blca_obj,
    features = gene.selected, split.by = "Condition",
    cols = mycolor$Condition,
    pt.size = 0, ncol = 2
)
dev.off()

pdf("Expression_Selected/Dot.Expression.Growth_factor.pdf", 8, 6)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = FALSE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Growth Factors")
dev.off()

pdf("Expression_Selected/Dot.Expression.Growth_factor.scaled.pdf", 8, 6)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = TRUE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Scaled Expression of Growth Factors")
dev.off()

plot.data <- DotPlot(blca_obj,
    features = gene.selected, scale = FALSE,
    group.by = "Cell_Type_2_Condition"
)$data

plot.data$Condition <- gsub("^.+?-", "", plot.data$id)
plot.data$Condition <- factor(plot.data$Condition, levels = names(mycolor$Condition))
plot.data$Cell_Type_2 <- gsub("-.+?*$", "", plot.data$id)
plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(blca_obj$Cell_Type_2))
plot.data <- plot.data[!is.na(plot.data$Cell_Type_2), ]
colnames(plot.data)
plot.data[plot.data$avg.exp > 6, ]$avg.exp <- 6

pdf("Expression_Selected/Dot.Expression.Growth_factor.split.pdf", 15, 15)
ggplot(plot.data, aes(x = features.plot, y = Condition)) +
    geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
    scale_fill_viridis_c() +
    facet_wrap(~Cell_Type_2, scales = "free_x", ncol = 3) +
    theme_classic() +
    labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Growth Factors")
dev.off()
rm(plot.data, gene.selected)

# cytokines and chemokines
gene.selected <- c(
    "Egf", "Igf1", "Igfbp5", "Igfbp6", "Igfbp2",
    "Ccl11", "Cxcl13", "Ccl19",
    "Cxcl5",
    "Cd40", "Selp", "Tnfrsf18",
    "Il20", "Cd70",
    "Vegfa", "Cxcl2", "Igfbp3", "Ccl24", "Mmp9"
)

pdf("Expression_Selected/Violin.Expression.Cytokine.pdf", 24, 16)
VlnPlot(blca_obj,
    features = gene.selected, cols = mycolor$Cell_Type_2,
    group.by = "Cell_Type_2",
    pt.size = 0, ncol = 4
)
dev.off()

pdf("Expression_Selected/Violin.Expression.Cytokine.split.pdf", 28, 36)
VlnPlot(blca_obj,
    features = gene.selected, split.by = "Condition",
    cols = mycolor$Condition,
    pt.size = 0, ncol = 2
)
dev.off()

pdf("Expression_Selected/Dot.Expression.Cytokine.pdf", 8, 6)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = FALSE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Cytokines")
dev.off()

pdf("Expression_Selected/Dot.Expression.Cytokine.scaled.pdf", 8, 6)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = TRUE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Scaled Expression of Cytokines")
dev.off()

plot.data <- DotPlot(blca_obj,
    features = gene.selected, scale = FALSE,
    group.by = "Cell_Type_2_Condition"
)$data

plot.data$Condition <- gsub("^.+?-", "", plot.data$id)
plot.data$Condition <- factor(plot.data$Condition, levels = names(mycolor$Condition))
plot.data$Cell_Type_2 <- gsub("-.+?*$", "", plot.data$id)
plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(blca_obj$Cell_Type_2))
plot.data <- plot.data[!is.na(plot.data$Cell_Type_2), ]
colnames(plot.data)
plot.data[plot.data$avg.exp > 6, ]$avg.exp <- 6

pdf("Expression_Selected/Dot.Expression.Cytokine.split.pdf", 15, 15)
ggplot(plot.data, aes(x = features.plot, y = Condition)) +
    geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
    scale_fill_viridis_c() +
    facet_wrap(~Cell_Type_2, scales = "free_x", ncol = 3) +
    theme_classic() +
    labs(x = "Gene", y = "Condition") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Cytokines")
dev.off()
rm(plot.data, gene.selected)

# additianal genes 
gene.selected <- c(
    "Calca", "Calcb", "Crcp",
    "Tac1", "Tacr1", "Tacr2",
    "Vip", "Vipr1", "Vipr2",
    "Adcyap1", "Adcyap1r1",
    "Gria1", "Gria2", "Gria3", "Gria4",
    "Grin1", "Grin2a", "Grin2b", "Grin2c", "Grin2d",
    "Grik1", "Grik2", "Grik3", "Grik4", "Grik5",
    "Grm1", "Grm2", "Grm3", "Grm4", "Grm5", "Grm6", "Grm7", "Grm8",
    "P2rx3", "P2rx7",
    "P2ry1", "P2ry2", "P2ry4", "P2ry6", "P2ry12", "P2ry13", "P2ry14"
)

pdf("Expression_Selected/Dot.Expression.Selected2.pdf", 8, 9)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = FALSE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Cytokines")
dev.off()

pdf("Expression_Selected/Dot.Expression.Selected2.scaled.pdf", 8, 9)
DotPlot(blca_obj,
    features = gene.selected %>% rev(), scale = TRUE,
    group.by = "Cell_Type_2",
) +
    coord_flip() +
    theme_classic() +
    scale_color_viridis_c() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Scaled Expression of Cytokines")
dev.off()

# 2. Joint CellChat analysis ----
## 2.1. prepare CellChat object ----
table(blca_obj$Cell_Type_2) %>% sort()
cell_types <- levels(blca_obj$Cell_Type_2) %>%
    grep("Cycling", ., invert = TRUE, value = TRUE) %>%
    setdiff(c("Adipocyte", "Schwann", "Smooth_Muscle", "Doublet"))
table(blca_obj$Cell_Type_2)[cell_types] %>% sort()

blca_sub <- subset(blca_obj, downsample = 2000)
blca_sub <- subset(blca_sub, Cell_Type_2 %in% cell_types)
table(blca_sub$Cell_Type_2) %>% sort()

Idents(blca_sub) <- blca_sub$Cell_Type_2

cellchat <- createCellChat(object = blca_sub, group.by = "ident", assay = "RNA")

saveRDS(blca_obj, file = "blca_obj.rds")
rm(blca_sub, blca_obj)

# primary cellchat pipeline
showDatabaseCategory(CellChatDB.mouse)
cellchat@DB <- CellChatDB.mouse

# filter genes & LR pairs
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # 3220 variable L-R pairs

## 2.2. infer communication network ----
# communication probability
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)

# pathway-level communication
cellchat <- computeCommunProbPathway(cellchat)

# cc network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

pdf("Circle.Network.All.Count.pdf", 6, 7)
netVisual_circle(cellchat@net$count,
    color.use = mycolor$Cell_Type_2[cell_types],
    vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE,
    title.name = "Number of interactions"
)
dev.off()

pdf("Circle.Network.All.Weight.pdf", 6, 7)
netVisual_circle(cellchat@net$weight,
    color.use = mycolor$Cell_Type_2[cell_types],
    vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE,
    title.name = "Interaction weights/strength"
)
dev.off()

# split by cell type
mat <- cellchat@net$weight
pdf("Circle.Network.Split.Weigh.Out.pdf", 20, 16)
par(mfrow = c(4, 5), xpd = TRUE)
for (i in seq_along(rownames(mat))) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2,
        color.use = mycolor$Cell_Type_2[cell_types],
        vertex.weight = groupSize,
        weight.scale = TRUE,
        label.edge = TRUE,
        edge.weight.max = max(mat),
        title.name = rownames(mat)[i]
    )
}
dev.off()

pdf("Circle.Network.Split.Weigh.In.pdf", 20, 16)
par(mfrow = c(4, 5), xpd = TRUE)
for (i in seq_along(rownames(mat))) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[, i] <- mat[, i]
    netVisual_circle(mat2,
        color.use = mycolor$Cell_Type_2[cell_types],
        vertex.weight = groupSize,
        weight.scale = TRUE,
        label.edge = TRUE,
        edge.weight.max = max(mat),
        title.name = rownames(mat)[i]
    )
}
dev.off()
rm(i, mat, mat2)

## 2.3. analysis of whole network ----
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf("Scatter.In_out.pdf", 5, 4)
netAnalysis_signalingRole_scatter(cellchat,
    color.use = mycolor$Cell_Type_2[cell_types],
    dot.alpha = 0.4
)
dev.off()

pdf("Heatmap.In_out.pdf", 10, 15)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat,
    color.use = mycolor$Cell_Type_2[cell_types], pattern = "outgoing", height = 30
)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,
    color.use = mycolor$Cell_Type_2[cell_types], pattern = "incoming", height = 30
)
(ht1 + ht2)
dev.off()
rm(ht1, ht2)

pdf("Heatmap.Interaction.pdf", 6, 5.5)
netVisual_heatmap(cellchat,
    measure = "weight",
    color.use = mycolor$Cell_Type_2[cell_types],
    color.heatmap = "Reds"
)
dev.off()

# global communication patterns
dir.create("Global_Communication_Patterns")
library(NMF)
library(ggalluvial)

cellchat@netP$pathways %>% grep("GF", ., value = TRUE)
cellchat@netP$pathways %>% grep("IL", ., value = TRUE)
cellchat@netP$pathways %>% grep("IFN", ., value = TRUE)

pathways.selected <- c(
    "TGFb", "VEGF", "EGF", "IGF", "PDGF", "FGF",
    "CXCL", "CCL", "IL1", "IL2", "IL6", "FASLG", "MIF",
    "MHC-I", "MHC-II", "PD-L1", "PDL2",
    "COLLAGEN", "LAMININ"
)

## Outgoing
pdf("Global_Communication_Patterns/Elbow.Outgoing.pdf", 6, 2.5)
selectK(cellchat, pattern = "outgoing")
dev.off()

nPatterns <- 6

pdf("Global_Communication_Patterns/Heatmap.Pattern.Out.pdf", 8, 4)
cellchat <- identifyCommunicationPatterns(cellchat,
    color.use = mycolor$Cell_Type_2[cell_types], pattern = "outgoing", k = nPatterns
)
dev.off()

pdf("Global_Communication_Patterns/River.Pattern.Out.pdf", 8, 4)
netAnalysis_river(cellchat, color.use = mycolor$Cell_Type_2[cell_types], pattern = "outgoing")
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.Out.pdf", 4.5, 15)
netAnalysis_dot(cellchat, color.use = mycolor$Cell_Type_2[cell_types], pattern = "outgoing") + coord_flip()
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.Out.Selected.pdf", 5.5, 4)
netAnalysis_dot(cellchat,
    pattern = "outgoing",
    color.use = mycolor$Cell_Type_2[cell_types],
    pathway.show = pathways.selected %>% rev()
) + coord_flip()
dev.off()

## incoming
pdf("Global_Communication_Patterns/Elbow.Incoming.pdf", 6, 2.5)
selectK(cellchat, pattern = "incoming")
dev.off()

nPatterns <- 6

pdf("Global_Communication_Patterns/Heatmap.Pattern.In.pdf", 8, 4)
cellchat <- identifyCommunicationPatterns(cellchat,
    color.use = mycolor$Cell_Type_2[cell_types], pattern = "incoming", k = nPatterns
)
dev.off()

pdf("Global_Communication_Patterns/River.Pattern.In.pdf", 8, 4)
netAnalysis_river(cellchat, color.use = mycolor$Cell_Type_2[cell_types], pattern = "incoming")
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.In.pdf", 4.5, 15)
netAnalysis_dot(cellchat, color.use = mycolor$Cell_Type_2[cell_types], pattern = "incoming") + coord_flip()
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.In.Selected.pdf", 5.5, 4)
netAnalysis_dot(cellchat,
    pattern = "incoming",
    color.use = mycolor$Cell_Type_2[cell_types],
    pathway.show = pathways.selected %>% rev()
) + coord_flip()
dev.off()
rm(nPatterns)

## 2.4. pathway level analysis ----
dir.create("Pathway_Level")

table(pathways.selected %in% cellchat@netP$pathways)

lapply(pathways.selected %>% setdiff(c("FGF")), function(pathway) {
    pdf(paste0("Pathway_Level/Circle.Network.Pathway.", pathway, ".pdf"), 6.5, 7)
    netVisual_aggregate(cellchat,
        color.use = mycolor$Cell_Type_2[cell_types],
        signaling = pathway, layout = "circle",
    )
    dev.off()

    pdf(paste0("Pathway_Level/Heatmap.Pathway.", pathway, ".pdf"), 6, 5.5)
    p <- netVisual_heatmap(cellchat,
        signaling = pathway,
        color.use = mycolor$Cell_Type_2[cell_types],
        color.heatmap = "Reds"
    )
    plot(p)
    dev.off()

    pdf(paste0("Pathway_Level/Chord.Pathway.", pathway, ".pdf"), 8, 9)
    netVisual_aggregate(cellchat,
        signaling = pathway,
        color.use = mycolor$Cell_Type_2[cell_types],
        layout = "chord"
    )
    dev.off()

    pdf(paste0("Pathway_Level/Heatmap_signalingRole.", pathway, ".pdf"), 8, 3)
    netAnalysis_signalingRole_network(cellchat,
        signaling = pathway,
        color.use = mycolor$Cell_Type_2[cell_types],
        width = 8, height = 2.5, font.size = 10
    )
    dev.off()

    pdf(paste0("Pathway_Level/Bar.LR_Contribution.", pathway, ".pdf"), 5, 4)
    p <- netAnalysis_contribution(cellchat,
        signaling = pathway,
    )
    plot(p)
    dev.off()
})
rm(p)

## 2.5. L-R pair level analysis ----
dir.create("LR_Pairs_Level")
dir.create("LR_Pairs_Level/All")
dir.create("LR_Pairs_Level/Selected")

cell_types

lapply(cell_types, function(cell_type) {
    # all
    p <- netVisual_bubble(cellchat,
        sources.use = cell_type,
        remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE
    ) +
        facet_grid(
            cols = vars(source),
            rows = vars(pathway_name), scales = "free", space = "free"
        ) +
        theme(
            panel.grid.major = element_line(color = "gray90", size = 0.25),
            strip.background = element_rect(fill = "gray80"),
            strip.text = element_text(colour = "black"),
            strip.text.y = element_text(angle = 0)
        )
    height <- length(unique(p$data$interaction_name)) / 6 + 1
    pdf(paste0("LR_Pairs_Level/All/Bubble.", cell_type, "_Outgoing.pdf"), 7, height)
    plot(p)
    dev.off()

    p <- netVisual_bubble(cellchat,
        targets.use = cell_type,
        remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE
    ) +
        facet_grid(
            cols = vars(target),
            rows = vars(pathway_name), scales = "free", space = "free"
        ) +
        theme(
            panel.grid.major = element_line(color = "gray90", size = 0.25),
            strip.background = element_rect(fill = "gray80"),
            strip.text = element_text(colour = "black"),
            strip.text.y = element_text(angle = 0)
        )
    height <- length(unique(p$data$interaction_name)) / 6 + 1
    pdf(paste0("LR_Pairs_Level/All/Bubble.", cell_type, "_Incoming.pdf"), 7, height)
    plot(p)
    dev.off()

    # selected
    p <- netVisual_bubble(cellchat,
        sources.use = cell_type,
        signaling = pathways.selected %>% setdiff(c("COLLAGEN", "LAMININ", "MHC-I", "MHC-II")),
        remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE
    ) +
        facet_grid(
            cols = vars(source),
            rows = vars(pathway_name), scales = "free", space = "free"
        ) +
        theme(
            panel.grid.major = element_line(color = "gray90", size = 0.25),
            strip.background = element_rect(fill = "gray80"),
            strip.text = element_text(colour = "black"),
            strip.text.y = element_text(angle = 0)
        )
    height <- length(unique(p$data$interaction_name)) / 5.5 + 2
    pdf(paste0("LR_Pairs_Level/Selected/Bubble.", cell_type, "_Outgoing.pdf"), 6, height)
    plot(p)
    dev.off()

    p <- netVisual_bubble(cellchat,
        targets.use = cell_type,
        signaling = pathways.selected %>% setdiff(c("COLLAGEN", "LAMININ", "MHC-I", "MHC-II")),
        remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE
    ) +
        facet_grid(
            cols = vars(target),
            rows = vars(pathway_name), scales = "free", space = "free"
        ) +
        theme(
            panel.grid.major = element_line(color = "gray90", size = 0.25),
            strip.background = element_rect(fill = "gray80"),
            strip.text = element_text(colour = "black"),
            strip.text.y = element_text(angle = 0)
        )
    height <- length(unique(p$data$interaction_name)) / 5.5 + 2
    pdf(paste0("LR_Pairs_Level/Selected/Bubble.", cell_type, "_Incoming.pdf"), 6, height)
    plot(p)
    dev.off()
})
rm(p, cell_type)

# 3. Compare signal across conditions ----
blca_obj <- readRDS("blca_obj.rds")
dir.create("Compare_Conditions")
table(blca_obj$Cell_Type_2, blca_obj$Condition)

## 3.1. prepare CellChat objects ----
blca_sub <- subset(blca_obj, downsample = 20000)
rm(blca_obj)

blca_sub <- subset(blca_sub, Cell_Type_2 %in% cell_types)
table(blca_sub$Cell_Type_2) %>% sort()
table(blca_sub$Condition)

temp <- subset(blca_sub, Condition %in% c("CTRL"))
Idents(temp) <- temp$Cell_Type_2
cellchat.ctrl <- createCellChat(object = temp, group.by = "ident", assay = "RNA")
cellchat.ctrl@DB <- CellChatDB.mouse

temp <- subset(blca_sub, Condition %in% c("Capsaicin"))
Idents(temp) <- temp$Cell_Type_2
cellchat.capsaicin <- createCellChat(object = temp, group.by = "ident", assay = "RNA")
cellchat.capsaicin@DB <- CellChatDB.mouse

temp <- subset(blca_sub, Condition %in% c("LC_6OHDA"))
Idents(temp) <- temp$Cell_Type_2
cellchat.lc_6ohda <- createCellChat(object = temp, group.by = "ident", assay = "RNA")
cellchat.lc_6ohda@DB <- CellChatDB.mouse

temp <- subset(blca_sub, Condition %in% c("anti_PD1"))
Idents(temp) <- temp$Cell_Type_2
cellchat.anti_pd1 <- createCellChat(object = temp, group.by = "ident", assay = "RNA")
cellchat.anti_pd1@DB <- CellChatDB.mouse

temp <- subset(blca_sub, Condition %in% c("Rime"))
Idents(temp) <- temp$Cell_Type_2
cellchat.rime <- createCellChat(object = temp, group.by = "ident", assay = "RNA")
cellchat.rime@DB <- CellChatDB.mouse

rm(blca_sub, temp)

# processing
cellchat.ctrl <- cellchat.ctrl %>%
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")

cellchat.capsaicin <- cellchat.capsaicin %>%
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")

cellchat.lc_6ohda <- cellchat.lc_6ohda %>%
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")

cellchat.anti_pd1 <- cellchat.anti_pd1 %>%
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")

cellchat.rime <- cellchat.rime %>%
    subsetData() %>%
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    computeCommunProb(type = "triMean") %>%
    filterCommunication(min.cells = 10) %>%
    computeCommunProbPathway() %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")
object.list <- list(
    "CTRL" = cellchat.ctrl,
    "Capsaicin" = cellchat.capsaicin,
    "LC_6OHDA" = cellchat.lc_6ohda,
    "anti_PD1" = cellchat.anti_pd1,
    "Rime" = cellchat.rime
)

cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(object.list, file = "object.list.rds")
rm(cellchat.ctrl, cellchat.capsaicin, cellchat.lc_6ohda, cellchat.anti_pd1, cellchat.rime)

object.list <- readRDS("object.list.rds")

pdf("Compare_Conditions/Circle.Network.All.Weight.Conditions.pdf", 15, 10)
par(mfrow = c(2, 3), xpd = TRUE)
lapply(names(object.list), function(x) {
    groupSize <- as.numeric(table(object.list[[x]]@idents))
    p <- netVisual_circle(object.list[[x]]@net$weight,
        color.use = mycolor$Cell_Type_2[cell_types],
        vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE,
        title.name = paste0("Interaction weights/strength: ", x),
    )
})
dev.off()

pdf("Compare_Conditions/Circle.Network.All.Count.Conditions.pdf", 15, 10)
par(mfrow = c(2, 3), xpd = TRUE)
lapply(names(object.list), function(x) {
    groupSize <- as.numeric(table(object.list[[x]]@idents))
    p <- netVisual_circle(object.list[[x]]@net$count,
        color.use = mycolor$Cell_Type_2[cell_types],
        vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE,
        title.name = paste0("Interaction counts: ", x),
    )
})
dev.off()

## 3.2. compare overall network ----
p1 <- compareInteractions(cellchat.merge,
    show.legend = FALSE, color.use = mycolor$Condition %>% unname(),
    group = 1:5
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- compareInteractions(cellchat.merge,
    show.legend = FALSE, color.use = mycolor$Condition %>% unname(),
    group = 1:5, measure = "weight"
) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
pdf("Compare_Conditions/Bar.Compare_Interactions.pdf", 5, 3)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()
rm(p1, p2)

pdf("Compare_Conditions/Circle.Network.Capsaicin_vs_Ctrl.Count.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 2), weight.scale = TRUE
)
dev.off()

pdf("Compare_Conditions/Circle.Network.LC_6OHDA_vs_Ctrl.Count.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 3), weight.scale = TRUE
)
dev.off()

pdf("Compare_Conditions/Circle.Network.anti_PD1_vs_Ctrl.Count.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 4), weight.scale = TRUE
)
dev.off()

pdf("Compare_Conditions/Circle.Network.Rime_vs_Ctrl.Count.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 5), weight.scale = TRUE
)
dev.off()

pdf("Compare_Conditions/Circle.Network.Capsaicin_vs_Ctrl.Weight.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 2), weight.scale = TRUE, measure = "weight"
)
dev.off()

pdf("Compare_Conditions/Circle.Network.LC_6OHDA_vs_Ctrl.Weight.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 3), weight.scale = TRUE, measure = "weight"
)
dev.off()

pdf("Compare_Conditions/Circle.Network.anti_PD1_vs_Ctrl.Weight.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 4), weight.scale = TRUE, measure = "weight"
)
dev.off()

pdf("Compare_Conditions/Circle.Network.Rime_vs_Ctrl.Weight.pdf", 8, 9)
netVisual_diffInteraction(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types], comparison = c(1, 5), weight.scale = TRUE, measure = "weight"
)
dev.off()

# pdf("Compare_Conditions/Circle.Network.Capsaicin.Weight.CAF.pdf", 8, 9)
# netVisual_diffInteraction(cellchat.merge,
#     color.use = mycolor$Cell_Type_2[cell_types],
#     weight.scale = TRUE, measure = "weight", comparison = c(1, 2),
#     sources.use = c("iCAF", "myCAF", "apCAF")
# )
# dev.off()

# pdf("Compare_Conditions/Circle.Network.LC_6OHDA.Weight.CAF.pdf", 8, 9)
# netVisual_diffInteraction(cellchat.merge,
#     color.use = mycolor$Cell_Type_2[cell_types],
#     weight.scale = TRUE, measure = "weight", comparison = c(1, 3),
#     sources.use = c("iCAF", "myCAF", "apCAF")
# )
# dev.off()

pdf("Compare_Conditions/Heatmap.Capsaicin_vs_Ctrl.Weight.pdf", 6, 5.5)
netVisual_heatmap(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types],
    comparison = c(1, 2), measure = "weight"
)
dev.off()

pdf("Compare_Conditions/Heatmap.LC_6OHDA_vs_Ctrl.Weight.pdf", 6, 5.5)
netVisual_heatmap(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types],
    comparison = c(1, 3), measure = "weight"
)
dev.off()

pdf("Compare_Conditions/Heatmap.anti_PD1_vs_Ctrl.Weight.pdf", 6, 5.5)
netVisual_heatmap(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types],
    comparison = c(1, 4), measure = "weight"
)
dev.off()

pdf("Compare_Conditions/Heatmap.Rime_vs_Ctrl.Weight.pdf", 6, 5.5)
netVisual_heatmap(cellchat.merge,
    color.use = mycolor$Cell_Type_2[cell_types],
    comparison = c(1, 5), measure = "weight"
)
dev.off()

num.link <- sapply(object.list, function(x) {
    rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
p.list <- list()
for (i in seq_along(object.list)) {
    p.list[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
        color.use = mycolor$Cell_Type_2[cell_types],
        title = names(object.list)[i], weight.MinMax = weight.MinMax
    ) + lims(x = c(0, 45), y = c(0, 25))
}

pdf("Compare_Conditions/Scatter.In_out.Compare.pdf", 12, 7)
patchwork::wrap_plots(plots = p.list, ncol = 3)
dev.off()
rm(p.list, num.link, weight.MinMax, i)

## 3.3. identify altered pathways ----
# scatter altered pathways
p.list <- lapply(cell_types, function(cell_type) {
    netAnalysis_signalingChanges_scatter(cellchat.merge,
        comparison = c(1, 2), idents.use = cell_type
    )
})
pdf("Compare_Conditions/Scatter.Signaling_Changes.Capsaicin.pdf", 20, 15)
patchwork::wrap_plots(plots = p.list, ncol = 4)
dev.off()

p.list <- lapply(cell_types, function(cell_type) {
    netAnalysis_signalingChanges_scatter(cellchat.merge,
        signaling.exclude = c("COLLAGEN", "LAMININ"),
        comparison = c(1, 2), idents.use = cell_type
    )
})
pdf("Compare_Conditions/Scatter.Signaling_Changes.noECM.Capsaicin.pdf", 20, 15)
patchwork::wrap_plots(plots = p.list, ncol = 4)
dev.off()

p.list <- lapply(cell_types, function(cell_type) {
    netAnalysis_signalingChanges_scatter(cellchat.merge,
        comparison = c(1, 3), idents.use = cell_type
    )
})
pdf("Compare_Conditions/Scatter.Signaling_Changes.LC_6OHDA.pdf", 20, 15)
patchwork::wrap_plots(plots = p.list, ncol = 4)
dev.off()
rm(p.list)

p.list <- lapply(cell_types, function(cell_type) {
    netAnalysis_signalingChanges_scatter(cellchat.merge,
        signaling.exclude = c("COLLAGEN", "LAMININ"),
        comparison = c(1, 3), idents.use = cell_type
    )
})
pdf("Compare_Conditions/Scatter.Signaling_Changes.noECM.LC_6OHDA.pdf", 20, 15)
patchwork::wrap_plots(plots = p.list, ncol = 4)
dev.off()
rm(p.list)

# bar plot of ranked pathways
p1 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 2),
    targets.use = NULL, stacked = TRUE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
p2 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 2),
    targets.use = NULL, stacked = FALSE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
pdf("Compare_Conditions/Bar.Ranked.Pathways.Capsaicin.pdf", 7.5, 15)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()

p1 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 3),
    targets.use = NULL, stacked = TRUE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
p2 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 3),
    targets.use = NULL, stacked = FALSE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
pdf("Compare_Conditions/Bar.Ranked.Pathways.LC_6OHDA.pdf", 7.5, 15)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()
rm(p1, p2)

p1 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 4),
    targets.use = NULL, stacked = TRUE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
p2 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 4),
    targets.use = NULL, stacked = FALSE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
pdf("Compare_Conditions/Bar.Ranked.Pathways.anti_PD1.pdf", 7.5, 15)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()
rm(p1, p2)

p1 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 5),
    targets.use = NULL, stacked = TRUE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
p2 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, comparison = c(1, 5),
    targets.use = NULL, stacked = FALSE, do.stat = TRUE
) + scale_fill_manual(values = mycolor$Condition)
pdf("Compare_Conditions/Bar.Ranked.Pathways.Rime.pdf", 7.5, 15)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()
rm(p1, p2)

dir.create("Compare_Conditions/Pathway_Level")

lapply(pathways.selected %>% setdiff(c("FGF", "IL6")), function(pathways.show) {
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

    pdf(paste0("Compare_Conditions/Pathway_Level/Circle.Network.Pathway.", pathways.show, ".pdf"), 15, 10)
    par(mfrow = c(2, 3), xpd = TRUE)
    for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]],
            color.use = mycolor$Cell_Type_2[cell_types],
            signaling = pathways.show,
            layout = "circle", edge.weight.max = weight.max[1],
            edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i])
        )
    }
    dev.off()

    par(mfrow = c(2, 3), xpd = TRUE)
    ht <- list()
    for (i in seq_along(object.list)) {
        ht[[i]] <- netVisual_heatmap(object.list[[i]],
            color.use = mycolor$Cell_Type_2[cell_types],
            signaling = pathways.show, color.heatmap = "Reds",
            title.name = paste(pathways.show, "signaling ", names(object.list)[i])
        )
    }

    pdf(paste0("Compare_Conditions/Pathway_Level/Heatmap.Pathway.", pathways.show, ".pdf"), 15, 10)
    ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm"))
    dev.off()

    pdf(paste0("Compare_Conditions/Pathway_Level/Chord.Pathway.", pathways.show, ".pdf"), 15, 10)
    par(mfrow = c(2, 3), xpd = TRUE)
    for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]],
            color.use = mycolor$Cell_Type_2[cell_types],
            signaling = pathways.show,
            layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i])
        )
    }
    dev.off()
})
rm(pathways.show, ht, i, weight.max)

# ## 3.4. altered interactions (Capsaicin) ----
# dir.create("Compare_Conditions/LR_Pairs_Capsaicin")

# # identify DEG
# pos.dataset <- "Capsaicin"
# features.name <- paste0(pos.dataset, ".merged")

# cellchat.merge <- identifyOverExpressedGenes(cellchat.merge,
#     group.dataset = "datasets", pos.dataset = pos.dataset,
#     features.name = features.name,
#     only.pos = FALSE,
#     thresh.pc = 0.1, thresh.fc = 0.25, thresh.p = 0.05, group.DE.combined = FALSE
# )
# net.Capsaicin <- netMappingDEG(cellchat.merge, features.name = features.name, variable.all = TRUE)
# net.Capsaicin.up <- subsetCommunication(cellchat.merge,
#     net = net.Capsaicin,
#     datasets = "Capsaicin", ligand.logFC = 0.25, receptor.logFC = NULL
# )
# net.Capsaicin.down <- subsetCommunication(cellchat.merge,
#     net = net.Capsaicin,
#     datasets = "CTRL", ligand.logFC = -0.25, receptor.logFC = NULL
# )

# write.table(net.Capsaicin.up,
#     file = "Compare_Conditions/LR_Pairs_Capsaicin/DEG.Up.LR.Pairs.Capsaicin.txt",
#     sep = "\t", quote = FALSE, row.names = FALSE
# )
# write.table(net.Capsaicin.down,
#     file = "Compare_Conditions/LR_Pairs_Capsaicin/DEG.Down.LR.Pairs.Capsaicin.txt",
#     sep = "\t", quote = FALSE, row.names = FALSE
# )
# rm(features.name, net.Capsaicin, pos.dataset)

# # plot significatnt L-R pairs
# lapply(cell_types, function(cell_type) {
#     # Increased incoming
#     pairLR.use <- net.Capsaicin.up %>%
#         # filter(pathway_name %in% setdiff(pathways.selected, c("COLLAGEN", "LAMININ")), ) %>%
#         filter(target == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_types,
#             targets.use = cell_type,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 2), max.dataset = 2, angle.x = 45,
#             title.name = paste0("Increased ", cell_type, " incoming signaling in Capsaicin"),
#         ) +
#             facet_grid(
#                 cols = vars(source),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_Capsaicin/Bubble.", cell_type, "_Increased_Incoming.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }

#     # Decreased incoming
#     pairLR.use <- net.Capsaicin.down %>%
#         # filter(pathway_name %in% setdiff(pathways.selected, c("COLLAGEN", "LAMININ")), ) %>%
#         filter(target == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_types,
#             targets.use = cell_type,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 2), max.dataset = 1, angle.x = 45,
#             title.name = paste0("Decreased ", cell_type, " incoming signaling in Capsaicin"),
#         ) +
#             facet_grid(
#                 cols = vars(source),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_Capsaicin/Bubble.", cell_type, "_Decreased_Incoming.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }

#     # Increased outgoing
#     pairLR.use <- net.Capsaicin.up %>%
#         # filter(pathway_name %in% setdiff(pathways.selected, c("COLLAGEN", "LAMININ")), ) %>%
#         filter(source == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_type,
#             targets.use = cell_types,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 2), max.dataset = 2, angle.x = 45,
#             title.name = paste0("Increased ", cell_type, " outgoing signaling in Capsaicin"),
#         ) +
#             facet_grid(
#                 cols = vars(target),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_Capsaicin/Bubble.", cell_type, "_Increased_Outgoing.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }

#     # Decreased outgoing
#     pairLR.use <- net.Capsaicin.down %>%
#         # filter(pathway_name %in% setdiff(pathways.selected, c("COLLAGEN", "LAMININ")), ) %>%
#         filter(source == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_type,
#             targets.use = cell_types,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 2), max.dataset = 1, angle.x = 45,
#             title.name = paste0("Decreased ", cell_type, " outgoing signaling in Capsaicin"),
#         ) +
#             facet_grid(
#                 cols = vars(target),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_Capsaicin/Bubble.", cell_type, "_Decreased_Outgoing.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }
# })

# rm(p, pairLR.use, cell_type, height)

# ## 3.5. altered interactions (LC_6OHDA) ----
# dir.create("Compare_Conditions/LR_Pairs_LC_6OHDA")

# # identify DEG
# pos.dataset <- "LC_6OHDA"
# features.name <- paste0(pos.dataset, ".merged")

# cellchat.merge <- identifyOverExpressedGenes(cellchat.merge,
#     group.dataset = "datasets", pos.dataset = pos.dataset,
#     features.name = features.name,
#     only.pos = FALSE,
#     thresh.pc = 0.1, thresh.fc = 0.25, thresh.p = 0.05, group.DE.combined = FALSE
# )
# net.LC_6OHDA <- netMappingDEG(cellchat.merge, features.name = features.name, variable.all = TRUE)
# net.LC_6OHDA.up <- subsetCommunication(cellchat.merge,
#     net = net.LC_6OHDA,
#     datasets = "LC_6OHDA", ligand.logFC = 0.25, receptor.logFC = NULL
# )
# net.LC_6OHDA.down <- subsetCommunication(cellchat.merge,
#     net = net.LC_6OHDA,
#     datasets = "CTRL", ligand.logFC = -0.25, receptor.logFC = NULL
# )

# write.table(net.LC_6OHDA.up,
#     file = "Compare_Conditions/LR_Pairs_LC_6OHDA/DEG.Up.LR.Pairs.LC_6OHDA.txt",
#     sep = "\t", quote = FALSE, row.names = FALSE
# )
# write.table(net.LC_6OHDA.down,
#     file = "Compare_Conditions/LR_Pairs_LC_6OHDA/DEG.Down.LR.Pairs.LC_6OHDA.txt",
#     sep = "\t", quote = FALSE, row.names = FALSE
# )
# rm(features.name, net.LC_6OHDA, pos.dataset)

# # plot significatnt L-R pairs
# lapply(cell_types, function(cell_type) {
#     # Increased incoming
#     pairLR.use <- net.LC_6OHDA.up %>%
#         filter(target == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_types,
#             targets.use = cell_type,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 3), max.dataset = 3, angle.x = 45,
#             title.name = paste0("Increased ", cell_type, " incoming signaling in LC_6OHDA"),
#         ) +
#             facet_grid(
#                 cols = vars(source),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_LC_6OHDA/Bubble.", cell_type, "_Increased_Incoming.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }

#     # Decreased incoming
#     pairLR.use <- net.LC_6OHDA.down %>%
#         filter(target == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_types,
#             targets.use = cell_type,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 3), max.dataset = 1, angle.x = 45,
#             title.name = paste0("Decreased ", cell_type, " incoming signaling in LC_6OHDA"),
#         ) +
#             facet_grid(
#                 cols = vars(source),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_LC_6OHDA/Bubble.", cell_type, "_Decreased_Incoming.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }

#     # Increased outgoing
#     pairLR.use <- net.LC_6OHDA.up %>%
#         filter(source == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_type,
#             targets.use = cell_types,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 3), max.dataset = 3, angle.x = 45,
#             title.name = paste0("Increased ", cell_type, " outgoing signaling in LC_6OHDA"),
#         ) +
#             facet_grid(
#                 cols = vars(target),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_LC_6OHDA/Bubble.", cell_type, "_Increased_Outgoing.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }

#     # Decreased outgoing
#     pairLR.use <- net.LC_6OHDA.down %>%
#         filter(source == cell_type) %>%
#         dplyr::select(interaction_name)
#     if (nrow(pairLR.use) > 0) {
#         p <- netVisual_bubble(cellchat.merge,
#             sources.use = cell_type,
#             targets.use = cell_types,
#             pairLR.use = pairLR.use,
#             remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#             comparison = c(1, 3), max.dataset = 1, angle.x = 45,
#             title.name = paste0("Decreased ", cell_type, " outgoing signaling in LC_6OHDA"),
#         ) +
#             facet_grid(
#                 cols = vars(target),
#                 rows = vars(pathway_name), scales = "free", space = "free"
#             ) +
#             theme(
#                 panel.grid.major = element_line(color = "gray90", size = 0.25),
#                 strip.background = element_rect(fill = "gray80"),
#                 strip.text = element_text(colour = "black"),
#                 strip.text.y = element_text(angle = 0)
#             )
#         height <- length(unique(p$data$interaction_name)) / 6.5 + 2.5
#         pdf(paste0("Compare_Conditions/LR_Pairs_LC_6OHDA/Bubble.", cell_type, "_Decreased_Outgoing.pdf"), 10, height)
#         plot(p)
#         dev.off()
#     }
# })

# rm(p, pairLR.use, cell_type, height)

# ## 3.6. zoom: iCAF to immune cells (Capsaicin) ----
# dir.create("Zoom_CAF_to_Immune")
# table(blca_obj$Cell_Type_2, blca_obj$Condition)

# immune.types <- c(
#     "T_CD8", "T_CD4", "Treg", "NK", "B",
#     "Macrophage", "Monocyte", "cDC", "pDC", "Neutrophil"
# )
# # compare pathway intensity
# for (cell.type.select in immune.types) {
#     for (source in c("iCAF", "apCAF", "myCAF")) {
#         # merged intensity
#         p1 <- rankNet(cellchat,
#             mode = "single", measure = "weight", sources.use = source,
#             targets.use = cell.type.select, stacked = TRUE, do.stat = TRUE
#         )
#         p2 <- rankNet(cellchat,
#             mode = "single", measure = "weight", sources.use = source,
#             targets.use = cell.type.select, stacked = FALSE, do.stat = TRUE
#         )

#         pdf(paste0("Zoom_CAF_to_Immune/Bar.Ranked.Single.Capsaicin.", source, "_to_", cell.type.select, ".pdf"), 7, 3.5)
#         patchwork::wrap_plots(p1, p2, ncol = 2) %>% plot()
#         dev.off()

#         # compare intensity
#         p1 <- rankNet(cellchat.merge,
#             mode = "comparison", measure = "weight", sources.use = source, comparison = c(1, 2),
#             targets.use = cell.type.select, stacked = TRUE, do.stat = TRUE,
#         ) + scale_fill_manual(values = mycolor$Condition)
#         p2 <- rankNet(cellchat.merge,
#             mode = "comparison", measure = "weight", sources.use = source, comparison = c(1, 2),
#             targets.use = cell.type.select, stacked = FALSE, do.stat = TRUE
#         ) + scale_fill_manual(values = mycolor$Condition)
#         pdf(paste0("Zoom_CAF_to_Immune/Bar.Ranked.CMP.Capsaicin.", source, "_to_", cell.type.select, ".pdf"), 7, 4)
#         patchwork::wrap_plots(p1, p2, ncol = 2) %>% plot()
#         dev.off()
#     }
# }
# rm(cell.type.select, source, p1, p2)

# # L-R interactions
# pdf("Zoom_CAF_to_Immune/Bubble.CAF_to_T.pdf", 6, 12)
# netVisual_bubble(cellchat,
#     signaling = cellchat@netP$pathways %>% setdiff(c("COLLAGEN", "LAMININ", "MHC-I", "MHC-II")),
#     sources.use = c("iCAF", "apCAF", "myCAF"),
#     targets.use = c("T_CD8", "T_CD4", "Treg"),
#     remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
# ) +
#     facet_grid(
#         cols = vars(source),
#         rows = vars(pathway_name), scales = "free", space = "free"
#     ) +
#     theme(
#         panel.grid.major = element_line(color = "gray90", size = 0.25),
#         strip.background = element_rect(fill = "gray80"),
#         strip.text = element_text(colour = "black"),
#         strip.text.y = element_text(angle = 0)
#     )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Bubble.CAF_to_Immune.pdf", 9, 20)
# netVisual_bubble(cellchat,
#     signaling = cellchat@netP$pathways %>% setdiff(c("COLLAGEN", "LAMININ", "MHC-I", "MHC-II")),
#     sources.use = c("iCAF", "apCAF", "myCAF"),
#     targets.use = immune.types,
#     remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
# ) +
#     facet_grid(
#         cols = vars(source),
#         rows = vars(pathway_name), scales = "free", space = "free"
#     ) +
#     theme(
#         panel.grid.major = element_line(color = "gray90", size = 0.25),
#         strip.background = element_rect(fill = "gray80"),
#         strip.text = element_text(colour = "black"),
#         strip.text.y = element_text(angle = 0)
#     )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Chord.gene.iCAF_to_T.pdf", 8, 8)
# netVisual_chord_gene(cellchat,
#     sources.use = c("iCAF"),
#     targets.use = c("T_CD8", "T_CD4", "Treg"),
#     lab.cex = 0.5, legend.pos.y = 30,
#     # color.use = mycolor$Cell_Type_2[immune.types]
# )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Chord.gene.iCAF_to_Immune.pdf", 8, 8)
# netVisual_chord_gene(cellchat,
#     sources.use = c("iCAF"),
#     targets.use = immune.types,
#     lab.cex = 0.5, legend.pos.y = 30,
#     # color.use = mycolor$Cell_Type_2[immune.types]
# )
# dev.off()

# # expression of selected genes
# plot.data <- netVisual_bubble(cellchat,
#     signaling = cellchat@netP$pathways %>% setdiff(c("COLLAGEN", "LAMININ", "MHC-I", "MHC-II")),
#     sources.use = c("iCAF", "apCAF", "myCAF"),
#     targets.use = immune.types,
#     remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
# )$data

# gene.selected <- list(
#     "ligand" = unique(plot.data$ligand) %>% intersect(rownames(blca_obj)),
#     "receptor" = unique(plot.data$receptor) %>% intersect(rownames(blca_obj))
# )

# pdf("Zoom_CAF_to_Immune/Violin.Expression.Ligands.pdf", 40, 27)
# VlnPlot(blca_obj,
#     features = gene.selected$ligand, cols = mycolor$Cell_Type_2,
#     group.by = "Cell_Type_2",
#     pt.size = 0, ncol = 8
# )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Violin.Expression.Receptors.pdf", 40, 24)
# VlnPlot(blca_obj,
#     features = gene.selected$receptor, cols = mycolor$Cell_Type_2,
#     group.by = "Cell_Type_2",
#     pt.size = 0, ncol = 8
# )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Dot.Expression.Ligands.pdf", 9, 15)
# DotPlot(blca_obj,
#     features = gene.selected$ligand %>% rev(), scale = FALSE,
#     group.by = "Cell_Type_2",
# ) +
#     coord_flip() +
#     theme_classic() +
#     scale_color_viridis_c() +
#     theme(
#         axis.text.x = element_text(angle = 45, hjust = 1)
#     )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Dot.Expression.Receptors.pdf", 9, 12)
# DotPlot(blca_obj,
#     features = gene.selected$receptor %>% rev(), scale = FALSE,
#     group.by = "Cell_Type_2",
# ) +
#     coord_flip() +
#     theme_classic() +
#     scale_color_viridis_c() +
#     theme(
#         axis.text.x = element_text(angle = 45, hjust = 1)
#     )
# dev.off()

# rm(plot.data, gene.selected)

# # differential interactions
# pairLR.use <- net.Capsaicin.up %>%
#     # filter(!pathway_name %in% c("COLLAGEN", "LAMININ", "MHC-I", "MHC-II")) %>%
#     filter(target %in% c("T_CD8", "T_CD4", "Treg", "Neutrophil")) %>%
#     filter(source %in% c("iCAF", "apCAF", "myCAF")) %>%
#     dplyr::select(interaction_name)

# pdf("Zoom_CAF_to_Immune/Bubble.Increase.CAF_to_Tcells.pdf", 7, 6)
# netVisual_bubble(cellchat.merge,
#     sources.use = c("iCAF", "apCAF", "myCAF"),
#     targets.use = c("T_CD8", "T_CD4", "Treg", "Neutrophil"),
#     pairLR.use = pairLR.use,
#     remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#     comparison = c(1, 2), max.dataset = 2, angle.x = 45,
#     title.name = paste0("Increased CAF to T cell signal in Capsaicin"),
# ) +
#     facet_grid(
#         cols = vars(source),
#         rows = vars(pathway_name), scales = "free", space = "free"
#     ) +
#     theme(
#         panel.grid.major = element_line(color = "gray90", size = 0.25),
#         strip.background = element_rect(fill = "gray80"),
#         strip.text = element_text(colour = "black"),
#         strip.text.y = element_text(angle = 0)
#     )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Bubble.Increase.CAF_to_Tcells.All.pdf", 7, 6)
# netVisual_bubble(cellchat.merge,
#     sources.use = c("iCAF", "apCAF", "myCAF"),
#     targets.use = c("T_CD8", "T_CD4", "Treg", "Neutrophil"),
#     pairLR.use = pairLR.use,
#     remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#     comparison = c(1, 2), angle.x = 45,
#     title.name = paste0("Increased CAF to T cell signal in Capsaicin"),
# ) +
#     facet_grid(
#         cols = vars(source),
#         rows = vars(pathway_name), scales = "free", space = "free"
#     ) +
#     theme(
#         panel.grid.major = element_line(color = "gray90", size = 0.25),
#         strip.background = element_rect(fill = "gray80"),
#         strip.text = element_text(colour = "black"),
#         strip.text.y = element_text(angle = 0)
#     )
# dev.off()

# pairLR.use <- net.Capsaicin.up %>%
#     # filter(!pathway_name %in% c("COLLAGEN", "LAMININ", "MHC-I", "MHC-II")) %>%
#     filter(target %in% immune.types) %>%
#     filter(source %in% c("iCAF", "apCAF", "myCAF")) %>%
#     dplyr::select(interaction_name)

# pdf("Zoom_CAF_to_Immune/Bubble.Increase.CAF_to_Immune.pdf", 12, 8)
# netVisual_bubble(cellchat.merge,
#     sources.use = c("iCAF", "apCAF", "myCAF"),
#     targets.use = immune.types,
#     pairLR.use = pairLR.use,
#     remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#     comparison = c(1, 2), max.dataset = 2, angle.x = 45,
#     title.name = paste0("Increased CAF to T cell signal in Capsaicin"),
# ) +
#     facet_grid(
#         cols = vars(source),
#         rows = vars(pathway_name), scales = "free", space = "free"
#     ) +
#     theme(
#         panel.grid.major = element_line(color = "gray90", size = 0.25),
#         strip.background = element_rect(fill = "gray80"),
#         strip.text = element_text(colour = "black"),
#         strip.text.y = element_text(angle = 0)
#     )
# dev.off()

# pdf("Zoom_CAF_to_Immune/Bubble.Increase.CAF_to_Immune.All.pdf", 12, 8)
# netVisual_bubble(cellchat.merge,
#     sources.use = c("iCAF", "apCAF", "myCAF"),
#     targets.use = immune.types,
#     pairLR.use = pairLR.use,
#     remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
#     comparison = c(1, 2), angle.x = 45,
#     title.name = paste0("Increased CAF to T cell signal in Capsaicin"),
# ) +
#     facet_grid(
#         cols = vars(source),
#         rows = vars(pathway_name), scales = "free", space = "free"
#     ) +
#     theme(
#         panel.grid.major = element_line(color = "gray90", size = 0.25),
#         strip.background = element_rect(fill = "gray80"),
#         strip.text = element_text(colour = "black"),
#         strip.text.y = element_text(angle = 0)
#     )
# dev.off()
# rm(pairLR.use)

# save data ----
gc()
saveRDS(mycolor, file = "mycolor.rds")
saveRDS(blca_obj, file = "blca_obj.rds")
rm(blca_obj, object.list)
save.image("04.Communication.RData")
