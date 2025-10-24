setwd("/thinker/tb269store/user/liuzhenyu/Project/BLCA_Neuro/BLCA_Neuro/results/04.Conmmunication")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(CellChat)

load("../01.Cell_Type/01.Cell_Type.RData")
rm(DEG.cancer.conditions, DEG.cancer.clusters, DEG.cancer.states, data.dir, marker.genes)

load("04.Communication.RData")

# 1. Integrate cell types ----
## 1.1. merge cell types from sub-objects ----
blca.obj$Cell_Type_2 <- "none"

immune.obj <- readRDS("../01.Cell_Type/sub1_Immune/immune.obj.rds")
blca.obj@meta.data[immune.obj$Cell, "Cell_Type_2"] <- immune.obj$Cell_Type_2 %>% as.character()

cancer.obj <- readRDS("../01.Cell_Type/sub2_Cancer/cancer.obj.rds")
blca.obj@meta.data[cancer.obj$Cell, "Cell_Type_2"] <- cancer.obj$Cell_Type_1 %>% as.character()

stromal.obj <- readRDS("../01.Cell_Type/sub3_Stromal/stromal.obj.rds")
blca.obj@meta.data[stromal.obj$Cell, "Cell_Type_2"] <- stromal.obj$Cell_Type_2 %>% as.character()

blca.obj@meta.data %>% filter(Cell_Type_2 == "none") %>% pull(Cell_Type_1) %>% table()
blca.obj@meta.data[blca.obj$Cell_Type_2 == "none", "Cell_Type_2"] <- "Bladder_Epithelial"
table(blca.obj$Cell_Type_2, blca.obj$Cell_Type_1)

mycolor$Cell_Type_2["Bladder_Epithelial"] <- mycolor$Cell_Type_1["Bladder_epithelial"]
mycolor$Cell_Type_2["Cancer"] <- mycolor$Cell_Type_1["Cancer"]

table(blca.obj$Cell_Type_2 %in% names(mycolor$Cell_Type_2))

blca.obj$Cell_Type_2 <- factor(
    blca.obj$Cell_Type_2,
    levels = intersect(names(mycolor$Cell_Type_2), unique(blca.obj$Cell_Type_2))
)

pdf("UMAP.Cell_Type_2.pdf", 9, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_2", label = TRUE) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

pdf("UMAP.Cell_Type_2.NoLabel.pdf", 9, 6)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_2", label = FALSE) +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

pdf("UMAP.Cell_Type_2.split.Condition.pdf", 12, 4)
DimPlot(blca.obj, reduction = "umap", group.by = "Cell_Type_2", label = FALSE, split.by = "Condition") +
    scale_color_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    labs(title = "UMAP plot colored by Cell_Type_2")
dev.off()

p <- DimPlot(blca.obj,
    reduction = "umap", group.by = "Cell_Type_2", split.by = "Condition",
    label = FALSE, raster = TRUE, pt.size = NULL,
    cols = mycolor$Cell_Type_2
) +
    theme_classic()

pdf("UMAP.Cell_Type_2.split.density.pdf", 12, 4)
plot(p + stat_density_2d(aes_string(x = "umap_1", y = "umap_2", fill = "after_stat(level)"),
    linewidth = 0, geom = "density_2d_filled",
    n = 150, h = c(1.2, 1.2)
) +
    scale_fill_viridis_c(option = "D") +
    guides(fill = guide_colorbar(title = "Density")) +
    theme(legend.position = "right"))
dev.off()
rm(p)

rm(immune.obj, cancer.obj, stromal.obj)

## 1.2. Cell types cross conditions ----
plot.data <- as.data.frame(table(blca.obj$Cell_Type_2, blca.obj$Condition))
colnames(plot.data) <- c("Cell_Type_2", "Condition", "Count")

pdf("Barplot.Cell_Type_2.Condition.pdf", 9, 5)
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
plot.data <- as.data.frame(table(blca.obj$Cell_Type_2, blca.obj$Sample))
colnames(plot.data) <- c("Cell_Type_2", "Sample", "Count")
plot.data$Condition <- gsub("_00.+?*$", "", plot.data$Sample)
table(plot.data$Condition)
table(plot.data$Sample)

temp <- table(blca.obj$Sample)
plot.data$Fraction <- plot.data$Count / temp[plot.data$Sample] * 100
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC_6OHDA"))

pdf("Boxplot.Fraction.Cell_Type_2.pdf", 8, 12)
ggplot(plot.data, aes(x = Condition, y = Fraction)) +
    geom_boxplot(aes(fill = Cell_Type_2), outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 1) +
    ggpubr::stat_compare_means(
        label = "p.signif", label.y.npc = 1, size = 2.5,
        comparisons = list(c("CTRL", "Capsaicin"), c("Capsaicin", "LC_6OHDA"), c("CTRL", "LC_6OHDA"))
    ) +
    facet_wrap(~Cell_Type_2, ncol = 5, scales = "free") +
    labs(title = "Cell Type composition across Samples") +
    ylab("Fraction of all cells (%)") +
    scale_fill_manual(values = mycolor$Cell_Type_2) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NoLegend()
dev.off()
rm(plot.data, temp)

## 1.3. expression of selected genes ----
gene.selected <- c("Ngf", "Nrg1", "Bdnf", "Ntf3", "Gdnf", "Cntf", "Lif")

plot.data <- DotPlot(blca.obj,
    features = gene.selected, scale = FALSE,
    group.by = "Condition", split.by = "Cell_Type_2", cols = rep(RColorBrewer::brewer.pal(9, "YlGnBu"), 5)
)$data

plot.data$Condition <- gsub("_.+?$", "", plot.data$id)
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC"))
plot.data$Cell_Type_2 <- gsub("Capsaicin_|CTRL_|LC_6OHDA_", "", plot.data$id)
plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(blca.obj$Cell_Type_2))
plot.data <- plot.data[!is.na(plot.data$Cell_Type_2), ]
colnames(plot.data)

pdf("Dot.Expression.Neuro_factor.pdf", 10, 8)
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

plot.data <- DotPlot(blca.obj,
    features = gene.selected[1], scale = FALSE,
    group.by = "Condition", split.by = "Cell_Type_2", cols = rep(RColorBrewer::brewer.pal(9, "YlGnBu"), 5)
)$data

plot.data$Condition <- gsub("_.+?$", "", plot.data$id)
plot.data$Condition <- factor(plot.data$Condition, levels = c("CTRL", "Capsaicin", "LC"))
plot.data$Cell_Type_2 <- gsub("Capsaicin_|CTRL_|LC_6OHDA_", "", plot.data$id)
plot.data$Cell_Type_2 <- factor(plot.data$Cell_Type_2, levels = levels(blca.obj$Cell_Type_2))
plot.data <- plot.data[!is.na(plot.data$Cell_Type_2), ]
colnames(plot.data)

pdf("Dot.Expression.Ngf.pdf", 10, 7)
ggplot(plot.data, aes(x = Condition, y = features.plot)) +
    geom_point(aes(size = pct.exp, fill = avg.exp), shape = 21) +
    scale_fill_viridis_c() +
    facet_wrap(~Cell_Type_2, scales = "free_x", ncol = 5) +
    theme_classic() +
    labs(x = "Condition", y = "Gene") +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Expression of Neuronal Factors")
dev.off()
rm(plot.data, gene.selected)

# 2. Joint CellChat analysis ----
## 2.1. prepare CellChat object ----
cell_types <- levels(blca.obj$Cell_Type_2) %>%
    grep("Cycling", ., invert = TRUE, value = TRUE) %>%
    setdiff(c("Adipocyte", "Smooth_Muscle"))
table(blca.obj$Cell_Type_2)[cell_types] %>% sort()

blca.sub <- subset(blca.obj, Cell_Type_2 %in% cell_types)
Idents(blca.sub) <- blca.sub$Cell_Type_2

cellchat <- createCellChat(object = blca.sub, group.by = "ident", assay = "RNA")
rm(blca.sub)

showDatabaseCategory(CellChatDB.mouse)
cellchat@DB <- CellChatDB.mouse

# filter genes & LR pairs
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) # 3212 variable L-R pairs

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
    vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE,
     title.name = "Number of interactions"
)
dev.off()

pdf("Circle.Network.All.Weight.pdf", 6, 7)
netVisual_circle(cellchat@net$weight,
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
netAnalysis_signalingRole_scatter(cellchat)
dev.off()

pdf("Heatmap.In_out.pdf", 10, 15)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 30)
(ht1 + ht2)
dev.off()
rm(ht1, ht2)

# global communication patterns
dir.create("Global_Communication_Patterns")
library(NMF)
library(ggalluvial)

pathways.selected <- c(
    "TGFb", "VEGF", "EGF", "CXCL", "IL1", "IL2",
    "MHC-I", "MHC-II", "PD-L1", "PDL2", "COLLAGEN", "LAMININ"
)

## Outgoing
pdf("Global_Communication_Patterns/Elbow.Outgoing.pdf", 6, 2.5)
selectK(cellchat, pattern = "outgoing")
dev.off()

nPatterns <- 5

pdf("Global_Communication_Patterns/Heatmap.Pattern.Out.pdf", 8, 4)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
dev.off()

pdf("Global_Communication_Patterns/River.Pattern.Out.pdf", 8, 4)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.Out.pdf", 4.5, 15)
netAnalysis_dot(cellchat, pattern = "outgoing") + coord_flip()
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.Out.Selected.pdf", 7, 3)
netAnalysis_dot(cellchat,
    pattern = "outgoing",
    pathway.show = union(c("COLLAGEN", "LAMININ", "APP", "MHC-I", "FN1", "THBS", "CCL", "SPP1", "ADGRE", "MHC-II",
    "CypA", "CD80", "PTN", "CD86", "MIF", "GALECTIN", "CD45", "JAM", "PECAM1", "TENASCIN", "SEMA3"),
     pathways.selected) %>%
    rev()
)
dev.off()

## incoming
pdf("Global_Communication_Patterns/Elbow.Incoming.pdf", 6, 2.5)
selectK(cellchat, pattern = "incoming")
dev.off()

nPatterns <- 6

pdf("Global_Communication_Patterns/Heatmap.Pattern.In.pdf", 8, 4)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
dev.off()

pdf("Global_Communication_Patterns/River.Pattern.In.pdf", 8, 4)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.In.pdf", 4.5, 15)
netAnalysis_dot(cellchat, pattern = "incoming") + coord_flip()
dev.off()

pdf("Global_Communication_Patterns/Dot.Pattern.In.Selected.pdf", 7, 3)
netAnalysis_dot(cellchat,
    pattern = "incoming",
    pathway.show = union(
        c(
            "COLLAGEN", "LAMININ", "APP", "MHC-I", "FN1", "THBS", "CCL", "SPP1", "ADGRE", "MHC-II",
            "CypA", "CD80", "PTN", "CD86", "MIF", "GALECTIN", "CD45", "JAM", "PECAM1", "TENASCIN", "SEMA3"
        ),
        pathways.selected
    ) %>%
        rev()
)
dev.off()
rm(nPatterns)

## 2.4. pathway level analysis ----
dir.create("Pathway_Level")
cellchat@netP$pathways %>%
    grep("GF", ., value = TRUE)

table(pathways.selected %in% cellchat@netP$pathways)

lapply(pathways.selected, function(pathway) {
    pdf(paste0("Pathway_Level/Circle.Network.Pathway.", pathway, ".pdf"), 6, 7)
    netVisual_aggregate(cellchat,
        signaling = pathway, layout = "circle",
    )
    dev.off()

    pdf(paste0("Pathway_Level/Heatmap.Pathway.", pathway, ".pdf"), 6, 5.5)
    p <- netVisual_heatmap(cellchat,
        signaling = pathway,
        color.heatmap = "Reds"
    )
    plot(p)
    dev.off()

    pdf(paste0("Pathway_Level/Chord.Pathway.", pathway, ".pdf"), 8, 9)
    netVisual_aggregate(cellchat,
        signaling = pathway,
        layout = "chord"
    )
    dev.off()

    pdf(paste0("Pathway_Level/Heatmap_signalingRole.", pathway, ".pdf"), 8, 3)
    netAnalysis_signalingRole_network(cellchat,
        signaling = pathway,
         width = 8, height = 2.5, font.size = 10
    )
    dev.off()
})
rm(p)

## 2.5. cell type level analysis ----
dir.create("Cell_Type_Level")
cell_types

lapply(cell_types, function(cell_type) {
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
    pdf(paste0("Cell_Type_Level/Bubble.", cell_type, "_Outgoing.pdf"), 7, height)
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
    pdf(paste0("Cell_Type_Level/Bubble.", cell_type, "_Incoming.pdf"), 7, height)
    plot(p)
    dev.off()
})
rm(p, cell_type)

# 3. Compare signal across conditions ----
blca.obj <- readRDS("blca.obj.rds")
dir.create("Compare_Conditions")

## 3.1. prepare CellChat objects ----
blca.sub <- subset(blca.obj, Cell_Type_2 %in% cell_types & Condition %in% c("CTRL"))
Idents(blca.sub) <- blca.sub$Cell_Type_2
cellchat.ctrl <- createCellChat(object = blca.sub, group.by = "ident", assay = "RNA")
cellchat.ctrl@DB <- CellChatDB.mouse

blca.sub <- subset(blca.obj, Cell_Type_2 %in% cell_types & Condition %in% c("Capsaicin"))
Idents(blca.sub) <- blca.sub$Cell_Type_2
cellchat.capsaicin <- createCellChat(object = blca.sub, group.by = "ident", assay = "RNA")
cellchat.capsaicin@DB <- CellChatDB.mouse
rm(blca.sub)
rm(blca.obj)

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
    aggregateNet()

object.list <- list(
    "CTRL" = cellchat.ctrl,
    "Capsaicin" = cellchat.capsaicin
)

cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(object.list, file = "CellChat.object.list.rds")
rm(cellchat.ctrl, cellchat.capsaicin)

## 3.2. compare overall network ----
p1 <- compareInteractions(cellchat.merge, show.legend = FALSE, group = c(1, 2))
p2 <- compareInteractions(cellchat.merge, show.legend = FALSE, group = c(1, 2), measure = "weight")
pdf("Compare_Conditions/Bar.Compare_Interactions.pdf", 5, 2.5)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()
rm(p1, p2)

pdf("Compare_Conditions/Circle.Network.Compare.Count.pdf", 6, 7)
netVisual_diffInteraction(cellchat.merge, weight.scale = TRUE)
dev.off()

pdf("Compare_Conditions/Circle.Network.Compare.Weight.pdf", 6, 7)
netVisual_diffInteraction(cellchat.merge, weight.scale = TRUE, measure = "weight")
dev.off()

pdf("Compare_Conditions/Circle.Network.Compare.Weight.CAF.pdf", 6, 7)
netVisual_diffInteraction(cellchat.merge,
    weight.scale = TRUE, measure = "weight",
    sources.use = c("iCAF", "myCAF"),
)
dev.off()

pdf("Compare_Conditions/Heatmap.Compare.Weight.pdf", 6, 5.5)
netVisual_heatmap(cellchat.merge, measure = "weight")
dev.off()

num.link <- sapply(object.list, function(x) {
    rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
p.list <- list()
for (i in seq_along(object.list)) {
    p.list[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
        title = names(object.list)[i], weight.MinMax = weight.MinMax
    )
}

pdf("Compare_Conditions/Scatter.In_out.Compare.pdf", 9, 4)
patchwork::wrap_plots(plots = p.list)
dev.off()

## 3.3. identify altered pathways ----
p1 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "iCAF")
p2 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "Cancer")
p3 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "T_CD8")
p4 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "T_CD4")
pdf("Compare_Conditions/Scatter.Signaling_Changes.iCAF.T_CD8.pdf", 11, 6.5)
patchwork::wrap_plots(plots = list(p1, p2, p3, p4), ncol = 2)
dev.off()
rm(p1, p2, p3, p4, p.list, num.link, weight.MinMax, i)

p1 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL,
    stacked = TRUE, do.stat = TRUE
)
p2 <- rankNet(cellchat.merge,
    mode = "comparison", measure = "weight", sources.use = NULL,
    targets.use = NULL, stacked = FALSE, do.stat = TRUE
)
pdf("Compare_Conditions/Bar.Ranked.Pathways.Compare.pdf", 7.5, 15)
patchwork::wrap_plots(p1, p2, ncol = 2)
dev.off()
rm(p1, p2)

dir.create("Compare_Conditions/Pathway_Level")

lapply(pathways.selected, function(pathways.show) {
    weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)

    pdf(paste0("Compare_Conditions/Pathway_Level/Circle.Network.Pathway.", pathways.show, ".pdf"), 12, 7)
    par(mfrow = c(1, 2), xpd = TRUE)
    for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]],
            signaling = pathways.show,
            layout = "circle", edge.weight.max = weight.max[1],
            edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i])
        )
    }
    dev.off()

    par(mfrow = c(1, 2), xpd = TRUE)
    ht <- list()
    for (i in seq_along(object.list)) {
        ht[[i]] <- netVisual_heatmap(object.list[[i]],
            signaling = pathways.show, color.heatmap = "Reds",
            title.name = paste(pathways.show, "signaling ", names(object.list)[i])
        )
    }

    pdf(paste0("Compare_Conditions/Pathway_Level/Heatmap.Pathway.", pathways.show, ".pdf"), 10, 5.5)
    ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
    dev.off()

    pdf(paste0("Compare_Conditions/Pathway_Level/Chord.Pathway.", pathways.show, ".pdf"), 14, 9)
    par(mfrow = c(1, 2), xpd = TRUE)
    for (i in seq_along(object.list)) {
        netVisual_aggregate(object.list[[i]],
            signaling = pathways.show,
            layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i])
        )
    }
    dev.off()
})
rm(pathways.show, ht, i, weight.max)

## 3.4. identify altered interactions ----
dir.create("Compare_Conditions/Cell_Type_Level")

# identify DEG
pos.dataset <- "Capsaicin"
features.name <- paste0(pos.dataset, ".merged")

cellchat.merge <- identifyOverExpressedGenes(cellchat.merge,
    group.dataset = "datasets", pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 0.05, group.DE.combined = FALSE
)
net <- netMappingDEG(cellchat.merge, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat.merge,
    net = net,
    datasets = "Capsaicin", ligand.logFC = 0.05, receptor.logFC = NULL
)
net.down <- subsetCommunication(cellchat.merge,
    net = net,
    datasets = "CTRL", ligand.logFC = -0.05, receptor.logFC = NULL
)

write.table(net.up,
    file = "Compare_Conditions/DEG.Up.LR.Pairs.Capsaicin.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(net.down,
    file = "Compare_Conditions/DEG.Down.LR.Pairs.Capsaicin.txt",
    sep = "\t", quote = FALSE, row.names = FALSE
)
rm(features.name, net, pos.dataset)

# iCAF outgoing
cell_type.selected <- c("T_CD4", "T_CD8", "NK", "Neutrophil", "Cancer", "iCAF")
pairLR.use.up <- net.up %>%
    filter(source == "iCAF" & target %in% cell_type.selected) %>%
    select(interaction_name)
p <- netVisual_bubble(cellchat.merge,
    sources.use = "iCAF",
    targets.use = cell_type.selected,
    pairLR.use = pairLR.use.up,
    remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
    comparison = c(1, 2), max.dataset = 2, angle.x = 45,
    title.name = "Increased iCAF outgoing signaling in Capsaicin",
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
height <- length(unique(p$data$interaction_name)) / 6.5 + 2
pdf(paste0("Compare_Conditions/Cell_Type_Level/Bubble.iCAF_Increased_Outgoing.pdf"), 6, height)
plot(p)
dev.off()

pairLR.use.up <- net.up %>%
    filter(source %in% cell_type.selected & target == "iCAF") %>%
    select(interaction_name)
p <- netVisual_bubble(cellchat.merge,
    sources.use = cell_type.selected,
    targets.use = "iCAF",
    pairLR.use = pairLR.use.up,
    remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
    comparison = c(1, 2), max.dataset = 2, angle.x = 45,
    title.name = "Increased iCAF incoming signaling in Capsaicin",
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
height <- length(unique(p$data$interaction_name)) / 6.5 + 2
pdf(paste0("Compare_Conditions/Cell_Type_Level/Bubble.iCAF_Increased_Incoming.pdf"), 6, height)
plot(p)
dev.off()

pairLR.use.down <- net.down %>%
    filter(source == "iCAF" & target %in% cell_type.selected) %>%
    select(interaction_name)
p <- netVisual_bubble(cellchat.merge,
    sources.use = "iCAF",
    targets.use = cell_type.selected,
    pairLR.use = pairLR.use.down,
    remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
    comparison = c(1, 2), max.dataset = 1, angle.x = 45,
    title.name = "Decreased iCAF outgoing signaling in Capsaicin",
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
height <- length(unique(p$data$interaction_name)) / 6.5 + 2
pdf(paste0("Compare_Conditions/Cell_Type_Level/Bubble.iCAF_Decreased_Outgoing.pdf"), 6, height)
plot(p)
dev.off()

# T_CD8 incoming
pairLR.use.up <- net.up %>%
    filter(source %in% cell_type.selected & target == "T_CD8") %>%
    select(interaction_name)
p <- netVisual_bubble(cellchat.merge,
    sources.use = cell_type.selected,
    targets.use = "T_CD8",
    pairLR.use = pairLR.use.up,
    remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
    comparison = c(1, 2), max.dataset = 2, angle.x = 45,
    title.name = "Increased T_CD8 incoming signaling in Capsaicin",
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
height <- length(unique(p$data$interaction_name)) / 6.5 + 2
pdf(paste0("Compare_Conditions/Cell_Type_Level/Bubble.T_CD8_Increased_Incoming.pdf"), 6, height)
plot(p)
dev.off()

pairLR.use.down <- net.down %>%
    filter(source %in% cell_type.selected & target == "T_CD8") %>%
    select(interaction_name)
p <- netVisual_bubble(cellchat.merge,
    sources.use = cell_type.selected,
    targets.use = "T_CD8",
    pairLR.use = pairLR.use.down,
    remove.isolate = TRUE, thresh = 0.05, grid.on = FALSE,
    comparison = c(1, 2), max.dataset = 1, angle.x = 45,
    title.name = "Decreased T_CD8 incoming signaling in Capsaicin",
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
height <- length(unique(p$data$interaction_name)) / 6.5 + 2
pdf(paste0("Compare_Conditions/Cell_Type_Level/Bubble.T_CD8_Decreased_Incoming.pdf"), 6, height)
plot(p)
dev.off()

rm(p, pairLR.use.up, pairLR.use.down, cell_type.selected, height)

# save data ----
gc()
saveRDS(blca.obj, file = "blca.obj.rds")
rm(blca.obj, object.list)
save.image("04.Communication.RData")
