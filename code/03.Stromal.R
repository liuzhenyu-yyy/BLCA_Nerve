## 5.4. Differential analysis ----
dir.create("sub3_Stromal/DEG")

unique(sample.info$Condition)
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

lapply(names(DEG.stromal.Capsaicin), function(x) {
    write.csv(DEG.stromal.Capsaicin[[x]],
        file = paste0("sub3_Stromal/DEG/DEG_Capsaicin_vs_CTRL_in_", x, ".csv"),
        row.names = FALSE
    )
})
lapply(names(DEG.stromal.LC_6OHDA), function(x) {
    write.csv(DEG.stromal.LC_6OHDA[[x]],
        file = paste0("sub3_Stromal/DEG/DEG_LC_6OHDA_vs_CTRL_in_", x, ".csv"),
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
pdf("sub3_Stromal/DEG/Volcano.Capsaicin_vs_CTRL.pdf", 7.5, 4.5)
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

pdf("sub3_Stromal/DEG/Volcano.LC_6OHDA_vs_CTRL.pdf", 7.5, 4.5)
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