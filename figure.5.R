# Please change the path here before you run the script.
setwd('./')

library(dplyr)
library(ggplot2)
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(ggpubr)
library(msigdbr)
library(enrichplot)
library(aplot)
library(GSVA)
library(AUCell)
library(monocle)
library(ClusterGVis)
library(stringr)
library(tidyr)
library(rlang)

fig_path <- '../results/'
# dir.create(fig_path, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 0. Load Data
# ------------------------------------------------------------------------------
cat("Loading pre-computed RData...\n")
load('../data/PDAC_Drug_Prediction_Plotting.RData')

# ------------------------------------------------------------------------------
# 1. Global UMAP & Sample Barplots
# ------------------------------------------------------------------------------
cat("Generating Global UMAP and Sample Barplots...\n")

png(paste0(fig_path, 'pdac_treatment_umap.png'), units = 'in', res = 300, width = 6, height = 5)
DimPlot(obj.int, group.by = 'treatment', pt.size = 1) + 
  theme(text = element_text(size = 24)) + 
  labs(title = NULL, x = "UMAP1", y = "UMAP2") + 
  scale_color_manual(values = c(my_colors[1], my_colors[3], my_colors[4]))
dev.off()

png(paste0(fig_path, 'pal_dmso_umap.png'), units = 'in', res = 300, width = 6, height = 5)
DimPlot(pal_dmso, group.by = 'treatment', pt.size = 1) + 
  theme(text = element_text(size = 24)) + 
  labs(title = NULL, x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c(my_colors[1], my_colors[3]))
dev.off()

png(paste0(fig_path, 'tram_dmso_umap.png'), units = 'in', res = 300, width = 6, height = 5)
DimPlot(tram_dmso, group.by = 'treatment', pt.size = 1) + 
  theme(text = element_text(size = 24)) + 
  labs(title = NULL, x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c(my_colors[1], my_colors[4]))
dev.off()

df_bar <- obj.int@meta.data %>%
  as.data.frame() %>%
  group_by(sample, treatment) %>%
  summarise(count = n(), .groups = "drop")

png(paste0(fig_path, 'pdac_sample_barplot.png'), units = 'in', res = 300, width = 5, height = 5, bg = "transparent")
ggplot(df_bar, aes(x = sample, y = count, fill = treatment)) +
  geom_bar(stat = "identity", color = 'black') +
  labs(x = NULL, y = "Cell counts", fill = NULL) +
  theme_classic() +
  theme(
    text = element_text(size = 24),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) + 
  scale_fill_manual(values = c(my_colors[1], my_colors[3], my_colors[4]))
dev.off()

# ------------------------------------------------------------------------------
# 2. Model Performance Plots
# ------------------------------------------------------------------------------
cat("Generating Model Performance Barplots...\n")

plot_single_drug <- function(data, target_drug, my_colors) {
  tool_order <- data %>%
    filter(drug == target_drug) %>%
    group_by(Tool) %>%
    summarize(mean_auc = mean(AUC, na.rm = TRUE)) %>%
    arrange(desc(mean_auc)) %>%
    pull(Tool)
  
  plot_df <- data %>%
    filter(drug == target_drug) %>%
    mutate(Tool = factor(Tool, levels = tool_order)) %>%
    pivot_longer(
      cols = c(AUC, ACC, F1, Recall, Precision), 
      names_to = "Metric", 
      values_to = "Value"
    ) %>%
    group_by(Metric) %>%
    arrange(Metric, desc(Value)) %>%
    mutate(rank = row_number()) %>% 
    ungroup()
  
  n_tools <- length(unique(plot_df$Tool))
  gap <- 2 
  plot_df <- plot_df %>%
    mutate(
      metric_idx = as.numeric(factor(Metric, levels = c("AUC", "ACC", "F1", "Recall", "Precision"))),
      x_pos = (metric_idx - 1) * (n_tools + gap) + rank
    )
  
  axis_labels <- plot_df %>%
    group_by(Metric) %>%
    summarize(center_pos = mean(x_pos))
  
  ggplot(plot_df, aes(x = x_pos, y = Value, fill = Tool)) +
    geom_col(width = 0.9) + 
    scale_x_continuous(breaks = axis_labels$center_pos, labels = axis_labels$Metric) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1), add = c(0, 0))) +
    scale_fill_manual(values = my_colors) +
    theme_bw() +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme(
      text = element_text(size = 28),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right",
      axis.text.x = element_text(margin = margin(t = 20, b = 0), vjust = 0.5, hjust = 0.5),
      legend.text = element_text(margin = margin(l = 5, r = 15))
    )
}

png(paste0(fig_path, 'pal_perform_all.png'), res = 300, units = 'in', width = 12, height = 4)
plot_single_drug(markus_all_eval_df, "paclitaxel", my_colors_by_tools)
dev.off()

png(paste0(fig_path, 'tram_perform_all.png'), res = 300, units = 'in', width = 12, height = 4)
plot_single_drug(markus_all_eval_df, "trametinib", my_colors_by_tools)
dev.off()

# ------------------------------------------------------------------------------
# 3. scDEAL Predictions (UMAP, Density, Violin)
# ------------------------------------------------------------------------------
cat("Generating scDEAL Prediction Plots...\n")

# -- Paclitaxel --
pal_cells_ordered <- rownames(pal_dmso@meta.data)[order(pal_dmso$scDEAL_pred, decreasing = TRUE)]
pal_min_val <- min(pal_dmso$scDEAL_pred, na.rm = TRUE)
pal_max_val <- max(pal_dmso$scDEAL_pred, na.rm = TRUE)

png(paste0(fig_path, 'pal_scDEAL_umap.png'), units = 'in', res = 300, width = 6, height = 6)
FeaturePlot(pal_dmso, features = 'scDEAL_pred', cells = pal_cells_ordered, cols = c(my_colors[3], my_colors[1]), pt.size = 1) + 
  theme(text = element_text(size = 24),
        axis.text = element_text(size = 20, color = "black"),
        legend.position = 'top',
        legend.justification = "center",
        legend.box.just = "center",
        legend.box = "horizontal",
        legend.text = element_text(margin = margin(t = 5))) + 
  labs(title = NULL, x = "UMAP1", y = "UMAP2") + 
  scale_color_gradientn(colors = c(my_colors[3], my_colors[1]), breaks = c(pal_min_val, pal_max_val), labels = c('Resistant', 'Sensitive')) +
  guides(color = guide_colorbar(barwidth = unit(20, "lines"), barheight = unit(1, "lines"), title = NULL, ticks = FALSE, label.position = "bottom", label.hjust = c(0, 1))) + 
  annotate("text", x = Inf, y = Inf, label = 'Paclitaxel', hjust = 1.1, vjust = 1.5, size = 10, fontface = "bold", family = "sans")
dev.off()

# -- Trametinib --
tram_cells_ordered <- rownames(tram_dmso@meta.data)[order(tram_dmso$scDEAL_pred, decreasing = TRUE)]
tram_min_val <- min(tram_dmso$scDEAL_pred, na.rm = TRUE)
tram_max_val <- max(tram_dmso$scDEAL_pred, na.rm = TRUE)

png(paste0(fig_path, 'tram_scDEAL_umap.png'), units = 'in', res = 300, width = 6, height = 7)
FeaturePlot(tram_dmso, features = 'scDEAL_pred', cells = tram_cells_ordered, cols = c(my_colors[3], my_colors[1]), pt.size = 1) + 
  theme(text = element_text(size = 24),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 22),
        legend.position = 'top',
        legend.justification = "center",
        legend.box.just = "center",
        legend.box = "horizontal",
        legend.text = element_text(margin = margin(t = 5))) + 
  labs(title = NULL, x = "UMAP1", y = "UMAP2") + 
  scale_color_gradientn(colors = c(my_colors[4], my_colors[1]), breaks = c(tram_min_val, tram_max_val), labels = c('Resistant', 'Sensitive')) +
  guides(color = guide_colorbar(barwidth = unit(20, "lines"), barheight = unit(1, "lines"), title = NULL, ticks = FALSE, label.position = "bottom", label.hjust = c(0, 1))) + 
  annotate("text", x = Inf, y = Inf, label = 'Trametinib', hjust = 1.1, vjust = 1.5, size = 10, fontface = "bold", family = "sans")
dev.off()

# -- Density Plots --
png(paste0(fig_path, 'pal_scDEAL_pred_density.png'), units = 'in', res = 300, width = 6, height = 7)
ggplot(pal_dmso@meta.data, aes(x = scDEAL_pred, color = condition)) +
  geom_density(linewidth = 2.2, key_glyph = "path") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_color_manual(values = c(my_colors[1], my_colors[3]), breaks = c("sensitive", "resistant"), labels = c("Non-treated", "Treated")) + 
  theme_classic() +
  theme(text = element_text(size = 24), axis.text = element_text(size = 20, color = "black"), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 28), legend.key = element_blank()) +
  labs(title = NULL, x = NULL, y = "Density")
dev.off()

png(paste0(fig_path, 'tram_scDEAL_pred_density.png'), units = 'in', res = 300, width = 6, height = 7)
ggplot(tram_dmso@meta.data, aes(x = scDEAL_pred, color = condition)) +
  geom_density(linewidth = 2.2, key_glyph = "path") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  scale_color_manual(values = c(my_colors[1], my_colors[4]), breaks = c("sensitive", "resistant"), labels = c("Non-treated", "Treated")) + 
  theme_classic() +
  theme(text = element_text(size = 24), axis.text = element_text(size = 20, color = "black"), legend.position = "top", legend.title = element_blank(), legend.text = element_text(size = 28), legend.key = element_blank()) +
  labs(title = NULL, x = NULL, y = "Density")
dev.off()

# -- Violin Plots --
png(paste0(fig_path, 'pal_scDEAL_pred_violin.png'), units = 'in', res = 300, width = 4, height = 5)
ggplot(pal_dmso@meta.data, aes(x = condition, y = scDEAL_pred, fill = condition)) +
  geom_violin(trim = TRUE, alpha = 1, color = NA) + 
  geom_boxplot(width = 0.03, fill = "white", color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c(my_colors[3], my_colors[1])) +
  theme_classic() +
  labs(x = NULL, y = "Sensitive probability") +
  annotate("segment", x = 1, xend = 2, y = 1.05, yend = 1.05, colour = "black", size = 0.8) +
  annotate("text", x = 1.5, y = 1.1, label = "p < 2.2 %*% 10^-16", parse = TRUE, size = 7) +
  annotate("segment", x = 1, xend = 1, y = 1.05, yend = 1.03, colour = "black", size = 0.8) +
  annotate("segment", x = 2, xend = 2, y = 1.05, yend = 1.03, colour = "black", size = 0.8) +
  ylim(NA, 1.1) + 
  theme(legend.position = "none", text = element_text(size = 28), axis.title.y = element_text(size = 20)) 
dev.off()

png(paste0(fig_path, 'tram_scDEAL_pred_violin.png'), units = 'in', res = 300, width = 4, height = 5)
ggplot(tram_dmso@meta.data, aes(x = condition, y = scDEAL_pred, fill = condition)) +
  geom_violin(trim = TRUE, alpha = 1, color = NA) + 
  geom_boxplot(width = 0.03, fill = "white", color = "black", outlier.shape = NA) + 
  scale_fill_manual(values = c(my_colors[4], my_colors[1])) +
  theme_classic() +
  labs(x = NULL, y = "Sensitive probability") +
  annotate("segment", x = 1, xend = 2, y = 1.05, yend = 1.05, colour = "black", size = 0.8) +
  annotate("text", x = 1.5, y = 1.1, label = "p < 2.2 %*% 10^-16", parse = TRUE, size = 7) +
  annotate("segment", x = 1, xend = 1, y = 1.05, yend = 1.03, colour = "black", size = 0.8) +
  annotate("segment", x = 2, xend = 2, y = 1.05, yend = 1.03, colour = "black", size = 0.8) +
  ylim(NA, 1.1) + 
  theme(legend.position = "none", text = element_text(size = 28), axis.title.y = element_text(size = 20)) 
dev.off()


# ------------------------------------------------------------------------------
# 4. Bioinformatics (Gene Expr & ssGSEA Violin Plots)
# ------------------------------------------------------------------------------
cat("Generating Bioinformatics Violin Plots...\n")

Idents(pal_dmso) <- "scDEAL_pred_cond"
Idents(tram_dmso) <- "scDEAL_pred_cond"

# -- Target Genes Expression --
png(paste0(fig_path, 'tram_tar_MAP2K1.png'), units = 'in', res = 300, height = 6, width = 5)
VlnPlot(tram_dmso, features = 'MAP2K1') + 
  scale_fill_manual(values = c('skyblue','tomato')) +
  annotate("segment", x = 1, xend = 2, y = 2.3, yend = 2.3, colour = "black") +
  annotate("text", x = 1.5, y = 2.4, label = "p = 0.035", size = 5) + 
  ylim(0, 2.5) + ggtitle("MAP2K1") + 
  theme(text = element_text(size = 18), axis.text = element_text(size = 18), axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

png(paste0(fig_path, 'tram_tar_MAP2K2.png'), units = 'in', res = 300, height = 6, width = 5)
VlnPlot(tram_dmso, features = 'MAP2K2') + 
  scale_fill_manual(values = c('skyblue','tomato')) +
  annotate("segment", x = 1, xend = 2, y = 2.8, yend = 2.8, colour = "black") +
  annotate("text", x = 1.5, y = 2.9, label = "p = 0.0015", size = 5) + 
  ylim(0, 3) + ggtitle("MAP2K2") + 
  theme(text = element_text(size = 18), axis.text = element_text(size = 18), axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

png(paste0(fig_path, 'pal_tar_TUBB1.png'), units = 'in', res = 300, height = 6, width = 5)
VlnPlot(pal_dmso, features = 'TUBB1') + 
  scale_fill_manual(values = c('skyblue','tomato')) +
  annotate("segment", x = 1, xend = 2, y = 1.3, yend = 1.3, colour = "black") +
  annotate("text", x = 1.5, y = 1.4, label = "NS", size = 5) + 
  ylim(0, 1.5) + ggtitle("TUBB1") + 
  theme(text = element_text(size = 18), axis.text = element_text(size = 18), axis.title.x = element_blank(), axis.text.x = element_text(angle = 0, hjust = 0.5)) + NoLegend()
dev.off()

# -- ssGSEA Enrichment Scores --
pal_plot_df <- pal_dmso@meta.data
tram_plot_df <- tram_dmso@meta.data

png(paste0(fig_path, 'pal_ssGSEA_apoptosis.png'), units = 'in', res = 300, height = 6, width = 5)
ggplot(pal_plot_df, aes(x = scDEAL_pred_cond, y = HALLMARK_APOPTOSIS, fill = scDEAL_pred_cond)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c('tomato', 'skyblue')) + theme_classic() +
  labs(title = "APOPTOSIS (ssGSEA)", y = "Enrichment Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_blank(), axis.text = element_text(size = 18), text = element_text(size = 18), legend.position = "none") + 
  annotate("segment", x = 1, xend = 2, y = 0.4, yend = 0.4, colour = "black") +
  annotate("text", x = 1.5, y = 0.42, label = "p = 0.071", size = 5) + ylim(NA, 0.45)
dev.off()

png(paste0(fig_path, 'pal_ssGSEA_g2m.png'), units = 'in', res = 300, height = 6, width = 5)
ggplot(pal_plot_df, aes(x = scDEAL_pred_cond, y = HALLMARK_G2M_CHECKPOINT, fill = scDEAL_pred_cond)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c('tomato', 'skyblue')) + theme_classic() +
  labs(title = "G2M CHECKPOINT (ssGSEA)", y = "Enrichment Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_blank(), axis.text = element_text(size = 18), text = element_text(size = 18), legend.position = "none") + 
  annotate("segment", x = 1, xend = 2, y = 0.55, yend = 0.55, colour = "black") +
  annotate("text", x = 1.5, y = 0.58, label = "p == 2.2 %*% 10^-16", parse = TRUE, size = 5) + ylim(NA, 0.6)
dev.off()

png(paste0(fig_path, 'pal_ssGSEA_Microtubule_Stability.png'), units = 'in', res = 300, height = 6, width = 5)
ggplot(pal_plot_df, aes(x = scDEAL_pred_cond, y = Microtubule_Stability_UCell, fill = scDEAL_pred_cond)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c('tomato', 'skyblue')) + theme_classic() +
  labs(title = "Microtubule Stability", y = "Enrichment Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), text = element_text(size = 18), axis.title.x = element_blank(), axis.text = element_text(size = 18), legend.position = "none") + 
  annotate("segment", x = 1, xend = 2, y = 0.65, yend = 0.65, colour = "black") +
  annotate("text", x = 1.5, y = 0.68, label = "p == 7.8 %*% 10^-12", parse = TRUE, size = 5) + ylim(NA, 0.7)
dev.off()

png(paste0(fig_path, 'pal_ssGSEA_Paclitaxel_Resistance_UCell.png'), units = 'in', res = 300, height = 6, width = 5)
ggplot(subset(pal_plot_df, Paclitaxel_Resistance_UCell <= 0.05), aes(x = scDEAL_pred_cond, y = Paclitaxel_Resistance_UCell, fill = scDEAL_pred_cond)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c('tomato', 'skyblue')) + theme_classic() +
  labs(title = "Paclitaxel Resistance", y = "Enrichment Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_blank(), text = element_text(size = 18), axis.text = element_text(size = 18), legend.position = "none") + 
  annotate("segment", x = 1, xend = 2, y = 0.055, yend = 0.055, colour = "black") +
  annotate("text", x = 1.5, y = 0.058, label = "p = 0.018", size = 5) + ylim(NA, 0.06)
dev.off()

png(paste0(fig_path, 'tram_ssGSEA_KRAS_UP.png'), units = 'in', res = 300, height = 6, width = 5)
ggplot(tram_plot_df, aes(x = scDEAL_pred_cond, y = HALLMARK_KRAS_SIGNALING_UP, fill = scDEAL_pred_cond)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c('tomato', 'skyblue')) + theme_classic() +
  labs(title = "KRAS SIGNALING (ssGSEA)", y = "Enrichment Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_blank(), axis.text = element_text(size = 18), text = element_text(size = 18), legend.position = "none") + 
  annotate("segment", x = 1, xend = 2, y = 0.23, yend = 0.23, colour = "black") +
  annotate("text", x = 1.5, y = 0.25, label = "p == 2.2 %*% 10^-16", parse = TRUE, size = 5) + ylim(NA, 0.26)
dev.off()

png(paste0(fig_path, 'tram_ssGSEA_EMT.png'), units = 'in', res = 300, height = 6, width = 5)
ggplot(tram_plot_df, aes(x = scDEAL_pred_cond, y = HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, fill = scDEAL_pred_cond)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c('tomato', 'skyblue')) + theme_classic() +
  labs(title = "EMT (ssGSEA)", y = "Enrichment Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_blank(), axis.text = element_text(size = 18), text = element_text(size = 18), legend.position = "none") + 
  annotate("segment", x = 1, xend = 2, y = 0.27, yend = 0.27, colour = "black") +
  annotate("text", x = 1.5, y = 0.29, label = "p == 2.2 %*% 10^-16", parse = TRUE, size = 5) + ylim(NA, 0.3)
dev.off()

png(paste0(fig_path, 'tram_ssGSEA_MAPK_Feedback_Targets_UCell.png'), units = 'in', res = 300, height = 6, width = 5)
ggplot(tram_plot_df, aes(x = scDEAL_pred_cond, y = MAPK_Feedback_Targets_UCell, fill = scDEAL_pred_cond)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c('tomato', 'skyblue')) + theme_classic() +
  labs(title = "MAPK Feedback Targets", y = "Enrichment Score") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title.x = element_blank(), axis.text = element_text(size = 18), text = element_text(size = 18), legend.position = "none") + 
  annotate("segment", x = 1, xend = 2, y = 0.37, yend = 0.37, colour = "black") +
  annotate("text", x = 1.5, y = 0.39, label = "p == 8.6 %*% 10^-7", parse = TRUE, size = 5) + ylim(NA, 0.4)
dev.off()


# ------------------------------------------------------------------------------
# 5. Monocle2 Trajectories Visualization
# ------------------------------------------------------------------------------
cat("Generating Monocle2 Trajectory Plots...\n")

# =========================================================================
# Modernized plot_cell_trajectory
# (Uses explicit namespaces to resolve dplyr deprecation and masking issues)
# =========================================================================
plot_cell_trajectory_modern <- function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, 
                                         show_backbone = TRUE, backbone_color = "black", markers = NULL, 
                                         use_color_gradient = FALSE, markers_linear = FALSE, show_cell_names = FALSE, 
                                         show_state_number = FALSE, cell_size = 1.5, cell_link_size = 0.75, 
                                         cell_name_size = 2, state_number_size = 2.9, show_branch_points = TRUE, 
                                         theta = 0, ...) 
{
  requireNamespace("igraph")
  gene_short_name <- NA
  sample_name <- NA
  sample_state <- pData(cds)$State
  data_dim_1 <- NA
  data_dim_2 <- NA
  lib_info_with_pseudo <- pData(cds)
  
  if (is.null(cds@dim_reduce_type)) {
    stop("Error: dimensionality not yet reduced. Please call reduceDimension() before calling this function.")
  }
  
  if (cds@dim_reduce_type == "ICA") {
    reduced_dim_coords <- reducedDimS(cds)
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree")) {
    reduced_dim_coords <- reducedDimK(cds)
  } else {
    stop("Error: unrecognized dimensionality reduction method.")
  }
  
  # Explicitly use dplyr:: and tibble:: to prevent namespace conflicts
  ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
    dplyr::select(prin_graph_dim_1 = dplyr::all_of(x), prin_graph_dim_2 = dplyr::all_of(y)) %>% 
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
  
  dp_mst <- minSpanningTree(cds)
  if (is.null(dp_mst)) {
    stop("You must first call orderCells() before using this function")
  }
  
  edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
    dplyr::select(source = from, target = to) %>% 
    dplyr::left_join(ica_space_df %>% dplyr::select(source = sample_name, 
                                                    source_prin_graph_dim_1 = prin_graph_dim_1, 
                                                    source_prin_graph_dim_2 = prin_graph_dim_2), 
                     by = "source") %>% 
    dplyr::left_join(ica_space_df %>% dplyr::select(target = sample_name, 
                                                    target_prin_graph_dim_1 = prin_graph_dim_1, 
                                                    target_prin_graph_dim_2 = prin_graph_dim_2), 
                     by = "target")
  
  data_df <- t(monocle::reducedDimS(cds)) %>% as.data.frame() %>% 
    dplyr::select(data_dim_1 = dplyr::all_of(x), data_dim_2 = dplyr::all_of(y)) %>% 
    tibble::rownames_to_column("sample_name") %>% 
    dplyr::mutate(sample_state) %>% 
    dplyr::left_join(lib_info_with_pseudo %>% tibble::rownames_to_column("sample_name"), by = "sample_name")
  
  return_rotation_mat <- function(theta) {
    theta <- theta/180 * pi
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  
  rot_mat <- return_rotation_mat(theta)
  cn1 <- c("data_dim_1", "data_dim_2")
  cn2 <- c("source_prin_graph_dim_1", "source_prin_graph_dim_2")
  cn3 <- c("target_prin_graph_dim_1", "target_prin_graph_dim_2")
  data_df[, cn1] <- as.matrix(data_df[, cn1]) %*% t(rot_mat)
  edge_df[, cn2] <- as.matrix(edge_df[, cn2]) %*% t(rot_mat)
  edge_df[, cn3] <- as.matrix(edge_df[, cn3]) %*% t(rot_mat)
  markers_exprs <- NULL
  
  if (is.null(markers) == FALSE) {
    markers_fData <- subset(fData(cds), gene_short_name %in% markers)
    if (nrow(markers_fData) >= 1) {
      markers_exprs <- reshape2::melt(as.matrix(exprs(cds[row.names(markers_fData), ])))
      colnames(markers_exprs)[1:2] <- c("feature_id", "cell_id")
      markers_exprs <- merge(markers_exprs, markers_fData, by.x = "feature_id", by.y = "row.names")
      markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
      markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
    }
  }
  
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
    data_df <- merge(data_df, markers_exprs, by.x = "sample_name", by.y = "cell_id")
    if (use_color_gradient) {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
          geom_point(aes(color = value), size = I(cell_size), na.rm = TRUE) + 
          scale_color_viridis(name = paste0("value"), ...) + facet_wrap(~feature_label)
      } else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2)) + 
          geom_point(aes(color = log10(value + 0.1)), size = I(cell_size), na.rm = TRUE) + 
          scale_color_viridis(name = paste0("log10(value + 0.1)"), ...) + facet_wrap(~feature_label)
      }
    } else {
      if (markers_linear) {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, size = (value * 0.1))) + 
          facet_wrap(~feature_label)
      } else {
        g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2, size = log10(value + 0.1))) + 
          facet_wrap(~feature_label)
      }
    }
  } else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
  }
  
  if (show_tree) {
    g <- g + geom_segment(aes(x = source_prin_graph_dim_1, 
                              y = source_prin_graph_dim_2, 
                              xend = target_prin_graph_dim_1, 
                              yend = target_prin_graph_dim_2), 
                          size = cell_link_size, linetype = "solid", na.rm = TRUE, data = edge_df)
  }
  
  if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 0) {
    if (use_color_gradient) {
    } else {
      g <- g + geom_point(aes(color = .data[[color_by]]), na.rm = TRUE)
    }
  } else {
    if (use_color_gradient) {
    } else {
      g <- g + geom_point(aes(color = .data[[color_by]]), size = I(cell_size), na.rm = TRUE)
    }
  }
  
  if (show_branch_points && cds@dim_reduce_type == "DDRTree") {
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    branch_point_df <- ica_space_df %>% 
      dplyr::slice(match(mst_branch_nodes, sample_name)) %>% 
      dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
    
    g <- g + geom_point(aes(x = prin_graph_dim_1, y = prin_graph_dim_2), 
                        size = 5, na.rm = TRUE, data = branch_point_df) + 
      geom_text(aes(x = prin_graph_dim_1, y = prin_graph_dim_2, label = branch_point_idx), 
                size = 4, color = "white", na.rm = TRUE, data = branch_point_df)
  }
  
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_state_number) {
    g <- g + geom_text(aes(label = sample_state), size = state_number_size)
  }
  
  # Fallback to standard ggplot2 theme to avoid missing private monocle functions
  g <- g + theme_classic() + 
    xlab(paste("Component", x)) + ylab(paste("Component", y)) + 
    theme(legend.position = "top", legend.key.height = grid::unit(0.35, "in")) + 
    theme(legend.key = element_blank()) + 
    theme(panel.background = element_rect(fill = "white"))
  
  g
}
# =========================================================================

# -- Paclitaxel Trajectory --
png(paste0(fig_path, 'pal_monocle_state.png'), units = 'in', res = 300, width = 6, height = 6)
plot_cell_trajectory_modern(pal_cds, color_by = "State", cell_size = 3) + 
  scale_color_manual(values = my_colors[1:5]) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 28))
dev.off()

png(paste0(fig_path, 'pal_monocle_pseudotime.png'), units = 'in', res = 300, width = 6, height = 6)
plot_cell_trajectory_modern(pal_cds, color_by = "Pseudotime", cell_size = 3) + 
  scale_color_gradient(high = my_colors[3], low = my_colors[2], breaks = c(0, 5, 10)) + 
  theme(text = element_text(size = 22))
dev.off()

png(paste0(fig_path, 'pal_monocle_scdeal_prob.png'), units = 'in', res = 300, width = 6, height = 6)
plot_cell_trajectory_modern(pal_cds, color_by = "scDEAL_pred", cell_size = 3) + 
  scale_color_gradient(high = 'tomato', low = 'skyblue', breaks = c(0, 0.5, 1), limits = c(0, 1), name = 'Sensitivity prob') + 
  theme(text = element_text(size = 22))
dev.off()

# -- Trametinib Trajectory --
png(paste0(fig_path, 'tram_monocle_state.png'), units = 'in', res = 300, width = 6, height = 6)
plot_cell_trajectory_modern(tram_cds_subset, color_by = "State", cell_size = 3) + 
  scale_color_manual(values = my_colors[1:6]) + 
  theme(text = element_text(size = 22), legend.text = element_text(size = 28))
dev.off()

png(paste0(fig_path, 'tram_monocle_pseudotime.png'), units = 'in', res = 300, width = 6, height = 6)
plot_cell_trajectory_modern(tram_cds_subset, color_by = "Pseudotime", cell_size = 3) + 
  scale_color_gradient(high = my_colors[3], low = my_colors[2], breaks = c(0, 5, 10)) + 
  theme(text = element_text(size = 22))
dev.off()

png(paste0(fig_path, 'tram_monocle_scdeal_prob.png'), units = 'in', res = 300, width = 6, height = 6)
plot_cell_trajectory_modern(tram_cds_subset, color_by = "scDEAL_pred", cell_size = 3) + 
  scale_color_gradient(high = 'tomato', low = 'skyblue', breaks = c(0, 0.5, 1), limits = c(0, 1), name = 'Sensitivity prob') + 
  theme(text = element_text(size = 22))
dev.off()

# ------------------------------------------------------------------------------
# 6. ClusterGVis Heatmaps
# ------------------------------------------------------------------------------
cat("Generating ClusterGVis Heatmaps...\n")

text_colors <- c("blue4", "red4", "darkgreen", "purple4", "darkorange3", "grey30")

# ==============================================================================
# IMPORTANT: PASTE YOUR FULL `my_visCluster` FUNCTION HERE
# Replace the placeholder function below with your complete custom function.
# ==============================================================================
my_visCluster<-function (object = NULL, htColList = list(col_range = c(-2, 
                                                                       0, 2), col_color = c("#08519C", "white", "#A50F15")), border = TRUE, 
                         plotType = c("line", "heatmap", "both"), msCol = c("#0099CC", 
                                                                            "grey90", "#CC3333"), lineSize = 0.1, lineCol = "grey90", 
                         addMline = TRUE, mlineSize = 2, mlineCol = "#CC3333", ncol = 4, 
                         ctAnnoCol = NULL, setMd = "median", textboxPos = c(0.5, 
                                                                            0.8), textboxSize = 8, panelArg = c(2, 0.25, 4, "grey90", 
                                                                                                                NA), ggplotPanelArg = c(2, 0.25, 4, "grey90", NA), annoTermData = NULL, 
                         annoTermMside = "right", termAnnoArg = c("grey95", "grey50"), 
                         addBar = FALSE, barWidth = 8, textbarPos = c(0.8, 0.8), 
                         goCol = NULL, goSize = NULL, byGo = "anno_link", annoKeggData = NULL, 
                         annoKeggMside = "right", keggAnnoArg = c("grey95", "grey50"), 
                         addKeggBar = FALSE, keggCol = NULL, keggSize = NULL, byKegg = "anno_link", 
                         wordWrap = TRUE, addNewLine = TRUE, addBox = FALSE, boxCol = NULL, 
                         boxArg = c(0.1, "grey50"), addPoint = FALSE, pointArg = c(19, 
                                                                                   "orange", "orange", 1), addLine = TRUE, lineSide = "right", 
                         markGenes = NULL, markGenesSide = "right", genesGp = c("italic", 
                                                                                10, NA), termTextLimit = c(10, 18), mulGroup = NULL, 
                         lgdLabel = NULL, showRowNames = FALSE, subgroupAnno = NULL, 
                         annnoblockText = TRUE, annnoblockGp = c("white", 8), addSampleAnno = TRUE, 
                         sampleGroup = NULL, sampleCol = NULL, sampleOrder = NULL, 
                         clusterOrder = NULL, sampleCellOrder = NULL, heatmapAnnotation = NULL, 
                         columnSplit = NULL, clusterColumns = FALSE, pseudotimeCol = NULL, 
                         gglist = NULL, rowAnnotationObj = NULL, ...) 
{
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. Please install it.")
  }
  ComplexHeatmap::ht_opt(message = FALSE)
  if (is.null(htColList[["col_range"]])) {
    col_range <- c(-2, 0, 2)
  }
  else {
    col_range <- htColList[["col_range"]]
  }
  if (is.null(htColList[["col_color"]])) {
    col_color <- c("#08519C", "white", "#A50F15")
  }
  else {
    col_color <- htColList[["col_color"]]
  }
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' is required. Please install it.")
  }
  col_fun <- circlize::colorRamp2(col_range, col_color)
  plotType <- match.arg(plotType)
  if (plotType == "line") {
    if (object$type %in% c("scRNAdata", "monocle", "wgcna")) {
      data <- dplyr::arrange(data.frame(object$long.res), 
                             .data[["cluster"]])
    }
    else {
      if (object$type %in% c("mfuzz", "TCseq")) {
        data <- dplyr::arrange(data.frame(object$long.res), 
                               .data[["cluster"]], .data[["membership"]])
      }
      else {
        data <- dplyr::arrange(data.frame(object$long.res), 
                               .data[["cluster"]])
      }
    }
    data$gene <- factor(data$gene, levels = unique(data$gene))
    if (!is.null(sampleOrder)) {
      data$cell_type <- factor(data$cell_type, levels = sampleOrder)
    }
    line <- ggplot2::ggplot(data, ggplot2::aes(x = cell_type, 
                                               y = norm_value))
    if (object$type %in% c("mfuzz", "TCseq")) {
      line <- line + ggplot2::geom_line(ggplot2::aes(color = membership, 
                                                     group = gene), linewidth = lineSize) + ggplot2::scale_color_gradient2(low = msCol[1], 
                                                                                                                           mid = msCol[2], high = msCol[3], midpoint = 0.5)
    }
    else {
      line <- line + ggplot2::geom_line(ggplot2::aes(group = gene), 
                                        color = lineCol, linewidth = lineSize)
    }
    if (addMline == TRUE) {
      if (object$type == "wgcna") {
        linec <- unique(data$modulecol)
        names(linec) <- linec
        line <- line + ggplot2::geom_line(stat = "summary", 
                                          fun = "median", linewidth = mlineSize, ggplot2::aes(group = 1, 
                                                                                              color = modulecol)) + ggplot2::scale_color_manual(values = linec)
      }
      else {
        line <- line + ggplot2::geom_line(stat = "summary", 
                                          fun = "median", colour = mlineCol, linewidth = mlineSize, 
                                          ggplot2::aes(group = 1))
      }
    }
    else {
      line <- line
    }
    line1 <- line + ggplot2::theme_classic(base_size = 14) + 
      ggplot2::ylab("Normalized expression") + ggplot2::xlab("") + 
      ggplot2::theme(axis.ticks.length = ggplot2::unit(0.1, 
                                                       "cm"), axis.text.x = ggplot2::element_text(angle = 45, 
                                                                                                  hjust = 1, color = "black"), strip.background = ggplot2::element_blank()) + 
      ggplot2::facet_wrap(~cluster_name, ncol = ncol, 
                          scales = "free")
    return(line1)
  }
  else {
    data <- dplyr::arrange(data.frame(object$wide.res, check.names = FALSE), 
                           as.numeric(as.character(cluster)))
    if (object$type %in% c("mfuzz", "TCseq")) {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster, -membership)
    }
    else if (object$type == "wgcna") {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster, -modulecol)
    }
    else if (object$type == "scRNAdata") {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster)
    }
    else if (object$type == "monocle") {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster)
    }
    else {
      mat <- dplyr::select(dplyr::arrange(data, as.numeric(as.character(cluster))), 
                           -gene, -cluster)
    }
    rownames(mat) <- data$gene
    if (object$geneMode == "all" | ncol(mat) > 20) {
      use_raster <- TRUE
    }
    else {
      use_raster <- FALSE
    }
    if (!is.null(sampleOrder)) {
      mat <- mat[, sampleOrder]
    }
    cl.info <- dplyr::arrange(dplyr::mutate(data.frame(table(data$cluster)), 
                                            Var1 = as.numeric(as.character(Var1))), Var1)
    cluster.num <- nrow(cl.info)
    subgroup <- unlist(lapply(seq_len(nrow(cl.info)), function(x) {
      nm <- rep(as.character(cl.info$Var1[x]), cl.info$Freq[x])
      paste("C", nm, sep = "")
    }))
    if (!is.null(clusterOrder)) {
      subgroup <- factor(subgroup, levels = paste("C", 
                                                  clusterOrder, sep = ""))
      cluster_row_slices <- FALSE
    }
    else {
      cluster_row_slices <- TRUE
    }
    if (object$geneMode == "all" & object$type == "scRNAdata") {
      celltype <- vapply(strsplit(colnames(mat), split = "\\|"), 
                         function(x) {
                           x[2]
                         }, character(1))
      cell.num.info <- table(celltype)[unique(celltype)]
      if (is.null(sampleCellOrder)) {
        column_split <- factor(rep(names(cell.num.info), 
                                   cell.num.info), levels = unique(celltype))
      }
      else {
        column_split <- factor(rep(names(cell.num.info), 
                                   cell.num.info), levels = sampleCellOrder)
      }
      if (is.null(sampleCol)) {
        block.col <- seq_len(length(cell.num.info))
      }
      else {
        block.col <- sampleCol
      }
    }
    else {
      if (is.null(sampleGroup)) {
        sample.info <- colnames(mat)
        if (is.null(heatmapAnnotation)) {
          if (object$geneType == "branched") {
            if (ncol(mat) == 200) {
              column_split <- rep(c("branch1", "branch2"), 
                                  each = 100)
            }
            else {
              column_split <- rep(levels(object$pseudotime), 
                                  rep(100, ncol(mat)/100))
            }
          }
          else {
            column_split <- NULL
          }
        }
        else {
          column_split <- columnSplit
        }
      }
      else {
        sample.info <- sampleGroup
        column_split <- factor(sampleGroup, levels = unique(sampleGroup))
      }
      if (object$type != "monocle") {
        sample.info <- factor(sample.info, levels = unique(sample.info))
        if (is.null(sampleCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }
          scol <- circlize::rand_color(n = length(sample.info))
          names(scol) <- sample.info
        }
        else {
          scol <- sampleCol
          names(scol) <- sample.info
        }
      }
      else {
        sample.info <- factor(object$pseudotime, levels = unique(object$pseudotime))
        if (is.null(pseudotimeCol)) {
          if (object$geneType == "branched") {
            if (length(unique(object$pseudotime)) == 
                3) {
              pseudotimeCol <- c("red", "grey80", "blue")
            }
            else {
              if (!requireNamespace("circlize", quietly = TRUE)) {
                stop("Package 'circlize' is required. Please install it.")
              }
              pseudotimeCol <- circlize::rand_color(n = length(unique(object$pseudotime)))
            }
          }
          else {
            pseudotimeCol <- c("blue", "red")
          }
        }
        else {
          pseudotimeCol <- pseudotimeCol
        }
        if (is.null(sampleCol)) {
          if (object$type != "monocle") {
            if (!requireNamespace("circlize", quietly = TRUE)) {
              stop("Package 'circlize' is required. Please install it.")
            }
            scol <- circlize::rand_color(n = length(sample.info))
            names(scol) <- sample.info
          }
          else {
            if (object$geneType == "branched") {
              if (length(unique(object$pseudotime)) == 
                  3) {
                scol <- rep(pseudotimeCol, table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
              else {
                scol <- rep(pseudotimeCol, table(object$pseudotime)[unique(object$pseudotime)])
                names(scol) <- sample.info
              }
            }
            else {
              scol <- (grDevices::colorRampPalette(pseudotimeCol))(100)
              names(scol) <- sample.info
            }
          }
        }
        else {
          scol <- sampleCol
          names(scol) <- sample.info
        }
      }
    }
    if (addSampleAnno == TRUE) {
      if (object$geneMode == "all" & object$type == "scRNAdata") {
        topanno <- ComplexHeatmap::HeatmapAnnotation(cluster = ComplexHeatmap::anno_block(gp = grid::gpar(fill = block.col), 
                                                                                          labels = NULL), show_annotation_name = FALSE)
      }
      else {
        if (is.null(heatmapAnnotation)) {
          topanno <- ComplexHeatmap::HeatmapAnnotation(sample = sample.info, 
                                                       col = list(sample = scol), gp = grid::gpar(col = ifelse(object$type == 
                                                                                                                 "monocle", NA, "white")), show_legend = ifelse(object$type == 
                                                                                                                                                                  "monocle", FALSE, TRUE), show_annotation_name = FALSE)
        }
        else {
          topanno <- heatmapAnnotation
        }
      }
    }
    else {
      topanno <- NULL
    }
    if (is.null(ctAnnoCol)) {
      if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required. Please install it.")
      }
      colanno <- circlize::rand_color(n = cluster.num)
    }
    else {
      colanno <- ctAnnoCol
    }
    names(colanno) <- seq_len(cluster.num)
    align_to <- split(seq_len(nrow(mat)), subgroup)
    anno.block <- ComplexHeatmap::anno_block(align_to = align_to, 
                                             panel_fun = function(index, nm) {
                                               npos <- as.numeric(unlist(strsplit(nm, split = "C"))[2])
                                               grid::grid.rect(gp = grid::gpar(fill = colanno[npos], 
                                                                               col = NA))
                                               if (annnoblockText == TRUE) {
                                                 grid::grid.text(label = paste("Num: ", length(index), 
                                                                               sep = ""), rot = 90, gp = grid::gpar(col = annnoblockGp[1], 
                                                                                                                    fontsize = as.numeric(annnoblockGp[2])))
                                               }
                                             }, which = "row")
    if (!is.null(markGenes)) {
      rowGene <- rownames(mat)
      annoGene <- markGenes
      gene.col <- dplyr::filter(dplyr::select(data, gene, 
                                              cluster), gene %in% annoGene)
      gene.col <- purrr::map_df(seq_len(cluster.num), 
                                function(x) {
                                  tmp <- dplyr::mutate(dplyr::filter(gene.col, 
                                                                     as.numeric(cluster) == x), col = colanno[x])
                                })
      gene.col <- gene.col[match(annoGene, gene.col$gene), 
      ]
      if (is.na(genesGp[3])) {
        gcol <- gene.col$col
      }
      else {
        gcol <- genesGp[3]
      }
      index <- match(annoGene, rowGene)
      geneMark <- ComplexHeatmap::anno_mark(at = index, 
                                            labels = annoGene, which = "row", side = markGenesSide, 
                                            labels_gp = grid::gpar(fontface = genesGp[1], 
                                                                   fontsize = as.numeric(genesGp[2]), col = gcol))
    }
    else {
      geneMark <- NULL
    }
    right_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                      cluster = anno.block)
    if (object$type == "monocle" | object$geneMode == "all" | 
        ncol(mat) > 20) {
      show_column_names <- FALSE
    }
    else {
      show_column_names <- TRUE
    }
    if (object$geneType == "non-branched") {
      rg <- range(as.numeric(as.character(sample.info)))
      if (!requireNamespace("circlize", quietly = TRUE)) {
        stop("Package 'circlize' is required. Please install it.")
      }
      col_fun2 <- circlize::colorRamp2(c(rg[1], rg[2]), 
                                       pseudotimeCol)
      lgd <- ComplexHeatmap::Legend(col_fun = col_fun2, 
                                    title = "pseudotime")
      lgd_list <- list(lgd)
    }
    else if (object$geneType == "branched") {
      if (length(levels(sample.info)) == 3) {
        lgd <- ComplexHeatmap::Legend(labels = levels(sample.info), 
                                      legend_gp = grid::gpar(fill = pseudotimeCol), 
                                      title = "branch")
      }
      else {
        lgd <- ComplexHeatmap::Legend(labels = levels(sample.info), 
                                      legend_gp = grid::gpar(fill = pseudotimeCol), 
                                      title = "branch")
      }
      lgd_list <- list(lgd)
    }
    else {
      lgd_list <- NULL
    }
    if (plotType == "heatmap") {
      if (!is.null(rowAnnotationObj)) {
        left_annotation_ht <- rowAnnotationObj
      }
      else {
        left_annotation_ht <- NULL
      }
      htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                                     cluster_columns = clusterColumns, show_row_names = showRowNames, 
                                     border = border, column_split = column_split, 
                                     row_split = subgroup, cluster_row_slices = cluster_row_slices, 
                                     column_names_side = "top", show_column_names = show_column_names, 
                                     top_annotation = topanno, left_annotation = left_annotation_ht, 
                                     right_annotation = right_annotation, col = col_fun, 
                                     use_raster = use_raster, ...)
      ComplexHeatmap::draw(htf, merge_legend = TRUE, annotation_legend_list = lgd_list)
    }
    else {
      rg <- range(mat)
      if (!is.null(gglist)) {
        anno_ggplot2 <- ComplexHeatmap::anno_zoom(align_to = align_to, 
                                                  which = "row", panel_fun = function(index, 
                                                                                      nm) {
                                                    g <- gglist[[nm]]
                                                    g <- grid::grid.grabExpr(grid::grid.draw(g))
                                                    grid::pushViewport(grid::viewport())
                                                    grid::grid.rect()
                                                    grid::grid.draw(g)
                                                    grid::popViewport()
                                                  }, size = grid::unit(as.numeric(ggplotPanelArg[1]), 
                                                                       "cm"), gap = grid::unit(as.numeric(ggplotPanelArg[2]), 
                                                                                               "cm"), width = grid::unit(as.numeric(ggplotPanelArg[3]), 
                                                                                                                         "cm"), side = "right", link_gp = grid::gpar(fill = ggplotPanelArg[4], 
                                                                                                                                                                     col = ggplotPanelArg[5]))
      }
      else {
        anno_ggplot2 <- NULL
      }
      panel_fun <- function(index, nm) {
        if (addBox == TRUE & addLine != TRUE) {
          xscale <- c(-0.1, 1.1)
        }
        else {
          xscale <- c(-0.1, 1.1)
          panel_scale <- c(0.1, 0.9)
        }
        grid::pushViewport(grid::viewport(xscale = xscale, 
                                          yscale = c(0, 1)))
        grid::grid.rect()
        if (object$geneMode == "all" & object$type == 
            "scRNAdata") {
          mulGroup <- cell.num.info
          grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                                                2), gp = grid::gpar(col = "black", lty = "dashed"))
          cu <- cumsum(mulGroup)
          seqn <- data.frame(st = c(1, cu[seq(1, (length(cu) - 
                                                    1), 1)] + 1), sp = c(cu[1], cu[seq(2, length(cu), 
                                                                                       1)]))
        }
        else {
          if (is.null(mulGroup)) {
            mulGroup <- ncol(mat)
            seqn <- data.frame(st = 1, sp = ncol(mat))
          }
          else {
            mulGroup <- mulGroup
            grid::grid.lines(x = c(0, 1), y = rep(0.5, 
                                                  2), gp = grid::gpar(col = "black", lty = "dashed"))
            cu <- cumsum(mulGroup)
            seqn <- data.frame(st = c(1, cu[seq(1, length(cu) - 
                                                  1, 1)] + 1), sp = c(cu[1], cu[seq(2, length(cu), 
                                                                                    1)]))
          }
        }
        if (object$geneMode == "all" && object$type == 
            "scRNAdata") {
          cell.ave <- purrr::map_dfr(seq_len(nrow(seqn)), 
                                     function(x) {
                                       tmp <- seqn[x, ]
                                       tmpmat <- mat[index, seq(tmp$st, tmp$sp, 
                                                                1)]
                                       rg <- base::range(mat[index, ])
                                       if (setMd == "mean") {
                                         mdia <- base::mean(base::rowMeans(tmpmat))
                                       }
                                       else if (setMd == "median") {
                                         mdia <- stats::median(base::apply(tmpmat, 
                                                                           1, stats::median))
                                       }
                                       else {
                                         message("supply mean/median !")
                                       }
                                       res <- data.frame(x = x, val = mdia)
                                       return(res)
                                     })
          if (addLine == TRUE) {
            grid::grid.lines(x = scales::rescale(cell.ave$x, 
                                                 to = c(0, 1)), y = scales::rescale(cell.ave$val, 
                                                                                    to = c(0.1, 0.9)), gp = grid::gpar(lwd = 3, 
                                                                                                                       col = mlineCol))
          }
        }
        else {
          base::lapply(seq_len(nrow(seqn)), function(x) {
            tmp <- seqn[x, ]
            tmpmat <- mat[index, seq(tmp$st, tmp$sp, 
                                     1)]
            if (setMd == "mean") {
              mdia <- base::colMeans(tmpmat)
            }
            else if (setMd == "median") {
              mdia <- base::apply(tmpmat, 2, stats::median)
            }
            else {
              message("supply mean/median !")
            }
            pos <- scales::rescale(seq_len(ncol(tmpmat)), 
                                   to = c(0, 1))
            if (is.null(boxCol)) {
              boxCol <- rep("grey90", ncol(tmpmat))
            }
            else {
              boxCol <- boxCol
            }
            if (addBox == TRUE) {
              lapply(seq_len(ncol(tmpmat)), function(x) {
                ComplexHeatmap::grid.boxplot(scales::rescale(tmpmat[, 
                                                                    x], to = c(0, 1), from = c(rg[1] - 
                                                                                                 0.5, rg[2] + 0.5)), pos = pos[x], 
                                             direction = "vertical", box_width = as.numeric(boxArg[1]), 
                                             outline = FALSE, gp = grid::gpar(col = boxArg[2], 
                                                                              fill = boxCol[x]))
              })
            }
            if (addPoint == TRUE) {
              grid::grid.points(x = scales::rescale(seq_len(ncol(tmpmat)), 
                                                    to = c(0.1, 0.9)), y = scales::rescale(mdia, 
                                                                                           to = c(0, 1), from = c(rg[1] - 0.5, 
                                                                                                                  rg[2] + 0.5)), pch = as.numeric(pointArg[1]), 
                                gp = grid::gpar(fill = pointArg[2], 
                                                col = pointArg[3]), size = grid::unit(as.numeric(pointArg[4]), 
                                                                                      "char"))
            }
            if (addLine == TRUE) {
              grid::grid.lines(x = scales::rescale(seq_len(ncol(tmpmat)), 
                                                   to = c(0.1, 0.9)), y = scales::rescale(mdia, 
                                                                                          to = c(0, 1), from = c(rg[1] - 0.5, 
                                                                                                                 rg[2] + 0.5)), gp = grid::gpar(lwd = 3, 
                                                                                                                                                col = mlineCol[x]))
            }
          })
        }
        text <- paste("Gene size:", nrow(mat[index, ]), sep = " ")
        grid::grid.text(text, x = textboxPos[1], y = textboxPos[2], 
                        gp = grid::gpar(fontsize = textboxSize, fontface = "italic"))
        grid::popViewport()
      }
      if (!is.null(subgroupAnno)) {
        align_to <- split(seq_len(nrow(mat)), subgroup)
        align_to <- align_to[subgroupAnno]
      }
      else {
        align_to <- subgroup
      }
      anno <- ComplexHeatmap::anno_link(align_to = align_to, 
                                        which = "row", panel_fun = panel_fun, size = grid::unit(as.numeric(panelArg[1]), 
                                                                                                "cm"), gap = grid::unit(as.numeric(panelArg[2]), 
                                                                                                                        "cm"), width = grid::unit(as.numeric(panelArg[3]), 
                                                                                                                                                  "cm"), side = lineSide, link_gp = grid::gpar(fill = panelArg[4], 
                                                                                                                                                                                               col = panelArg[5]))
      if (!is.null(annoTermData)) {
        termanno <- annoTermData
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        }
        else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        }
        else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", 
                                  "ratio")
        }
        else {
          message("No more than 4 columns!")
        }
        if (is.null(goCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }
          gocol <- circlize::rand_color(n = nrow(termanno))
        }
        else {
          gocol <- goCol
        }
        if (is.null(goSize)) {
          gosize <- rep(12, nrow(termanno))
        }
        else {
          if (goSize == "pval") {
            termanno.tmp <- purrr::map_df(unique(termanno$id), 
                                          function(x) {
                                            tmp <- dplyr::mutate(dplyr::filter(termanno, 
                                                                               id == x), size = scales::rescale(-log10(pval), 
                                                                                                                to = termTextLimit))
                                          })
            gosize <- termanno.tmp$size
          }
          else {
            gosize <- goSize
          }
        }
        termanno <- dplyr::mutate(dplyr::ungroup(termanno), 
                                  col = gocol, fontsize = gosize)
        term.list <- lapply(seq_len(length(unique(termanno$id))), 
                            function(x) {
                              tmp <- termanno[which(termanno$id == unique(termanno$id)[x]), 
                              ]
                              df <- data.frame(text = tmp$term, col = tmp$col, 
                                               fontsize = tmp$fontsize)
                              return(df)
                            })
        names(term.list) <- unique(termanno$id)
        if (!is.null(subgroupAnno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroupAnno]
          term.list <- term.list[subgroupAnno]
        }
        else {
          align_to2 <- subgroup
          term.list <- term.list
        }
        textbox <- ComplexHeatmap::anno_textbox(align_to2, 
                                                term.list, word_wrap = wordWrap, add_new_line = addNewLine, 
                                                side = annoTermMside, background_gp = grid::gpar(fill = termAnnoArg[1], 
                                                                                                 col = termAnnoArg[2]), by = byGo)
        if (ncol(termanno) - 2 > 2) {
          anno_gobar <- function(data = NULL, barWidth = 0.1, 
                                 align_to = NULL, panelArg = panelArg, ...) {
            if (ncol(data) - 2 == 3) {
              data <- dplyr::mutate(data, bary = -log10(pval))
            }
            else {
              data <- dplyr::mutate(data, bary = ratio)
            }
            ComplexHeatmap::anno_zoom(align_to = align_to, 
                                      which = "row", panel_fun = function(index, 
                                                                          nm) {
                                        grid::pushViewport(grid::viewport(xscale = c(0, 
                                                                                     1), yscale = c(0, 1)))
                                        grid::grid.rect()
                                        tmp <- dplyr::filter(data, id == nm)
                                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                                                            x1 = scales::rescale(rev(tmp$bary), 
                                                                                 to = c(0.1, 0.9)), y0 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                         to = c(0.1, 0.9)), y1 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                                                                 to = c(0.1, 0.9)), gp = grid::gpar(lwd = barWidth, 
                                                                                                                                                                                                    col = rev(tmp$col), lineend = "butt"))
                                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                                                                                "ComplexHeatmap")
                                        text <- nm
                                        grid.textbox(text, x = textbarPos[1], 
                                                     y = textbarPos[2], gp = grid::gpar(fontsize = textboxSize, 
                                                                                        fontface = "italic", col = unique(tmp$col), 
                                                                                        ...))
                                        grid::popViewport()
                                      }, size = grid::unit(as.numeric(panelArg[1]), 
                                                           "cm"), gap = grid::unit(as.numeric(panelArg[2]), 
                                                                                   "cm"), width = grid::unit(as.numeric(panelArg[3]), 
                                                                                                             "cm"), side = "right", link_gp = grid::gpar(fill = termAnnoArg[1], 
                                                                                                                                                         col = termAnnoArg[2]), ...)
          }
          baranno <- anno_gobar(data = termanno, align_to = align_to2, 
                                panelArg = panelArg, barWidth = barWidth)
        }
        if (addBar == TRUE) {
          baranno
        }
        else {
          baranno <- NULL
        }
      }
      else {
        textbox <- NULL
        baranno <- NULL
      }
      if (!is.null(annoKeggData)) {
        termanno <- annoKeggData
        if (ncol(termanno) == 2) {
          colnames(termanno) <- c("id", "term")
        }
        else if (ncol(termanno) == 3) {
          colnames(termanno) <- c("id", "term", "pval")
        }
        else if (ncol(termanno) == 4) {
          colnames(termanno) <- c("id", "term", "pval", 
                                  "ratio")
        }
        else {
          message("No more than 4 columns!")
        }
        if (is.null(keggCol)) {
          if (!requireNamespace("circlize", quietly = TRUE)) {
            stop("Package 'circlize' is required. Please install it.")
          }
          gocol <- circlize::rand_color(n = nrow(termanno))
        }
        else {
          gocol <- keggCol
        }
        if (is.null(keggSize)) {
          gosize <- rep(12, nrow(termanno))
        }
        else {
          if (keggSize == "pval") {
            termanno.tmp <- purrr::map_df(unique(termanno$id), 
                                          function(x) {
                                            tmp <- dplyr::mutate(dplyr::filter(termanno, 
                                                                               id == x), size = scales::rescale(-log10(pval), 
                                                                                                                to = termTextLimit))
                                          })
            gosize <- termanno.tmp$size
          }
          else {
            gosize <- keggSize
          }
        }
        termanno <- dplyr::mutate(dplyr::ungroup(termanno), 
                                  col = gocol, fontsize = gosize)
        term.list <- lapply(seq_len(length(unique(termanno$id))), 
                            function(x) {
                              tmp <- termanno[which(termanno$id == unique(termanno$id)[x]), 
                              ]
                              df <- data.frame(text = tmp$term, col = tmp$col, 
                                               fontsize = tmp$fontsize)
                              return(df)
                            })
        names(term.list) <- unique(termanno$id)
        if (!is.null(subgroupAnno)) {
          align_to2 <- split(seq_along(subgroup), subgroup)
          align_to2 <- align_to2[subgroupAnno]
          term.list <- term.list[subgroupAnno]
        }
        else {
          align_to2 <- subgroup
          term.list <- term.list
        }
        textbox.kegg <- ComplexHeatmap::anno_textbox(align_to2, 
                                                     term.list, word_wrap = wordWrap, add_new_line = addNewLine, 
                                                     side = annoKeggMside, background_gp = grid::gpar(fill = keggAnnoArg[1], 
                                                                                                      col = keggAnnoArg[2]), by = byKegg)
        if (ncol(termanno) - 2 > 2) {
          anno_keggbar <- function(data = NULL, barWidth = 0.1, 
                                   align_to = NULL, panelArg = panelArg, ...) {
            if (ncol(data) - 2 == 3) {
              data <- dplyr::mutate(data, bary = -log10(pval))
            }
            else {
              data <- dplyr::mutate(data, bary = ratio)
            }
            ComplexHeatmap::anno_zoom(align_to = align_to, 
                                      which = "row", panel_fun = function(index, 
                                                                          nm) {
                                        grid::pushViewport(grid::viewport(xscale = c(0, 
                                                                                     1), yscale = c(0, 1)))
                                        grid::grid.rect()
                                        tmp <- dplyr::filter(data, id == nm)
                                        grid::grid.segments(x0 = rep(0, nrow(tmp)), 
                                                            x1 = scales::rescale(rev(tmp$bary), 
                                                                                 to = c(0.1, 0.9)), y0 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                         to = c(0.1, 0.9)), y1 = scales::rescale(seq_len(nrow(tmp)), 
                                                                                                                                                                 to = c(0.1, 0.9)), gp = grid::gpar(lwd = barWidth, 
                                                                                                                                                                                                    col = rev(tmp$col), lineend = "butt"))
                                        grid.textbox <- utils::getFromNamespace("grid.textbox", 
                                                                                "ComplexHeatmap")
                                        text <- nm
                                        grid.textbox(text, x = textbarPos[1], 
                                                     y = textbarPos[2], gp = grid::gpar(fontsize = textboxSize, 
                                                                                        fontface = "italic", col = unique(tmp$col), 
                                                                                        ...))
                                        grid::popViewport()
                                      }, size = grid::unit(as.numeric(panelArg[1]), 
                                                           "cm"), gap = grid::unit(as.numeric(panelArg[2]), 
                                                                                   "cm"), width = grid::unit(as.numeric(panelArg[3]), 
                                                                                                             "cm"), side = "right", link_gp = grid::gpar(fill = keggAnnoArg[1], 
                                                                                                                                                         col = keggAnnoArg[2]), ...)
          }
          baranno.kegg <- anno_keggbar(data = termanno, 
                                       align_to = align_to2, panelArg = panelArg, 
                                       barWidth = barWidth)
        }
        if (addKeggBar == TRUE) {
          baranno.kegg
        }
        else {
          baranno.kegg <- NULL
        }
      }
      else {
        textbox.kegg <- NULL
        baranno.kegg <- NULL
      }
      if (lineSide == "right") {
        if (markGenesSide == "right") {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                             cluster = anno.block, line = anno, anno_ggplot2 = anno_ggplot2, 
                                                             textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- rowAnnotationObj
          }
          else {
            left_annotation <- NULL
          }
        }
        else {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                                                             line = anno, anno_ggplot2 = anno_ggplot2, 
                                                             textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                                       modifyList(list(gene = geneMark), rowAnnotationObj))
          }
          else {
            left_annotation <- ComplexHeatmap::rowAnnotation(gene = geneMark)
          }
        }
      }
      else {
        if (markGenesSide == "right") {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(gene = geneMark, 
                                                             cluster = anno.block, anno_ggplot2 = anno_ggplot2, 
                                                             textbox = textbox, bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                                       modifyList(list(line = anno), rowAnnotationObj))
          }
          else {
            left_annotation <- ComplexHeatmap::rowAnnotation(line = anno)
          }
        }
        else {
          right_annotation2 <- ComplexHeatmap::rowAnnotation(cluster = anno.block, 
                                                             anno_ggplot2 = anno_ggplot2, textbox = textbox, 
                                                             bar = baranno, textbox.kegg = textbox.kegg, 
                                                             baranno.kegg = baranno.kegg)
          if (!is.null(rowAnnotationObj)) {
            left_annotation <- do.call(ComplexHeatmap::rowAnnotation, 
                                       modifyList(list(gene = geneMark, line = anno), 
                                                  rowAnnotationObj))
          }
          else {
            left_annotation <- ComplexHeatmap::rowAnnotation(line = anno, 
                                                             gene = geneMark)
          }
        }
      }
      if (object$type == "monocle" | object$geneMode == 
          "all" | ncol(mat) > 20) {
        show_column_names <- FALSE
      }
      else {
        show_column_names <- TRUE
      }
      htf <- ComplexHeatmap::Heatmap(as.matrix(mat), name = "Z-score", 
                                     cluster_columns = clusterColumns, show_row_names = showRowNames, 
                                     border = border, column_split = column_split, 
                                     top_annotation = topanno, right_annotation = right_annotation2, 
                                     left_annotation = left_annotation, column_names_side = "top", 
                                     show_column_names = show_column_names, row_split = subgroup, 
                                     cluster_row_slices = cluster_row_slices, col = col_fun, 
                                     use_raster = use_raster, ...)
      if (is.null(mulGroup)) {
        ComplexHeatmap::draw(htf, merge_legend = TRUE, 
                             annotation_legend_list = lgd_list)
      }
      else {
        if (is.null(lgdLabel)) {
          lgdLabel <- paste("group", seq_len(length(mulGroup)), 
                            sep = "")
        }
        else {
          lgdLabel <- lgdLabel
        }
        lgd_list2 <- ComplexHeatmap::Legend(labels = lgdLabel, 
                                            type = "lines", legend_gp = grid::gpar(col = mlineCol, 
                                                                                   lty = 1))
        if (!is.null(lgd_list)) {
          lgd_list_com <- ComplexHeatmap::packLegend(lgd_list, 
                                                     lgd_list2)
        }
        else {
          lgd_list_com <- lgd_list2
        }
        ComplexHeatmap::draw(htf, annotation_legend_list = lgd_list_com, 
                             merge_legend = TRUE)
      }
    }
  }
}

# -- Paclitaxel ClusterGVis --
pal_heatmap_res <- plot_pseudotime_heatmap2(pal_cds[pal_sig_genes_htmap,], num_clusters = 3, cores = 8)
pal_clusterGVis_gene <- sample(pal_heatmap_res$wide.res$gene, 20, replace = FALSE)

pal_clusterGVis_enrich <- enrichCluster(object = pal_heatmap_res, OrgDb = org.Hs.eg.db, type = "BP", organism = "hsa", pvalueCutoff = 0.05, topn = 10)

pal_c3_filtered <- pal_clusterGVis_enrich[pal_clusterGVis_enrich$group == "C3" & !grepl("digestion|digestive", pal_clusterGVis_enrich$Description, ignore.case = TRUE), ]
pal_c3_filtered <- pal_c3_filtered[-2, ] 

pal_c1_df <- pal_clusterGVis_enrich[pal_clusterGVis_enrich$group == "C1", ]
pal_c2_df <- pal_clusterGVis_enrich[pal_clusterGVis_enrich$group == "C2", ]

pal_clusterGVis_enrich <- bind_rows(
  pal_c1_df %>% head(5),
  pal_c2_df %>% head(5),
  pal_c3_filtered %>% head(5)
)

png(paste0(fig_path, 'pal_monocle_clusterGVis_compact_text.png'), units = 'in', res = 300, width = 12, height = 8)
my_visCluster(object = pal_heatmap_res, plotType = "both", column_names_rot = 45, show_row_dend = FALSE,
              markGenes = pal_clusterGVis_gene, markGenesSide = "left", pseudotimeCol = c(my_colors[2], my_colors[3]),
              annoTermData = pal_clusterGVis_enrich, ctAnnoCol = text_colors[1:3], goCol = rep(text_colors[1:3], each = 5),
              genesGp = c("italic", 15, NA), goSize = 16, addBar = FALSE, textboxSize = 12, textboxPos = c(0.5, 0.88), lineSide = "left")
dev.off()

# -- Trametinib ClusterGVis --
tram_heatmap_res <- plot_pseudotime_heatmap2(tram_cds_subset[tram_sig_genes_htmap,], num_clusters = 3, cores = 8)
tram_clusterGVis_gene <- sample(tram_heatmap_res$wide.res$gene, 20, replace = FALSE)

tram_clusterGVis_enrich <- enrichCluster(object = tram_heatmap_res, OrgDb = org.Hs.eg.db, type = "BP", organism = "hsa", pvalueCutoff = 0.05, topn = 5)

tram_palette <- c("Blues2", "Purples2", "Purples3", "Greens2")
tram_clusterGVis_enrich_data <- tram_clusterGVis_enrich
if(!"ratio" %in% colnames(tram_clusterGVis_enrich_data)){
  tram_clusterGVis_enrich_data$ratio <- apply(tram_clusterGVis_enrich_data, 1, function(x) {
    gr <- as.numeric(unlist(strsplit(x["GeneRatio"], "/")))
    return(gr[1]/gr[2] * 100)
  })
}

png(paste0(fig_path, 'tram_monocle_clusterGVis_compact_text.png'), units = 'in', res = 300, width = 12, height = 8)
my_visCluster(object = tram_heatmap_res, plotType = "both", column_names_rot = 45, show_row_dend = FALSE,
              markGenes = tram_clusterGVis_gene, markGenesSide = "left", pseudotimeCol = c(my_colors[2], my_colors[3]),
              annoTermData = tram_clusterGVis_enrich, ctAnnoCol = text_colors[1:3], goCol = rep(text_colors[1:3], each = 5),
              genesGp = c("italic", 15, NA), goSize = 16, textboxSize = 12, textboxPos = c(0.5, 0.88), addBar = FALSE, lineSide = "left")
dev.off()

cat("Visualization processing finished successfully.\n")