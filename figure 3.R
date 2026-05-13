# Please change the path here before you run the script.
setwd('./')

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ComplexHeatmap)
library(circlize)
library(ggbreak)
library(gridExtra)

# ==============================================================================
# 1. Dataset preparation
# ==============================================================================

fig_path <- '../results/'
# if(!dir.exists(fig_path)) dir.create(fig_path)

color_ref <- c(
  "scDEAL"      = "#E7A673",
  "CaDRReS-sc"  = "#69C6AC",
  "SCAD"        = "#DFABC8",
  "DREEP"       = "#61A7CF",
  "DrugFormer"  = "#B87FED",
  "Precily"     = "#C5C5C5",
  "scDr"        = "#F1EE9A",
  "Beyondcell"  = "#AAD3E3",
  "scIDUC"      = "#F3D48E"
)

tool_order <- c("scDEAL", "CaDRReS-sc", "SCAD", "DREEP", "DrugFormer", "Precily", "scDr", "Beyondcell", "scIDUC")
ratio_levels <- c("1:3", "1:10", "1:30", "1:100")

merged_df <- read.csv("../data/unb_merged_performance_df.csv")
merged_df$Ratio <- factor(merged_df$Ratio, levels = ratio_levels)
balanced_rd <- read.csv("../data/bal_final_merged_data.csv")
unb_stats_df <- read.csv("../data/unb_cell_counts.csv")
gse108397_scores <- read.csv("../data/GSE108397_pred_scores.csv")
gse163836_scores <- read.csv("../data/GSE163836_pred_scores.csv")


# ==============================================================================
# 2. Dataset Proportion Bar Plot
# ==============================================================================
unb_stats_plot_df <- unb_stats_df %>%
  mutate(
    total = Sensitive + Resistant,
    S_percent = Sensitive / total * 100,
    R_percent = Resistant / total * 100
  ) %>%
  select(Dataset, Ratio, S_percent, R_percent) %>%
  pivot_longer(cols = c(S_percent, R_percent), names_to = "Condition", values_to = "Percent") %>%
  mutate(
    Condition = factor(Condition, levels = c("S_percent", "R_percent"), labels = c("Sensitive", "Resistant")),
    Ratio = factor(gsub("_", ":", Ratio), levels = ratio_levels)
  )

# Clean dataset names
unb_stats_plot_df$Dataset <- gsub("_Vem", "", unb_stats_plot_df$Dataset)
unb_stats_plot_df$Dataset <- gsub("_PacBlood", "_Blood", unb_stats_plot_df$Dataset)
unb_stats_plot_df$Dataset <- gsub("_PacTissue", "_Tissue", unb_stats_plot_df$Dataset)

png(paste0(fig_path, 'unb_proportion_plot.png'), units = 'in', res = 300, width = 47, height = 13)
ggplot(unb_stats_plot_df, aes(x = Ratio, y = Percent, fill = Condition)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ Dataset, nrow = 3, scales = "free_y") + 
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  scale_fill_manual(values = c("Sensitive" = "#4DBBD5", "Resistant" = "#E64B35")) +
  theme_bw(base_size = 14) +
  labs(x = "", y = "Proportion", fill = NULL) + 
  theme(
    panel.grid = element_blank(),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(size = 30, face = "bold", margin = margin(t = 15, b = 15)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.text.y = element_text(size = 25),
    axis.title.y = element_text(size = 30),
    legend.text = element_text(size = 30),
    panel.spacing = unit(2, "lines")
  )
dev.off()

# ==============================================================================
# 3. Performance Boxplots & Trends (AUROC, ACC, NormAUPRC)
# ==============================================================================

plot_performance_metrics <- function(df, metric_col, plot_prefix, y_label) {
  
  # Prepare data
  df_metric <- df %>%
    group_by(Ratio, Tool) %>%
    mutate(Mean_Val = mean(!!sym(metric_col), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Tool_sorted = fct_reorder(Tool, -Mean_Val, .fun = mean))
  
  # Trend plot data summary
  df_summary <- df_metric %>%
    group_by(Ratio, Tool_sorted) %>%
    summarise(
      Mean = mean(!!sym(metric_col), na.rm = TRUE),
      se = sd(!!sym(metric_col), na.rm = TRUE)/sqrt(n()),
      ymin = Mean - 1.96 * se,
      ymax = Mean + 1.96 * se,
      .groups = "drop"
    )
  
  # Boxplot
  png(paste0(fig_path, plot_prefix, '_boxplot.png'), units = 'in', res = 300, width = 10, height = 5)
  p_box <- ggplot(df_metric, aes_string(x="Ratio", y=metric_col, fill="Tool_sorted")) + 
    geom_boxplot() + 
    theme_bw() + 
    theme(text = element_text(size = 18), panel.grid = element_blank()) + 
    scale_fill_manual(values = color_ref, name = "Tool", breaks = tool_order) +
    ylab(y_label)
  print(p_box)
  dev.off()
  
  # Trend plot
  png(paste0(fig_path, plot_prefix, '_trend_plot.png'), units = 'in', res = 300, width = 8, height = 5)
  p_trend <- ggplot(df_summary, aes(x = Ratio, y = Mean, color = Tool_sorted, fill = Tool_sorted, group = Tool_sorted)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.2, color = NA) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    theme_bw() +
    theme(text = element_text(size = 18), legend.position = "bottom") +
    scale_color_manual(values = color_ref, name = "Tool", breaks = tool_order) +
    scale_fill_manual(values = color_ref, name = "Tool", breaks = tool_order) + 
    ylab(paste("Mean", y_label))
  print(p_trend)
  dev.off()
}

plot_performance_metrics(merged_df, "AUROC", "unb_auc", "AUROC")
plot_performance_metrics(merged_df, "ACC", "unb_acc", "Accuracy")

# NormAUPRC Boxplot (With break)
merged_df_nAP <- merged_df %>%
  group_by(Ratio, Tool) %>%
  mutate(NormAP_mean = mean(NormAUPRC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Tool_sorted = fct_reorder(Tool, -NormAP_mean, .fun = mean))

png(paste0(fig_path, 'unb_normAP_boxplot.png'), units = 'in', res = 300, width = 10, height = 10)
ggplot(merged_df_nAP, aes(x=Ratio, y=NormAUPRC, fill=Tool_sorted)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(
    text = element_text(size = 28),
    panel.grid = element_blank(),
    legend.title = element_blank()
  ) + 
  scale_y_break(c(-5, -7)) + 
  scale_fill_manual(values = color_ref, name = "Tool", breaks = tool_order) + 
  coord_cartesian(ylim = c(-1, 1)) 
dev.off()

# ==============================================================================
# 4. Balanced vs Unbalanced AUROC Comparison
# ==============================================================================
df_unbalanced <- merged_df %>% select(Dataset, AUROC, Tool) %>% mutate(Scenario = 'Unbalanced')
df_balanced <- balanced_rd %>% select(Dataset, AUROC = AUC, Tool) %>% mutate(Scenario = 'Balanced')
df_compare <- bind_rows(df_unbalanced, df_balanced)

p_df <- df_compare %>%
  group_by(Tool) %>%
  summarise(p = wilcox.test(AUROC ~ Scenario)$p.value) %>%
  mutate(
    label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    )
  )

png(paste0(fig_path, 'unb_vs_bal_auc.png'), units = 'in', res = 300, width = 12, height = 6)
ggplot(df_compare, aes(x=Tool, y=AUROC)) + 
  geom_violin(aes(fill = Scenario), position = position_dodge(0.7), trim = T, adjust = 1) + 
  geom_boxplot(aes(color=Scenario, group = interaction(Tool, Scenario)), width = 0.1, fill = "white", position = position_dodge(0.7)) +
  geom_text(
    data = p_df,
    aes(x = Tool, y = 1.12, label = label, vjust = ifelse(grepl("\\*", label), 0.78, 0.5)),
    size = 10, fontface = "bold", inherit.aes = FALSE
  ) +
  theme_bw() + 
  theme(
    text = element_text(size = 28),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = 'top'
  ) + 
  scale_color_manual(values = c('black','black')) + 
  scale_fill_manual(values = c('skyblue', 'tomato')) + 
  xlab('') + ylim(0, 1.25)
dev.off()

# ==============================================================================
# 5. Dataset Majority F1 Boxplot (S_majority vs R_majority)
# ==============================================================================
unb_dataset_majority <- unb_stats_plot_df %>%
  group_by(Dataset) %>%
  summarise(
    S_total = sum((Condition == "Sensitive") * Percent, na.rm = TRUE),
    R_total = sum((Condition == "Resistant") * Percent, na.rm = TRUE)
  ) %>%
  mutate(Majority = ifelse(S_total >= R_total, "S_majority", "R_majority"))

merged_df_maj <- merged_df %>%
  mutate(
    Clean_Dataset = gsub("_DMSO_rm|_Vem", "", Dataset),
    Clean_Dataset = gsub("_PacBlood", "_Blood", Clean_Dataset),
    Clean_Dataset = gsub("_PacTissue", "_Tissue", Clean_Dataset)
  ) %>%
  left_join(
    unb_dataset_majority %>% select(Dataset, Majority), 
    by = c("Clean_Dataset" = "Dataset") 
  ) %>%
  mutate(
    Tool = factor(Tool, levels = sort(unique(as.character(Tool)))),
    Ratio = factor(Ratio, levels = c("1:3", "1:10", "1:30", "1:100"))
  )

tool_levels <- levels(merged_df_maj$Tool)

background_data <- data.frame(
  Tool = tool_levels,
  xmin = seq_along(tool_levels) - 0.5,
  xmax = seq_along(tool_levels) + 0.5,
  fill = ifelse(seq_along(tool_levels) %% 2 == 0, "grey60", "white")
)

background_data_full <- merge(background_data, data.frame(Ratio = c("1:3", "1:10", "1:30", "1:100")))
background_data_full$Ratio <- factor(background_data_full$Ratio, levels = c("1:3", "1:10", "1:30", "1:100"))

png(paste0(fig_path, 'dataset_majority_F1.png'), units = 'in', res = 300, width = 9, height = 8)
p_maj <- ggplot(merged_df_maj, aes(x = Tool, y = F1, fill = Majority)) + 
  geom_rect(
    data = background_data_full, 
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = NULL),
    fill = background_data_full$fill, alpha = 0.3, inherit.aes = FALSE
  ) +
  geom_boxplot(outlier.size = 0.5) + 
  facet_wrap(~ Ratio, nrow = 2) + 
  scale_fill_manual(values = c("S_majority"="#4DBBD5", "R_majority"="#E64B35")) + 
  theme_bw(base_size = 14) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(size = 22, face = "bold", margin = margin(t = 6, r = 4, b = 6, l = 4)),
    panel.grid = element_blank(),
    axis.text = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.title = element_blank(),
    legend.position = 'bottom'
  ) + 
  labs(x = "", y = "F1 Score", fill = "Dataset Majority") +
  scale_x_discrete(limits = tool_levels)

print(p_maj)
dev.off()

# ==============================================================================
# 6. F1, Precision, Recall Heatmaps 
# ==============================================================================
df_agg <- merged_df %>%
  group_by(Dataset, Tool) %>%
  summarise(
    F1 = mean(F1, na.rm = TRUE),
    Precision = mean(Precision, na.rm = TRUE),
    Recall = mean(Recall, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(F1, Precision, Recall), ~ ifelse(is.nan(.), NA_real_, .)))

df_full <- expand_grid(Dataset = unique(df_agg$Dataset), Tool = unique(df_agg$Tool)) %>%
  left_join(df_agg, by = c("Dataset", "Tool")) %>%
  mutate(Dataset = gsub("_DMSO_rm|_Vem", "", Dataset)) %>%
  mutate(Dataset = gsub("_PacBlood", "_Blood", Dataset)) %>%
  mutate(Dataset = gsub("_PacTissue", "_Tissue", Dataset))

make_metric_matrix <- function(df, metric_col) {
  df_wide <- df %>%
    select(Dataset, Tool, value = {{ metric_col }}) %>%
    pivot_wider(names_from = Tool, values_from = value) %>%
    arrange(Dataset)
  mat <- as.matrix(df_wide[, -1, drop = FALSE])
  rownames(mat) <- df_wide$Dataset
  return(mat)
}

col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
font_size_main <- 20; font_size_title <- 24; font_size_X <- 18; font_size_legend_title <- 18; font_size_legend_label <- 15

create_heatmap <- function(mat, title) {
  Heatmap(
    mat, name = title, col = col_fun, na_col = "white",
    cluster_rows = FALSE, cluster_columns = FALSE,
    column_title = title, column_title_gp = gpar(fontsize = font_size_title, fontface = "bold"),
    row_names_gp = gpar(fontsize = font_size_main), column_names_gp = gpar(fontsize = font_size_main),
    heatmap_legend_param = list(direction = "horizontal", title_gp = gpar(fontsize = font_size_legend_title, fontface = "bold"), labels_gp = gpar(fontsize = font_size_legend_label)),
    cell_fun = function(j, i, x, y, w, h, col) {
      if (is.na(mat[i, j])) grid.text("X", x = x, y = y, gp = gpar(fontsize = font_size_X, col = "black"))
    }
  )
}

ht_list <- create_heatmap(make_metric_matrix(df_full, F1), "F1") + 
  create_heatmap(make_metric_matrix(df_full, Precision), "Precision") + 
  create_heatmap(make_metric_matrix(df_full, Recall), "Recall")

png(paste0(fig_path, 'unb_f1_pr_recall_heatmap.png'), units = 'in', res = 300, width = 11, height = 13)
draw(ht_list, heatmap_legend_side = "bottom", padding = unit(c(0.5, 0.5, 0.5, 2), "cm"))
dev.off()

# ==============================================================================
# 7. Detailed F1/Precision/Recall Heatmaps per Ratio (Figure S1)
# ==============================================================================
make_ratio_heatmap <- function(df, metric, ratio_value) {
  
  df_wide <- df %>%
    filter(Ratio == ratio_value) %>%
    group_by(Dataset, Tool) %>%
    summarise(value = mean(!!sym(metric), na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Dataset = gsub("_DMSO_rm|_Vem", "", Dataset),
      Dataset = gsub("_PacBlood", "_Blood", Dataset),
      Dataset = gsub("_PacTissue", "_Tissue", Dataset)
    ) %>%
    tidyr::pivot_wider(names_from = Tool, values_from = value) %>%
    arrange(Dataset)
  
  mat <- as.matrix(df_wide[, -1, drop = FALSE])
  rownames(mat) <- df_wide$Dataset

  
  if (metric == "NormAUPRC") {
    col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  } else {
    col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  }
  
  Heatmap(
    mat,
    name = paste0(metric, "_", ratio_value),
    col = col_fun,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = paste0(metric, " (", ratio_value, ")"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cell_fun = function(j, i, x, y, w, h, col) {
      if (is.na(mat[i, j])) {
        grid.text("X", x = x, y = y, gp = gpar(fontsize = 10, col = "black"))
      }
    }
  )
}

metrics <- c("F1", "Precision", "Recall", "AUROC", "NormAUPRC")
heatmap_list <- list()

for (metric in metrics) {
  for (rt in ratio_levels) { 
    heatmap_name <- paste(metric, rt, sep = "_")
    heatmap_list[[heatmap_name]] <- make_ratio_heatmap(merged_df, metric, rt)
  }
}

grob_list <- lapply(heatmap_list, function(ht) grid.grabExpr(draw(ht)))

png(paste0(fig_path, 'unb_all_metrics_per_ratio_heatmap.png'), units = 'in', res = 300, width = 22, height = 33)
grid.arrange(grobs = grob_list, nrow = 5, ncol = 4)
dev.off()

# ==============================================================================
# 8. Density Plots of Prediction Scores 
# ==============================================================================

plot_tool_density <- function(score_df, tool_name, dataset_name) {
  
  df_filtered <- score_df %>% 
    filter(Tool == tool_name) %>%
    mutate(
      ratio = factor(ratio, levels = c('1:3', '1:30', '1:100')),
      condition = factor(condition, levels = c('sensitive', 'resistant'))
    )
  
  if(nrow(df_filtered) == 0) return(NULL)
  
  p <- ggplot(df_filtered, aes(x = score, color = condition)) +
    geom_density(size = 1) +
    geom_rug(alpha = 0.2) +
    facet_wrap(~ ratio, nrow = 1) +
    labs(title = tool_name) + 
    geom_text(
      aes(label = ratio), x = Inf, y = Inf, hjust = 1.2, vjust = 1.5, 
      color = "black", size = 6, fontface = "bold", check_overlap = TRUE
    ) +
    theme_bw(base_size = 14) + 
    scale_color_manual(values = c('sensitive' = 'skyblue', 'resistant' = 'tomato')) + 
    xlab('Score / Prob') +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_blank(),
      strip.text = element_blank(), 
      panel.grid = element_blank(),
      axis.text = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    )
  
  file_name <- paste0(fig_path, tolower(tool_name), "_", dataset_name, "_example_pred_score_density.png")
  png(file_name, units = 'in', res = 300, width = 9, height = 4)
  print(p)
  dev.off()
}
tools_gse108397 <- c("CaDRReS-sc", "SCAD", "Beyondcell", "Precily", "DrugFormer", "scDEAL")

for (tool in tools_gse108397) {
  plot_tool_density(gse108397_scores, tool, "gse108397")
}

tools_gse163836 <- c("scDr", "DREEP", "scDEAL")

for (tool in tools_gse163836) {
  plot_tool_density(gse163836_scores, tool, "gse163836")
}