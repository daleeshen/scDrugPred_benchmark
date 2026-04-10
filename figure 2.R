setwd('D:/scDrugPredict_benchmark/')

library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)
library(tidyr)

INPUT_DIR  <- "./"
OUTPUT_DIR <- "./plots"


# ==============================================================================
# 1. Data Preprocessing & Merging
# ==============================================================================
final_merged_df<-read.csv('./bal_final_merged_data.csv', header = T)
meta <- read.csv(file.path(INPUT_DIR, 'drmref_dataInfo.csv'), header = TRUE)
meta_brief <- read.csv(file.path(INPUT_DIR, 'dataset_drug_single.csv'), header = TRUE)

my_colors <- c(
  "scDEAL"     = "#D55E00", "DREEP"      = "#0072B2", "CaDRReS-sc" = "#009E73",
  "SCAD"       = "#CC79A7", "scDr"       = "#F0E442", "Beyondcell" = "#56B4E9",
  "scIDUC"     = "#E69F00", "Precily"    = "#999999", "DrugFormer" = "#8A2BE2"
)
tool_order <- c("scDEAL", "CaDRReS-sc", "SCAD", "DREEP", "DrugFormer", "Precily", "scDr", "Beyondcell", "scIDUC")

tissue_target_order <- c("Cell line", "PBMC", "Tumor tissue", "BMA")
drug_target_order   <- c("Chemotherapy", "Targeted therapy", "Immunotherapy")

# ==============================================================================
# 2. Performance vs. Cell Count
# ==============================================================================

# Prepare factor ordering for Boxplots
df_box_auc <- final_merged_df %>%
  group_by(cells) %>%
  mutate(Tool = fct_reorder(Tool, AUC, .fun = median, .desc = TRUE)) %>%
  ungroup()

# Figure 2A (Up)
png(file.path(OUTPUT_DIR, 'balanced_auc_all_boxplot.png'), units = 'in', res = 300, width = 13, height = 4)
ggplot(df_box_auc, aes(x = as.factor(cells), y = AUC, fill = Tool)) +
  geom_boxplot(alpha = 0.8) +
  labs(x = "Cells", y = "AUROC") +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 22), legend.position = "right") +
  scale_fill_manual(values = my_colors)
dev.off()

df_box_acc <- final_merged_df %>%
  group_by(cells) %>%
  mutate(Tool = fct_reorder(Tool, ACC, .fun = median, .desc = TRUE)) %>%
  ungroup()

# Figure 2A (Down)
png(file.path(OUTPUT_DIR, 'balanced_acc_all_boxplot.png'), units = 'in', res = 300, width = 13, height = 4)
ggplot(df_box_acc, aes(x = as.factor(cells), y = ACC, fill = Tool)) +
  geom_boxplot(alpha = 0.8) +
  labs(x = "Cells", y = "ACC") +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 22), legend.position = "right") +
  scale_fill_manual(values = my_colors)
dev.off()

plot_df_auc <- final_merged_df %>%
  group_by(Tool, cells) %>%
  summarise(mean_auc = mean(AUC, na.rm = TRUE), se_auc = sd(AUC, na.rm = TRUE) / sqrt(n()), .groups = "drop")

# Figure 2B (Left)
png(file.path(OUTPUT_DIR, 'balanced_auc_all_line.png'), units = 'in', res = 300, width = 9, height = 5)
ggplot(plot_df_auc, aes(x = cells, y = mean_auc, color = Tool, group = Tool)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  geom_ribbon(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc, fill = Tool), alpha = 0.2, color = NA) +
  theme_bw() + labs(x = "Cells", y = "AUROC") +
  theme(panel.grid = element_blank(), text = element_text(size = 22), legend.title = element_blank()) +
  scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
dev.off()

plot_df_acc <- final_merged_df %>%
  group_by(Tool, cells) %>%
  summarise(mean_acc = mean(ACC, na.rm = TRUE), se_acc = sd(ACC, na.rm = TRUE) / sqrt(n()), .groups = "drop")

# Figure 2B (right)
png(file.path(OUTPUT_DIR, 'balanced_acc_all_line.png'), units = 'in', res = 300, width = 9, height = 5)
ggplot(plot_df_acc, aes(x = cells, y = mean_acc, color = Tool, group = Tool)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  geom_ribbon(aes(ymin = mean_acc - se_acc, ymax = mean_acc + se_acc, fill = Tool), alpha = 0.2, color = NA) +
  theme_bw() + labs(x = "Cells", y = "ACC") +
  theme(panel.grid = element_blank(), text = element_text(size = 22), legend.title = element_blank()) +
  scale_color_manual(values = my_colors) + scale_fill_manual(values = my_colors)
dev.off()

# ==============================================================================
# 3. Performance by Tissue & Drug Type
# ==============================================================================

generate_bar_plot <- function(data, group_col, metric, y_label, order_levels, filename) {
  
  # 1. Base Statistics
  df_summary <- data %>%
    group_by(!!sym(group_col), Tool) %>%
    summarise(
      mean_val = mean(!!sym(metric), na.rm = TRUE),
      se_val   = sd(!!sym(metric), na.rm = TRUE) / sqrt(n()),
      .groups  = "drop"
    ) %>%
    filter(!is.na(!!sym(group_col)))
  
  # 2. Core Logic: Calculate Absolute X Coordinates for Each Bar
  bar_width <- 0.7 
  dodge_width <- 0.8
  
  df_summary <- df_summary %>%
    mutate(!!sym(group_col) := factor(!!sym(group_col), levels = order_levels)) %>%
    arrange(!!sym(group_col), desc(mean_val)) %>%
    group_by(!!sym(group_col)) %>%
    mutate(
      n_tools = n(), 
      tool_idx = row_number(),
      group_center = as.numeric(!!sym(group_col)),
      x_offset = (tool_idx - (n_tools + 1) / 2) * (dodge_width / max(n_tools)),
      x_pos = group_center + x_offset
    ) %>%
    ungroup()
  
  # Dynamically calculate physical bar width to keep all bars uniform
  max_tools_in_a_group <- max(df_summary$n_tools)
  actual_bar_width <- bar_width / max_tools_in_a_group
  
  png(file.path(OUTPUT_DIR, filename), units = 'in', res = 300, width = 11, height = 6)
  
  # 3. Plotting: Continuous X-axis with Discrete Labels
  p <- ggplot(df_summary, aes(x = x_pos, y = mean_val, fill = Tool)) +
    geom_bar(
      stat = "identity", 
      width = actual_bar_width, 
      na.rm = TRUE
    ) +
    geom_errorbar(
      aes(ymin = mean_val - se_val, ymax = mean_val + se_val), 
      width = actual_bar_width * 0.4, 
      color = "black", 
      linewidth = 0.8, 
      na.rm = TRUE
    ) +
    scale_fill_manual(values = my_colors) +
    scale_x_continuous(
      breaks = 1:length(order_levels), 
      labels = order_levels,
      limits = c(0.5, length(order_levels) + 0.5) 
    ) +
    theme_bw() +
    labs(x = "", y = y_label) +
    theme(
      text = element_text(size = 22),
      axis.text.x = element_text(angle = 0),
      legend.title = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  print(p)
  dev.off()
}

generate_bar_plot(final_merged_df, "Tissue", "AUC", "Mean AUROC", tissue_target_order, "balanced_auc_by_tissue.png")
generate_bar_plot(final_merged_df, "Tissue", "ACC", "Mean ACC", tissue_target_order, "balanced_acc_by_tissue.png")
generate_bar_plot(final_merged_df, "Drug.type", "AUC", "Mean AUROC", drug_target_order, "balanced_auc_by_drugType.png")
generate_bar_plot(final_merged_df, "Drug.type", "ACC", "Mean ACC", drug_target_order, "balanced_acc_by_drugType.png")


# ==============================================================================
# 4. Global Comparison: Cell line vs Tumor / Chemo vs Targeted
# ==============================================================================

# general tissue comparison
final_merged_df_tissue_sub <- subset(final_merged_df, Tissue %in% c('Cell line', 'Tumor tissue'))

png(file.path(OUTPUT_DIR, 'auc_comp_by_tissue_overall.png'), units = 'in', res = 300, width = 7, height = 5)
ggplot(final_merged_df_tissue_sub, aes(x = Tissue, y = AUC, fill = Tissue)) +
  geom_boxplot(width = 0.45, color = "black", outlier.shape = NA, alpha = 0.9, size = 0.8) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1.8, color = "grey50") +
  scale_fill_manual(values = c("Cell line" = "#F4A9A8", "Tumor tissue" = "#A7D3E8")) +
  annotate("segment", x = 1, xend = 2, y = 1.06, yend = 1.06, size = 1) +
  annotate("text", x = 1.5, y = 1.1, label = expression(p == 2.8 %*% 10^-11), size = 8) + 
  labs(x = NULL,y = "AUROC") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(), text = element_text(size = 25))
dev.off()

final_merged_df_drug_sub <- subset(final_merged_df, Drug.type %in% c('Chemotherapy', 'Targeted therapy'))

png(file.path(OUTPUT_DIR, 'auc_comp_by_drugtype_overall.png'), units = 'in', res = 300, width = 7, height = 5)
ggplot(final_merged_df_drug_sub, aes(x = Drug.type, y = AUC, fill = Drug.type)) +
  geom_boxplot(width = 0.45, color = "black", outlier.shape = NA, alpha = 0.9, size = 0.8) +
  geom_jitter(width = 0.12, alpha = 0.4, size = 1.8, color = "grey50") +
  scale_fill_manual(values = c("Chemotherapy" = "#F4A9A8", "Targeted therapy" = "#009E73")) +
  annotate("segment", x = 1, xend = 2, y = 1.06, yend = 1.06, size = 1) +
  annotate("text", x = 1.5, y = 1.1, label = expression(p == 3 %*% 10^-7), size = 8) +
  labs(x = NULL, y = "AUROC") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(), text = element_text(size = 25))
dev.off()

# ==============================================================================
# 5. Tool-specific Significance (Facetted / Grouped)
# ==============================================================================

final_merged_df_tissue_sub <- final_merged_df_tissue_sub %>%
  mutate(Tool = factor(Tool, levels = tool_order))

p_val_tissue_by_tool <- final_merged_df_tissue_sub %>%
  group_by(Tool) %>%
  summarise(
    p_value = if(n_distinct(Tissue) == 2) wilcox.test(AUC ~ Tissue)$p.value else NA,
    y_pos  =1.05,
    .groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value <= 0.001 ~ "***", p_value <= 0.01 ~ "**", 
    p_value <= 0.05 ~ "*", TRUE ~ "ns"
  ))

p_val_tissue_by_tool <- p_val_tissue_by_tool %>%
  mutate(Tool = factor(Tool, levels = tool_order))

png(file.path(OUTPUT_DIR, 'auc_comp_by_tissue_per_tool.png'), units = 'in', res = 300, width = 12, height = 7)
ggplot(final_merged_df_tissue_sub, aes(x = Tool, y = AUC, fill = Tissue)) +
  geom_boxplot(width = 0.45,
               color = "black",
               outlier.shape = NA,
               alpha = 0.9, size = 0.8) +
  scale_fill_manual(values = c(
    "Cell line" = "#F4A9A8",
    "Tumor tissue"  = "#A7D3E8"
  ), name = 'Sample source') +
  
  geom_text(
    data = p_val_tissue_by_tool,
    aes(x = Tool, y = y_pos, label = signif),
    inherit.aes = FALSE,
    size = 8
  ) + 
  labs(x = "",
       y = "AUROC") +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    text = element_text(size = 25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) + ylim(0,1.1)
dev.off()

# ==============================================================================
# 6. Tool-specific Significance: Drug Type (Chemo vs Targeted)
# ==============================================================================

final_merged_df_drug_sub <- final_merged_df_drug_sub %>%
  mutate(Tool = factor(Tool, levels = tool_order))

p_val_dt_by_tool <- final_merged_df_drug_sub %>%
  group_by(Tool) %>%
  summarise(
    p_value = if(n_distinct(Drug.type) == 2) wilcox.test(AUC ~ Drug.type)$p.value else NA,
    y_pos   = 1.05,
    .groups = "drop"
  ) %>%
  mutate(signif = case_when(
    p_value <= 0.0001 ~ "****",
    p_value <= 0.001  ~ "***", 
    p_value <= 0.01   ~ "**",  
    p_value <= 0.05   ~ "*",   
    TRUE              ~ "ns"
  ))

p_val_dt_by_tool <- p_val_dt_by_tool %>%
  mutate(Tool = factor(Tool, levels = tool_order))

png(file.path(OUTPUT_DIR, 'auc_comp_by_drugtype_per_tool.png'), units = 'in', res = 300, width = 12, height = 7)
ggplot(final_merged_df_drug_sub, aes(x = Tool, y = AUC, fill = Drug.type)) +
  geom_boxplot(width = 0.45,
               color = "black",
               outlier.shape = NA,
               alpha = 0.9, size = 0.8) +
  scale_fill_manual(values = c(
    "Chemotherapy" = "#F4A9A8",
    "Targeted therapy"  = "#009E73"
  ), name = 'Drug Type') +
  
  geom_text(
    data = p_val_dt_by_tool,
    aes(x = Tool, y = y_pos, label = signif),
    inherit.aes = FALSE,
    size = 8
  ) + 
  labs(x = "",
       y = "AUROC") + 
  theme_bw(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    text = element_text(size = 25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) + 
  ylim(0, 1.1)
dev.off()

