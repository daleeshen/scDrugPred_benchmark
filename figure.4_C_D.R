# Please change the path here before you run the script.
setwd('D:/scDrugPredict_benchmark/final_submission_files/new_files/codes/')

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidytext)

data_dir <- "./data" 
fig_dir <- "./plots"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

my_colors <- c(
  "scDEAL"     = "#D55E00", 
  "DREEP"      = "#0072B2", 
  "CaDRReS-sc" = "#009E73", 
  "SCAD"       = "#CC79A7", 
  "scDr"       = "#F0E442", 
  "Beyondcell" = "#56B4E9", 
  "scIDUC"     = "#E69F00", 
  "Precily"    = "#999999", 
  "DrugFormer" = "#8A2BE2"  
)

final_df <- read.csv(file.path(data_dir, "preS_vs_preR_perform.csv"), header = TRUE)

final_df$Tool[final_df$Tool == "beyondcell"] <- "Beyondcell"
final_df$Tool[final_df$Tool == "precily"]    <- "Precily"
final_df$Tool[final_df$Tool == "Cadrres-sc"] <- "CaDRReS-sc"

df_sum <- final_df %>%
  pivot_longer(
    cols = c(AUC, ACC, F1, Recall, Precision),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  group_by(Metric, Tool) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se   = sd(Value, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(Method_sorted = reorder_within(Tool, -mean, Metric))

auc_order <- df_sum %>% filter(Metric == "AUC") %>% arrange(desc(mean)) %>% pull(Tool)
df_sum$Tool <- factor(df_sum$Tool, levels = auc_order)

png(file.path(fig_dir, "shuyu_data_perform_preS_vs_preR.png"), units = "in", res = 300, width = 13, height = 6)
ggplot(df_sum, aes(x = Method_sorted, y = mean, fill = Tool)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_x_reordered() + 
  facet_wrap(~Metric, scales = "free_x", ncol = 5) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(size = 20),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  scale_fill_manual(values = my_colors) +
  labs(y = "Mean Value", x = NULL)
dev.off()

load(file.path(data_dir, "figure.4C_D.RData"))

my_umap_theme <- theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    text = element_text(size = 22),
    legend.position = "none",
    legend.title = element_text(margin = margin(r = 10, t = -15))
  )

p_true <- DimPlot(ola_pres_prer_obj, group.by = "drugSens", pt.size = 1.5) + 
  scale_color_manual(values = c("skyblue", "tomato")) +  
  my_umap_theme + 
  labs(title = NULL, x = "UMAP1", y = "UMAP2") +
  annotate("text", x = -Inf, y = Inf, label = "True label", 
           hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold")


png(file.path(fig_dir, "preS_vs_preR_trueLabel_umap.png"), units = "in", res = 300, width = 6, height = 6)
print(p_true + theme(legend.position = "top", legend.text = element_text(size = 26)))
dev.off()

png(file.path(fig_dir, "preS_vs_preR_trueLabel_umap_without_legend.png"), units = "in", res = 300, width = 6, height = 6)
print(p_true)
dev.off()


plot_configs <- list(
  "scDr"       = list(feature = "scdr_pred_score", converse = FALSE, legend_name = "Score"),
  "DREEP"      = list(feature = "dreep_pred_score", converse = FALSE, legend_name = "Score"),
  "scIDUC"     = list(feature = "sciduc_pred_score", converse = FALSE, legend_name = "Score"),
  "CaDRReS-sc" = list(feature = "cadrres_pred_score", converse = FALSE, legend_name = "Score"),
  "Precily"    = list(feature = "precily_pred_score", converse = FALSE, legend_name = "Score"),
  "Beyondcell" = list(feature = "beyondcell_pred_score", converse = FALSE, legend_name = "Score"),
  "SCAD"       = list(feature = "scad_pred_prob", converse = FALSE, legend_name = "Score"),
  "scDEAL"     = list(feature = "scdeal_pred_prob", converse = TRUE, legend_name = "sensitive probability"),
  "DrugFormer" = list(feature = "drugformer_pred_prob", converse = TRUE, legend_name = "Prediction scaled")
)


for (tool in names(plot_configs)) {
  cfg <- plot_configs[[tool]]
  
  # Handle color inversion
  high_col <- ifelse(cfg$converse, "skyblue", "tomato")
  low_col  <- ifelse(cfg$converse, "tomato", "skyblue")
  
  # Format filename securely
  file_name <- paste0("preS_vs_preR_", gsub("-", "_", tool), "_umap.png")
  
  png(file.path(fig_dir, file_name), units = "in", res = 300, width = 6, height = 6)
  
  p <- FeaturePlot(ola_pres_prer_obj, features = cfg$feature, pt.size = 1.5, order = FALSE) + 
    scale_color_gradient(high = high_col, low = low_col, name = cfg$legend_name) + 
    my_umap_theme + 
    labs(title = NULL, x = "UMAP1", y = "UMAP2") +
    annotate("text", x = -Inf, y = Inf, label = tool, hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold") + 
    guides(color = guide_colorbar(barwidth = 12, barheight = 1)) 
  
  print(p)
  dev.off()
}

message("All preS vs preR figures generated successfully.")