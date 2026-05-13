# Please change the path here before you run the script.
setwd('./')

library(Seurat)
library(tidyverse)
library(tidytext)

data_dir <- "../data"
fig_dir <- '../results/'
# dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(data_dir, "preS_vs_R_perform.csv"), header = TRUE)

# 统一命名规范
df$Method[df$Method == "beyondcell"] <- "Beyondcell"
df$Method[df$Method == "precily"]    <- "Precily"
df$Method[df$Method == "Cadrres-sc"] <- "CaDRReS-sc"

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


df_sum <- df %>%
  pivot_longer(
    cols = c("AUC", "ACC", "F1", "Precision", "Recall"), 
    names_to = "Metric", 
    values_to = "value"
  ) %>%
  group_by(Method, Metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n), 
    .groups = "drop"
  ) %>%
  mutate(Method_sorted = reorder_within(Method, -mean, Metric))

auc_order <- df_sum %>% filter(Metric == "AUC") %>% arrange(desc(mean)) %>% pull(Method)
df_sum$Method <- factor(df_sum$Method, levels = auc_order)


png(file.path(fig_dir, "shuyu_data_perform.png"), units = "in", res = 300, width = 13, height = 6)
ggplot(df_sum, aes(x = Method_sorted, y = mean, fill = Method)) +
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
  scale_fill_manual(values = my_colors) + labs(y = "Mean Value", x = NULL)
dev.off()


load(file.path(data_dir, "/figure.4A_B.RData"))

my_umap_theme <- theme_bw() + 
  theme(
    panel.grid = element_blank(), 
    text = element_text(size = 22),
    legend.position = "none" 
  )

png(file.path(fig_dir, "ola_preS_vs_R_umap_TrueLabel.png"), units = "in", res = 300, width = 6, height = 6)
DimPlot(ola_obj, group.by = "drugSens", pt.size = 1.5) + 
  my_umap_theme + 
  theme(legend.position = "top", legend.text = element_text(size = 26)) +
  scale_color_manual(values = c("tomato", "skyblue")) +  
  labs(title = NULL, x = "UMAP1", y = "UMAP2") +
  annotate("text", x = -Inf, y = Inf, label = "True label", hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold")
dev.off()

methods_to_plot <- c(
  "scDEAL"     = "scdeal_pred_prob",
  "Beyondcell" = "beyondcell_pred_score",
  "DREEP"      = "dreep_pred_score",
  "CaDRReS-sc" = "cadrres_pred_score",
  "Precily"    = "precily_pred_score",
  "SCAD"       = "scad_pred_prob",
  "scIDUC"     = "sciduc_pred_score",
  "scDr"       = "scdr_pred_score",
  "DrugFormer" = "drugformer_pred_prob"
)

for (method_name in names(methods_to_plot)) {
  feature_col <- methods_to_plot[[method_name]]
  
  file_name <- paste0("ola_preS_vs_R_umap_", tolower(gsub("-", "", method_name)), ".png")
  
  png(file.path(fig_dir, file_name), units = "in", res = 300, width = 6, height = 6)
  
  p <- FeaturePlot(ola_obj, features = feature_col, pt.size = 1.5, order = FALSE) + 
    scale_color_gradient(high = "tomato", low = "skyblue") + 
    my_umap_theme + 
    labs(title = NULL, x = "UMAP1", y = "UMAP2") +
    annotate("text", x = -Inf, y = Inf, label = method_name, hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold")
  
  print(p)
  dev.off()
}


