# Please change the path here before you run the script.
setwd('./')

library(tidyverse)
library(tidytext)
library(patchwork)
library(scales)
library(ggnewscale)
library(ggrain)

my_colors <- c('#0DA8BB','#14A97A','#833839','#DC5805','#FF9AB3','#711C9A','#74FCD4','lightgrey','#046488')

out_dir <- '../results/'

# ==============================================================================
# Part 1: Dataset Statistics Barplot
# ==============================================================================
# Note: Ensure your data files are placed in the './data/' directory
df <- read.csv('../data/drmref_info_final.csv')

# Clean and format data
df <- df %>%
  mutate(
    Cancer = case_when(
      Cancer.type.level1 == 'Acute lymphoblastic leukemia' ~ 'ALL',
      Cancer.type.level1 == 'Acute myeloid leukemia' ~ 'AML',
      Cancer.type.level1 == 'Chronic lymphocytic leukemia' ~ 'CLL',
      TRUE ~ Cancer.type.level1
    ),
    Origin = ifelse(Tissue %in% c('Bone marrow aspirate', 'Tumor tissue', 'Peripheral blood mononuclear cells'), 'Patient', 'Cell line'),
    Tissue = case_when(
      Tissue == 'Peripheral blood mononuclear cells' ~ 'PBMC',
      Tissue == 'Bone marrow aspirate' ~ 'BMA',
      TRUE ~ Tissue
    )
  )

# Calculate counts for plotting
df_counts <- df %>%
  select(Tissue, Cancer, Drug.type) %>%
  pivot_longer(cols = everything(), names_to = "feature", values_to = "level") %>%
  count(feature, level, name = "count") %>%
  mutate(level = ifelse(level %in% c("Melanoma", "Multiple myeloma"), "Melanoma", level)) %>%
  group_by(feature, level) %>%
  summarise(count = sum(count), .groups = "drop")

# Coordinates for segment lines
line_df <- df_counts %>%
  group_by(feature) %>%
  summarise(
    y_line = max(count) * 1.05,
    x_min = 0.5,
    x_max = n_distinct(level) + 0.5,
    .groups = "drop"
  )

p_dataset_stats <- ggplot(df_counts, aes(x = reorder_within(level, count, feature), y = count, fill = feature)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_grid(rows = vars(feature), scales = "free_y", space = "free_y", 
             labeller = as_labeller(c("Cancer" = "Cancer Type", "Drug.type" = 'Drug Type', 'Tissue' = 'Tissue'))) +
  scale_x_reordered() +
  scale_fill_manual(values = my_colors[1:3]) +
  geom_segment(data = line_df, aes(x = x_min, xend = x_max, y = 16, yend = 16), 
               inherit.aes = FALSE, linewidth = 0.8, color = "black") +
  labs(x = NULL, y = "Number of Datasets") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = -90, size = 24),
    panel.grid = element_blank(),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    panel.border = element_blank(),
    strip.background = element_blank()
  )

png(paste0(out_dir, 'Dataset_stats_barplot.png'), units = 'in', res = 300, width = 10, height = 10.5)
print(p_dataset_stats)
dev.off()

# ==============================================================================
# Part 2: Cell Count Barplot
# ==============================================================================
cell_count_df <- read.csv('../data/cell_count.csv', header = TRUE) %>%
  mutate(Dataset = gsub('_DMSO_rm', '', Dataset))

df_agg <- cell_count_df %>%
  left_join(df %>% select(Dataset, Tissue, Drug.type, Cancer, Origin), by = "Dataset") %>%
  pivot_longer(cols = c(nCount_resistant, nCount_sensitive), names_to = "Condition", values_to = "Cell_Count") %>%
  mutate(Condition = recode(Condition, nCount_resistant = "Resistant", nCount_sensitive = "Sensitive")) %>%
  select(Cancer, Drug.type, Tissue, Condition, Cell_Count) %>%
  pivot_longer(cols = c(Cancer, Drug.type, Tissue), names_to = "Block", values_to = "Label") %>%
  mutate(
    Block = factor(Block, levels = c("Cancer", "Drug.type", "Tissue")),
    Label = ifelse(Label %in% c("Melanoma", "Multiple myeloma"), "Melanoma", as.character(Label))
  ) %>%
  group_by(Block, Label, Condition) %>%
  summarise(Cell_Count = sum(Cell_Count, na.rm = TRUE), .groups = "drop")

# Order labels by total count
order_tbl <- df_agg %>% group_by(Block, Label) %>% summarise(Total = sum(Cell_Count), .groups = "drop")
df_agg <- df_agg %>%
  left_join(order_tbl, by = c("Block", "Label")) %>%
  mutate(Label = reorder(factor(Label), Total))

line_df_2 <- df_agg %>%
  group_by(Block) %>%
  summarise(y_min = 0.5, y_max = n_distinct(Label) + 0.5, .groups = "drop")

p_cell_count <- ggplot(df_agg, aes(x = Cell_Count, y = reorder(Label, Cell_Count), fill = Condition)) +
  geom_col() +
  facet_grid(rows = vars(Block), scales = "free_y", space = "free_y", 
             labeller = as_labeller(c("Cancer" = "Cancer Type", "Drug.type" = 'Drug Type', 'Tissue' = 'Tissue'))) +
  scale_fill_manual(values = c("Resistant" = "tomato", "Sensitive" = "skyblue")) +
  scale_x_continuous(labels = function(x) ifelse(x == 0, "0", paste0(x / 1000, "k"))) +  
  geom_segment(data = line_df_2, aes(x = 470000, xend = 470000, y = y_min, yend = y_max), 
               inherit.aes = FALSE, linewidth = 0.8, color = "black") +
  labs(x = "Number of Cells", y = NULL, fill = "Condition") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = -90, size = 24),
    panel.grid = element_blank(),
    axis.text = element_text(size = 24),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 24),
    panel.border = element_blank(),
    strip.background = element_blank(),
    legend.position = 'bottom'
  )

png(paste0(out_dir, 'cell_count_bar.png'), units = 'in', res = 300, width = 10, height = 12)
print(p_cell_count)
dev.off()

# ==============================================================================
# Part 3: Average Cosine Similarity & PCA
# ==============================================================================
avg_cos_df <- read.csv('../data/avg_cosine_similarity_SvsR.csv') %>%
  left_join(df, by = 'Dataset') %>%
  mutate(Origin = ifelse(Tissue %in% c('BMA', 'Tumor tissue', 'PBMC'), 'Patient', 'Cell line'))

p_avg_cos <- ggplot(avg_cos_df, aes(x = Origin, y = Mean, fill = Origin)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA) +
  geom_jitter(aes(color = Tissue), width = 0.15, size = 1.8, alpha = 0.8) +
  scale_fill_manual(values = c("#C0392B", "#2980B9")) + 
  scale_color_manual(values = c("green", "skyblue", "red", "orange")) + 
  geom_segment(aes(x = 1, xend = 2, y = 1.05, yend = 1.05), color = "black", linewidth = 0.6)  +
  geom_segment(aes(x = 1, xend = 1, y = 1.05, yend = 1.03), color = "black", linewidth = 0.6) + 
  geom_segment(aes(x = 2, xend = 2, y = 1.05, yend = 1.03), color = "black", linewidth = 0.6) +
  annotate("text", x = 1.5, y = 1.09, label = "p = 0.00053", size = 5.5) +
  labs(x = "", y = "Mean Cosine Similarity") +
  theme_classic() +
  theme(text = element_text(size = 20), legend.position = "right")

png(paste0(out_dir, 'avg_cos.png'), units = 'in', res = 300, width = 5, height = 5)
print(p_avg_cos)
dev.off()

# PCA
pca_input <- avg_cos_df[, c("Mean", "Median", "SD", "Min", "Max", "Q1", "Q3")]
pca_res <- prcomp(pca_input, scale. = TRUE)
pca_df <- data.frame(
  Origin = avg_cos_df$Origin,
  Tissue = avg_cos_df$Tissue,
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2]
)

p_avg_cos_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = Origin), geom = "polygon", alpha = 0.15, color = NA) +
  scale_color_manual(values = c("green", "skyblue", "red", "orange")) + 
  scale_fill_manual(values = c("#C0392B", "#2980B9")) + 
  labs(x = "PC1", y = "PC2") + 
  theme_classic() +
  theme(text = element_text(size = 20))

png(paste0(out_dir, 'avg_cos_pca.png'), units = 'in', res = 300, width = 5, height = 5)
print(p_avg_cos_pca)
dev.off()

# ==============================================================================
# Part 4: Lineage Tracing Benchmarking (Shuyu Data)
# ==============================================================================
meta_rd_ola <- read.csv('../data/meta_rd_ola.csv', row.names = 1)

# Plot 4: Stacked Barplot
df_bar <- meta_rd_ola %>%
  select(drugSens, time_point) %>%
  mutate(
    drugSens = recode(drugSens, "pre-resistant" = "Pre-resistant", "pre-sensitive" = "Pre-sensitive", "resistant" = "Resistant"),
    time_point = factor(recode(time_point, "AfterTreatment" = "After treatment", "BeforeTreatment" = "Before treatment"),
                        levels = c("Before treatment", "After treatment"))
  )

p_shuyu_stats <- ggplot(df_bar, aes(x = time_point, fill = drugSens)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("Pre-resistant" = my_colors[4], "Pre-sensitive" = my_colors[2], "Resistant" = my_colors[5])) +
  labs(x = "Time Point", y = "Number of Cells") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    legend.title = element_blank()
  )

png(paste0(out_dir, 'shuyu_data_stacked_barplot.png'), units = 'in', res = 300, width = 5, height = 6)
print(p_shuyu_stats)
dev.off()

# ==============================================================================
# Part 5: Sister vs Random Cells Cosine Similarity Comparison
# ==============================================================================
cosine_sim <- function(x, y) {
  sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

emb <- readRDS('../data/ola_emb_matrix.rds')
meta <- read.csv('../data/meta_rd_ola.csv', row.names = 1)
sister_groups <- split(rownames(meta), meta$sisters)

# Calculate similarity for sister pairs
sister_pairs <- lapply(sister_groups, function(cells) {
  if (length(cells) < 2) return(NULL)
  apply(combn(cells, 2), 2, function(p) cosine_sim(emb[p[1], ], emb[p[2], ]))
})
sister_sim <- unlist(sister_pairs)

# Calculate similarity for random pairs
set.seed(123)
all_cells <- rownames(meta)
n_random <- length(sister_sim)
random_sim <- numeric(n_random)

i <- 1
while (i <= n_random) {
  c1 <- sample(all_cells, 1)
  c2 <- sample(all_cells, 1)
  s1 <- meta[c1, "sisters"]
  s2 <- meta[c2, "sisters"]
  
  if (!is.na(s1) && !is.na(s2) && s1 != s2) {
    random_sim[i] <- cosine_sim(emb[c1,], emb[c2,])
    i <- i + 1
  }
}

df_sim <- data.frame(
  similarity = c(sister_sim, random_sim),
  type = c(rep("Sister pairs", length(sister_sim)), rep("Random pairs", length(random_sim)))
)

p_sisters <- ggplot(df_sim, aes(x = type, y = similarity, fill = type)) +
  geom_rain(point.args = list(color = "grey60"), boxplot.args = list(fill = "white"), violin.args = list(alpha = 0.5)) +
  annotate("segment", x = c(1, 1, 2), xend = c(2, 1, 2), y = c(0.985, 0.98, 0.98), yend = 0.985, linewidth = 0.8) +
  annotate("text", x = 1.5, y = 0.995, label = "p == 2.2 %*% 10^-16", size = 7, parse = TRUE) +
  scale_fill_manual(values = c(my_colors[9], my_colors[5])) +
  labs(x = "", y = "Cosine Similarity") +
  theme_classic() +
  theme(legend.position = "none", text = element_text(size = 24))

png(paste0(out_dir, 'sisters_compare.png'), units = 'in', res = 300, width = 5, height = 6)
print(p_sisters)
dev.off()