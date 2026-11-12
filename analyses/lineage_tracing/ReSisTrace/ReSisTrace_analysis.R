# ============================================================
# ReSisTrace lineage-tracing analysis
#
# Includes:
#   1. Benchmark performance
#   2. True-label and prediction-score UMAPs
#   3. Pre-resistant vs pre-sensitive DEG analysis
#   4. Hallmark / KEGG / GO-MF GSEA
#   5. Intrinsic vs pre/post-treatment comparison
#
# Please change the paths below before running.
# ============================================================

library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(tidytext)
library(clusterProfiler)
library(GSVA)
library(ggrepel)
library(MAST)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("./")

data_dir <- "../data"
fig_dir <- "../results"
gmt_dir <- "/media/shenli/Expansion/GSEA_gmt_file"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "bioinformatics_results"), recursive = TRUE, showWarnings = FALSE)
bioinfo_dir <- file.path(fig_dir, "bioinformatics_results")

# ============================================================
# 1. Benchmark performance
# ============================================================

my_colors <- c(
  "scDEAL" = "#D55E00",
  "DREEP" = "#0072B2",
  "CaDRReS-sc" = "#009E73",
  "SCAD" = "#CC79A7",
  "scDr" = "#F0E442",
  "Beyondcell" = "#56B4E9",
  "scIDUC" = "#E69F00",
  "Precily" = "#999999",
  "DrugFormer" = "#8A2BE2"
)

final_df <- read.csv(file.path(data_dir, "preS_vs_preR_perform.csv"), header = TRUE)

final_df$Tool[final_df$Tool == "beyondcell"] <- "Beyondcell"
final_df$Tool[final_df$Tool == "precily"] <- "Precily"
final_df$Tool[final_df$Tool == "Cadrres-sc"] <- "CaDRReS-sc"

df_sum <- final_df %>%
  pivot_longer(cols = c(AUC, ACC, F1, Recall, Precision), names_to = "Metric", values_to = "Value") %>%
  group_by(Metric, Tool) %>%
  summarise(mean = mean(Value, na.rm = TRUE), se = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))), .groups = "drop") %>%
  mutate(Method_sorted = reorder_within(Tool, -mean, Metric))

auc_order <- df_sum %>% filter(Metric == "AUC") %>% arrange(desc(mean)) %>% pull(Tool)
df_sum$Tool <- factor(df_sum$Tool, levels = auc_order)

p_performance <- ggplot(df_sum, aes(x = Method_sorted, y = mean, fill = Tool)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  scale_x_reordered() +
  facet_wrap(~Metric, scales = "free_x", ncol = 5) +
  scale_fill_manual(values = my_colors) +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(size = 20),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(y = "Mean Value", x = NULL)

png(file.path(fig_dir, "ReSisTrace_preS_vs_preR_performance.png"), units = "in", res = 300, width = 13, height = 6)
print(p_performance)
dev.off()

# ============================================================
# 2. Load pre-sensitive / pre-resistant Seurat object
# ============================================================

load(file.path(data_dir, "ReSisTrace_pres_prer.RData"))

rna <- ola_pres_prer_obj[["RNA"]]
gene_names <- unname(rownames(rna@counts))
rownames(rna@counts) <- gene_names
rownames(rna@data) <- gene_names
rownames(rna@meta.features) <- gene_names
ola_pres_prer_obj@assays$RNA <- rna

stopifnot(
  identical(rownames(ola_pres_prer_obj[["RNA"]]@counts), rownames(ola_pres_prer_obj[["RNA"]]@data)),
  identical(rownames(ola_pres_prer_obj[["RNA"]]@data), rownames(ola_pres_prer_obj[["RNA"]]@meta.features))
)

ola_pres_prer_obj$condition <- ola_pres_prer_obj$drugSens
ola_pres_prer_obj$condition <- factor(ola_pres_prer_obj$condition, levels = c("pre-sensitive", "pre-resistant"))

# ============================================================
# 3. True-label UMAP
# ============================================================

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
  annotate("text", x = -Inf, y = Inf, label = "True label", hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold")

png(file.path(fig_dir, "preS_vs_preR_trueLabel_umap.png"), units = "in", res = 300, width = 6, height = 6)
print(p_true + theme(legend.position = "top", legend.text = element_text(size = 26)))
dev.off()

png(file.path(fig_dir, "preS_vs_preR_trueLabel_umap_without_legend.png"), units = "in", res = 300, width = 6, height = 6)
print(p_true)
dev.off()

# ============================================================
# 4. Prediction-score UMAPs
# ============================================================

plot_configs <- list(
  "scDr" = list(feature = "scdr_pred_score", converse = FALSE, legend_name = "Score"),
  "DREEP" = list(feature = "dreep_pred_score", converse = FALSE, legend_name = "Score"),
  "scIDUC" = list(feature = "sciduc_pred_score", converse = FALSE, legend_name = "Score"),
  "CaDRReS-sc" = list(feature = "cadrres_pred_score", converse = FALSE, legend_name = "Score"),
  "Precily" = list(feature = "precily_pred_score", converse = FALSE, legend_name = "Score"),
  "Beyondcell" = list(feature = "beyondcell_pred_score", converse = FALSE, legend_name = "Score"),
  "SCAD" = list(feature = "scad_pred_prob", converse = FALSE, legend_name = "Score"),
  "scDEAL" = list(feature = "scdeal_pred_prob", converse = TRUE, legend_name = "Sensitive probability"),
  "DrugFormer" = list(feature = "drugformer_pred_prob", converse = TRUE, legend_name = "Prediction scaled")
)

for (tool in names(plot_configs)) {
  cfg <- plot_configs[[tool]]
  if (!cfg$feature %in% colnames(ola_pres_prer_obj@meta.data)) {
    warning("Prediction column not found: ", cfg$feature)
    next
  }
  
  high_col <- ifelse(cfg$converse, "skyblue", "tomato")
  low_col <- ifelse(cfg$converse, "tomato", "skyblue")
  file_name <- paste0("preS_vs_preR_", gsub("-", "_", tool), "_umap.png")
  
  p <- FeaturePlot(ola_pres_prer_obj, features = cfg$feature, pt.size = 1.5, order = FALSE) +
    scale_color_gradient(high = high_col, low = low_col, name = cfg$legend_name) +
    my_umap_theme +
    labs(title = NULL, x = "UMAP1", y = "UMAP2") +
    annotate("text", x = -Inf, y = Inf, label = tool, hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold") +
    guides(color = guide_colorbar(barwidth = 12, barheight = 1))
  
  png(file.path(fig_dir, file_name), units = "in", res = 300, width = 6, height = 6)
  print(p)
  dev.off()
}

# ============================================================
# 5. Cell-cycle-adjusted pre-resistant vs pre-sensitive DEG
# ============================================================

s_genes <- intersect(cc.genes.updated.2019$s.genes, rownames(ola_pres_prer_obj))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(ola_pres_prer_obj))

ola_pres_prer_obj <- CellCycleScoring(
  ola_pres_prer_obj,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)

Idents(ola_pres_prer_obj) <- "condition"

deg_intrinsic <- FindMarkers(
  ola_pres_prer_obj,
  ident.1 = "pre-resistant",
  ident.2 = "pre-sensitive",
  test.use = "MAST",
  logfc.threshold = 0,
  min.pct = 0,
  latent.vars = c("S.Score", "G2M.Score")
)

deg_intrinsic <- deg_intrinsic %>% rownames_to_column("gene")
write.csv(deg_intrinsic, file.path(bioinfo_dir, "olaparib_preR_vs_preS_DEG.csv"), row.names = FALSE)

fc_col <- grep("^avg_log", colnames(deg_intrinsic), value = TRUE)[1]

cat("\nOlaparib pre-resistant vs pre-sensitive\n")
print(c(
  Up = sum(deg_intrinsic$p_val < 0.05 & deg_intrinsic[[fc_col]] > 0),
  Down = sum(deg_intrinsic$p_val < 0.05 & deg_intrinsic[[fc_col]] < 0),
  FC01_Up = sum(deg_intrinsic$p_val < 0.05 & deg_intrinsic[[fc_col]] > 0.1),
  FC01_Down = sum(deg_intrinsic$p_val < 0.05 & deg_intrinsic[[fc_col]] < -0.1),
  FDR = sum(deg_intrinsic$p_val_adj < 0.05)
))

# ============================================================
# 6. Volcano plot
# ============================================================

plot_volcano <- function(deg) {
  fc_col <- grep("^avg_log", colnames(deg), value = TRUE)[1]
  
  deg <- deg %>% mutate(
    status = case_when(
      p_val < 0.05 & .data[[fc_col]] > 0.1 ~ "Up (FC > 0.1)",
      p_val < 0.05 & .data[[fc_col]] > 0 & .data[[fc_col]] <= 0.1 ~ "Up (0 < FC ≤ 0.1)",
      p_val < 0.05 & .data[[fc_col]] < -0.1 ~ "Down (FC < -0.1)",
      p_val < 0.05 & .data[[fc_col]] < 0 & .data[[fc_col]] >= -0.1 ~ "Down (-0.1 ≤ FC < 0)",
      TRUE ~ "NS"
    ),
    label = ifelse(p_val < 0.05 & abs(.data[[fc_col]]) > 0.1, gene, NA)
  )
  
  ggplot(deg, aes(x = .data[[fc_col]], y = -log10(p_val), color = status)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-0.1, 0, 0.1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_text_repel(aes(label = label), max.overlaps = 20, size = 6, show.legend = FALSE) +
    theme_bw() +
    labs(x = "Average log2FC", y = "-log10(P)") +
    theme(
      text = element_text(size = 22),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
    ) +
    scale_color_manual(
      values = c(
        "Up (FC > 0.1)" = "#B2182B",
        "Up (0 < FC ≤ 0.1)" = "#EF8A8A",
        "NS" = "grey80",
        "Down (-0.1 ≤ FC < 0)" = "#91BFDB",
        "Down (FC < -0.1)" = "#2166AC"
      ),
      breaks = c(
        "Up (FC > 0.1)",
        "Up (0 < FC ≤ 0.1)",
        "NS",
        "Down (-0.1 ≤ FC < 0)",
        "Down (FC < -0.1)"
      )
    )
}

p_volcano <- plot_volcano(deg_intrinsic)

png(file.path(bioinfo_dir, "olaparib_preR_vs_preS_volcano.png"), units = "in", width = 6, height = 5, res = 300)
print(p_volcano)
dev.off()

# ============================================================
# 7. GSEA
# ============================================================

hallmark <- read.gmt(file.path(gmt_dir, "h.all.v2026.1.Hs.symbols.gmt"))
c2 <- read.gmt(file.path(gmt_dir, "c2.all.v2026.1.Hs.symbols.gmt"))
c5 <- read.gmt(file.path(gmt_dir, "c5.all.v2026.1.Hs.symbols.gmt"))

colnames(hallmark) <- c("pathway", "gene")
colnames(c2) <- c("pathway", "gene")
colnames(c5) <- c("pathway", "gene")

kegg <- c2 %>% filter(grepl("^KEGG_", pathway), !grepl("MEDICUS", pathway))
go_mf <- c5 %>% filter(grepl("^GOMF_", pathway))

run_gsea <- function(deg, term2gene) {
  fc_col <- grep("^avg_log", colnames(deg), value = TRUE)[1]
  gene_rank <- sign(deg[[fc_col]]) * -log10(deg$p_val + 1e-300)
  names(gene_rank) <- deg$gene
  gene_rank <- gene_rank[!is.na(gene_rank) & is.finite(gene_rank)]
  gene_rank <- gene_rank[!duplicated(names(gene_rank))]
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  as.data.frame(GSEA(geneList = gene_rank, TERM2GENE = term2gene, pvalueCutoff = 1, verbose = FALSE))
}

gsea_hallmark <- run_gsea(deg_intrinsic, hallmark)
gsea_kegg <- run_gsea(deg_intrinsic, kegg)
gsea_go_mf <- run_gsea(deg_intrinsic, go_mf)

write.csv(gsea_hallmark, file.path(bioinfo_dir, "olaparib_preR_vs_preS_Hallmark_GSEA.csv"), row.names = FALSE)
write.csv(gsea_kegg, file.path(bioinfo_dir, "olaparib_preR_vs_preS_KEGG_GSEA.csv"), row.names = FALSE)
write.csv(gsea_go_mf, file.path(bioinfo_dir, "olaparib_preR_vs_preS_GO_MF_GSEA.csv"), row.names = FALSE)

clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("^GOMF_", "", x)
  gsub("_", " ", x)
}

plot_gsea <- function(gsea, top_n = 15, title = "") {
  plot_df <- gsea %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway = clean_pathway(ID),
      direction = ifelse(NES > 0, "Pre-resistant", "Pre-sensitive"),
      pathway = reorder(pathway, NES)
    )
  
  ggplot(plot_df, aes(x = NES, y = pathway, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("Pre-resistant" = "#F4A6A6", "Pre-sensitive" = "#9EC5E6")) +
    theme_bw() +
    labs(x = "NES", y = NULL, fill = NULL, title = title) +
    theme(
      text = element_text(size = 18),
      panel.grid = element_blank(),
      plot.title = element_text(size = 22, face = "bold")
    )
}

p_gsea_hallmark <- plot_gsea(gsea_hallmark, top_n = 15, title = "Hallmark")
p_gsea_kegg <- plot_gsea(gsea_kegg, top_n = 15, title = "KEGG")
p_gsea_go_mf <- plot_gsea(gsea_go_mf, top_n = 20, title = "GO Molecular Function")

png(file.path(bioinfo_dir, "olaparib_preR_preS_gsea_hallmark.png"), units = "in", width = 10, height = 5, res = 300)
print(p_gsea_hallmark)
dev.off()

png(file.path(bioinfo_dir, "olaparib_preR_preS_gsea_kegg.png"), units = "in", width = 10, height = 5, res = 300)
print(p_gsea_kegg)
dev.off()

png(file.path(bioinfo_dir, "olaparib_preR_preS_gsea_GOMF.png"), units = "in", width = 12, height = 5, res = 300)
print(p_gsea_go_mf)
dev.off()

# ============================================================
# 8. Selected pathway summary plot
# ============================================================

hallmark_selected <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_UV_RESPONSE_UP"
)

kegg_selected <- c(
  "KEGG_NON_HOMOLOGOUS_END_JOINING",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_RIBOSOME",
  "KEGG_PROTEASOME",
  "KEGG_SPLICEOSOME"
)

go_mf_selected <- c(
  "GOMF_DNA_POLYMERASE_BINDING",
  "GOMF_FOUR_WAY_JUNCTION_DNA_BINDING",
  "GOMF_NADH_DEHYDROGENASE_UBIQUINONE_ACTIVITY",
  "GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME",
  "GOMF_ELECTRON_TRANSFER_ACTIVITY"
)

gsea_barplot_df <- bind_rows(
  gsea_hallmark %>% filter(ID %in% hallmark_selected) %>% mutate(database = "Hallmark"),
  gsea_kegg %>% filter(ID %in% kegg_selected) %>% mutate(database = "KEGG"),
  gsea_go_mf %>% filter(ID %in% go_mf_selected) %>% mutate(database = "GO_MF")
) %>%
  mutate(
    pathway = clean_pathway(ID),
    pathway_label = paste0(database, ": ", pathway),
    group = ifelse(NES > 0, "Pre-resistant", "Pre-sensitive"),
    score = abs(NES),
    neglogP = -log10(pvalue)
  ) %>%
  filter(pvalue < 0.05)

gsea_barplot_df$group <- factor(gsea_barplot_df$group, levels = c("Pre-resistant", "Pre-sensitive"))

gsea_barplot_df_preR <- gsea_barplot_df %>% filter(group == "Pre-resistant") %>% arrange(score)
gsea_barplot_df_preS <- gsea_barplot_df %>% filter(group == "Pre-sensitive") %>% arrange(score)
pathway_order <- c(gsea_barplot_df_preR$pathway_label, gsea_barplot_df_preS$pathway_label)
gsea_barplot_df$pathway_label <- factor(gsea_barplot_df$pathway_label, levels = pathway_order)

p_gsea_barplot <- ggplot(gsea_barplot_df, aes(x = score, y = pathway_label, fill = neglogP)) +
  geom_col(width = 0.8) +
  geom_text(aes(x = 0.03, label = pathway_label), hjust = 0, size = 5.2, lineheight = 0.9) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#F3E79B", high = "#C53E3E", name = expression(-LogP)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  labs(x = "Absolute NES", y = NULL) +
  theme(
    text = element_text(size = 18),
    panel.grid = element_blank(),
    strip.text = element_text(size = 22, face = "bold"),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

png(file.path(bioinfo_dir, "olaparib_preR_preS_gsea_all.png"), units = "in", width = 8, height = 9, res = 300)
print(p_gsea_barplot)
dev.off()

# ============================================================
# 9. Intrinsic transcriptional signal summary
# ============================================================

get_signal_summary <- function(deg, comparison) {
  fc_col <- grep("^avg_log", colnames(deg), value = TRUE)[1]
  data.frame(
    comparison = comparison,
    nominal_DEG = sum(deg$p_val < 0.05),
    FDR_DEG = sum(deg$p_val_adj < 0.05),
    FC01_DEG = sum(deg$p_val < 0.05 & abs(deg[[fc_col]]) > 0.1),
    median_abs_FC = median(abs(deg[[fc_col]]), na.rm = TRUE),
    mean_abs_FC = mean(abs(deg[[fc_col]]), na.rm = TRUE)
  )
}

intrinsic_summary <- get_signal_summary(deg_intrinsic, "Intrinsic")
write.csv(intrinsic_summary, file.path(bioinfo_dir, "olaparib_intrinsic_signal_summary.csv"), row.names = FALSE)

# ============================================================
# 10. Pre/post-treatment comparison
# ============================================================

load(file.path(data_dir, "figure.4A_B.RData"))
ola_obj$condition <- ola_obj$drugSens

gene_info <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(ola_obj),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
)

s_ensg <- gene_info$ENSEMBL[gene_info$SYMBOL %in% cc.genes.updated.2019$s.genes]
g2m_ensg <- gene_info$ENSEMBL[gene_info$SYMBOL %in% cc.genes.updated.2019$g2m.genes]
s_ensg <- intersect(s_ensg, rownames(ola_obj))
g2m_ensg <- intersect(g2m_ensg, rownames(ola_obj))

ola_obj <- CellCycleScoring(
  ola_obj,
  s.features = s_ensg,
  g2m.features = g2m_ensg,
  set.ident = FALSE
)

Idents(ola_obj) <- "condition"

deg_traditional <- FindMarkers(
  ola_obj,
  ident.1 = "resistant",
  ident.2 = "sensitive",
  test.use = "MAST",
  logfc.threshold = 0,
  min.pct = 0,
  latent.vars = c("S.Score", "G2M.Score")
)

deg_traditional <- deg_traditional %>% rownames_to_column("gene")
write.csv(deg_traditional, file.path(bioinfo_dir, "olaparib_pre_post_DEG.csv"), row.names = FALSE)

traditional_summary <- get_signal_summary(deg_traditional, "Pre/Post-treatment")
signal_summary <- bind_rows(intrinsic_summary, traditional_summary)
write.csv(signal_summary, file.path(bioinfo_dir, "intrinsic_vs_pre_post_summary.csv"), row.names = FALSE)

# ============================================================
# 11. Intrinsic vs pre/post DEG summary
# ============================================================

signal_plot_df <- signal_summary %>%
  select(comparison, FC01_DEG, median_abs_FC) %>%
  pivot_longer(-comparison, names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(
    metric,
    FC01_DEG = "Number of DEGs",
    median_abs_FC = "Median |log2FC|"
  ))

p_signal <- ggplot(signal_plot_df, aes(x = comparison, y = value, fill = comparison)) +
  geom_col(width = 0.7) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c("Intrinsic" = "#4DBBD5", "Pre/Post-treatment" = "#E64B35"),
    name = NULL
  ) +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20, face = "bold"),
    strip.background = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 14)
  )

png(file.path(bioinfo_dir, "intrinsic_vs_pre_post_bar.png"), units = "in", res = 300, width = 6, height = 4)
print(p_signal)
dev.off()

# ============================================================
# 12. Intrinsic vs pre/post |log2FC| distribution
# ============================================================

get_fc_distribution <- function(deg, comparison) {
  fc_col <- grep("^avg_log", colnames(deg), value = TRUE)[1]
  data.frame(comparison = comparison, abs_log2FC = abs(deg[[fc_col]]))
}

fc_distribution <- bind_rows(
  get_fc_distribution(deg_intrinsic, "Intrinsic"),
  get_fc_distribution(deg_traditional, "Pre/Post")
)

fc_distribution$comparison <- factor(fc_distribution$comparison, levels = c("Intrinsic", "Pre/Post"))

fc_plot <- fc_distribution %>%
  group_by(comparison) %>%
  mutate(
    cutoff = quantile(abs_log2FC, 0.9, na.rm = TRUE),
    abs_log2FC_plot = pmin(abs_log2FC, cutoff)
  ) %>%
  ungroup()

stat_test <- wilcox.test(abs_log2FC ~ comparison, data = fc_distribution)
p_value <- stat_test$p.value

sig_label <- case_when(
  p_value < 0.0001 ~ "****",
  p_value < 0.001 ~ "***",
  p_value < 0.01 ~ "**",
  p_value < 0.05 ~ "*",
  TRUE ~ "ns"
)

y_max <- max(fc_plot$abs_log2FC_plot, na.rm = TRUE)
y_sig <- y_max * 1.08
h_sig <- y_max * 0.025

p_fc_distribution <- ggplot(fc_plot, aes(x = comparison, y = abs_log2FC_plot, fill = comparison)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.85, alpha = 0.9, linewidth = 0.4) +
  geom_boxplot(width = 0.10, outlier.shape = NA, fill = "white", linewidth = 0.45) +
  scale_fill_manual(values = c("Intrinsic" = "#4DBBD5", "Pre/Post" = "#E64B35")) +
  scale_y_continuous(limits = c(0, y_max * 1.18), expand = expansion(mult = c(0, 0))) +
  annotate("segment", x = 1.10, xend = 1.90, y = y_sig, yend = y_sig, linewidth = 0.5) +
  annotate("segment", x = 1.10, xend = 1.10, y = y_sig, yend = y_sig - h_sig, linewidth = 0.5) +
  annotate("segment", x = 1.90, xend = 1.90, y = y_sig, yend = y_sig - h_sig, linewidth = 0.5) +
  annotate("text", x = 1.5, y = y_sig + y_max * 0.035, label = sig_label, size = 5.5) +
  theme_classic() +
  labs(x = NULL, y = "|log2FC|") +
  theme(
    text = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 17, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = 0.6),
    legend.position = "none"
  )

png(file.path(bioinfo_dir, "intrinsic_vs_pre_post_fc_violin.png"), units = "in", res = 300, width = 6, height = 4)
print(p_fc_distribution)
dev.off()

message("ReSisTrace analysis completed successfully.")