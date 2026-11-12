library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(clusterProfiler)
library(GSVA)
library(ggrepel)
library(MAST)
library(harmony)

setwd('/media/shenli/Expansion/DRMREF/lineage_tracing/GSE150949')

pc9_count <- read.csv("GSE150949_pc9_count_matrix.csv", row.names = 1, check.names = FALSE)
pc9_meta <- read.delim("GSE150949_metaData_with_lineage.txt", header = TRUE, row.names = 1, check.names = FALSE)
pc9_meta <- pc9_meta[colnames(pc9_count), , drop = FALSE]
pc9_count <- Matrix::Matrix(as.matrix(pc9_count), sparse = TRUE)

watermelon <- CreateSeuratObject(counts = pc9_count, meta.data = pc9_meta, project = "Watermelon_PC9", min.cells = 0, min.features = 0)
watermelon$time_label <- factor(watermelon$time_point, levels = c(0, 3, 7, 14), labels = c("D0", "D3", "D7", "D14"))

build_watermelon_strict_all <- function(obj) {
  meta <- obj@meta.data %>% rownames_to_column("cell_id") %>% mutate(n_barcode = ifelse(is.na(lineage_barcode) | lineage_barcode == "", 0, str_count(lineage_barcode, ",") + 1))
  
  single_meta <- meta %>% filter(n_barcode == 1)
  
  barcode_history <- single_meta %>%
    group_by(lineage_barcode) %>%
    summarise(
      D0 = any(time_point == 0),
      D3 = any(time_point == 3),
      D7 = any(time_point == 7),
      D14 = any(time_point == 14),
      .groups = "drop"
    )
  
  preR_barcodes <- barcode_history %>% filter(D0 & D3 & D7 & D14) %>% pull(lineage_barcode)
  preS_barcodes <- barcode_history %>% filter(D0 & !D3 & !D7 & !D14) %>% pull(lineage_barcode)
  
  meta <- meta %>% mutate(condition = case_when(
    n_barcode != 1 ~ "uncertain",
    time_point == 0 & lineage_barcode %in% preR_barcodes ~ "pre-resistant",
    time_point == 0 & lineage_barcode %in% preS_barcodes ~ "pre-sensitive",
    time_point %in% c(3, 7, 14) & lineage_barcode %in% preR_barcodes ~ "resistant",
    TRUE ~ "uncertain"
  ))
  
  rownames(meta) <- meta$cell_id
  obj$condition <- meta[colnames(obj), "condition"]
  obj$condition <- factor(obj$condition, levels = c("pre-sensitive", "pre-resistant", "resistant", "uncertain"))
  return(obj)
}

watermelon <- build_watermelon_strict_all(watermelon)

watermelon_D0 <- subset(watermelon, subset = time_point == 0)
watermelon_D3 <- subset(watermelon, subset = time_point == 3)
watermelon_D7 <- subset(watermelon, subset = time_point == 7)
watermelon_D14 <- subset(watermelon, subset = time_point == 14)

table(watermelon$time_point, watermelon$condition)

watermelon <- NormalizeData(watermelon, verbose = FALSE)
watermelon_D0 <- NormalizeData(watermelon_D0, verbose = FALSE)
watermelon_D3 <- NormalizeData(watermelon_D3, verbose = FALSE)
watermelon_D7 <- NormalizeData(watermelon_D7, verbose = FALSE)
watermelon_D14 <- NormalizeData(watermelon_D14, verbose = FALSE)

s_genes <- intersect(cc.genes.updated.2019$s.genes, rownames(watermelon))
g2m_genes <- intersect(cc.genes.updated.2019$g2m.genes, rownames(watermelon))

watermelon <- CellCycleScoring(watermelon, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
watermelon_D0 <- CellCycleScoring(watermelon_D0, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
watermelon_D3 <- CellCycleScoring(watermelon_D3, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
watermelon_D7 <- CellCycleScoring(watermelon_D7, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
watermelon_D14 <- CellCycleScoring(watermelon_D14, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

watermelon_D3$batch <- ifelse(grepl("rep1", watermelon_D3$sample_name), "rep1", "rep2")
watermelon_D7$batch <- ifelse(grepl("rep1", watermelon_D7$sample_name), "rep1", "rep2")
watermelon_D14$batch <- ifelse(grepl("rep1", watermelon_D14$sample_name), "rep1", "rep2")

dir.create("bioinformatics_results", showWarnings = FALSE)

# Step 0: UMAP analysis
run_umap <- function(obj) {
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  return(obj)
}

run_umap_harmony <- function(obj) {
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "batch", reduction.use = "pca", dims.use = 1:20, verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
  return(obj)
}

watermelon <- run_umap(watermelon)
watermelon_D0 <- run_umap(watermelon_D0)
watermelon_D3 <- run_umap_harmony(watermelon_D3)
watermelon_D7 <- run_umap_harmony(watermelon_D7)
watermelon_D14 <- run_umap_harmony(watermelon_D14)

p_umap_time <- DimPlot(watermelon, reduction = "umap", group.by = "time_label", pt.size = 0.3) +
  scale_color_manual(values = c("D0" = "#4DBBD5", "D3" = "#00A087", "D7" = "#E64B35", "D14" = "#7E6148")) +
  theme_bw() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(text = element_text(size = 22),
        panel.grid = element_blank(),
        plot.title = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0))

condition_colors_D0 <- c("pre-sensitive" = "#4DBBD5", "pre-resistant" = "#E64B35", "uncertain" = "grey75")
condition_colors_post <- c("resistant" = "#B2182B", "uncertain" = "grey75")

p_umap_D0 <- DimPlot(subset(watermelon_D0, subset = condition != 'uncertain'), reduction = "umap", group.by = "condition", pt.size = 1) +
  scale_color_manual(values = condition_colors_D0, drop = TRUE) +
  annotate("text", x = -Inf, y = -Inf, label = "D0", hjust = -0.2, vjust = -0.8, size = 8, fontface = "bold") +
  theme_bw() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(text = element_text(size = 22),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.title = element_blank())

plot_post_umap <- function(obj, title) {
  DimPlot(obj, reduction = "umap", group.by = "condition", pt.size = 1) +
    scale_color_manual(values = condition_colors_post, drop = TRUE) +
    annotate("text", x = -Inf, y = -Inf, label = title, hjust = -0.2, vjust = -0.8, size = 8, fontface = "bold") +
    theme_bw() +
    labs(x = "UMAP1", y = "UMAP2") +
    theme(text = element_text(size = 22),
          panel.grid = element_blank(),
          legend.position = "none",
          plot.title = element_blank())
}

p_umap_D3 <- plot_post_umap(subset(watermelon_D3, subset = condition != 'uncertain'), "D3")
p_umap_D7 <- plot_post_umap(subset(watermelon_D7, subset = condition != 'uncertain'), "D7")
p_umap_D14 <- plot_post_umap(subset(watermelon_D14, subset = condition != 'uncertain'), "D14")

png('./bioinformatics_results/UMAP_time_point.png', units = 'in', width = 5, height = 5, res = 300)
p_umap_time
dev.off()

png('./bioinformatics_results/UMAP_condition_D0.png', units = 'in', width = 5, height = 5, res = 300)
p_umap_D0
dev.off()

png('./bioinformatics_results/UMAP_condition_D3.png', units = 'in', width = 5, height = 5, res = 300)
p_umap_D3
dev.off()

png('./bioinformatics_results/UMAP_condition_D7.png', units = 'in', width = 5, height = 5, res = 300)
p_umap_D7
dev.off()

png('./bioinformatics_results/UMAP_condition_D14.png', units = 'in', width = 5, height = 5, res = 300)
p_umap_D14
dev.off()

# Step 1: cell-cycle-adjusted pre-resistant vs pre-sensitive DEG
run_deg <- function(obj, group1, group2) {
  Idents(obj) <- "condition"
  deg <- FindMarkers(obj, ident.1 = group1, ident.2 = group2, test.use = "MAST", logfc.threshold = 0, min.pct = 0, latent.vars = c("S.Score", "G2M.Score"))
  deg <- deg %>% rownames_to_column("gene")
  return(deg)
}

deg_intrinsic <- run_deg(watermelon_D0, "pre-resistant", "pre-sensitive")

write.csv(deg_intrinsic, "bioinformatics_results/preR_vs_preS_DEG.csv", row.names = FALSE)

fc_col <- grep("^avg_log", colnames(deg_intrinsic), value = TRUE)[1]

cat("\npre-resistant vs pre-sensitive\n")
print(c(
  Up = sum(deg_intrinsic$p_val < 0.05 & deg_intrinsic[[fc_col]] > 0),
  Down = sum(deg_intrinsic$p_val < 0.05 & deg_intrinsic[[fc_col]] < 0),
  FDR = sum(deg_intrinsic$p_val_adj < 0.05)
))

plot_volcano <- function(deg, title) {
  fc_col <- grep("^avg_log", colnames(deg), value = TRUE)[1]
  deg <- deg %>% mutate(
    status = case_when(
      p_val < 0.05 & .data[[fc_col]] > 0.1 ~ "Up (FC > 0.1)",
      p_val < 0.05 & .data[[fc_col]] > 0 & .data[[fc_col]] <= 0.1 ~ "Up (0 < FC ≤ 0.1)",
      p_val < 0.05 & .data[[fc_col]] < -0.1 ~ "Down (FC < -0.1)",
      p_val < 0.05 & .data[[fc_col]] < 0 & .data[[fc_col]] >= -0.1 ~ "Down (-0.1 ≤ FC < 0)",
      TRUE ~ "NS"
    ),
    label = ifelse(p_val_adj < 0.05, gene, NA)
  )
  
  ggplot(deg, aes(x = .data[[fc_col]], y = -log10(p_val), color = status)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(xintercept = c(-0.1, 0, 0.1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_text_repel(aes(label = label), max.overlaps = 20, size = 6, show.legend = FALSE) +
    #annotate("text", x = -Inf, y = -Inf, label = title, hjust = -0.2, vjust = -0.8, size = 8, fontface = "bold") +
    theme_bw() +
    labs(x = "Average log2FC", y = "-log10(P)") +
    theme(text = element_text(size = 22),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.98, 0.98),
          legend.justification = c(1, 1),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
    ) +
    scale_color_manual(values = c(
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
    ))
}

p_intrinsic <- plot_volcano(deg_intrinsic, "Pre-resistant vs pre-sensitive")

png('./bioinformatics_results/preR_vs_preS_volcano.png', units = 'in', width = 6, height = 5, res = 300)
p_intrinsic
dev.off()

# Step 2: cell-cycle-adjusted GSEA analysis
gmt_dir <- "/media/shenli/Expansion/GSEA_gmt_file"

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
  
  gsea <- GSEA(geneList = gene_rank, TERM2GENE = term2gene, pvalueCutoff = 1, verbose = FALSE)
  gsea <- as.data.frame(gsea)
  return(gsea)
}

gsea_hallmark <- run_gsea(deg_intrinsic, hallmark)
gsea_kegg <- run_gsea(deg_intrinsic, kegg)
gsea_go_mf <- run_gsea(deg_intrinsic, go_mf)

write.csv(gsea_hallmark, "bioinformatics_results/preR_vs_preS_Hallmark_GSEA.csv", row.names = FALSE)
write.csv(gsea_kegg, "bioinformatics_results/preR_vs_preS_KEGG_GSEA.csv", row.names = FALSE)
write.csv(gsea_go_mf, "bioinformatics_results/preR_vs_preS_GO_MF_GSEA.csv", row.names = FALSE)

clean_pathway <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^KEGG_", "", x)
  x <- gsub("^GOMF_", "", x)
  x <- gsub("_", " ", x)
  return(x)
}

plot_gsea <- function(gsea, top_n = 15, title = "") {
  plot_df <- gsea %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = top_n) %>%
    mutate(pathway = clean_pathway(ID),
           direction = ifelse(NES > 0, "Pre-resistant", "Pre-sensitive"),
           pathway = reorder(pathway, NES))
  
  ggplot(plot_df, aes(x = NES, y = pathway, fill = direction)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_fill_manual(values = c("Pre-resistant" = "#F4A6A6", "Pre-sensitive" = "#9EC5E6")) +
    theme_bw() +
    labs(x = "NES", y = NULL, fill = NULL, title = title) +
    theme(text = element_text(size = 18),
          panel.grid = element_blank(),
          plot.title = element_text(size = 22, face = "bold"))
}

p_gsea_hallmark <- plot_gsea(gsea_hallmark, top_n = 15, title = "Hallmark")
p_gsea_kegg <- plot_gsea(gsea_kegg, top_n = 15, title = "KEGG")
p_gsea_go_mf <- plot_gsea(gsea_go_mf, top_n = 20, title = "GO Molecular Function")

png('./bioinformatics_results/preR_preS_gsea_hallmark.png', units = 'in', width = 10, height = 5, res = 300)
p_gsea_hallmark
dev.off()
png('./bioinformatics_results/preR_preS_gsea_kegg.png', units = 'in', width = 10, height = 5, res = 300)
p_gsea_kegg
dev.off()
png('./bioinformatics_results/preR_preS_gsea_GOMF.png', units = 'in', width = 10, height = 5, res = 300)
p_gsea_go_mf
dev.off()

# Selected pathways for one-sided barplot
hallmark_selected <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_HYPOXIA",
  "HALLMARK_HEME_METABOLISM",
  "HALLMARK_COMPLEMENT"
)

kegg_selected <- c(
  "KEGG_HOMOLOGOUS_RECOMBINATION",
  "KEGG_OXIDATIVE_PHOSPHORYLATION",
  "KEGG_SPLICEOSOME",
  "KEGG_VEGF_SIGNALING_PATHWAY",
  "KEGG_RIBOSOME",
  "KEGG_LYSOSOME",
  "KEGG_GLUTATHIONE_METABOLISM",
  "KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450"
)

go_mf_selected <- c(
  "GOMF_UBIQUITIN_LIKE_PROTEIN_CONJUGATING_ENZYME_ACTIVITY",
  "GOMF_HISTONE_READER_ACTIVITY",
  "GOMF_CHROMATIN_PROTEIN_ADAPTOR_ACTIVITY",
  "GOMF_NUCLEOSIDE_DIPHOSPHATE_KINASE_ACTIVITY",
  "GOMF_IRON_ION_BINDING",
  "GOMF_LIPOPROTEIN_PARTICLE_RECEPTOR_ACTIVITY"
)

gsea_barplot_df <- bind_rows(
  gsea_hallmark %>% filter(ID %in% hallmark_selected) %>% mutate(database = "Hallmark"),
  gsea_kegg %>% filter(ID %in% kegg_selected) %>% mutate(database = "KEGG"),
  gsea_go_mf %>% filter(ID %in% go_mf_selected) %>% mutate(database = "GO_MF")
) %>%
  mutate(
    pathway = clean_pathway(ID),
    pathway_label = paste0(database, ": ", pathway),
    #pathway_label = str_wrap(pathway_label, width = 45),
    group = ifelse(NES > 0, "Pre-resistant", "Pre-sensitive"),
    score = abs(NES),
    neglogP = -log10(pvalue)
  ) %>%
  filter(pvalue < 0.05)

gsea_barplot_df$group <- factor(gsea_barplot_df$group, levels = c("Pre-resistant", "Pre-sensitive"))

gsea_barplot_df_preR <- gsea_barplot_df %>%
  filter(group == "Pre-resistant") %>%
  arrange(score)

gsea_barplot_df_preS <- gsea_barplot_df %>%
  filter(group == "Pre-sensitive") %>%
  arrange(score)

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

png('./bioinformatics_results/preR_preS_gsea_all.png', units = 'in', width = 8, height = 7, res = 300)
p_gsea_barplot
dev.off()

# Step 3: intrinsic vs treatment-induced signal strength
run_treatment_deg <- function(obj, endpoint) {
  obj_sub <- subset(obj, subset = (time_point == 0 & condition == "pre-resistant") | (time_point == endpoint & condition == "resistant"))
  obj_sub$comparison <- ifelse(obj_sub$time_point == 0, "pre-resistant", "resistant")
  Idents(obj_sub) <- "comparison"
  deg <- FindMarkers(obj_sub, ident.1 = "resistant", ident.2 = "pre-resistant", test.use = "MAST", logfc.threshold = 0, min.pct = 0, latent.vars = c("S.Score", "G2M.Score"))
  deg <- deg %>% rownames_to_column("gene")
  return(deg)
}

deg_treatment_D3 <- run_treatment_deg(watermelon, 3)
deg_treatment_D7 <- run_treatment_deg(watermelon, 7)
deg_treatment_D14 <- run_treatment_deg(watermelon, 14)

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

signal_summary <- bind_rows(
  get_signal_summary(deg_intrinsic, "Intrinsic"),
  get_signal_summary(deg_treatment_D3, "D3 treatment"),
  get_signal_summary(deg_treatment_D7, "D7 treatment"),
  get_signal_summary(deg_treatment_D14, "D14 treatment")
)

signal_summary$comparison <- factor(signal_summary$comparison, levels = c("Intrinsic", "D3 treatment", "D7 treatment", "D14 treatment"))

print(signal_summary)

signal_plot_df <- signal_summary %>%
  select(comparison, FC01_DEG, median_abs_FC) %>%
  pivot_longer(-comparison, names_to = "metric", values_to = "value") %>%
  mutate(
    comparison = recode(comparison,
                        "Intrinsic" = "Intrinsic",
                        "D3 treatment" = "D3",
                        "D7 treatment" = "D7",
                        "D14 treatment" = "D14"),
    metric = recode(metric,
                    FC01_DEG = "Number of DEGs",
                    median_abs_FC = "Median |log2FC|")
  )

p_signal <- ggplot(signal_plot_df, aes(x = comparison, y = value, fill = comparison)) +
  geom_col(width = 0.7) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c("Intrinsic" = "#4DBBD5", "D3" = "#F39B7F", "D7" = "#E64B35", "D14" = "#B2182B"),
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

png('./bioinformatics_results/intrinsic_vs_pre_post_bar.png', units = 'in', res = 300, width = 6, height = 4)
p_signal
dev.off()

get_fc_distribution <- function(deg, comparison) {
  fc_col <- grep("^avg_log", colnames(deg), value = TRUE)[1]
  data.frame(
    comparison = comparison,
    abs_log2FC = abs(deg[[fc_col]])
  )
}

fc_distribution <- bind_rows(
  get_fc_distribution(deg_intrinsic, "Intrinsic"),
  get_fc_distribution(deg_treatment_D3, "D3 treatment"),
  get_fc_distribution(deg_treatment_D7, "D7 treatment"),
  get_fc_distribution(deg_treatment_D14, "D14 treatment")
)

fc_distribution$comparison <- factor(fc_distribution$comparison, levels = c("Intrinsic", "D3 treatment", "D7 treatment", "D14 treatment"))

fc_ylim <- quantile(fc_distribution$abs_log2FC, 0.95, na.rm = TRUE)

p_fc_distribution <- ggplot(fc_distribution, aes(x = comparison, y = abs_log2FC, fill = comparison)) +
  geom_violin(trim = TRUE, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("Intrinsic" = "#4DBBD5", "D3 treatment" = "#F39B7F", "D7 treatment" = "#E64B35", "D14 treatment" = "#B2182B")) +
  coord_cartesian(ylim = c(0, fc_ylim)) +
  theme_classic() +
  labs(x = NULL, y = "|log2FC|") +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )

png('./bioinformatics_results/intrinsic_vs_pre_post_FC_violin.png', units = 'in', res = 300, width = 6, height = 4)
p_fc_distribution
dev.off()

save.image('./watermelon_analysis.RData')
