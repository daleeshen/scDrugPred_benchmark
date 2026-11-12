setwd('/media/shenli/Expansion/DRMREF/lineage_tracing/GSE150949/')

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(pROC)
library(Seurat)
library(tidytext)

watermelon_D0<-readRDS('./watermelon_D0.rds')
watermelon_D0<-subset(watermelon_D0, subset = condition != 'uncertain')
root <- "/media/shenli/Expansion/DRMREF/lineage_tracing/GSE150949/preR_vs_preS"
method_order <- c("BeyondCell", "DREEP", "DrugFormer", "Precily", "SCAD", "scDEAL", "scIDUC")

norm_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  case_when(
    x %in% c("1", "s", "sens", "sensitive", "pre-sensitive", "pre_sensitive", "pres", "pre-s") ~ "pre-sensitive",
    x %in% c("0", "r", "res", "resistant", "pre-resistant", "pre_resistant", "prer", "pre-r") ~ "pre-resistant",
    TRUE ~ NA_character_
  )
}

truth <- data.frame(cells = colnames(watermelon_D0), condition = norm_label(watermelon_D0$condition), stringsAsFactors = FALSE)
if (anyNA(truth$condition)) stop("watermelon_D0$condition contains labels other than pre-sensitive/pre-resistant.")
if (anyDuplicated(truth$cells)) stop("Duplicated cell IDs found in watermelon_D0.")

require_cols <- function(df, cols, file = "") {
  miss <- setdiff(cols, colnames(df))
  if (length(miss) > 0) stop("Missing columns in ", file, ": ", paste(miss, collapse = ", "))
}

guess_cell_ids <- function(df, file = "") {
  overlap <- sapply(df, function(x) sum(as.character(x) %in% truth$cells, na.rm = TRUE))
  if (length(overlap) > 0 && max(overlap) > 0) return(as.character(df[[which.max(overlap)]]))
  rn_overlap <- sum(rownames(df) %in% truth$cells)
  if (rn_overlap > 0) return(rownames(df))
  stop("Cannot identify cell-ID column in: ", file)
}

make_pred <- function(method, round, cells, score_raw, pred_condition, higher_resistant, config = NA_character_, source = NA_character_) {
  cells <- as.character(cells)
  score_raw <- as.numeric(score_raw)
  pred_condition <- norm_label(pred_condition)
  if (anyDuplicated(cells)) stop(method, " round ", round, ": duplicated cell IDs.")
  idx <- match(cells, truth$cells)
  if (anyNA(idx)) stop(method, " round ", round, ": ", sum(is.na(idx)), " cells cannot be found in watermelon_D0. Example: ", paste(head(cells[is.na(idx)], 5), collapse = ", "))
  data.frame(
    method = method, round = round, config = config, source = source,
    cells = cells, condition = truth$condition[idx],
    score_raw = score_raw,
    score_res = if (higher_resistant) score_raw else -score_raw,
    pred_condition = pred_condition,
    stringsAsFactors = FALSE
  )
}

safe_div <- function(x, y) ifelse(y == 0, NA_real_, x / y)

calc_metrics <- function(x) {
  meta_method <- x$method[1]
  meta_round <- x$round[1]
  meta_config <- x$config[1]
  auc_df <- x[!is.na(x$condition) & !is.na(x$score_res), ]
  cls_df <- x[!is.na(x$condition) & !is.na(x$pred_condition), ]
  if (length(unique(auc_df$condition)) == 2) {
    y_auc <- as.integer(auc_df$condition == "pre-resistant")
    roc_obj <- pROC::roc(y_auc, auc_df$score_res, levels = c(0, 1), direction = "<", quiet = TRUE)
    auroc <- as.numeric(pROC::auc(roc_obj))
  } else {
    auroc <- NA_real_
  }
  y <- cls_df$condition == "pre-resistant"
  p <- cls_df$pred_condition == "pre-resistant"
  tp <- sum(y & p)
  tn <- sum(!y & !p)
  fp <- sum(!y & p)
  fn <- sum(y & !p)
  precision <- safe_div(tp, tp + fp)
  recall <- safe_div(tp, tp + fn)
  f1 <- ifelse(is.na(precision) | is.na(recall) | precision + recall == 0, NA_real_, 2 * precision * recall / (precision + recall))
  data.frame(
    method = meta_method, round = meta_round, config = meta_config,
    n = nrow(cls_df), n_preS = sum(cls_df$condition == "pre-sensitive"), n_preR = sum(cls_df$condition == "pre-resistant"),
    AUROC = auroc, ACC = safe_div(tp + tn, tp + tn + fp + fn), F1 = f1, Precision = precision, Recall = recall,
    stringsAsFactors = FALSE
  )
}

selected_pred_list <- list()
selected_metric_list <- list()
scad_candidate_metric_list <- list()
scdeal_candidate_metric_list <- list()


#### BeyondCell: higher score = sensitive; use provided predicted condition ####
bc_root <- file.path(root, "BeyondCell")
for (r in 1:3) {
  f <- file.path(bc_root, paste0("watermelon_round", r, "_pred_res.csv"))
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  require_cols(d, c("cells", "beyondcell_osimertinib", "beyondcell_pred_condition"), f)
  p <- make_pred("BeyondCell", r, d$cells, d$beyondcell_osimertinib, d$beyondcell_pred_condition, higher_resistant = FALSE, source = f)
  selected_pred_list[[length(selected_pred_list) + 1]] <- p
  selected_metric_list[[length(selected_metric_list) + 1]] <- calc_metrics(p)
}

#### DREEP: >0 resistant, <=0 sensitive; higher score = resistant ####
dreep_root <- file.path(root, "DREEP")
for (r in 1:3) {
  f <- file.path(dreep_root, paste0("watermelon_round", r, "_pred_res.csv"))
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  score_col <- "GDSC_OSIMERTINIB:GDSC2:1919"
  require_cols(d, c("cells", score_col), f)
  score <- d[[score_col]]
  pred <- ifelse(score > 0, "pre-resistant", "pre-sensitive")
  p <- make_pred("DREEP", r, d$cells, score, pred, higher_resistant = TRUE, source = f)
  selected_pred_list[[length(selected_pred_list) + 1]] <- p
  selected_metric_list[[length(selected_metric_list) + 1]] <- calc_metrics(p)
}

#### DrugFormer: concatenate epoch 4 fold 0-4 OOF predictions; median threshold; higher score = resistant ####
drugformer_root <- file.path(root, "DrugFormer")
for (r in 1:3) {
  dir_candidates <- c(
    file.path(drugformer_root, paste0("round", r, "_pred_res")),
    file.path(drugformer_root, paste0("round", r), "pred_res"),
    file.path(drugformer_root, paste0("round", r))
  )
  round_dir <- dir_candidates[dir.exists(dir_candidates)][1]
  if (is.na(round_dir)) stop("Cannot find DrugFormer round ", r, " directory.")
  fs <- list.files(round_dir, pattern = "^epoch_4_fold_[0-4]_res\\.csv$", full.names = TRUE)
  if (length(fs) != 5) stop("DrugFormer round ", r, ": expected 5 epoch-4 fold files, found ", length(fs))
  d <- bind_rows(lapply(fs, function(f) {
    x <- fread(f, data.table = FALSE, check.names = FALSE)
    require_cols(x, c("sample_id", "pred_score"), f)
    x
  }))
  if (anyDuplicated(d$sample_id)) stop("DrugFormer round ", r, ": duplicated sample_id after concatenating 5 test folds.")
  med <- median(d$pred_score, na.rm = TRUE)
  pred <- ifelse(d$pred_score > med, "pre-resistant", "pre-sensitive")
  p <- make_pred("DrugFormer", r, d$sample_id, d$pred_score, pred, higher_resistant = TRUE, config = "epoch4_OOF", source = paste(fs, collapse = ";"))
  selected_pred_list[[length(selected_pred_list) + 1]] <- p
  selected_metric_list[[length(selected_metric_list) + 1]] <- calc_metrics(p)
}

#### Precily: higher score = resistant; use provided predicted condition ####
precily_root <- file.path(root, "Precily")
for (r in 1:3) {
  f <- file.path(precily_root, paste0("watermelon_round", r, "_pred_res.csv"))
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  require_cols(d, c("cells", "precily_osimertinib", "precily_pred_condition"), f)
  p <- make_pred("Precily", r, d$cells, d$precily_osimertinib, d$precily_pred_condition, higher_resistant = TRUE, source = f)
  selected_pred_list[[length(selected_pred_list) + 1]] <- p
  selected_metric_list[[length(selected_metric_list) + 1]] <- calc_metrics(p)
}

#### SCAD: norm/ppi/tp4k candidates; median threshold; higher score = sensitive; choose highest-AUROC candidate per round ####
scad_root <- file.path(root, "SCAD")
for (r in 1:3) {
  cand_pred <- list()
  cand_metric <- list()
  for (setting in c("norm", "ppi", "tp4k")) {
    search_root <- file.path(scad_root, paste0("round", r), paste0("best_param_results_", setting))
    fs <- list.files(search_root, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
    fs <- fs[grepl(paste0("Watermelon_round", r, "_Osimertinib"), fs, fixed = TRUE)]
    if (length(fs) != 1) stop("SCAD round ", r, " / ", setting, ": expected exactly one TSV under Watermelon_round", r, "_Osimertinib, found ", length(fs))
    f <- fs[1]
    d <- fread(f, header = FALSE, sep = "\t", data.table = FALSE)
    if (ncol(d) < 4) stop("SCAD file has fewer than 4 columns: ", f)
    cells <- as.character(d[[2]])
    score <- as.numeric(d[[4]])
    if (anyDuplicated(cells)) stop("SCAD round ", r, " / ", setting, ": duplicated cell IDs found.")
    med <- median(score, na.rm = TRUE)
    pred <- ifelse(score >= med, "pre-sensitive", "pre-resistant")
    p <- make_pred("SCAD", r, cells, score, pred, higher_resistant = FALSE, config = setting, source = f)
    m <- calc_metrics(p)
    cand_pred[[setting]] <- p
    cand_metric[[setting]] <- m
    scad_candidate_metric_list[[length(scad_candidate_metric_list) + 1]] <- m
  }
  cand_metric_df <- bind_rows(cand_metric)
  best_setting <- cand_metric_df$config[which.max(replace(cand_metric_df$AUROC, is.na(cand_metric_df$AUROC), -Inf))]
  p <- cand_pred[[best_setting]]
  selected_pred_list[[length(selected_pred_list) + 1]] <- p
  selected_metric_list[[length(selected_metric_list) + 1]] <- calc_metrics(p)
}

#### scDEAL: use sens_preds + sens_label; higher score = sensitive; choose highest-AUROC parameter combination per round ####
scdeal_root <- file.path(root, "scDEAL", "prediction_results")
for (r in 1:3) {
  fs <- list.files(scdeal_root, pattern = paste0("^Watermelon_round", r, ".*\\.csv$"), full.names = TRUE, ignore.case = TRUE)
  if (length(fs) == 0) stop("No scDEAL files found for round ", r)
  cand_pred <- list()
  cand_metric <- list()
  for (f in fs) {
    d <- fread(f, data.table = FALSE, check.names = FALSE)
    require_cols(d, "sens_preds", f)
    label_col <- intersect(c("sens_label", "sens_label_pret"), colnames(d))[1]
    if (is.na(label_col)) stop("Cannot find sens_label in: ", f)
    cells <- guess_cell_ids(d, f)
    score <- d$sens_preds
    pred <- ifelse(as.numeric(d[[label_col]]) == 1, "pre-sensitive", "pre-resistant")
    config <- tools::file_path_sans_ext(basename(f))
    p <- make_pred("scDEAL", r, cells, score, pred, higher_resistant = FALSE, config = config, source = f)
    m <- calc_metrics(p)
    cand_pred[[config]] <- p
    cand_metric[[config]] <- m
    scdeal_candidate_metric_list[[length(scdeal_candidate_metric_list) + 1]] <- m
  }
  cand_metric_df <- bind_rows(cand_metric)
  best_config <- cand_metric_df$config[which.max(replace(cand_metric_df$AUROC, is.na(cand_metric_df$AUROC), -Inf))]
  p <- cand_pred[[best_config]]
  selected_pred_list[[length(selected_pred_list) + 1]] <- p
  selected_metric_list[[length(selected_metric_list) + 1]] <- calc_metrics(p)
}

#### scIDUC: median threshold; higher score = resistant ####
sciduc_root <- file.path(root, "scIDUC")
for (r in 1:3) {
  f <- file.path(sciduc_root, paste0("scIDUC_res_round", r, ".csv"))
  d <- fread(f, data.table = FALSE, check.names = FALSE)
  require_cols(d, c("Cell_Name", "Osimertinib"), f)
  score <- d$Osimertinib
  med <- median(score, na.rm = TRUE)
  pred <- ifelse(score >= med, "pre-resistant", "pre-sensitive")
  p <- make_pred("scIDUC", r, d$Cell_Name, score, pred, higher_resistant = TRUE, source = f)
  selected_pred_list[[length(selected_pred_list) + 1]] <- p
  selected_metric_list[[length(selected_metric_list) + 1]] <- calc_metrics(p)
}

metrics_wide <- bind_rows(selected_metric_list) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  arrange(method, round)

metrics_long <- metrics_wide %>%
  pivot_longer(cols = c(AUROC, ACC, F1, Precision, Recall), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("AUROC", "ACC", "F1", "Precision", "Recall"))) %>%
  arrange(method, round, metric)

selected_predictions <- bind_rows(selected_pred_list) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  arrange(method, round, cells)

#### barplot #####
my_colors <- c("scDEAL"="#D55E00","DREEP"="#0072B2","CaDRReS-sc"="#009E73","SCAD"="#CC79A7","scDr"="#F0E442","Beyondcell"="#56B4E9","scIDUC"="#E69F00","Precily"="#999999","DrugFormer"="#8A2BE2")

df_sum <- metrics_long %>% mutate(Tool=as.character(method),Metric=as.character(metric),Tool=recode(Tool,"BeyondCell"="Beyondcell"),Metric=recode(Metric,"AUROC"="AUC")) %>% group_by(Metric,Tool) %>% summarise(mean=mean(value,na.rm=TRUE),se=sd(value,na.rm=TRUE)/sqrt(n()),.groups="drop") %>% mutate(Method_sorted=reorder_within(Tool,-mean,Metric))

auc_order <- df_sum %>% filter(Metric=="AUC") %>% arrange(desc(mean)) %>% pull(Tool)
df_sum$Tool <- factor(df_sum$Tool,levels=auc_order)
df_sum$Metric <- factor(df_sum$Metric,levels=c("ACC","AUC","F1","Recall","Precision"))

png(file.path("./preR_vs_preS/perform_figs/watermelon_perform_preS_vs_preR.png"),units="in",res=300,width=13,height=6)
ggplot(df_sum,aes(x=Method_sorted,y=mean,fill=Tool)) + 
  geom_col(width=0.7) + 
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=0.2) + 
  scale_x_reordered() + 
  facet_wrap(~Metric,scales="free_x",ncol=5) + 
  theme_bw(base_size=14) + 
  theme(text=element_text(size=20),
        strip.background=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1),
        legend.position="top",
        legend.title=element_blank()) + 
  scale_fill_manual(values=my_colors) + labs(y="Mean Value",x=NULL)
dev.off()

#### umap ####
fig_dir <- './preR_vs_preS/perform_figs/'

# For later UMAP: repeated pre-sensitive cells are averaged across rounds;
# each pre-resistant cell occurs in only one round and therefore keeps its own score.
umap_score_table <- selected_predictions %>%
  group_by(method, cells, condition) %>%
  summarise(
    score_raw = mean(score_raw, na.rm = TRUE),
    score_res = mean(score_res, na.rm = TRUE),
    n_rounds = n_distinct(round),
    .groups = "drop"
  ) %>%
  arrange(method, cells)

scad_candidate_metrics <- bind_rows(scad_candidate_metric_list)
scdeal_candidate_metrics <- bind_rows(scdeal_candidate_metric_list)

# keep only cells that were actually used in the 3 benchmark rounds
benchmark_cells <- unique(selected_predictions$cells)
watermelon_pres_prer_obj <- subset(watermelon_D0, cells = benchmark_cells)
watermelon_pres_prer_obj$drugSens <- factor(norm_label(watermelon_pres_prer_obj$condition), levels = c("pre-sensitive","pre-resistant"))

score_wide <- umap_score_table %>%
  mutate(meta_col = recode(method,
                           "BeyondCell" = "beyondcell_pred_score",
                           "DREEP" = "dreep_pred_score",
                           "DrugFormer" = "drugformer_pred_prob",
                           "Precily" = "precily_pred_score",
                           "SCAD" = "scad_pred_prob",
                           "scDEAL" = "scdeal_pred_prob",
                           "scIDUC" = "sciduc_pred_score")) %>%
  dplyr::select(cells, meta_col, score_raw) %>%
  pivot_wider(names_from = meta_col, values_from = score_raw)

meta <- watermelon_pres_prer_obj@meta.data
for (col in setdiff(colnames(score_wide), "cells")) {
  tmp <- setNames(score_wide[[col]], score_wide$cells)
  meta[[col]] <- tmp[rownames(meta)]
}
watermelon_pres_prer_obj@meta.data <- meta

my_umap_theme <- theme_bw() + theme(panel.grid = element_blank(), text = element_text(size = 22), legend.position = "none", legend.title = element_text(margin = margin(r = 10, t = -15)))

p_true <- DimPlot(watermelon_pres_prer_obj, group.by = "drugSens", pt.size = 1.5) +
  scale_color_manual(values = c("pre-sensitive" = "skyblue", "pre-resistant" = "tomato")) +
  my_umap_theme +
  labs(title = NULL, x = "UMAP1", y = "UMAP2") +
  annotate("text", x = -Inf, y = Inf, label = "True label", hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold")

png(file.path(fig_dir, "watermelon_preS_vs_preR_trueLabel_umap.png"), units = "in", res = 300, width = 6, height = 6)
print(p_true + theme(legend.position = "top", legend.text = element_text(size = 26)))
dev.off()

png(file.path(fig_dir, "watermelon_preS_vs_preR_trueLabel_umap_without_legend.png"), units = "in", res = 300, width = 6, height = 6)
print(p_true)
dev.off()

plot_configs <- list(
  "DREEP"      = list(feature = "dreep_pred_score", converse = FALSE, legend_name = "Score"),
  "scIDUC"     = list(feature = "sciduc_pred_score", converse = FALSE, legend_name = "Score"),
  "Precily"    = list(feature = "precily_pred_score", converse = FALSE, legend_name = "Score"),
  "Beyondcell" = list(feature = "beyondcell_pred_score", converse = FALSE,  legend_name = "Score"),
  "SCAD"       = list(feature = "scad_pred_prob",      converse = TRUE,  legend_name = "Score"),
  "scDEAL"     = list(feature = "scdeal_pred_prob",    converse = TRUE,  legend_name = "Sensitive probability"),
  "DrugFormer" = list(feature = "drugformer_pred_prob",converse = FALSE,legend_name = "Prediction scaled")
)

for (tool in names(plot_configs)) {
  cfg <- plot_configs[[tool]]
  high_col <- ifelse(cfg$converse, "skyblue", "tomato")
  low_col  <- ifelse(cfg$converse, "tomato", "skyblue")
  file_name <- paste0("watermelon_preS_vs_preR_", gsub("-", "_", tool), "_umap.png")
  png(file.path(fig_dir, file_name), units = "in", res = 300, width = 6, height = 6)
  p <- FeaturePlot(watermelon_pres_prer_obj, features = cfg$feature, pt.size = 1.5, order = FALSE) +
    scale_color_gradient(high = high_col, low = low_col, name = cfg$legend_name, na.value = "grey90") +
    my_umap_theme +
    labs(title = NULL, x = "UMAP1", y = "UMAP2") +
    annotate("text", x = -Inf, y = Inf, label = tool, hjust = -0.1, vjust = 1.5, size = 12, fontface = "bold") +
    guides(color = guide_colorbar(barwidth = 12, barheight = 1))
  print(p)
  dev.off()
}

message("All Watermelon preS vs preR UMAP figures generated successfully.")