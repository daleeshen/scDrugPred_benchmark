library(tidyverse)
library(caret)
library(pROC)
library(ggpubr)
library(patchwork)

# ============================================================
# Paths
# ============================================================

out_dir <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5C/summary_metrics"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# scDEAL original benchmark results
scdeal_bal_original_file <-
  "/media/shenli/Expansion/DRMREF/1_round_1000cells_res/1000cells_scDEAL.csv"

scdeal_imb_original_file <-
  "/media/shenli/Expansion/DRMREF/unb_scDEAL_best_performance.csv"

# scDEAL replaced-label predictions
scdeal_bal_replaced_root <-
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5C/scDEAL_balanced"

scdeal_imb_replaced_root <-
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5C/scDEAL_imbalanced"

# SCAD results
scad_root <-
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5C/SCAD_fixed_original_params_results"

# SCAD dataset-drug mapping
scad_pair_file <-
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"


# ============================================================
# Helper functions
# ============================================================

pick_col <- function(df, candidates, required = TRUE) {
  nm <- names(df)
  hit <- which(tolower(nm) %in% tolower(candidates))
  if (length(hit) == 0) {
    if (required) stop("Cannot find column: ", paste(candidates, collapse = "/"))
    return(NULL)
  }
  nm[hit[1]]
}

norm_round <- function(x) {
  as.integer(stringr::str_extract(as.character(x), "\\d+"))
}

norm_ratio <- function(x) {
  x <- gsub(":", "_", as.character(x))
  x <- gsub("^1_", "1_", x)
  case_when(
    str_detect(x, "1_3$")  ~ "1_3",
    str_detect(x, "1_30$") ~ "1_30",
    TRUE ~ x
  )
}

norm_sampling <- function(x) {
  case_when(
    tolower(x) == "smote"        ~ "SMOTE",
    tolower(x) == "upsampling"   ~ "upsampling",
    tolower(x) == "downsampling" ~ "downsampling",
    TRUE ~ as.character(x)
  )
}

true_to_num <- function(x) {
  xx <- tolower(trimws(as.character(x)))
  case_when(
    xx %in% c("sensitive", "s", "1") ~ 1L,
    xx %in% c("resistant", "r", "0") ~ 0L,
    TRUE ~ suppressWarnings(as.integer(xx))
  )
}

# sklearn-compatible average precision
average_precision_sklearn <- function(y, score) {
  keep <- complete.cases(y, score)
  y <- as.integer(y[keep])
  score <- as.numeric(score[keep])
  
  P <- sum(y == 1)
  if (P == 0) return(NA_real_)
  
  ord <- order(score, decreasing = TRUE)
  y <- y[ord]
  score <- score[ord]
  
  group_end <- c(which(diff(score) != 0), length(score))
  
  tp <- cumsum(y == 1)[group_end]
  fp <- group_end - tp
  
  recall <- tp / P
  precision <- tp / (tp + fp)
  
  delta_recall <- c(recall[1], diff(recall))
  sum(delta_recall * precision)
}

weighted_f1 <- function(true, pred) {
  lv <- c("0", "1")
  true <- factor(as.character(true), levels = lv)
  pred <- factor(as.character(pred), levels = lv)
  
  tab <- table(true, pred)
  support <- rowSums(tab)
  
  f1 <- sapply(seq_along(lv), function(i) {
    tp <- tab[i, i]
    fp <- sum(tab[-i, i])
    fn <- sum(tab[i, -i])
    
    precision <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
    recall <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
    
    if (precision + recall == 0) 0 else 2 * precision * recall / (precision + recall)
  })
  
  sum(f1 * support) / sum(support)
}

calc_metrics <- function(true, score, pred = NULL, calc_norm_auprc = TRUE) {
  keep <- complete.cases(true, score)
  if (!is.null(pred)) keep <- keep & !is.na(pred)
  
  true <- as.integer(true[keep])
  score <- as.numeric(score[keep])
  
  if (is.null(pred)) {
    pred <- ifelse(score >= 0.5, 1L, 0L)
  } else {
    pred <- as.integer(pred[keep])
  }
  
  auc_value <- if (length(unique(true)) == 2) {
    as.numeric(
      pROC::auc(
        pROC::roc(
          response = true,
          predictor = score,
          levels = c(0, 1),
          direction = "<",
          quiet = TRUE
        )
      )
    )
  } else {
    NA_real_
  }
  
  cm <- caret::confusionMatrix(
    factor(pred, levels = c(0, 1)),
    factor(true, levels = c(0, 1))
  )
  
  acc_value <- unname(cm$overall["Accuracy"])
  f1_value <- weighted_f1(true, pred)
  ap_value <- average_precision_sklearn(true, score)
  
  baseline <- mean(true)
  norm_auprc <- if (
    calc_norm_auprc &&
    !is.na(ap_value) &&
    baseline < 1
  ) {
    (ap_value - baseline) / (1 - baseline)
  } else {
    NA_real_
  }
  
  tibble(
    AUROC = auc_value,
    ACC = acc_value,
    F1 = f1_value,
    AUPRC = ap_value,
    NormAUPRC = norm_auprc
  )
}

mean_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}


# ============================================================
# 1. scDEAL original benchmark results
# ============================================================

standardize_scdeal_original <- function(file, imbalance = FALSE) {
  d <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  
  dataset_col <- pick_col(d, c("Dataset", "dataset", "GEOID", "GEO_ID", "id"))
  drug_col <- pick_col(d, c("Drug", "drug"))
  round_col <- pick_col(d, c("round", "Round"))
  bottle_col <- pick_col(d, c("bottleneck", "Bottleneck"))
  sampling_col <- pick_col(d, c("sampling_method", "sampling", "Sampling"))
  auc_col <- pick_col(d, c("AUROC", "AUC"))
  acc_col <- pick_col(d, c("ACC", "Accuracy"))
  f1_col <- pick_col(d, c("F1", "f1_score", "F1_score"))
  norm_col <- pick_col(d, c("NormAUPRC", "normAUPRC", "Normalized_AUPRC"), required = FALSE)
  
  out <- tibble(
    dataset = as.character(d[[dataset_col]]),
    drug = as.character(d[[drug_col]]),
    round = norm_round(d[[round_col]]),
    bottleneck = as.integer(d[[bottle_col]]),
    sampling_method = norm_sampling(d[[sampling_col]]),
    AUROC = as.numeric(d[[auc_col]]),
    ACC = as.numeric(d[[acc_col]]),
    F1 = as.numeric(d[[f1_col]]),
    NormAUPRC = if (!is.null(norm_col)) as.numeric(d[[norm_col]]) else NA_real_
  )
  
  if (imbalance) {
    ratio_col <- pick_col(d, c("ratio", "Ratio"))
    out$ratio <- norm_ratio(d[[ratio_col]])
    out$scenario <- paste0("imbalance_", out$ratio)
  } else {
    out$ratio <- NA_character_
    out$scenario <- "balanced"
  }
  
  out
}

scdeal_bal_original <- standardize_scdeal_original(
  scdeal_bal_original_file,
  imbalance = FALSE
) %>%
  group_by(dataset, drug, round) %>%
  arrange(desc(AUROC), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

scdeal_imb_original <- standardize_scdeal_original(
  scdeal_imb_original_file,
  imbalance = TRUE
) %>%
  filter(ratio %in% c("1_3", "1_30")) %>%
  group_by(dataset, drug, round, ratio) %>%
  arrange(desc(NormAUPRC), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

scdeal_original_best <- bind_rows(
  scdeal_bal_original,
  scdeal_imb_original
) %>%
  mutate(
    Architecture = "scDEAL",
    Training_label = "scDEAL label"
  )


# ============================================================
# 2. scDEAL replaced-label predictions
# ============================================================

calc_scdeal_replaced_one <- function(dataset, drug, round, ratio,
                                     bottleneck, sampling_method,
                                     scenario) {
  
  pred_dir <- if (scenario == "balanced") {
    file.path(
      scdeal_bal_replaced_root,
      round,
      "prediction_results"
    )
  } else {
    file.path(
      scdeal_imb_replaced_root,
      ratio,
      round,
      "prediction_results"
    )
  }
  
  pred_file <- file.path(
    pred_dir,
    paste0(
      dataset, "_", drug,
      "_DAE_bottle", bottleneck,
      "_modnew_", sampling_method,
      "_prediction_results.csv"
    )
  )
  
  if (!file.exists(pred_file)) {
    warning("Missing scDEAL prediction: ", pred_file)
    return(
      tibble(
        dataset = dataset, drug = drug, round = round,
        ratio = ratio, scenario = scenario,
        bottleneck = bottleneck,
        sampling_method = sampling_method,
        AUROC = NA_real_, ACC = NA_real_,
        F1 = NA_real_, AUPRC = NA_real_,
        NormAUPRC = NA_real_,
        prediction_file = pred_file,
        file_exists = FALSE
      )
    )
  }
  
  d <- read.csv(
    pred_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  true <- true_to_num(d$condition)
  score <- as.numeric(d$sens_preds)
  pred <- true_to_num(d$sens_label)
  
  m <- calc_metrics(
    true = true,
    score = score,
    pred = pred,
    calc_norm_auprc = scenario != "balanced"
  )
  
  bind_cols(
    tibble(
      dataset = dataset,
      drug = drug,
      round = round,
      ratio = ratio,
      scenario = scenario,
      bottleneck = bottleneck,
      sampling_method = sampling_method
    ),
    m,
    tibble(
      prediction_file = pred_file,
      file_exists = TRUE
    )
  )
}

scdeal_replaced <- pmap_dfr(
  scdeal_original_best %>%
    select(
      dataset, drug, round, ratio,
      bottleneck, sampling_method, scenario
    ),
  calc_scdeal_replaced_one
) %>%
  mutate(
    Architecture = "scDEAL",
    Training_label = "SCAD label"
  )

scdeal_all <- bind_rows(
  scdeal_original_best %>%
    mutate(
      AUPRC = NA_real_,
      prediction_file = NA_character_,
      file_exists = TRUE
    ),
  scdeal_replaced
)

write.csv(
  scdeal_all,
  file.path(out_dir, "scDEAL_Figure5C_metrics.csv"),
  row.names = FALSE
)


# ============================================================
# 3. SCAD mapping
# ============================================================

scad_pairs <- read.csv(
  scad_pair_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(scad_pairs)[1:2] <- c("dataset", "drug")

scad_pairs <- scad_pairs %>%
  mutate(
    drug_scad = case_when(
      drug == "PLX-4720" ~ "PLX4720",
      drug == "Palbociclib" ~ "PD-0332991",
      TRUE ~ drug
    ),
    dataset_folder = paste0(dataset, "_", drug_scad)
  )


# ============================================================
# 4. SCAD prediction files
# ============================================================

scad_files <- list.files(
  scad_root,
  pattern = "_predicted_y\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

parse_scad_file <- function(f) {
  rel <- substring(f, nchar(scad_root) + 2)
  p <- strsplit(rel, "/", fixed = TRUE)[[1]]
  
  if (length(p) < 6) stop("Unexpected SCAD path: ", f)
  
  label_folder <- p[1]
  scenario <- p[2]
  round <- norm_round(p[3])
  dataset_folder <- p[4]
  feature <- p[5]
  param_id <- sub("_predicted_y\\.tsv$", "", basename(f))
  
  map_row <- scad_pairs %>%
    filter(dataset_folder == !!dataset_folder)
  
  if (nrow(map_row) != 1) {
    warning("Cannot map SCAD dataset folder: ", dataset_folder)
    dataset <- dataset_folder
    drug <- NA_character_
  } else {
    dataset <- map_row$dataset[1]
    drug <- map_row$drug[1]
  }
  
  d <- read.delim(
    f,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  d <- d %>%
    mutate(
      true_label = as.integer(true_label),
      prediction = as.numeric(prediction),
      pred_label = ifelse(prediction >= 0.5, 1L, 0L)
    )
  
  calc_norm <- scenario != "balanced"
  
  split_metrics <- d %>%
    group_by(split) %>%
    group_modify(~{
      calc_metrics(
        true = .x$true_label,
        score = .x$prediction,
        pred = .x$pred_label,
        calc_norm_auprc = calc_norm
      )
    }) %>%
    ungroup()
  
  summary_metrics <- split_metrics %>%
    summarise(
      AUROC = mean_na(AUROC),
      ACC = mean_na(ACC),
      F1 = mean_na(F1),
      AUPRC = mean_na(AUPRC),
      NormAUPRC = mean_na(NormAUPRC),
      n_splits = n()
    )
  
  ratio <- case_when(
    scenario == "imbalance_1_3" ~ "1_3",
    scenario == "imbalance_1_30" ~ "1_30",
    TRUE ~ NA_character_
  )
  
  bind_cols(
    tibble(
      Architecture = "SCAD",
      Training_label = ifelse(
        label_folder == "SCAD_label",
        "SCAD label",
        "scDEAL label"
      ),
      label_folder = label_folder,
      scenario = scenario,
      ratio = ratio,
      round = round,
      dataset = dataset,
      drug = drug,
      dataset_folder = dataset_folder,
      feature = feature,
      param_id = param_id,
      prediction_file = f
    ),
    summary_metrics
  )
}

scad_all_metrics <- map_dfr(
  scad_files,
  parse_scad_file
)

write.csv(
  scad_all_metrics,
  file.path(out_dir, "SCAD_all_configuration_metrics.csv"),
  row.names = FALSE
)


# ============================================================
# 5. Pair SCAD original/replaced labels using identical config
# ============================================================

scad_original <- scad_all_metrics %>%
  filter(Training_label == "SCAD label") %>%
  select(
    scenario, ratio, round, dataset, drug,
    dataset_folder, feature, param_id,
    AUROC, ACC, F1, AUPRC, NormAUPRC
  ) %>%
  rename_with(
    ~paste0(.x, "_SCADlabel"),
    c(AUROC, ACC, F1, AUPRC, NormAUPRC)
  )

scad_replaced <- scad_all_metrics %>%
  filter(Training_label == "scDEAL label") %>%
  select(
    scenario, ratio, round, dataset, drug,
    dataset_folder, feature, param_id,
    AUROC, ACC, F1, AUPRC, NormAUPRC
  ) %>%
  rename_with(
    ~paste0(.x, "_scDEALlabel"),
    c(AUROC, ACC, F1, AUPRC, NormAUPRC)
  )

scad_paired_all <- inner_join(
  scad_original,
  scad_replaced,
  by = c(
    "scenario", "ratio", "round",
    "dataset", "drug", "dataset_folder",
    "feature", "param_id"
  )
)

write.csv(
  scad_paired_all,
  file.path(out_dir, "SCAD_all_matched_configurations.csv"),
  row.names = FALSE
)


# ============================================================
# 6. One matched SCAD configuration per dataset
#
# balanced  : best SCAD-label AUROC
# imbalance : best SCAD-label NormAUPRC
# Then use the SAME configuration for replaced label.
# ============================================================

scad_best_matched <- scad_paired_all %>%
  mutate(
    selection_metric = ifelse(
      scenario == "balanced",
      AUROC_SCADlabel,
      NormAUPRC_SCADlabel
    )
  ) %>%
  group_by(
    scenario, ratio, round,
    dataset, drug
  ) %>%
  arrange(desc(selection_metric), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

write.csv(
  scad_best_matched,
  file.path(out_dir, "SCAD_best_matched_metrics.csv"),
  row.names = FALSE
)

scad_plot_data <- bind_rows(
  scad_best_matched %>%
    transmute(
      Architecture = "SCAD",
      Training_label = "SCAD label",
      scenario, ratio, round, dataset, drug,
      feature, param_id,
      AUROC = AUROC_SCADlabel,
      ACC = ACC_SCADlabel,
      F1 = F1_SCADlabel,
      NormAUPRC = NormAUPRC_SCADlabel
    ),
  
  scad_best_matched %>%
    transmute(
      Architecture = "SCAD",
      Training_label = "scDEAL label",
      scenario, ratio, round, dataset, drug,
      feature, param_id,
      AUROC = AUROC_scDEALlabel,
      ACC = ACC_scDEALlabel,
      F1 = F1_scDEALlabel,
      NormAUPRC = NormAUPRC_scDEALlabel
    )
)


# ============================================================
# 7. Prepare Figure 5C plotting data
# ============================================================

scdeal_plot_data <- scdeal_all %>%
  select(
    Architecture, Training_label,
    scenario, ratio, round, dataset, drug,
    AUROC, ACC, F1, NormAUPRC
  )

fig5c_data <- bind_rows(
  scdeal_plot_data,
  scad_plot_data %>%
    select(
      Architecture, Training_label,
      scenario, ratio, round, dataset, drug,
      AUROC, ACC, F1, NormAUPRC
    )
) %>%
  mutate(
    Scenario = case_when(
      scenario == "balanced" ~ "Balanced",
      ratio == "1_3" ~ "1:3",
      ratio == "1_30" ~ "1:30",
      TRUE ~ scenario
    ),
    Scenario = factor(
      Scenario,
      levels = c("Balanced", "1:3", "1:30")
    ),
    pair_id = paste(
      Architecture,
      dataset,
      drug,
      round,
      sep = "__"
    )
  ) %>%
  pivot_longer(
    cols = c(AUROC, ACC, F1, NormAUPRC),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  filter(
    !(Scenario == "Balanced" & Metric == "NormAUPRC")
  )

write.csv(
  fig5c_data,
  file.path(out_dir, "Figure5C_plotting_data_long.csv"),
  row.names = FALSE
)


# ============================================================
# 8. Visualization
# ============================================================
get_paired_pvalue <- function(data, architecture) {
  
  dat <- data %>%
    filter(
      Architecture == architecture,
      !is.na(Value)
    ) %>%
    select(
      Metric,
      Scenario,
      pair_id,
      Training_label,
      Value
    ) %>%
    distinct()
  
  wide <- dat %>%
    pivot_wider(
      names_from = Training_label,
      values_from = Value
    ) %>%
    filter(
      !is.na(`SCAD label`),
      !is.na(`scDEAL label`)
    )
  
  wide %>%
    group_by(Metric, Scenario) %>%
    summarise(
      p = wilcox.test(
        `SCAD label`,
        `scDEAL label`,
        paired = TRUE,
        exact = FALSE
      )$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p.signif = case_when(
        p <= 0.0001 ~ "∗∗∗∗",
        p <= 0.001  ~ "∗∗∗",
        p <= 0.01   ~ "∗∗",
        p <= 0.05   ~ "∗",
        TRUE        ~ "ns"
      )
    )
}

p_scdeal_stat <- get_paired_pvalue(
  fig5c_data,
  "scDEAL"
)

p_scad_stat <- get_paired_pvalue(
  fig5c_data,
  "SCAD"
)

p_scdeal <- fig5c_data %>%
  filter(Architecture == "scDEAL") %>%
  ggplot(
    aes(
      x = Training_label,
      y = Value,
      group = pair_id,
      fill = Training_label
    )
  ) +
  geom_line(alpha = 0.25) +
  geom_point(
    aes(shape = Training_label),
    size = 2,
    alpha = 0.75
  ) +
  geom_boxplot(
    aes(group = Training_label),
    width = 0.45,
    alpha = 0.15,
    outlier.shape = NA
  ) +
  
  # 固定显著性横线
  geom_segment(
    data = p_scdeal_stat,
    aes(
      x = 1,
      xend = 2,
      y = 1.02,
      yend = 1.02
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # 横线中间显示显著性
  geom_text(
    data = p_scdeal_stat,
    aes(
      x = 1.5,
      y = 1.02,
      label = p.signif,
      vjust = -0.05
    ),
    inherit.aes = FALSE,
    size = 6,
    vjust = -0.35
  ) +
  
  facet_grid(
    Metric ~ Scenario,
    scales = "free_y"
  ) + expand_limits(y = 1.2) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Performance"
  ) +
  theme(
    text = element_text(size = 18),
    legend.position = "none",
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

p_scad <- fig5c_data %>%
  filter(Architecture == "SCAD") %>%
  ggplot(
    aes(
      x = Training_label,
      y = Value,
      group = pair_id,
      fill = Training_label
    )
  ) +
  geom_line(alpha = 0.25) +
  geom_point(
    aes(shape = Training_label),
    size = 2,
    alpha = 0.75
  ) +
  geom_boxplot(
    aes(group = Training_label),
    width = 0.45,
    alpha = 0.15,
    outlier.shape = NA
  ) +
  
  geom_segment(
    data = p_scad_stat,
    aes(
      x = 1,
      xend = 2,
      y = 1.02,
      yend = 1.02
    ),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  geom_text(
    data = p_scad_stat,
    aes(
      x = 1.5,
      y = 1.02,
      label = p.signif,
      vjust = -0.05
    ),
    inherit.aes = FALSE,
    size = 6,
    vjust = -0.35
  ) +
  
  facet_grid(
    Metric ~ Scenario,
    scales = "free_y"
  ) + expand_limits(y = 1.2) +
  theme_classic(base_size = 12) +
  labs(
    x = NULL,
    y = "Performance"
  ) +
  theme(
    text = element_text(size = 18),
    legend.position = "none",
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

png(
  paste0(out_dir, "Figure5C_scDEAL.png"),
  units = "in",
  width = 8,
  height = 10,
  res = 300
)
print(p_scdeal)
dev.off()

png(
  paste0(out_dir, "Figure5C_SCAD.png"),
  units = "in",
  width = 8,
  height = 10,
  res = 300
)
print(p_scad)
dev.off()

# ============================================================
# Quick checks
# ============================================================

cat("\n===== scDEAL missing prediction files =====\n")
print(
  scdeal_replaced %>%
    filter(!file_exists) %>%
    select(dataset, drug, round, ratio, bottleneck, sampling_method)
)

cat("\n===== Number of SCAD prediction files =====\n")
print(table(scad_all_metrics$Training_label, scad_all_metrics$scenario))

cat("\n===== Figure5C final sample counts =====\n")
print(
  fig5c_data %>%
    distinct(
      Architecture, Training_label,
      Scenario, dataset, drug, round
    ) %>%
    count(Architecture, Training_label, Scenario)
)

cat("\nDone. Output directory:\n", out_dir, "\n")

# ============================================================
# Figure 5D: effect decomposition
# ============================================================

# 统一drug名称
norm_drug<-function(x){
  x<-toupper(x)
  x<-gsub("[^A-Z0-9]","",x)
  x<-case_when(
    x=="PD0332991" ~ "PALBOCICLIB",
    TRUE ~ x
  )
  x
}

# ============================================================
# 1. Prepare 2 × 2 factorial data
#
# A = SCAD architecture + SCAD label
# B = SCAD architecture + scDEAL label
# C = scDEAL architecture + SCAD label
# D = scDEAL architecture + scDEAL label
#
# 3 rounds first averaged within each dataset-drug
# ============================================================

effect_base<-fig5c_data %>%
  filter(!is.na(Value)) %>%
  mutate(drug_key=norm_drug(drug)) %>%
  group_by(dataset,drug_key,Scenario,Metric,Architecture,Training_label) %>%
  summarise(Value=mean(Value,na.rm=TRUE),.groups="drop") %>%
  mutate(config=case_when(
    Architecture=="SCAD" & Training_label=="SCAD label" ~ "A",
    Architecture=="SCAD" & Training_label=="scDEAL label" ~ "B",
    Architecture=="scDEAL" & Training_label=="SCAD label" ~ "C",
    Architecture=="scDEAL" & Training_label=="scDEAL label" ~ "D"
  )) %>%
  select(dataset,drug_key,Scenario,Metric,config,Value) %>%
  pivot_wider(names_from=config,values_from=Value) %>%
  drop_na(A,B,C,D)

# ============================================================
# 2. Calculate effects
#
# Label effect > 0:
#   scDEAL label better than SCAD label
#
# Architecture effect > 0:
#   scDEAL architecture better than SCAD architecture
#
# Interaction > 0:
#   advantage of scDEAL label is stronger under scDEAL architecture
# ============================================================

effect_data<-effect_base %>%
  mutate(
    `Label effect`=((B-A)+(D-C))/2,
    `Architecture effect`=((C-A)+(D-B))/2,
    `Label × architecture`=(D-C)-(B-A)
  ) %>%
  pivot_longer(
    c(`Label effect`,`Architecture effect`,`Label × architecture`),
    names_to="Effect",
    values_to="Effect_size"
  ) %>%
  mutate(
    Effect=factor(Effect,levels=c("Label × architecture","Architecture effect","Label effect")),
    Scenario=factor(Scenario,levels=c("Balanced","1:3","1:30")),
    Metric=factor(Metric,levels=c("ACC","AUROC","F1","NormAUPRC"))
  )

# Check sample size
print(effect_data %>% count(Metric,Scenario,Effect))

# ============================================================
# 3. Bootstrap 95% CI
# ============================================================

boot_ci<-function(x,B=5000){
  x<-x[is.finite(x)]
  n<-length(x)
  if(n<2) return(tibble(mean=mean(x),lower=NA_real_,upper=NA_real_,n=n))
  b<-replicate(B,mean(x[sample.int(n,n,replace=TRUE)]))
  tibble(
    mean=mean(x),
    lower=unname(quantile(b,0.025)),
    upper=unname(quantile(b,0.975)),
    n=n
  )
}

set.seed(42)

forest_data<-effect_data %>%
  group_by(Metric,Scenario,Effect) %>%
  group_modify(~boot_ci(.x$Effect_size)) %>%
  ungroup()

print(forest_data)

write.csv(effect_data,file.path(out_dir,"Figure5D_effect_data.csv"),row.names=FALSE)
write.csv(forest_data,file.path(out_dir,"Figure5D_forest_summary.csv"),row.names=FALSE)

# ============================================================
# 4. Combined forest plot
# ============================================================

library(patchwork)

make_forest_row<-function(metric_name,show_strip=TRUE){
  fd<-forest_data %>% filter(Metric==metric_name)
  xlim<-max(abs(c(fd$mean,fd$lower,fd$upper)),na.rm=TRUE)*1.15
  
  p<-ggplot(fd,aes(x=mean,y=Effect))+
    geom_vline(xintercept=0,linetype=2,linewidth=0.5,color="grey50")+
    geom_segment(aes(x=lower,xend=upper,y=Effect,yend=Effect),linewidth=0.8)+
    geom_point(size=3.5)+
    facet_grid(.~Scenario,drop=FALSE)+
    coord_cartesian(xlim=c(-xlim,xlim))+
    scale_x_continuous(breaks=pretty(c(-xlim,xlim),n=5))+
    theme_classic(base_size=14)+
    labs(x=NULL,y=NULL)+
    theme(
      text=element_text(size=17),
      strip.background=element_rect(fill="white",color="black"),
      strip.text=element_text(size=16),
      panel.spacing=unit(0.8,"lines"),
      axis.text.y=element_text(size=14),
      plot.margin=margin(2,5,2,5)
    )
  
  if(!show_strip){
    p<-p+theme(strip.text.x=element_blank(),strip.background.x=element_blank())
  }
  
  p
}

p_acc<-make_forest_row("ACC",TRUE)
p_auc<-make_forest_row("AUROC",FALSE)
p_f1<-make_forest_row("F1",FALSE)
p_norm<-make_forest_row("NormAUPRC",FALSE)

# 给每一行加Metric名称
p_acc<-p_acc+labs(tag="ACC")+theme(plot.tag=element_text(size=17,face="bold"))
p_auc<-p_auc+labs(tag="AUROC")+theme(plot.tag=element_text(size=17,face="bold"))
p_f1<-p_f1+labs(tag="F1")+theme(plot.tag=element_text(size=17,face="bold"))
p_norm<-p_norm+labs(tag="NormAUPRC", x="Effect size (Δ performance)")+theme(plot.tag=element_text(size=17,face="bold")) 

# 四行拼成一张图
p_forest_all<-(p_acc/p_auc/p_f1/p_norm)+
  plot_layout(heights=c(1,1,1,1))+
  plot_annotation(
    theme=theme(plot.margin=margin(5,5,5,5))
  )

# 最下面统一加x轴标题
p_forest_all<-p_forest_all & theme(axis.title.x=element_blank())

png(paste0(out_dir,"Figure5D_effect_forest_all.png"), units = 'in', res = 300, width = 10, height = 6)
print(p_forest_all)
dev.off()

# ============================================================
# Figure 5D supplementary: model-centered effects + interaction
# ============================================================

specific_effect_data<-effect_base %>%
  transmute(
    dataset,drug_key,Scenario,Metric,
    `scDEAL: Label effect`=D-C,
    `scDEAL: Architecture effect`=D-B,
    `SCAD: Label effect`=A-B,
    `SCAD: Architecture effect`=A-C,
    `Label × architecture`=(D-C)-(B-A)
  ) %>%
  pivot_longer(
    cols=c(`scDEAL: Label effect`,`scDEAL: Architecture effect`,
           `SCAD: Label effect`,`SCAD: Architecture effect`,
           `Label × architecture`),
    names_to="Effect",values_to="Effect_size"
  ) %>%
  mutate(
    Effect=factor(
      Effect,
      levels=c(
        "Label × architecture",
        "SCAD: Architecture effect",
        "SCAD: Label effect",
        "scDEAL: Architecture effect",
        "scDEAL: Label effect"
      )
    ),
    Scenario=factor(Scenario,levels=c("Balanced","1:3","1:30")),
    Metric=factor(Metric,levels=c("ACC","AUROC","F1","NormAUPRC"))
  )

set.seed(42)

specific_forest_data<-specific_effect_data %>%
  group_by(Metric,Scenario,Effect) %>%
  group_modify(~boot_ci(.x$Effect_size)) %>%
  ungroup()

write.csv(
  specific_effect_data,
  file.path(out_dir,"Figure5D_model_centered_effect_data.csv"),
  row.names=FALSE
)

write.csv(
  specific_forest_data,
  file.path(out_dir,"Figure5D_model_centered_forest_summary.csv"),
  row.names=FALSE
)


# ============================================================
# Plot function
# ============================================================

make_specific_forest_row<-function(metric_name,show_strip=TRUE){
  fd<-specific_forest_data %>% filter(Metric==metric_name)
  xlim<-max(abs(c(fd$mean,fd$lower,fd$upper)),na.rm=TRUE)*1.15
  
  p<-ggplot(fd,aes(x=mean,y=Effect))+
    geom_vline(xintercept=0,linetype=2,linewidth=0.5,color="grey50")+
    geom_segment(aes(x=lower,xend=upper,y=Effect,yend=Effect),linewidth=0.8)+
    geom_point(size=3.5)+
    facet_grid(.~Scenario,drop=FALSE)+
    coord_cartesian(xlim=c(-xlim,xlim))+
    scale_x_continuous(breaks=pretty(c(-xlim,xlim),n=5))+
    theme_classic(base_size=14)+
    labs(x=NULL,y=NULL)+
    theme(
      text=element_text(size=17),
      strip.background=element_rect(fill="white",color="black"),
      strip.text=element_text(size=16),
      panel.spacing=unit(0.8,"lines"),
      axis.text.y=element_text(size=13),
      plot.margin=margin(2,5,2,5)
    )
  
  if(!show_strip){
    p<-p+theme(
      strip.text.x=element_blank(),
      strip.background.x=element_blank()
    )
  }
  
  p
}


# ============================================================
# Four metrics
# ============================================================

p_specific_acc<-make_specific_forest_row("ACC",TRUE)+
  labs(tag="ACC")+
  theme(plot.tag=element_text(size=17,face="bold"), axis.text.y = element_text(size = 18))

p_specific_auc<-make_specific_forest_row("AUROC",FALSE)+
  labs(tag="AUROC")+
  theme(plot.tag=element_text(size=17,face="bold"), axis.text.y = element_text(size = 18))

p_specific_f1<-make_specific_forest_row("F1",FALSE)+
  labs(tag="F1")+
  theme(plot.tag=element_text(size=17,face="bold"), axis.text.y = element_text(size = 18))

p_specific_norm<-make_specific_forest_row("NormAUPRC",FALSE)+
  labs(tag="NormAUPRC",x="Effect size (Δ performance)")+
  theme(plot.tag=element_text(size=17,face="bold"), axis.text.y = element_text(size = 18))


# ============================================================
# Combine
# ============================================================

p_specific_all<-(p_specific_acc/p_specific_auc/p_specific_f1/p_specific_norm)+
  plot_layout(heights=c(1,1,1,1))+
  plot_annotation(
    theme=theme(plot.margin=margin(5,5,5,5))
  )

print(p_specific_all)

png(
  paste0(out_dir,"Figure5D_model_centered_specific_effects.png"),
  units="in",res=300,width=11,height=12
)
print(p_specific_all)
dev.off()

# ============================================================
# Heatmap
# ============================================================

heatmap_data<-fig5c_data %>%
  filter(!is.na(Value)) %>%
  group_by(dataset,drug,Scenario,Metric,Architecture,Training_label) %>%
  summarise(Value=mean(Value,na.rm=TRUE),.groups="drop") %>%
  pivot_wider(names_from=Training_label,values_from=Value) %>%
  drop_na(`SCAD label`,`scDEAL label`) %>%
  mutate(
    Delta=case_when(
      Architecture=="scDEAL" ~ `SCAD label`-`scDEAL label`,
      Architecture=="SCAD" ~ `scDEAL label`-`SCAD label`
    ),
    Dataset_drug=paste0(dataset," | ",drug)
  )

p_heatmap<-ggplot(heatmap_data,aes(x=Scenario,y=Dataset_drug,fill=Delta))+
  geom_tile(color="white",linewidth=0.3)+
  geom_text(aes(label=sprintf("%.2f",Delta)),size=3)+
  scale_fill_gradient2(low="#3B70B2",mid="white",high="#C84A42",midpoint=0)+
  facet_grid(Architecture~Metric)+
  theme_classic(base_size=13)+
  labs(x=NULL,y=NULL,fill="Δ performance\n(replaced − original)")+
  theme(
    axis.text.x=element_text(angle=45,hjust=1),
    panel.spacing=unit(0.8,"lines"),
    strip.background=element_blank()
  )

print(p_heatmap)

png(paste0(out_dir,"Figure5C_label_replacement_heatmap.png"), units = 'in', res = 300, width = 10, height = 6)
print(p_heatmap)
dev.off()
