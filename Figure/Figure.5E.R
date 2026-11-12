# ============================================================
# scDEAL Figure 5E ablation result summary
#
# Input:
#   Figure5E/scDEAL_ablation/
#
# Output:
#   1. scDEAL_ablation_performance.csv
#   2. scDEAL_ablation_counts.csv
#   3. scDEAL_ablation_missing_results.csv
#
# Main ablations:
#   no_cluster
#   no_dae
#   gdsc_only
#   ccle_only
#
# Additional:
#   no_transfer
#   -> extracted from no_cluster files using
#      sens_preds_pret / sens_label_pret
# ============================================================

library(tidyverse)
library(caret)
library(pROC)

# ============================================================
# 1. Paths
# ============================================================
root_dir <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_ablation"
mapping_csv <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"
param_table_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_ablation_best_params.tsv"
output_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_ablation_performance.csv"
count_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_ablation_counts.csv"
missing_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_ablation_missing_results.csv"

# ============================================================
# 2. Read dataset-drug mapping
# ============================================================
mapping_raw <- readr::read_csv(mapping_csv, show_col_types = FALSE)

# only use first two columns
mapping <- mapping_raw[, 1:2]
colnames(mapping) <- c("dataset", "drug")
mapping <- mapping %>%
  mutate(
    dataset = trimws(as.character(dataset)),
    drug = trimws(as.character(drug))
  ) %>%
  distinct()

print(mapping)

# ============================================================
# 3. Utility functions
# ============================================================
safe_name <- function(x) {
  x <- gsub("/", "_", x, fixed = TRUE)
  x <- gsub(" ", "_", x, fixed = TRUE)
  x
}

identify_dataset_drug <- function(file) {
  fname <- basename(file)
  prefixes <- paste0(safe_name(mapping$dataset), "_", safe_name(mapping$drug), "_")
  hit <- which(vapply(prefixes, function(x) startsWith(fname, x), logical(1)))
  
  if (length(hit) == 0) {
    stop(paste0("\nCannot identify dataset/drug from filename:\n", fname))
  }
  
  if (length(hit) > 1) {
    prefix_lengths <- nchar(prefixes[hit])
    hit <- hit[which.max(prefix_lengths)]
  }
  
  tibble(dataset = mapping$dataset[hit], drug = mapping$drug[hit])
}

parse_file_path <- function(file) {
  root_norm <- normalizePath(root_dir, winslash = "/")
  file_norm <- normalizePath(file, winslash = "/")
  relative_path <- substring(file_norm, nchar(root_norm) + 2)
  parts <- strsplit(relative_path, "/", fixed = TRUE)[[1]]
  
  ablation <- parts[1]
  scenario_type <- parts[2]
  
  if (scenario_type == "balanced") {
    setting <- "balanced"
    round_id <- as.integer(sub("round", "", parts[3]))
    ratio <- NA_character_
  } else if (scenario_type == "imbalanced") {
    setting <- parts[3]
    ratio <- parts[3]
    round_id <- as.integer(sub("round", "", parts[4]))
  } else {
    stop(paste("Unknown directory structure:", file))
  }
  
  tibble(
    source_ablation = ablation,
    scenario = scenario_type,
    setting = setting,
    ratio = ratio,
    round = round_id
  )
}

extract_bottleneck <- function(file) {
  fname <- basename(file)
  tmp <- stringr::str_match(fname, "_bottle([0-9]+)_")
  as.integer(tmp[, 2])
}

extract_sampling <- function(file) {
  fname <- basename(file)
  tmp <- stringr::str_match(fname, "_(upsampling|downsampling|SMOTE)_prediction_results\\.csv$")
  tmp[, 2]
}

# ============================================================
# 4. sklearn-like Average Precision
# ============================================================
average_precision_binary <- function(y_true, score) {
  y_true <- as.integer(y_true)
  score <- as.numeric(score)
  n_positive <- sum(y_true == 1)
  
  if (n_positive == 0) return(NA_real_)
  
  ord <- order(score, decreasing = TRUE)
  y <- y_true[ord]
  s <- score[ord]
  
  tp_all <- cumsum(y == 1)
  fp_all <- cumsum(y == 0)
  
  end_idx <- c(which(diff(s) != 0), length(s))
  tp <- tp_all[end_idx]
  fp <- fp_all[end_idx]
  
  precision <- tp / (tp + fp)
  recall <- tp / n_positive
  previous_recall <- c(0, head(recall, -1))
  ap <- sum((recall - previous_recall) * precision)
  
  as.numeric(ap)
}

# ============================================================
# 5. Weighted metrics from confusion matrix
# ============================================================
calculate_weighted_metrics <- function(cm_table) {
  tab <- as.matrix(cm_table)
  classes <- intersect(rownames(tab), colnames(tab))
  
  class_metrics <- lapply(classes, function(cl) {
    TP <- tab[cl, cl]
    FP <- sum(tab[cl, ]) - TP
    FN <- sum(tab[, cl]) - TP
    support <- sum(tab[, cl])
    
    precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
    recall <- ifelse(TP + FN == 0, 0, TP / (TP + FN))
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
    
    tibble(class = cl, support = support, precision = precision, recall = recall, f1 = f1)
  }) %>% bind_rows()
  
  total_support <- sum(class_metrics$support)
  weights <- class_metrics$support / total_support
  
  tibble(
    Precision = sum(class_metrics$precision * weights),
    Recall = sum(class_metrics$recall * weights),
    F1 = sum(class_metrics$f1 * weights)
  )
}

# ============================================================
# 6. Calculate metrics for one prediction file
# ============================================================
calculate_metrics <- function(df, score_col, label_col) {
  required_cols <- c("condition", score_col, label_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
  }
  
  condition <- as.character(df$condition)
  condition_levels <- sort(unique(condition[!is.na(condition)]))
  
  if (length(condition_levels) != 2) {
    stop(paste("condition must contain exactly 2 classes:", paste(condition_levels, collapse = ", ")))
  }
  
  negative_class <- condition_levels[1]
  positive_class <- condition_levels[2]
  y_true <- ifelse(condition == positive_class, 1L, 0L)
  score <- suppressWarnings(as.numeric(df[[score_col]]))
  pred_raw <- as.character(df[[label_col]])
  y_pred <- suppressWarnings(as.integer(as.numeric(pred_raw)))
  
  if (any(is.na(y_pred))) {
    y_pred <- ifelse(pred_raw == positive_class, 1L, ifelse(pred_raw == negative_class, 0L, NA_integer_))
  }
  
  keep <- complete.cases(y_true, score, y_pred)
  if (sum(!keep) > 0) {
    warning(paste(sum(!keep), "cells with missing values removed."))
  }
  
  y_true <- y_true[keep]
  score <- score[keep]
  y_pred <- y_pred[keep]
  
  if (length(unique(y_true)) != 2) {
    stop("Only one true class is present.")
  }
  
  roc_obj <- pROC::roc(response = y_true, predictor = score, levels = c(0, 1), direction = "<", quiet = TRUE)
  AUROC <- as.numeric(pROC::auc(roc_obj))
  AUPRC <- average_precision_binary(y_true, score)
  baseline <- mean(y_true)
  NormAUPRC <- (AUPRC - baseline) / (1 - baseline)
  
  truth_factor <- factor(y_true, levels = c(0, 1))
  prediction_factor <- factor(y_pred, levels = c(0, 1))
  cm <- caret::confusionMatrix(data = prediction_factor, reference = truth_factor, positive = "1", mode = "everything")
  
  ACC <- unname(cm$overall["Accuracy"])
  weighted <- calculate_weighted_metrics(cm$table)
  
  get_cm_value <- function(x, possible_names) {
    for (nm in possible_names) {
      if (nm %in% names(x)) return(as.numeric(x[nm]))
    }
    NA_real_
  }
  
  Precision_pos <- get_cm_value(cm$byClass, c("Precision", "Pos Pred Value"))
  Recall_pos <- get_cm_value(cm$byClass, c("Recall", "Sensitivity"))
  F1_pos <- get_cm_value(cm$byClass, "F1")
  Specificity <- get_cm_value(cm$byClass, "Specificity")
  
  tibble(
    N = length(y_true),
    Negative_N = sum(y_true == 0),
    Positive_N = sum(y_true == 1),
    Negative_class = negative_class,
    Positive_class = positive_class,
    Baseline = baseline,
    AUROC = AUROC,
    AUPRC = AUPRC,
    NormAUPRC = NormAUPRC,
    ACC = ACC,
    Precision = weighted$Precision,
    Recall = weighted$Recall,
    F1 = weighted$F1,
    Precision_pos = Precision_pos,
    Recall_pos = Recall_pos,
    F1_pos = F1_pos,
    Specificity = Specificity
  )
}

# ============================================================
# 7. Model configuration table
# ============================================================
model_config <- tribble(
  ~ablation,     ~bulk_source, ~encoder, ~mod,
  "no_cluster",  "integrate",  "DAE",    "ori",
  "no_dae",      "integrate",  "AE",     "new",
  "gdsc_only",   "GDSC",       "DAE",    "new",
  "ccle_only",   "CCLE",       "DAE",    "new",
  "no_transfer", "integrate",  "DAE",    NA_character_
)

# ============================================================
# 8. Find all prediction files
# ============================================================
prediction_files <- list.files(root_dir, pattern = "_prediction_results\\.csv$", recursive = TRUE, full.names = TRUE)
prediction_files <- prediction_files[grepl("/prediction_results/", gsub("\\\\", "/", prediction_files))]

cat("\nPrediction files found:", length(prediction_files), "\n\n")
if (length(prediction_files) == 0) stop("No prediction files were found.")

# ============================================================
# 9. Process one prediction CSV
# ============================================================
process_prediction_file <- function(file) {
  path_info <- parse_file_path(file)
  
  if (!path_info$source_ablation %in% c("no_cluster", "no_dae", "gdsc_only", "ccle_only")) return(NULL)
  
  pair_info <- identify_dataset_drug(file)
  bottleneck <- extract_bottleneck(file)
  sampling_method <- extract_sampling(file)
  df <- readr::read_csv(file, show_col_types = FALSE)
  
  metrics_main <- calculate_metrics(df = df, score_col = "sens_preds", label_col = "sens_label")
  result_main <- bind_cols(
    tibble(
      ablation = path_info$source_ablation,
      source_ablation = path_info$source_ablation,
      scenario = path_info$scenario,
      setting = path_info$setting,
      ratio = path_info$ratio,
      round = path_info$round,
      dataset = pair_info$dataset,
      drug = pair_info$drug,
      bottleneck = bottleneck,
      sampling_method = sampling_method,
      score_column = "sens_preds",
      label_column = "sens_label",
      source_file = file
    ),
    metrics_main
  )
  
  result_list <- list(result_main)
  
  if (path_info$source_ablation == "no_cluster") {
    if (all(c("sens_preds_pret", "sens_label_pret") %in% colnames(df))) {
      metrics_pret <- calculate_metrics(df = df, score_col = "sens_preds_pret", label_col = "sens_label_pret")
      result_pret <- bind_cols(
        tibble(
          ablation = "no_transfer",
          source_ablation = "no_cluster",
          scenario = path_info$scenario,
          setting = path_info$setting,
          ratio = path_info$ratio,
          round = path_info$round,
          dataset = pair_info$dataset,
          drug = pair_info$drug,
          bottleneck = bottleneck,
          sampling_method = sampling_method,
          score_column = "sens_preds_pret",
          label_column = "sens_label_pret",
          source_file = file
        ),
        metrics_pret
      )
      result_list[[2]] <- result_pret
    } else {
      warning(paste("No pret prediction columns found:", file))
    }
  }
  
  bind_rows(result_list)
}

# ============================================================
# 10. Process all files
# ============================================================
cat("Calculating metrics...\n")
performance <- purrr::map_dfr(prediction_files, process_prediction_file)
performance <- performance %>% left_join(model_config, by = "ablation")

# ============================================================
# 11. Check duplicates
# ============================================================
duplicate_check <- performance %>%
  count(ablation, setting, round, dataset, drug, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_check) > 0) {
  print(duplicate_check)
  stop(paste0(
    "\nDuplicate prediction results found.\n",
    "This usually means old/wrong prediction files still exist in the folders.\n",
    "Please clean them before analysis."
  ))
}

# ============================================================
# 12. Validate parameters against best-parameter table
# ============================================================
if (file.exists(param_table_file)) {
  expected_params <- readr::read_tsv(param_table_file, show_col_types = FALSE) %>%
    mutate(
      round = as.integer(round),
      bottleneck = as.integer(bottleneck),
      sampling_method = as.character(sampling_method)
    ) %>%
    rename(expected_bottleneck = bottleneck, expected_sampling = sampling_method)
  
  performance <- performance %>%
    left_join(expected_params, by = c("setting", "round", "dataset", "drug")) %>%
    mutate(
      parameter_match = case_when(
        is.na(expected_bottleneck) ~ NA,
        bottleneck == expected_bottleneck & sampling_method == expected_sampling ~ TRUE,
        TRUE ~ FALSE
      )
    )
  
  wrong_parameter <- performance %>% filter(parameter_match == FALSE)
  
  if (nrow(wrong_parameter) > 0) {
    print(wrong_parameter %>% select(ablation, setting, round, dataset, drug, bottleneck, sampling_method, expected_bottleneck, expected_sampling, source_file))
    stop("\nERROR: Some prediction files were generated using the wrong bottleneck/sampling parameters.")
  }
  cat("\nParameter validation passed.\n")
} else {
  warning("Best-parameter table was not found; parameter validation skipped.")
}

# ============================================================
# 13. Arrange result table
# ============================================================
performance <- performance %>%
  mutate(
    ablation = factor(ablation, levels = c("no_cluster", "no_transfer", "no_dae", "gdsc_only", "ccle_only")),
    setting = factor(setting, levels = c("balanced", "1_3", "1_30"))
  ) %>%
  arrange(ablation, setting, round, dataset, drug) %>%
  mutate(ablation = as.character(ablation), setting = as.character(setting))

# ============================================================
# 14. Count available results
# ============================================================
result_counts <- performance %>%
  count(ablation, setting, round, name = "N_results") %>%
  arrange(ablation, setting, round)

cat("\n========================================\nAvailable result counts\n========================================\n")
print(result_counts, n = Inf)

# ============================================================
# 15. Identify missing results
# ============================================================
missing_results <- tibble()

if (file.exists(param_table_file)) {
  expected_base <- readr::read_tsv(param_table_file, show_col_types = FALSE) %>%
    select(setting, round, dataset, drug)
  
  expected_all <- tidyr::crossing(
    ablation = c("no_cluster", "no_transfer", "no_dae", "gdsc_only", "ccle_only"),
    expected_base
  )
  
  observed <- performance %>% select(ablation, setting, round, dataset, drug) %>% distinct()
  missing_results <- expected_all %>% anti_join(observed, by = c("ablation", "setting", "round", "dataset", "drug"))
  
  cat("\n========================================\nMissing results\n========================================\n")
  if (nrow(missing_results) == 0) {
    cat("No missing results.\n")
  } else {
    print(missing_results, n = Inf)
  }
}

# ============================================================
# 16. Save results
# ============================================================
readr::write_csv(performance, output_file)
readr::write_csv(result_counts, count_file)
readr::write_csv(missing_results, missing_file)

cat("\n========================================\nFiles saved\n========================================\n")
cat("\nPerformance:\n", output_file, "\n")
cat("\nCounts:\n", count_file, "\n")
cat("\nMissing results:\n", missing_file, "\n")

# ============================================================
# 17. Sanity check
# ============================================================
sanity_check <- performance %>%
  filter(ablation == "no_cluster", setting == "1_3", round == 1, dataset == "GSE108397") %>%
  select(ablation, setting, round, dataset, drug, bottleneck, sampling_method, AUROC, AUPRC, NormAUPRC, ACC, Precision, Recall, F1)

cat("\n========================================\nSanity check: GSE108397 / 1_3 / round1\n========================================\n")
print(sanity_check)


# ============================================================
# 18. Add Full scDEAL performance from original benchmark
#
# Balanced:
# /media/shenli/Expansion/DRMREF/bal_final_merged_data.csv
#
# Imbalanced:
# /media/shenli/Expansion/DRMREF/unb_scDEAL_best_performance.csv
#
# Final output:
# scDEAL_Figure5E_full_and_ablation_performance.csv
# ============================================================

balanced_full_file <- "/media/shenli/Expansion/DRMREF/bal_final_merged_data.csv"
imbalanced_full_file <- "/media/shenli/Expansion/DRMREF/unb_scDEAL_best_performance.csv"
combined_output_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_Figure5E_full_and_ablation_performance.csv"

# ============================================================
# 19. Read the best-parameter table
# ============================================================
best_parameter_table <- readr::read_tsv(param_table_file, show_col_types = FALSE) %>%
  mutate(
    dataset = trimws(as.character(dataset)),
    drug = trimws(as.character(drug)),
    setting = trimws(as.character(setting)),
    round = as.integer(round),
    bottleneck = as.integer(bottleneck),
    sampling_method = trimws(as.character(sampling_method))
  )

# ============================================================
# 20. Full model: Balanced 1000 cells
# ============================================================
balanced_full_raw <- readr::read_csv(balanced_full_file, show_col_types = FALSE)

cat("\n========================================\nBalanced benchmark columns\n========================================\n")
print(colnames(balanced_full_raw))

# ------------------------------------------------------------
# Standardize important columns
# ------------------------------------------------------------
balanced_full <- balanced_full_raw %>%
  mutate(
    Dataset = trimws(as.character(Dataset)),
    Tool = trimws(as.character(Tool)),
    round = as.integer(round),
    cells = suppressWarnings(as.numeric(as.character(cells)))
  ) %>%
  # Figure 5E only: scDEAL, 1000 cells, round 1/2/3
  filter(tolower(Tool) == "scdeal", cells == 1000, round %in% c(1, 2, 3)) %>%
  # Keep only the 9 dataset-drug pairs used in Figure 5
  inner_join(mapping %>% select(dataset, drug), by = c("Dataset" = "dataset")) %>%
  transmute(
    ablation = "full",
    source_ablation = "full",
    scenario = "balanced",
    setting = "balanced",
    ratio = NA_character_,
    round = round,
    dataset = Dataset,
    drug = drug,
    bottleneck = NA_integer_,
    sampling_method = NA_character_,
    score_column = NA_character_,
    label_column = NA_character_,
    source_file = balanced_full_file,
    AUROC = as.numeric(AUC),
    AUPRC = NA_real_,
    NormAUPRC = NA_real_,
    ACC = as.numeric(ACC),
    Precision = as.numeric(Precision),
    Recall = as.numeric(Recall),
    F1 = as.numeric(F1)
  )

# ============================================================
# 21. Attach Full-model best parameters to balanced results
# ============================================================
balanced_parameters <- best_parameter_table %>%
  filter(setting == "balanced") %>%
  select(setting, round, dataset, drug, expected_bottleneck = bottleneck, expected_sampling_method = sampling_method)

balanced_full <- balanced_full %>%
  left_join(balanced_parameters, by = c("setting", "round", "dataset", "drug")) %>%
  mutate(bottleneck = expected_bottleneck, sampling_method = expected_sampling_method) %>%
  select(-expected_bottleneck, -expected_sampling_method)

# ============================================================
# 22. Check balanced Full results
# ============================================================
balanced_full_check <- balanced_full %>%
  count(setting, round, dataset, drug, name = "n") %>%
  filter(n != 1)

if (nrow(balanced_full_check) > 0) {
  print(balanced_full_check)
  stop("Duplicated or invalid balanced Full-model records found.")
}

cat("\n========================================\nFull model: balanced 1000-cell results\n========================================\n")
print(balanced_full %>% count(setting, round, name = "N_results"), n = Inf)

# ============================================================
# 23. Normalize imbalance ratio names
# ============================================================
normalize_ratio_R <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub(":", "_", x, fixed = TRUE)
  x <- gsub("/", "_", x, fixed = TRUE)
  x <- gsub(" ", "", x, fixed = TRUE)
  x <- sub("^1_0+3$", "1_3", x)
  x <- sub("^1_0+30$", "1_30", x)
  x
}

# ============================================================
# 24. Full model: Imbalanced
# ============================================================

imbalanced_full_raw <- readr::read_csv(imbalanced_full_file, show_col_types = FALSE)

# 先检查原始数据数量
cat("\nRaw imbalanced counts:\n")
print(
  imbalanced_full_raw %>%
    mutate(setting = normalize_ratio_R(ratio), round = as.integer(round)) %>%
    filter(setting %in% c("1_3", "1_30"), round %in% c(1, 2, 3)) %>%
    count(setting, round),
  n = Inf
)

# ------------------------------------------------------------
# 用 dataset 进行匹配
# drug 使用 mapping 中的标准名字
# ------------------------------------------------------------
mapping_for_join <- mapping %>%
  mutate(dataset_key = tolower(trimws(dataset))) %>%
  select(dataset_key, dataset_mapping = dataset, drug_mapping = drug)

imbalanced_full <- imbalanced_full_raw %>%
  mutate(
    dataset = trimws(as.character(dataset)),
    drug_original = trimws(as.character(drug)),
    dataset_key = tolower(dataset),
    round = as.integer(round),
    setting = normalize_ratio_R(ratio),
    bottleneck = as.integer(bottleneck),
    sampling_method = trimws(as.character(sampling_method))
  ) %>%
  filter(setting %in% c("1_3", "1_30"), round %in% c(1, 2, 3)) %>%
  # 只保留Figure 5的9个dataset
  inner_join(mapping_for_join, by = "dataset_key") %>%
  transmute(
    ablation = "full",
    source_ablation = "full",
    scenario = "imbalanced",
    setting = setting,
    ratio = setting,
    round = round,
    # 使用mapping中的标准名称
    dataset = dataset_mapping,
    drug = drug_mapping,
    bottleneck = bottleneck,
    sampling_method = sampling_method,
    score_column = NA_character_,
    label_column = NA_character_,
    source_file = imbalanced_full_file,
    AUROC = as.numeric(AUROC),
    AUPRC = as.numeric(AUPRC),
    NormAUPRC = as.numeric(NormAUPRC),
    ACC = as.numeric(ACC),
    Precision = as.numeric(Precision),
    Recall = as.numeric(Recall),
    F1 = as.numeric(F1)
  )

# ============================================================
# 25. Check imbalanced Full parameters against the parameter table
# ============================================================
imbalanced_parameter_check <- imbalanced_full %>%
  left_join(
    best_parameter_table %>%
      filter(setting %in% c("1_3", "1_30")) %>%
      select(setting, round, dataset, drug, expected_bottleneck = bottleneck, expected_sampling_method = sampling_method),
    by = c("setting", "round", "dataset", "drug")
  ) %>%
  mutate(parameter_match = bottleneck == expected_bottleneck & sampling_method == expected_sampling_method)

wrong_full_parameters <- imbalanced_parameter_check %>% filter(parameter_match == FALSE)

if (nrow(wrong_full_parameters) > 0) {
  cat("\nWARNING: Full-model imbalanced parameters do not match ablation parameter table:\n")
  print(wrong_full_parameters %>% select(setting, round, dataset, drug, bottleneck, sampling_method, expected_bottleneck, expected_sampling_method), n = Inf)
  stop("Full-model and ablation parameters are not matched.")
}

cat("\nImbalanced Full-model parameter validation passed.\n")

# ============================================================
# 26. Combine balanced + imbalanced Full model
# ============================================================
full_performance <- bind_rows(balanced_full, imbalanced_full)

# ============================================================
# 27. Check Full-model result counts
# ============================================================
full_counts <- full_performance %>%
  count(setting, round, name = "N_results") %>%
  arrange(factor(setting, levels = c("balanced", "1_3", "1_30")), round)

cat("\n========================================\nFull scDEAL result counts\n========================================\n")
print(full_counts, n = Inf)

# ============================================================
# 28. Combine Full + ablations
# ============================================================
full_performance <- full_performance %>%
  mutate(
    N = NA_integer_,
    Negative_N = NA_integer_,
    Positive_N = NA_integer_,
    Negative_class = NA_character_,
    Positive_class = NA_character_,
    Baseline = NA_real_,
    Precision_pos = NA_real_,
    Recall_pos = NA_real_,
    F1_pos = NA_real_,
    Specificity = NA_real_,
    bulk_source = "integrate",
    encoder = "DAE",
    mod = "new",
    parameter_match = TRUE
  )

# ============================================================
# 29. Combine
# ============================================================
performance_all <- bind_rows(full_performance, performance)

# ============================================================
# 30. Arrange model order
# ============================================================
performance_all <- performance_all %>%
  mutate(
    ablation = factor(ablation, levels = c("full", "no_transfer", "no_cluster", "no_dae", "gdsc_only", "ccle_only")),
    setting = factor(setting, levels = c("balanced", "1_3", "1_30"))
  ) %>%
  arrange(setting, round, dataset, drug, ablation) %>%
  mutate(ablation = as.character(ablation), setting = as.character(setting))

# ============================================================
# 31. Check duplicated Full/ablation records
# ============================================================
combined_duplicates <- performance_all %>%
  count(ablation, setting, round, dataset, drug, name = "n") %>%
  filter(n > 1)

if (nrow(combined_duplicates) > 0) {
  print(combined_duplicates, n = Inf)
  stop("Duplicate records found after combining Full and ablation results.")
}

# ============================================================
# 32. Count results after integration
# ============================================================
combined_counts <- performance_all %>%
  count(ablation, setting, round, name = "N_results") %>%
  arrange(
    factor(ablation, levels = c("full", "no_transfer", "no_cluster", "no_dae", "gdsc_only", "ccle_only")),
    factor(setting, levels = c("balanced", "1_3", "1_30")),
    round
  )

cat("\n========================================\nFull + ablation result counts\n========================================\n")
print(combined_counts, n = Inf)

# ============================================================
# 33. Save combined results
# ============================================================
readr::write_csv(performance_all, combined_output_file)

cat("\n========================================\nFull + ablation performance saved:\n", combined_output_file, "\n========================================\n")

# ============================================================
# 34. Quick preview
# ============================================================
performance_all %>%
  select(ablation, setting, round, dataset, drug, bottleneck, sampling_method, AUROC, AUPRC, NormAUPRC, ACC, Precision, Recall, F1) %>%
  print(n = 30)

# ============================================================
# 35. Calculate paired delta performance
# Delta = Full - Ablated
# ============================================================
fig5e_dir <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/plots"
dir.create(fig5e_dir, recursive = TRUE, showWarnings = FALSE)

key_cols <- c("setting","round","dataset","drug")
metrics_use <- c("AUROC","NormAUPRC","F1","ACC","Precision","Recall")

full_long <- performance_all %>%
  filter(ablation == "full") %>%
  select(all_of(key_cols), all_of(metrics_use)) %>%
  pivot_longer(cols = all_of(metrics_use), names_to = "metric", values_to = "full_value")

# ============================================================
# 36. Transfer / Cluster / DAE
# ============================================================
direct_ablation <- performance_all %>%
  filter(ablation %in% c("no_transfer","no_cluster","no_dae")) %>%
  select(ablation, all_of(key_cols), all_of(metrics_use)) %>%
  pivot_longer(cols = all_of(metrics_use), names_to = "metric", values_to = "ablated_value") %>%
  left_join(full_long, by = c(key_cols,"metric")) %>%
  mutate(delta = full_value - ablated_value,
         element = recode(ablation,
                          no_transfer = "Transfer learning",
                          no_cluster = "Cluster regularization",
                          no_dae = "DAE pretraining"))

# ============================================================
# 37. Source-data integration
# Full - best(GDSC only, CCLE only)
# ============================================================
single_source_long <- performance_all %>%
  filter(ablation %in% c("gdsc_only","ccle_only")) %>%
  select(ablation, all_of(key_cols), all_of(metrics_use)) %>%
  pivot_longer(cols = all_of(metrics_use), names_to = "metric", values_to = "value")

best_single_source <- single_source_long %>%
  group_by(across(all_of(key_cols)),metric) %>%
  summarise(
    ablated_value = if(all(is.na(value))) NA_real_ else max(value,na.rm=TRUE),
    best_single_source = if(all(is.na(value))) NA_character_ else ablation[which.max(replace(value,is.na(value),-Inf))][1],
    .groups="drop"
  )

integration_delta <- best_single_source %>%
  left_join(full_long,by=c(key_cols,"metric")) %>%
  mutate(delta=full_value-ablated_value,
         ablation="source_integration",
         element="Source-data integration")

# ============================================================
# 38. Final delta table
# ============================================================
delta_performance <- bind_rows(
  direct_ablation %>% select(ablation,element,all_of(key_cols),metric,full_value,ablated_value,delta),
  integration_delta %>% select(ablation,element,all_of(key_cols),metric,full_value,ablated_value,delta)
) %>%
  mutate(
    setting=factor(setting,levels=c("balanced","1_3","1_30")),
    element=factor(element,levels=c("Transfer learning","Cluster regularization","DAE pretraining","Source-data integration"))
  ) %>%
  arrange(setting,round,dataset,element,metric)

write_csv(delta_performance,
          "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_delta_performance.csv")

# ============================================================
# 39. Quick numerical summary
# ============================================================
delta_summary <- delta_performance %>%
  filter(!is.na(delta)) %>%
  group_by(element,setting,metric) %>%
  summarise(
    n=n(),
    mean_delta=mean(delta),
    median_delta=median(delta),
    sd_delta=sd(delta),
    se=sd_delta/sqrt(n),
    ci_low=mean_delta-qt(0.975,df=n-1)*se,
    ci_high=mean_delta+qt(0.975,df=n-1)*se,
    positive_fraction=mean(delta>0),
    .groups="drop"
  )

write_csv(delta_summary,
          "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_delta_summary.csv")
print(delta_summary,n=Inf)

# ============================================================
# 40. Balanced Figure 5E data
# Metrics:
# AUROC / F1 / ACC / Precision / Recall
# ============================================================
balanced_plot_data <- delta_performance %>%
  filter(setting=="balanced",
         metric %in% c("AUROC","F1","ACC","Precision","Recall"),
         !is.na(delta)) %>%
  mutate(metric=factor(metric,levels=c("AUROC","F1","ACC","Precision","Recall")))

balanced_plot_summary <- balanced_plot_data %>%
  group_by(metric,element) %>%
  summarise(
    n=n(),
    mean_delta=mean(delta),
    se=sd(delta)/sqrt(n),
    ci_low=mean_delta-qt(0.975,df=n-1)*se,
    ci_high=mean_delta+qt(0.975,df=n-1)*se,
    .groups="drop"
  )

# ============================================================
# 41. Figure 5E Balanced: ablation beeswarm
# Each point = one dataset × round
# Diamond = mean
# Horizontal line = 95% CI
# ============================================================
p_balanced <- ggplot(balanced_plot_data,aes(x=delta,y=element)) +
  geom_vline(xintercept=0,linetype=2,linewidth=0.5) +
  geom_jitter(aes(fill=element),shape=21,size=2.0,alpha=0.65,height=0.16,width=0,stroke=0.2) +
  geom_segment(data=balanced_plot_summary,
               aes(x=ci_low,xend=ci_high,y=element,yend=element),
               inherit.aes=FALSE,linewidth=0.8) +
  geom_point(data=balanced_plot_summary,
             aes(x=mean_delta,y=element),
             inherit.aes=FALSE,shape=23,size=3.3,fill="white",stroke=0.9) +
  facet_wrap(~metric,ncol=5,scales="free_x") +
  scale_x_continuous(
    breaks=scales::breaks_pretty(n=3),
    labels=scales::label_number(accuracy=0.01),
    expand=expansion(mult=c(0.05,0.08))
  ) +
  labs(x=expression(Delta*" performance (Full - Ablated)"),y=NULL) +
  theme_classic(base_size=11) +
  theme(
    text = element_text(size = 18),
    legend.position="none",
    strip.background=element_rect(fill="grey95",color=NA),
    strip.text=element_text(face="bold"),
    axis.text.x=element_text(size=18),
    axis.text.y=element_text(size=18),
    panel.spacing=grid::unit(0.8,"lines")
  )

png(file.path(fig5e_dir,"Figure5E_balanced_ablation_beeswarm.png"),
    width=13,height=4.8,units="in",res=300)
p_balanced
dev.off()

# ============================================================
# 42. Imbalanced Figure 5E data
# Metrics:
# NormAUPRC / AUROC / F1 / ACC / Precision / Recall
# Rows:
# 1:3 / 1:30
# ============================================================
imbalanced_plot_data<-delta_performance %>%
  filter(setting%in%c("1_3","1_30"),
         metric%in%c("NormAUPRC","AUROC","F1","ACC","Precision","Recall"),
         !is.na(delta)) %>%
  group_by(round,dataset,drug,element,metric) %>%
  summarise(delta=mean(delta),.groups="drop") %>%
  mutate(metric=factor(metric,levels=c("NormAUPRC","AUROC","F1","ACC","Precision","Recall")))

imbalanced_plot_summary<-imbalanced_plot_data %>%
  group_by(metric,element) %>%
  summarise(n=n(),mean_delta=mean(delta),se=sd(delta)/sqrt(n),
            ci_low=mean_delta-qt(0.975,df=n-1)*se,
            ci_high=mean_delta+qt(0.975,df=n-1)*se,.groups="drop")

# ============================================================
# 43. Figure 5E Imbalanced: ablation beeswarm
# Rows = imbalance ratios
# Columns = metrics
# ============================================================
p_imbalanced<-ggplot(imbalanced_plot_data,aes(x=delta,y=element))+
  geom_vline(xintercept=0,linetype=2,linewidth=0.5)+
  geom_jitter(aes(fill=element),shape=21,size=2,alpha=0.65,height=0.16,width=0,stroke=0.2)+
  geom_segment(data=imbalanced_plot_summary,
               aes(x=ci_low,xend=ci_high,y=element,yend=element),
               inherit.aes=FALSE,linewidth=0.8)+
  geom_point(data=imbalanced_plot_summary,
             aes(x=mean_delta,y=element),
             inherit.aes=FALSE,shape=23,size=3.3,fill="white",stroke=0.9)+
  facet_wrap(~metric,ncol=6,scales="free_x")+
  scale_x_continuous(breaks=scales::breaks_pretty(n=3),
                     labels=scales::label_number(accuracy=0.1),
                     expand=expansion(mult=c(0.05,0.08)))+
  labs(x=expression(Delta*" performance (Full - Ablated)"),y=NULL)+
  theme_classic(base_size=11)+
  theme(text=element_text(size=18),
        legend.position="none",
        strip.background=element_rect(fill="grey95",color=NA),
        strip.text=element_text(face="bold",size=18),
        axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        panel.spacing=grid::unit(0.8,"lines"))

png(file.path(fig5e_dir,"Figure5E_imbalanced_ablation_beeswarm.png"),
    width=13,height=7.5,units="in",res=300)
p_imbalanced
dev.off()

# ============================================================
# 44. Main metrics used for component importance ranking
#
# Balanced:
# AUROC / F1 / ACC
#
# Imbalanced:
# NormAUPRC / AUROC / F1
#
# Precision and Recall are shown in Figure 5E but are not used
# for overall ranking to avoid overweighting threshold metrics.
# ============================================================
delta_for_ranking <- delta_performance %>%
  filter(
    (setting=="balanced" & metric %in% c("AUROC","F1","ACC")) |
      (setting %in% c("1_3","1_30") & metric %in% c("NormAUPRC","AUROC","F1"))
  ) %>%
  filter(!is.na(delta))

# Future use:
# label_delta should contain:
# element, setting, round, dataset, drug, metric,
# full_value, ablated_value, delta
# Then:
# delta_for_ranking <- bind_rows(delta_for_ranking,label_delta)

# ============================================================
# 45. Mean delta for each element × condition
# ============================================================
importance_mean <- delta_for_ranking %>%
  group_by(element,setting,metric) %>%
  summarise(
    mean_delta=mean(delta),
    median_delta=median(delta),
    n=n(),
    .groups="drop"
  ) %>%
  mutate(
    panel=case_when(
      setting=="balanced" & metric=="AUROC" ~ "Balanced\nAUROC",
      setting=="balanced" & metric=="F1" ~ "Balanced\nF1",
      setting=="balanced" & metric=="ACC" ~ "Balanced\nACC",
      setting=="1_3" & metric=="NormAUPRC" ~ "1:3\nNormAUPRC",
      setting=="1_3" & metric=="AUROC" ~ "1:3\nAUROC",
      setting=="1_3" & metric=="F1" ~ "1:3\nF1",
      setting=="1_30" & metric=="NormAUPRC" ~ "1:30\nNormAUPRC",
      setting=="1_30" & metric=="AUROC" ~ "1:30\nAUROC",
      setting=="1_30" & metric=="F1" ~ "1:30\nF1"
    ),
    panel=factor(panel,levels=c(
      "Balanced\nAUROC","Balanced\nF1","Balanced\nACC",
      "1:3\nNormAUPRC","1:3\nAUROC","1:3\nF1",
      "1:30\nNormAUPRC","1:30\nAUROC","1:30\nF1"
    ))
  )

# ============================================================
# 46. Importance heatmap
# ============================================================
p_heatmap <- ggplot(importance_mean,aes(x=panel,y=element,fill=mean_delta)) +
  geom_tile(color="white",linewidth=0.8) +
  geom_text(aes(label=sprintf("%.3f",mean_delta)),color="black",size=6,fontface="bold") +
  scale_fill_gradient2(low="#3B6FB6",mid="white",high="#C94343",
                       midpoint=0,name=expression(mean~Delta)) +
  labs(x=NULL,y=NULL) +
  theme_minimal(base_size=14) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(face="bold"),
    axis.text.y=element_text(size=16,face="bold"),
    legend.position="right"
  )

png(file.path(fig5e_dir,"Figure5F_component_heatmap.png"),
    width=11,height=4.5,units="in",res=300)
p_heatmap
dev.off()

# ============================================================
# 47. Rank elements separately within each metric × setting
# Rank 1 = largest mean delta = most important
# ============================================================
importance_rank <- importance_mean %>%
  group_by(panel) %>%
  mutate(rank=rank(-mean_delta,ties.method="average")) %>%
  ungroup()

overall_rank <- importance_rank %>%
  group_by(element) %>%
  summarise(
    mean_rank=mean(rank),
    median_rank=median(rank),
    mean_delta_across_panels=mean(mean_delta),
    .groups="drop"
  ) %>%
  arrange(mean_rank) %>%
  mutate(overall_order=row_number())

print(importance_rank,n=Inf)
print(overall_rank,n=Inf)

write_csv(importance_rank,
          "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_element_rank_by_condition.csv")
write_csv(overall_rank,
          "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_element_overall_rank.csv")

# ============================================================
# 48. Overall mean-rank dot plot
# Smaller mean rank = more important
# ============================================================
overall_rank_plot <- overall_rank %>%
  arrange(mean_rank) %>%
  mutate(element=factor(element,levels=rev(as.character(element))))

rank_xmax <- max(overall_rank_plot$mean_rank)+0.8

p_rank <- ggplot(overall_rank_plot,aes(x=mean_rank,y=element)) +
  geom_segment(aes(x=1,xend=mean_rank,y=element,yend=element),
               linewidth=0.8,color="grey70") +
  geom_point(size=4,color="black") +
  geom_text(aes(label=sprintf("%.2f",mean_rank)),
            nudge_x=0.12,hjust=0,size=5,fontface="bold") +
  geom_text(aes(label=as.character(element)),
            x=0.95,nudge_y=0.30,hjust=0,size=4.8,fontface="bold") +
  scale_x_continuous(
    limits=c(0.8,rank_xmax),
    breaks=seq_len(ceiling(max(overall_rank_plot$mean_rank))),
    expand=expansion(mult=c(0,0))
  ) +
  scale_y_discrete(expand=expansion(add=c(0.45,0.75))) +
  coord_cartesian(clip="off") +
  labs(x="Mean rank across conditions",y=NULL) +
  theme_classic(base_size=20) +
  theme(
    axis.text.x=element_text(size=14,color="black"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_text(size=16,face="bold",margin=margin(t=10)),
    plot.margin=margin(t=20,r=30,b=10,l=5)
  )

png(file.path(fig5e_dir,"Figure5F_overall_mean_rank.png"),
    width=6,height=3.8,units="in",res=300)
p_rank
dev.off()

# ============================================================
# 49. Combined Figure 5F
# ============================================================
if(requireNamespace("patchwork",quietly=TRUE)){
  p5f_combined <- p_heatmap+p_rank+
    patchwork::plot_layout(widths=c(2.1,1))
  
  png(file.path(fig5e_dir,"Figure5F_heatmap_and_rank.png"),
      width=15,height=5,units="in",res=300)
  print(p5f_combined)
  dev.off()
}
# ============================================================
# 50. Fraction of observations with delta > 0
# ============================================================
positive_contribution <- delta_for_ranking %>%
  group_by(element,setting,metric) %>%
  summarise(
    n=n(),
    positive_n=sum(delta>0),
    positive_fraction=mean(delta>0),
    mean_delta=mean(delta),
    median_delta=median(delta),
    .groups="drop"
  ) %>%
  arrange(metric,setting,desc(mean_delta))

write_csv(positive_contribution,
          "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_positive_contribution_fraction.csv")
print(positive_contribution,n=Inf)

# ============================================================
# 51. Read Figure 5C data and calculate label-definition delta
# Delta = scDEAL label - SCAD label under scDEAL architecture
# ============================================================
fig5c_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5C/summary_metrics/Figure5C_plotting_data_long.csv"

fig5c_data <- readr::read_csv(fig5c_file,show_col_types=FALSE)

required_fig5c_cols <- c("Architecture","Training_label","Scenario","round",
                         "dataset","drug","Metric","Value")
missing_fig5c_cols <- setdiff(required_fig5c_cols,colnames(fig5c_data))

if(length(missing_fig5c_cols)>0){
  stop(paste("Missing columns in Figure5C_plotting_data_long.csv:",
             paste(missing_fig5c_cols,collapse=", ")))
}

mapping_importance <- mapping %>%
  transmute(
    dataset_key=tolower(trimws(as.character(dataset))),
    dataset_standard=trimws(as.character(dataset)),
    drug_standard=trimws(as.character(drug))
  ) %>%
  distinct()

label_long <- fig5c_data %>%
  filter(
    Architecture=="scDEAL",
    Training_label%in%c("SCAD label","scDEAL label"),
    !is.na(Value)
  ) %>%
  mutate(
    dataset_key=tolower(trimws(as.character(dataset))),
    original_drug=trimws(as.character(drug))
  ) %>%
  inner_join(mapping_importance,by="dataset_key") %>%
  transmute(
    setting=case_when(
      Scenario=="Balanced" ~ "balanced",
      Scenario=="1:3" ~ "1_3",
      Scenario=="1:30" ~ "1_30",
      TRUE ~ NA_character_
    ),
    round=as.integer(round),
    dataset=dataset_standard,
    drug=drug_standard,
    metric=as.character(Metric),
    Training_label=as.character(Training_label),
    value=as.numeric(Value)
  ) %>%
  filter(
    (setting=="balanced" & metric%in%c("AUROC","F1","ACC")) |
      (setting%in%c("1_3","1_30") & metric%in%c("NormAUPRC","AUROC","F1"))
  )

label_duplicate_check <- label_long %>%
  count(setting,round,dataset,drug,metric,Training_label,name="n") %>%
  filter(n>1)

if(nrow(label_duplicate_check)>0){
  print(label_duplicate_check,n=Inf)
  stop("Duplicated label-performance records were found.")
}

label_delta <- label_long %>%
  pivot_wider(names_from=Training_label,values_from=value) %>%
  filter(!is.na(`SCAD label`),!is.na(`scDEAL label`)) %>%
  transmute(
    element="Label definition",
    setting=setting,
    round=round,
    dataset=dataset,
    drug=drug,
    metric=metric,
    full_value=`scDEAL label`,
    ablated_value=`SCAD label`,
    delta=`scDEAL label`-`SCAD label`
  )

write_csv(
  label_delta,
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_label_definition_delta.csv"
)

print(label_delta,n=30)

# ============================================================
# 52. Combine label definition and architecture components
# Keep only observations available for all five elements
# ============================================================
component_delta_for_ranking <- delta_for_ranking %>%
  transmute(
    element=as.character(element),
    setting=as.character(setting),
    round=as.integer(round),
    dataset=trimws(as.character(dataset)),
    drug=trimws(as.character(drug)),
    metric=as.character(metric),
    full_value=as.numeric(full_value),
    ablated_value=as.numeric(ablated_value),
    delta=as.numeric(delta)
  )

delta_for_ranking_all_raw <- bind_rows(
  component_delta_for_ranking,
  label_delta
)

importance_key <- c("setting","round","dataset","drug","metric")

comprehensive_duplicate_check <- delta_for_ranking_all_raw %>%
  count(across(all_of(importance_key)),element,name="n") %>%
  filter(n>1)

if(nrow(comprehensive_duplicate_check)>0){
  print(comprehensive_duplicate_check,n=Inf)
  stop("Duplicated element records were found in the comprehensive importance data.")
}

complete_importance_keys <- delta_for_ranking_all_raw %>%
  distinct(across(all_of(importance_key)),element) %>%
  count(across(all_of(importance_key)),name="n_elements") %>%
  filter(n_elements==5) %>%
  select(all_of(importance_key))

incomplete_importance_keys <- delta_for_ranking_all_raw %>%
  distinct(across(all_of(importance_key)),element) %>%
  count(across(all_of(importance_key)),name="n_elements") %>%
  filter(n_elements<5)

if(nrow(complete_importance_keys)==0){
  stop("No observations were matched across all five importance elements.")
}

delta_for_ranking_all <- delta_for_ranking_all_raw %>%
  semi_join(complete_importance_keys,by=importance_key) %>%
  mutate(
    element=factor(
      element,
      levels=c(
        "Label definition",
        "DAE pretraining",
        "Cluster regularization",
        "Transfer learning",
        "Source-data integration"
      )
    ),
    setting=factor(setting,levels=c("balanced","1_3","1_30"))
  ) %>%
  arrange(setting,round,dataset,drug,metric,element)

write_csv(
  delta_for_ranking_all,
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_comprehensive_importance_delta.csv"
)

write_csv(
  incomplete_importance_keys,
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_comprehensive_incomplete_keys.csv"
)

cat("\nComplete matched observations:",nrow(complete_importance_keys),"\n")
cat("Incomplete observations removed:",nrow(incomplete_importance_keys),"\n")

# ============================================================
# 53. Mean delta for each element × condition
# ============================================================
importance_mean_all <- delta_for_ranking_all %>%
  group_by(element,setting,metric) %>%
  summarise(
    mean_delta=mean(delta),
    median_delta=median(delta),
    n=n(),
    .groups="drop"
  ) %>%
  mutate(
    panel=case_when(
      setting=="balanced" & metric=="AUROC" ~ "Balanced\nAUROC",
      setting=="balanced" & metric=="F1" ~ "Balanced\nF1",
      setting=="balanced" & metric=="ACC" ~ "Balanced\nACC",
      setting=="1_3" & metric=="NormAUPRC" ~ "1:3\nNormAUPRC",
      setting=="1_3" & metric=="AUROC" ~ "1:3\nAUROC",
      setting=="1_3" & metric=="F1" ~ "1:3\nF1",
      setting=="1_30" & metric=="NormAUPRC" ~ "1:30\nNormAUPRC",
      setting=="1_30" & metric=="AUROC" ~ "1:30\nAUROC",
      setting=="1_30" & metric=="F1" ~ "1:30\nF1"
    ),
    panel=factor(
      panel,
      levels=c(
        "Balanced\nAUROC","Balanced\nF1","Balanced\nACC",
        "1:3\nNormAUPRC","1:3\nAUROC","1:3\nF1",
        "1:30\nNormAUPRC","1:30\nAUROC","1:30\nF1"
      )
    )
  )

importance_condition_check <- importance_mean_all %>%
  count(panel,name="n_elements")

if(any(importance_condition_check$n_elements!=5)){
  print(importance_condition_check,n=Inf)
  stop("Not every condition contains all five importance elements.")
}

write_csv(
  importance_mean_all,
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_comprehensive_importance_mean.csv"
)

# ============================================================
# 54. Rank five elements within each condition
# Rank 1 = largest mean delta = most important
# ============================================================
importance_rank_all <- importance_mean_all %>%
  group_by(panel) %>%
  mutate(rank=rank(-mean_delta,ties.method="average")) %>%
  ungroup()

overall_rank_all <- importance_rank_all %>%
  group_by(element) %>%
  summarise(
    mean_rank=mean(rank),
    median_rank=median(rank),
    mean_delta_across_panels=mean(mean_delta),
    .groups="drop"
  ) %>%
  arrange(mean_rank) %>%
  mutate(overall_order=row_number())

element_order <- as.character(overall_rank_all$element)

importance_mean_all <- importance_mean_all %>%
  mutate(element=factor(as.character(element),levels=rev(element_order)))

overall_rank_all <- overall_rank_all %>%
  mutate(element=factor(as.character(element),levels=rev(element_order)))

write_csv(
  importance_rank_all,
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_comprehensive_rank_by_condition.csv"
)

write_csv(
  overall_rank_all,
  "/media/shenli/Expansion/DRMREF/New_fig5/Figure5E/scDEAL_comprehensive_overall_rank.csv"
)

print(importance_rank_all,n=Inf)
print(overall_rank_all,n=Inf)

# ============================================================
# 55. Comprehensive importance heatmap
# ============================================================
p_heatmap_all <- ggplot(importance_mean_all,aes(x=panel,y=element,fill=mean_delta)) +
  geom_tile(color="white",linewidth=0.8) +
  geom_text(aes(label=sprintf("%.3f",mean_delta)),
            color="black",size=5.5,fontface="bold") +
  scale_fill_gradient2(
    low="#3B6FB6",mid="white",high="#C94343",
    midpoint=0,name=expression(mean~Delta)
  ) +
  labs(x=NULL,y=NULL) +
  theme_minimal(base_size=14) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(size=12,face="bold"),
    axis.text.y=element_text(size=15,face="bold"),
    legend.position="right"
  )

png(file.path(fig5e_dir,"Figure5F_comprehensive_importance_heatmap.png"),
    width=11.5,height=5.2,units="in",res=300)
print(p_heatmap_all)
dev.off()

# ============================================================
# 56. Comprehensive overall mean-rank plot
# ============================================================
rank_xmax_all <- max(overall_rank_all$mean_rank)+0.8

p_rank_all <- ggplot(overall_rank_all,aes(x=mean_rank,y=element)) +
  geom_segment(
    aes(x=1,xend=mean_rank,y=element,yend=element),
    linewidth=0.8,color="grey70"
  ) +
  geom_point(size=4,color="black") +
  geom_text(
    aes(label=sprintf("%.2f",mean_rank)),
    nudge_x=0.12,hjust=0,size=5,fontface="bold"
  ) +
  geom_text(
    aes(label=as.character(element)),
    x=0.95,nudge_y=0.18,hjust=0,vjust=0,
    size=4.5,fontface="bold"
  ) +
  scale_x_continuous(
    limits=c(0.8,rank_xmax_all),
    breaks=seq_len(ceiling(max(overall_rank_all$mean_rank))),
    expand=expansion(mult=c(0,0))
  ) +
  scale_y_discrete(expand=expansion(add=c(0.45,0.85))) +
  coord_cartesian(clip="off") +
  labs(x="Mean rank across conditions",y=NULL) +
  theme_classic(base_size=20) +
  theme(
    axis.text.x=element_text(size=14,color="black"),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    axis.line.y=element_blank(),
    axis.title.x=element_text(size=16,face="bold",margin=margin(t=10)),
    plot.margin=margin(t=22,r=30,b=10,l=5)
  )

png(file.path(fig5e_dir,"Figure5F_comprehensive_mean_rank.png"),
    width=6.5,height=4.8,units="in",res=300)
print(p_rank_all)
dev.off()

# ============================================================
# 56.2 Rank matrix across nine metric × setting conditions
# ============================================================

rank_element_order <- overall_rank_all %>%
  arrange(mean_rank) %>%
  pull(element) %>%
  as.character()

element_label_table <- overall_rank_all %>%
  transmute(
    element=as.character(element),
    element_label=paste0(element,"  [mean rank = ",sprintf("%.2f",mean_rank),"]")
  )

element_label_levels <- element_label_table %>%
  mutate(element=factor(element,levels=rank_element_order)) %>%
  arrange(element) %>%
  pull(element_label)

rank_matrix_data <- importance_rank_all %>%
  mutate(
    element=as.character(element),
    setting=as.character(setting),
    metric=as.character(metric),
    setting_group=case_when(
      setting=="balanced" ~ "Balanced",
      setting=="1_3" ~ "1:3",
      setting=="1_30" ~ "1:30"
    ),
    setting_group=factor(setting_group,levels=c("Balanced","1:3","1:30")),
    condition=case_when(
      setting=="balanced" & metric=="AUROC" ~ "01_AUROC",
      setting=="balanced" & metric=="F1" ~ "02_F1",
      setting=="balanced" & metric=="ACC" ~ "03_ACC",
      setting%in%c("1_3","1_30") & metric=="NormAUPRC" ~ "01_NormAUPRC",
      setting%in%c("1_3","1_30") & metric=="AUROC" ~ "02_AUROC",
      setting%in%c("1_3","1_30") & metric=="F1" ~ "03_F1"
    ),
    condition=factor(
      condition,
      levels=c("01_AUROC","01_NormAUPRC","02_F1","02_AUROC","03_ACC","03_F1")
    ),
    rank_group=factor(paste0("Rank ",as.integer(rank)),levels=paste0("Rank ",1:5)),
    rank_label=ifelse(rank==round(rank),sprintf("%.0f",rank),sprintf("%.1f",rank)),
    text_color=ifelse(rank<=2,"white","black")
  ) %>%
  left_join(element_label_table,by="element") %>%
  mutate(element_label=factor(element_label,levels=rev(element_label_levels)))

rank_colors <- c(
  "Rank 1"="#244A73",
  "Rank 2"="#4F79A7",
  "Rank 3"="#8EADD0",
  "Rank 4"="#C9D7E6",
  "Rank 5"="#ECEFF2"
)

p_rank_matrix <- ggplot(rank_matrix_data,aes(x=condition,y=element_label)) +
  geom_point(
    aes(fill=rank_group),
    shape=22,size=13,color="white",stroke=0.8
  ) +
  geom_text(
    aes(label=rank_label,color=text_color),
    size=5,fontface="bold"
  ) +
  facet_grid(.~setting_group,scales="free_x",space="free_x") +
  scale_x_discrete(
    labels=function(x) sub("^\\d+_","",x),
    expand=expansion(add=0.55)
  ) +
  scale_y_discrete(expand=expansion(add=0.5)) +
  scale_fill_manual(values=rank_colors,drop=FALSE) +
  scale_color_identity() +
  labs(
    x=NULL,
    y=NULL
  ) +
  theme_minimal(base_size=15) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(size=12,face="bold",angle=25,hjust=1,color="black"),
    axis.text.y=element_text(size=14,face="bold",color="black"),
    strip.background=element_rect(fill="grey95",color="grey40",linewidth=0.6),
    strip.text=element_text(size=15,face="bold"),
    panel.spacing=grid::unit(0.7,"lines"),
    plot.subtitle=element_text(size=13,hjust=0.5,margin=margin(b=12)),
    plot.margin=margin(10,15,0,10),
    legend.position="none"
  )

rank_legend_data <- tibble(
  rank=1:5,
  rank_group=factor(paste0("Rank ",1:5),levels=paste0("Rank ",1:5))
)

p_rank_legend <- ggplot(rank_legend_data,aes(x=rank,y=0,fill=rank_group)) +
  geom_tile(width=1,height=0.42,color="grey30",linewidth=0.4) +
  geom_segment(
    aes(x=0.5,xend=5.5,y=0.43,yend=0.43),
    inherit.aes=FALSE,
    linewidth=0.5,
    arrow=grid::arrow(length=grid::unit(0.13,"cm"),type="closed")
  ) +
  annotate("text",x=0.42,y=0.43,label="1",hjust=1,vjust=0.5,size=4.5) +
  annotate("text",x=5.58,y=0.43,label="5",hjust=0,vjust=0.5,size=4.5) +
  annotate(
    "text",
    x=-0.1,y=0,
    label="Importance rank",
    hjust=1,vjust=0.5,
    size=4.8,fontface="bold"
  ) +
  scale_fill_manual(values=rank_colors,guide="none") +
  scale_x_continuous(limits=c(-1.8,5.9),expand=expansion(mult=c(0,0))) +
  coord_cartesian(ylim=c(-0.3,0.68),clip="off") +
  theme_void() +
  theme(plot.margin=margin(0,0,0,0))

# ============================================================
# 57. Combined comprehensive importance figure
# ============================================================
if(requireNamespace("patchwork",quietly=TRUE)){
  p_rank_matrix_final <- (p_rank_matrix/p_rank_legend) +
    patchwork::plot_layout(heights=c(1,0.18))
  
  png(file.path(fig5e_dir,"Figure5F_comprehensive_rank_matrix.png"),
      width=9,height=4.5,units="in",res=300)
  print(p_rank_matrix_final)
  dev.off()
}
