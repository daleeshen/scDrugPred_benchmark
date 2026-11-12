library(Seurat)
library(dplyr)
library(pROC)
library(caret)
library(stringr)

root_path <- "/set/imblanced/results/root/path/"
raw_data_path <- "/path/to/raw/data/"
dataset_drug_info <- read.csv("/media/shenli/Expansion/dataset_drug_single.csv", header = TRUE)
cell_list <- c(500,1000,2000,3000,5000)
round_list <- c(1,2,3)

load_meta <- function(path) {
  e <- new.env()
  obj_names <- load(path, envir = e)
  for (x in obj_names) {
    obj <- e[[x]]
    if (inherits(obj, "Seurat")) return(obj@meta.data)
  }
  stop("No Seurat object found: ", path)
}
mean_na <- function(x) if(all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
infer_recall <- function(precision, f1) {
  recall <- ifelse(is.na(precision) | is.na(f1) | (2 * precision - f1) == 0, NA, (f1 * precision)/(2 * precision - f1))
  ifelse(recall < 0 | recall > 1, NA, recall)
}
extract_metric <- function(txt, pattern) {
  x <- str_match(txt, pattern)[,2]
  if(length(x) == 0 || is.na(x[1])) NA_real_ else as.numeric(x[1])
}

#### Beyondcell ####
beyondcell_res_df <- data.frame()
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "beyondcell")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    for (csv_file in csv_files) {
      data_id <- sub("_beyondcell_res\\.csv$", "", basename(csv_file))
      pred_df <- read.csv(csv_file, header = TRUE)
      sc_data <- file.path(root_path, paste0(round, "_round_raw_", cells, "cells"), paste0("balanced_sampled_", data_id, "_seurat_afterAnno.RData"))
      meta_data <- load_meta(sc_data)
      rownames(pred_df) <- pred_df[,1]
      common_cells <- intersect(rownames(pred_df), rownames(meta_data))
      pred_df <- pred_df[common_cells,,drop = FALSE]
      pred_df$condition <- meta_data[common_cells,"condition"]
      sig_cols <- grep("^sig_", colnames(pred_df), value = TRUE)
      sig_cols <- sig_cols[!grepl("_pred_condition$", sig_cols)]
      for (sig_col in sig_cols) {
        pred_col <- paste0(sig_col, "_pred_condition")
        if (!(pred_col %in% colnames(pred_df))) next
        score <- pred_df[[sig_col]]
        pred_label <- factor(pred_df[[pred_col]], levels = c("resistant","sensitive"))
        true_label <- factor(pred_df$condition, levels = c("resistant","sensitive"))
        auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("resistant","sensitive"), quiet = TRUE)))
        cm <- caret::confusionMatrix(pred_label, true_label, positive = "sensitive")
        beyondcell_res_df <- rbind(beyondcell_res_df, data.frame(round = round, cells = cells, dataset = data_id, signature = sig_col, AUC = auc_value, ACC = as.numeric(cm$overall["Accuracy"]), Precision = as.numeric(cm$byClass["Precision"]), Recall = as.numeric(cm$byClass["Recall"]), F1 = as.numeric(cm$byClass["F1"])))
      }
    }
  }
}
beyondcell_best_df <- beyondcell_res_df %>%
  filter(!is.na(AUC)) %>%
  group_by(round, cells, dataset) %>%
  slice_max(order_by = AUC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### CaDRReS-sc ####
cadrres_sc_res_df <- data.frame()
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "CaDRReS-Sc")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    for (csv_file in csv_files) {
      data_id <- sub("_cadrres.*\\.csv$", "", basename(csv_file))
      pred_df <- read.csv(csv_file, header = TRUE, check.names = FALSE)
      sc_data <- file.path(root_path, paste0(round, "_round_raw_", cells, "cells"), paste0("balanced_sampled_", data_id, "_seurat_afterAnno.RData"))
      meta_data <- load_meta(sc_data)
      rownames(meta_data) <- gsub("\\.","-",rownames(meta_data))
      rownames(pred_df) <- gsub("\\.","-",pred_df[,1])
      common_cells <- intersect(rownames(pred_df), rownames(meta_data))
      pred_df <- pred_df[common_cells,,drop = FALSE]
      tar_drug_id <- dataset_drug_info[dataset_drug_info$Dataset == data_id,5]
      if(length(tar_drug_id) == 0 || is.na(tar_drug_id)) next
      tar_drug_id <- as.character(tar_drug_id)
      if (!(tar_drug_id %in% colnames(pred_df))) next
      pred_df <- pred_df[,c(1,which(colnames(pred_df) == tar_drug_id)),drop = FALSE]
      names(pred_df)[1] <- "cell_id"
      pred_df$condition <- meta_data[common_cells,"condition"]
      score <- pred_df[[tar_drug_id]]
      pred_df$pred_label <- ifelse(score > median(score, na.rm = TRUE), "resistant", "sensitive")
      true_label <- factor(pred_df$condition, levels = c("resistant","sensitive"))
      pred_label <- factor(pred_df$pred_label, levels = c("resistant","sensitive"))
      auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("resistant","sensitive"), quiet = TRUE)))
      cm <- caret::confusionMatrix(pred_label, true_label, positive = "sensitive")
      cadrres_sc_res_df <- rbind(cadrres_sc_res_df, data.frame(round = round, cells = cells, dataset = data_id, AUC = auc_value, ACC = as.numeric(cm$overall["Accuracy"]), Precision = as.numeric(cm$byClass["Precision"]), Recall = as.numeric(cm$byClass["Recall"]), F1 = as.numeric(cm$byClass["F1"])))
    }
  }
}
cadrres_sc_best_df <- cadrres_sc_res_df %>%
  group_by(round, cells, dataset) %>%
  slice_max(order_by = AUC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### DREEP ####
dreep_res_df <- data.frame()
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "DREEP")
    res_file <- file.path(res_folder_path, "dreep_auc.csv")
    if (!file.exists(res_file)) next
    dr_res <- read.csv(res_file, header = TRUE)
    dr_res$round <- round
    dr_res$cells <- cells
    dreep_res_df <- rbind(dreep_res_df, dr_res)
  }
}
colnames(dreep_res_df) <- tolower(colnames(dreep_res_df))
if (!("recall" %in% colnames(dreep_res_df))) dreep_res_df$recall <- infer_recall(dreep_res_df$precision, dreep_res_df$f1)
dreep_best_df <- dreep_res_df %>%
  group_by(round, cells, dataset) %>%
  slice_max(order_by = auc, n = 1, with_ties = FALSE) %>%
  ungroup()

#### DrugFormer ####
drugformer_res_df <- data.frame()
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "DrugFormer")
    log_files <- list.files(res_folder_path, pattern = "\\.txt$", full.names = TRUE)
    for (log_file in log_files) {
      data_id <- sub("_log\\.txt$", "", basename(log_file))
      test_lines <- grep("Testing", readLines(log_file), value = TRUE)
      if(length(test_lines) == 0) next
      get_metric <- function(pattern) as.numeric(sub(paste0(".*", pattern, ": ([0-9\\.eE+-]+).*"), "\\1", test_lines))
      acc <- get_metric("Accuracy")
      auc <- get_metric("AUROC")
      prec <- get_metric("Precision")
      f1 <- get_metric("F1")
      recall <- get_metric("Recall")
      avg <- colMeans(data.frame(acc,auc,prec,f1,recall), na.rm = TRUE)
      drugformer_res_df <- rbind(drugformer_res_df, data.frame(round = round, cells = cells, dataset = data_id, AUC = avg[2], ACC = avg[1], Precision = avg[3], F1 = avg[4], Recall = avg[5]))
    }
  }
}
drugformer_res_df <- na.omit(drugformer_res_df)

#### Precily ####
precily_sc_res_df <- data.frame()
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "precily", "pred_res")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    for (csv_file in csv_files) {
      data_id <- sub("_pred_res\\.csv$", "", basename(csv_file))
      pred_df <- read.csv(csv_file, header = TRUE, check.names = FALSE)
      rownames(pred_df) <- make.unique(pred_df[,1])
      pred_df <- as.data.frame(t(pred_df[,-1]))
      pred_df$cell_id <- rownames(pred_df)
      rownames(pred_df) <- gsub("\\.","-",rownames(pred_df))
      tar_drug_name <- dataset_drug_info[dataset_drug_info$Dataset == data_id,2]
      if(length(tar_drug_name) == 0 || is.na(tar_drug_name)) next
      sc_data <- file.path(root_path, paste0(round, "_round_raw_", cells, "cells"), paste0("balanced_sampled_", data_id, "_seurat_afterAnno.RData"))
      meta_data <- load_meta(sc_data)
      rownames(meta_data) <- gsub("\\.","-",rownames(meta_data))
      common_cells <- intersect(rownames(pred_df), rownames(meta_data))
      pred_df <- pred_df[common_cells,,drop = FALSE]
      if (!(tar_drug_name %in% colnames(pred_df))) next
      pred_df <- pred_df[,c("cell_id",tar_drug_name)]
      pred_df$condition <- meta_data[common_cells,"condition"]
      raw_score <- pred_df[[tar_drug_name]]
      score <- -raw_score
      pred_df$pred_label <- ifelse(raw_score > median(raw_score, na.rm = TRUE), "resistant", "sensitive")
      true_label <- factor(pred_df$condition, levels = c("resistant","sensitive"))
      pred_label <- factor(pred_df$pred_label, levels = c("resistant","sensitive"))
      auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("resistant","sensitive"), quiet = TRUE)))
      cm <- caret::confusionMatrix(pred_label, true_label, positive = "sensitive")
      precily_sc_res_df <- rbind(precily_sc_res_df, data.frame(round = round, cells = cells, dataset = data_id, AUC = auc_value, ACC = as.numeric(cm$overall["Accuracy"]), Precision = as.numeric(cm$byClass["Precision"]), Recall = as.numeric(cm$byClass["Recall"]), F1 = as.numeric(cm$byClass["F1"])))
    }
  }
}

#### SCAD: summary logs ####
scad_res_df <- data.frame()
sub_folder_list <- c("best_param_res_smote_norm","best_param_res_smote_ppi","best_param_res_smote_tp4k","best_param_res_weight_norm","best_param_res_weight_ppi","best_param_res_weight_tp4k")
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "SCAD")
    for (sub_folder in sub_folder_list) {
      mode_suffix <- sub("^best_param_res_", "", sub_folder)
      sub_folder_file_path <- file.path(res_folder_path, sub_folder)
      if (!dir.exists(sub_folder_file_path)) next
      log_files <- list.files(sub_folder_file_path, pattern = "\\.txt$", full.names = TRUE)
      for (log_file in log_files) {
        data_id_with_drug <- sub("_train_test_summary_5folds_5seeds\\.txt$", "", basename(log_file))
        data_id <- sub("_[^_]+$", "", data_id_with_drug)
        log_lines <- readLines(log_file, warn = FALSE)
        auc_line <- grep("Average Test AUC =", log_lines, value = TRUE)
        prec_line <- grep("Average Test Precision Recall =", log_lines, value = TRUE)
        if(length(auc_line) == 0) next
        auc <- as.numeric(sub(".*Average Test AUC = ([0-9.]+);.*", "\\1", auc_line[1]))
        precision <- if(length(prec_line) == 0) NA_real_ else as.numeric(sub(".*Average Test Precision Recall = ([0-9.]+);.*", "\\1", prec_line[1]))
        scad_res_df <- rbind(scad_res_df, data.frame(round = round, cells = cells, dataset = data_id, AUC = auc, Precision = precision, mode = paste0("SCAD_", mode_suffix)))
      }
    }
  }
}
scad_res_df_unique <- scad_res_df %>% distinct(dataset, round, cells, mode, .keep_all = TRUE)

#### SCAD: prediction tensors ####
parse_tensor <- function(txt) {
  txt <- gsub("([0-9])\\s*[eE]\\s*([+-]?\\d+)", "\\1e\\2", txt)
  as.numeric(unlist(str_extract_all(txt, "-?\\d*\\.?\\d+(?:[eE][+-]?\\d+)?")))
}

scad_pred_res_df <- data.frame()
log_file_list <- list.files(file.path(root_path, "balanced_SCAD_log_file"), pattern = "\\.out$", full.names = TRUE)
for (log_file in log_file_list) {
  log_lines <- readLines(log_file, warn = FALSE)
  dataset_starts <- grep("^开始处理数据集:", log_lines)
  if(length(dataset_starts) == 0) next
  dataset_ends <- c(dataset_starts[-1] - 1, length(log_lines))
  for (idx in seq_along(dataset_starts)) {
    block <- log_lines[dataset_starts[idx]:dataset_ends[idx]]
    dataset_name <- sub("^开始处理数据集:\\s*", "", block[1])
    dir_idx <- grep("Directory .*best_param_res_", block)
    eval_idx <- grep("^-- Evaluation --\\s*$", block)
    test_idx <- grep("^-- Test Results --\\s*$", block)
    if(length(dir_idx) == 0 || length(eval_idx) == 0) next
    for (k in seq_along(eval_idx)) {
      e_line <- eval_idx[k]
      pos <- max(which(dir_idx <= e_line), na.rm = TRUE)
      if(!is.finite(pos)) next
      dir_line <- block[dir_idx[pos]]
      mode <- str_match(dir_line, "best_param_res_([^/]+?)_(norm|ppi|tp4k)/")
      rc <- str_match(dir_line, "benchmark_data/([0-9]+)_round_raw_([0-9]+)cells")
      if(is.na(mode[1,1]) || is.na(rc[1,1])) next
      mode_name <- paste0("SCAD_", mode[1,2], "_", mode[1,3])
      round_num <- as.numeric(rc[1,2])
      cell_num <- as.numeric(rc[1,3])
      t_line <- test_idx[test_idx > e_line][1]
      if(is.na(t_line) || t_line <= e_line + 1) next
      section_text <- paste(block[(e_line + 1):(t_line - 1)], collapse = " ")
      tensor_pair <- str_match(section_text, "tensor\\((.*?)\\)\\s*and\\s*tensor\\((.*?)\\)")
      if(is.na(tensor_pair[1,1])) next
      y_true <- parse_tensor(tensor_pair[1,2])
      y_pred <- parse_tensor(tensor_pair[1,3])
      if(length(y_pred) - length(y_true) == 1) y_pred <- y_pred[-length(y_pred)]
      if(length(y_true) != length(y_pred) || length(y_true) == 0) next
      valid <- !is.na(y_true) & !is.na(y_pred)
      y_true <- y_true[valid]
      y_pred <- y_pred[valid]
      if(length(unique(y_true)) < 2) next
      y_pred_label <- ifelse(y_pred > median(y_pred, na.rm = TRUE), 1, 0)
      cm <- caret::confusionMatrix(factor(y_pred_label, levels = c(0,1)), factor(y_true, levels = c(0,1)), positive = "1")
      scad_pred_res_df <- rbind(scad_pred_res_df, data.frame(dataset = dataset_name, round = round_num, cells = cell_num, mode = mode_name, ACC = as.numeric(cm$overall["Accuracy"]), F1 = as.numeric(cm$byClass["F1"]), Recall = as.numeric(cm$byClass["Recall"])))
    }
  }
}

scad_pred_res_avg <- scad_pred_res_df %>%
  group_by(dataset, round, cells, mode) %>%
  summarise(ACC = mean_na(ACC), F1 = mean_na(F1), Recall = mean_na(Recall), .groups = "drop")

scad_final_df <- scad_res_df_unique %>%
  left_join(scad_pred_res_avg, by = c("dataset","round","cells","mode"))

scad_best_df <- scad_final_df %>%
  filter(!is.na(AUC)) %>%
  group_by(round, cells, dataset) %>%
  slice_max(order_by = AUC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### scDEAL ####
scdeal_res_df <- data.frame()
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "scDEAL")
    log_files <- list.files(res_folder_path, pattern = "\\.log$", full.names = TRUE)
    for (log_file in log_files) {
      log_name <- basename(log_file)
      dataset_id <- str_match(log_name, "^(.+?)_[^_]+_[^_]+_bottleneck\\d+_sampling[^_]+\\.log$")[,2]
      bottleneck <- str_match(log_name, "_bottleneck(\\d+)_")[,2]
      sampling <- str_match(log_name, "_sampling([^_\\.]+)")[,2]
      txt <- readChar(log_file, file.info(log_file)$size)
      auc <- extract_metric(txt, "auroc_score\\s+([0-9\\.eE+-]+)")
      acc <- extract_metric(txt, "Accuracy:\\s*([0-9\\.eE+-]+)")
      precision <- extract_metric(txt, "precision\\s*([0-9\\.eE+-]+)")
      recall <- extract_metric(txt, "recall\\s*([0-9\\.eE+-]+)")
      f1 <- extract_metric(txt, "f1score\\s*([0-9\\.eE+-]+)")
      if(is.na(recall) && !is.na(precision) && !is.na(f1)) recall <- infer_recall(precision, f1)
      scdeal_res_df <- rbind(scdeal_res_df, data.frame(round = round, cells = cells, bottleneck = bottleneck, sampling_method = sampling, dataset = dataset_id, AUC = auc, ACC = acc, Precision = precision, Recall = recall, F1 = f1))
    }
  }
}
scdeal_best_df <- scdeal_res_df %>%
  filter(!is.na(AUC)) %>%
  group_by(round, cells, dataset) %>%
  slice_max(order_by = AUC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### scDr ####
scdr_res_df <- data.frame()
for (cells in cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "scDr")
    res_file <- file.path(res_folder_path, "drmref_scDr_auc.csv")
    if (!file.exists(res_file)) next
    scdr_res <- read.csv(res_file, header = TRUE)
    scdr_res$round <- round
    scdr_res$cells <- cells
    scdr_res_df <- rbind(scdr_res_df, scdr_res)
  }
}
colnames(scdr_res_df) <- tolower(colnames(scdr_res_df))
if (!("recall" %in% colnames(scdr_res_df))) scdr_res_df$recall <- infer_recall(scdr_res_df$precision, scdr_res_df$f1)

#### scIDUC ####
sciduc_res_df <- data.frame()
sciduc_cell_list <- c(500,1000, 2000,3000, 5000)
for (cells in sciduc_cell_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_", cells, "cells_res"), "sciduc")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    pattern <- "^(.*?)_sciduc_res\\((\\w+)\\s*\\([^)]*\\)\\).*\\.csv$"
    for (csv_file in csv_files) {
      csv_name <- basename(csv_file)
      m <- str_match(csv_name, pattern)
      if (!is.na(m[1,1])) {
        data_id <- m[1,2]
      } else {
        m2 <- str_match(csv_name, "^(.*?)_sciduc_res\\((\\w+)\\).*\\.csv$")
        if(is.na(m2[1,1])) next
        data_id <- m2[1,2]
      }
      pred_df <- read.csv(csv_file, header = TRUE, row.names = 1, check.names = FALSE)
      rownames(pred_df) <- gsub("\\.","-",pred_df$Cell_Name)
      tar_drug <- names(pred_df)[2]
      sc_data <- file.path(root_path, paste0(round, "_round_raw_", cells, "cells"), paste0("balanced_sampled_", data_id, "_seurat_afterAnno.RData"))
      meta_data <- load_meta(sc_data)
      rownames(meta_data) <- gsub("\\.","-",rownames(meta_data))
      common_cells <- intersect(rownames(pred_df), rownames(meta_data))
      pred_df <- pred_df[common_cells,,drop = FALSE]
      pred_df$condition <- meta_data[common_cells,"condition"]
      score <- pred_df[[tar_drug]]
      pred_df$pred_label <- ifelse(score > median(score, na.rm = TRUE), "resistant", "sensitive")
      true_label <- factor(pred_df$condition, levels = c("resistant","sensitive"))
      pred_label <- factor(pred_df$pred_label, levels = c("resistant","sensitive"))
      auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("resistant","sensitive"), quiet = TRUE)))
      cm <- caret::confusionMatrix(pred_label, true_label, positive = "sensitive")
      sciduc_res_df <- rbind(sciduc_res_df, data.frame(round = round, cells = cells, dataset = data_id, AUC = auc_value, ACC = as.numeric(cm$overall["Accuracy"]), Precision = as.numeric(cm$byClass["Precision"]), Recall = as.numeric(cm$byClass["Recall"]), F1 = as.numeric(cm$byClass["F1"])))
    }
  }
}
sciduc_best_df <- sciduc_res_df %>%
  group_by(round, cells, dataset) %>%
  slice_max(order_by = AUC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### Merge ####
all_res_list <- list(beyondcell_best_df, cadrres_sc_best_df, dreep_best_df, drugformer_res_df, precily_sc_res_df, scad_best_df, scdeal_best_df, scdr_res_df, sciduc_best_df)
names(all_res_list) <- c("Beyondcell","CaDRReS-sc","DREEP","DrugFormer","Precily","SCAD","scDEAL","scDr","scIDUC")

target_cols <- c("dataset","round","cells","auc","acc","f1","precision","recall")
merged_df <- lapply(names(all_res_list), function(tool_name) {
  df <- all_res_list[[tool_name]]
  if(is.null(df) || nrow(df) == 0) return(NULL)
  colnames(df) <- tolower(colnames(df))
  keep_cols <- intersect(target_cols, colnames(df))
  df_sub <- df[,keep_cols,drop = FALSE]
  for(col in setdiff(target_cols, keep_cols)) df_sub[[col]] <- NA
  df_sub <- df_sub[,target_cols]
  df_sub$Tool <- tool_name
  df_sub
}) %>% bind_rows()

colnames(merged_df) <- c("Dataset","round","cells","AUC","ACC","F1","Precision","Recall","Tool")
merged_df <- merged_df[,c("Dataset","AUC","ACC","F1","Precision","Recall","Tool","round","cells")]

save(merged_df, beyondcell_best_df, cadrres_sc_best_df, dreep_best_df, drugformer_res_df, precily_sc_res_df, scad_best_df, scdeal_best_df, scdr_res_df, sciduc_best_df, file = file.path(root_path, "balanced_res_process.RData"))
write.csv(merged_df, file.path(root_path, "balanced_final_df.csv"), row.names = FALSE)