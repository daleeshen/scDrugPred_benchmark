library(Seurat)
library(dplyr)
library(PRROC)
library(pROC)
library(caret)
library(stringr)

root_path <- "/set/imblanced/results/root/path/"
raw_data_path <- "/path/to/raw/data/"
dataset_drug_info <- read.csv("~/scDrugBenchmark/evaluation/dataset_drug_single.csv", header = TRUE)
ratio_list <- c("1_3","1_10","1_30","1_100")
round_list <- c(1,2,3)

ratio_pred_label <- function(score, true_label, ratio, high_label, low_label) {
  ratio_n <- as.numeric(sub("^1_", "", ratio))
  high_majority <- sum(true_label == high_label, na.rm = TRUE) >= sum(true_label == low_label, na.rm = TRUE)
  high_prop <- ifelse(high_majority, ratio_n/(ratio_n + 1), 1/(ratio_n + 1))
  cutoff <- as.numeric(quantile(score, probs = 1 - high_prop, na.rm = TRUE, names = FALSE))
  ifelse(score > cutoff, high_label, low_label)
}

mean_na <- function(x) if(all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

#### Beyondcell ####
beyondcell_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_unbalanced_results"), ratio, "beyondcell")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    for (csv_file in csv_files) {
      data_id <- sub("_beyondcell_res\\.csv$", "", basename(csv_file))
      pred_df <- read.csv(csv_file, header = TRUE)
      sc_data <- file.path(raw_data_path, paste0(round, "_round_unbalanced_sampling"), ratio, paste0(data_id, ".RData"))
      loaded_objects <- load(sc_data)
      meta_data <- get(loaded_objects[1])@meta.data
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
        fg <- score[true_label == "sensitive"]
        bg <- score[true_label == "resistant"]
        AUPRC_value <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = FALSE)$auc.integral
        base_rate <- mean(true_label == "sensitive")
        norm_AUPRC <- (AUPRC_value - base_rate)/(1 - base_rate)
        if (is.nan(norm_AUPRC)) norm_AUPRC <- NA
        auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("resistant","sensitive"), quiet = TRUE)))
        cm <- caret::confusionMatrix(pred_label, true_label, positive = "sensitive")
        beyondcell_res_df <- rbind(beyondcell_res_df, data.frame(round = round, ratio = ratio, dataset = data_id, signature = sig_col, AUROC = auc_value, AUPRC = AUPRC_value, NormAUPRC = norm_AUPRC, ACC = as.numeric(cm$overall["Accuracy"]), Recall = as.numeric(cm$byClass["Recall"]), Precision = as.numeric(cm$byClass["Precision"]), F1 = as.numeric(cm$byClass["F1"])))
      }
    }
  }
}
beyondcell_best_df <- beyondcell_res_df %>%
  filter(!is.na(NormAUPRC)) %>%
  group_by(round, ratio, dataset) %>%
  slice_max(order_by = NormAUPRC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### CaDRReS-sc ####
cadrres_sc_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_unbalanced_results"), ratio, "CaDRres-sc")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    for (csv_file in csv_files) {
      data_id <- sub("_cadrres.*\\.csv$", "", basename(csv_file))
      pred_df <- read.csv(csv_file, header = TRUE, check.names = FALSE)
      sc_data <- file.path(raw_data_path, paste0(round, "_round_unbalanced_sampling"), ratio, paste0(data_id, ".RData"))
      loaded_objects <- load(sc_data)
      meta_data <- get(loaded_objects[1])@meta.data
      rownames(meta_data) <- gsub("\\.","-",rownames(meta_data))
      rownames(pred_df) <- gsub("\\.","-",pred_df[,1])
      common_cells <- intersect(rownames(pred_df), rownames(meta_data))
      pred_df <- pred_df[common_cells,,drop = FALSE]
      tar_drug_id <- dataset_drug_info[dataset_drug_info$Dataset == data_id,5]
      if (length(tar_drug_id) == 0 || is.na(tar_drug_id)) next
      tar_drug_id <- as.character(tar_drug_id)
      if (!(tar_drug_id %in% colnames(pred_df))) next
      pred_df <- pred_df[,c(1,which(colnames(pred_df) == tar_drug_id)),drop = FALSE]
      names(pred_df)[1] <- "cell_id"
      pred_df$condition <- meta_data[common_cells,"condition"]
      score <- pred_df[[tar_drug_id]]
      pred_df$pred_label <- ratio_pred_label(score, pred_df$condition, ratio, "resistant", "sensitive")
      true_label <- factor(pred_df$condition, levels = c("resistant","sensitive"))
      pred_label <- factor(pred_df$pred_label, levels = c("resistant","sensitive"))
      fg <- score[true_label == "sensitive"]
      bg <- score[true_label == "resistant"]
      AUPRC_value <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = FALSE)$auc.integral
      base_rate <- mean(true_label == "sensitive")
      norm_AUPRC <- (AUPRC_value - base_rate)/(1 - base_rate)
      if (is.nan(norm_AUPRC)) norm_AUPRC <- NA
      auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("resistant","sensitive"), quiet = TRUE)))
      cm <- caret::confusionMatrix(pred_label, true_label, positive = "sensitive")
      cadrres_sc_res_df <- rbind(cadrres_sc_res_df, data.frame(round = round, ratio = ratio, dataset = data_id, AUROC = auc_value, AUPRC = AUPRC_value, NormAUPRC = norm_AUPRC, ACC = as.numeric(cm$overall["Accuracy"]), Recall = as.numeric(cm$byClass["Recall"]), Precision = as.numeric(cm$byClass["Precision"]), F1 = as.numeric(cm$byClass["F1"])))
    }
  }
}
cadrres_sc_best_df <- cadrres_sc_res_df %>%
  group_by(round, ratio, dataset) %>%
  slice_max(order_by = NormAUPRC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### DREEP ####
dreep_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_unbalanced_results"), ratio, "DREEP")
    pred_res_list <- list.files(res_folder_path, pattern = "_pred_res\\.csv$", full.names = TRUE)
    for (pred_res in pred_res_list) {
      current_dataID <- sub("_pred_res\\.csv$", "", basename(pred_res))
      pred_df <- read.csv(pred_res, header = TRUE, check.names = FALSE)
      score_cols <- colnames(pred_df)[-(1:2)]
      true_label <- ifelse(pred_df[[2]] == "resistant", 1, 0)
      for (sc in score_cols) {
        scores <- pred_df[[sc]]
        valid_idx <- !is.na(scores)
        scores_clean <- scores[valid_idx]
        true_label_clean <- true_label[valid_idx]
        if (length(unique(true_label_clean)) < 2) {
          AUROC <- AUPRC <- NormAUPRC <- ACC <- Precision <- Recall <- F1 <- NA
        } else {
          AUROC <- as.numeric(pROC::roc(true_label_clean, scores_clean, quiet = TRUE)$auc)
          scores_pos <- scores_clean[true_label_clean == 1]
          scores_neg <- scores_clean[true_label_clean == 0]
          AUPRC <- if(length(scores_pos) == 0 || length(scores_neg) == 0) NA else PRROC::pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg, curve = FALSE)$auc.integral
          pos_rate <- mean(true_label_clean == 1)
          NormAUPRC <- ifelse(is.na(AUPRC), NA, (AUPRC - pos_rate)/(1 - pos_rate))
          pred_binary <- ifelse(scores_clean > 0, 1, 0)
          cm <- caret::confusionMatrix(factor(pred_binary, levels = c(0,1)), factor(true_label_clean, levels = c(0,1)), positive = "1")
          ACC <- as.numeric(cm$overall["Accuracy"])
          Precision <- as.numeric(cm$byClass["Precision"])
          Recall <- as.numeric(cm$byClass["Recall"])
          F1 <- as.numeric(cm$byClass["F1"])
        }
        dreep_res_df <- rbind(dreep_res_df, data.frame(Dataset = current_dataID, ratio = ratio, round = round, AUROC = AUROC, AUPRC = AUPRC, NormAUPRC = NormAUPRC, ACC = ACC, F1 = F1, Precision = Precision, Recall = Recall, Signature = sc))
      }
    }
  }
}
dreep_res_df <- dreep_res_df %>% filter(!grepl("\\(", Signature), !grepl("JQ12", Signature))
dreep_best_df <- dreep_res_df %>%
  group_by(round, ratio, Dataset) %>%
  slice_max(order_by = NormAUPRC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### DrugFormer ####
drugformer_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path("/media/shenli/Expansion/scDrug_benchmark/DrugFormer", paste0(round, "_cell_dt_", ratio, "_logs"))
    log_files <- list.files(res_folder_path, pattern = "\\.txt$", full.names = TRUE)
    for (log_file in log_files) {
      data_id <- sub("_log\\.txt$", "", basename(log_file))
      test_lines <- grep("Testing", readLines(log_file), value = TRUE)
      get_metric <- function(pattern) as.numeric(sub(paste0(".*",pattern,": ([0-9\\.eE+-]+).*"), "\\1", test_lines))
      ap <- get_metric("APscore")
      normAP <- get_metric("Normalized_AP")
      acc <- get_metric("Accuracy")
      auc <- get_metric("AUROC")
      prec <- get_metric("Precision")
      f1 <- get_metric("F1")
      recall <- get_metric("Recall")
      avg <- colMeans(data.frame(ap,normAP,acc,auc,prec,f1,recall), na.rm = TRUE)
      drugformer_res_df <- rbind(drugformer_res_df, data.frame(round = round, ratio = ratio, dataset = data_id, AUPRC = avg[1], NormAUPRC = avg[2], AUROC = avg[4], ACC = avg[3], Precision = avg[5], F1 = avg[6], Recall = avg[7]))
    }
  }
}

#### Precily ####
precily_sc_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_unbalanced_results"), ratio, "precily")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    for (csv_file in csv_files) {
      data_id <- sub("_predicted_results\\.csv$", "", basename(csv_file))
      pred_df <- read.csv(csv_file, header = TRUE, check.names = FALSE)
      rownames(pred_df) <- make.unique(pred_df[,1])
      pred_df <- as.data.frame(t(pred_df[,-1]))
      pred_df$cell_id <- rownames(pred_df)
      rownames(pred_df) <- gsub("\\.","-",rownames(pred_df))
      tar_drug_name <- dataset_drug_info[dataset_drug_info$Dataset == data_id,2]
      if (length(tar_drug_name) == 0 || is.na(tar_drug_name)) next
      sc_data <- file.path(raw_data_path, paste0(round, "_round_unbalanced_sampling"), ratio, paste0(data_id, ".RData"))
      loaded_objects <- load(sc_data)
      meta_data <- get(loaded_objects[1])@meta.data
      rownames(meta_data) <- gsub("\\.","-",rownames(meta_data))
      common_cells <- intersect(rownames(pred_df), rownames(meta_data))
      pred_df <- pred_df[common_cells,,drop = FALSE]
      if (!(tar_drug_name %in% colnames(pred_df))) next
      pred_df <- pred_df[,c("cell_id",tar_drug_name)]
      pred_df$condition <- meta_data[common_cells,"condition"]
      raw_score <- pred_df[[tar_drug_name]]
      score <- -raw_score
      pred_df$pred_label <- ratio_pred_label(raw_score, pred_df$condition, ratio, "resistant", "sensitive")
      true_label <- factor(pred_df$condition, levels = c("resistant","sensitive"))
      pred_label <- factor(pred_df$pred_label, levels = c("resistant","sensitive"))
      fg <- score[true_label == "sensitive"]
      bg <- score[true_label == "resistant"]
      AUPRC_value <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg, curve = FALSE)$auc.integral
      base_rate <- mean(true_label == "sensitive")
      norm_AUPRC <- (AUPRC_value - base_rate)/(1 - base_rate)
      if (is.nan(norm_AUPRC)) norm_AUPRC <- NA
      auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("resistant","sensitive"), quiet = TRUE)))
      cm <- caret::confusionMatrix(pred_label, true_label, positive = "sensitive")
      precily_sc_res_df <- rbind(precily_sc_res_df, data.frame(round = round, ratio = ratio, dataset = data_id, AUROC = auc_value, AUPRC = AUPRC_value, NormAUPRC = norm_AUPRC, ACC = as.numeric(cm$overall["Accuracy"]), Recall = as.numeric(cm$byClass["Recall"]), Precision = as.numeric(cm$byClass["Precision"]), F1 = as.numeric(cm$byClass["F1"])))
    }
  }
}

#### SCAD ####
scad_res_df <- data.frame()
mode_list <- c("norm","ppi","tp4k")
for (ratio in ratio_list) {
  for (round in round_list) {
    for (mode in mode_list) {
      scad_base_dir <- file.path(root_path, paste0(round, "_round_unbalanced_results"), ratio, "SCAD", "best_param_res_weight", paste0("best_param_res_weight_", mode))
      if (!dir.exists(scad_base_dir)) next
      scad_id_folders <- list.dirs(scad_base_dir, recursive = FALSE)
      for (id_dir in scad_id_folders) {
        pred_files <- list.files(id_dir, pattern = "_predicted_y\\.tsv$", recursive = TRUE, full.names = TRUE)
        if (length(pred_files) == 0) next
        id_name <- basename(id_dir)
        id_split <- strsplit(id_name, "_")[[1]]
        dataset_id <- id_split[1]
        drug <- ifelse(length(id_split) >= 2, id_split[2], NA)
        for (pf in pred_files) {
          scad_pred_df <- read.table(pf, header = TRUE)
          names(scad_pred_df) <- c("splits","cell_id","tr_label","predicted_val")
          truth <- scad_pred_df$tr_label
          pred_prob <- scad_pred_df$predicted_val
          valid_idx <- !is.na(pred_prob) & !is.na(truth)
          truth <- truth[valid_idx]
          pred_prob <- pred_prob[valid_idx]
          if (length(unique(truth)) < 2) next
          auroc <- as.numeric(pROC::roc(truth, pred_prob, quiet = TRUE)$auc)
          auprc <- PRROC::pr.curve(scores.class0 = pred_prob[truth == 1], scores.class1 = pred_prob[truth == 0], curve = FALSE)$auc.integral
          prevalence <- mean(truth == 1)
          norm_auprc <- (auprc - prevalence)/(1 - prevalence)
          pred_label <- ratio_pred_label(pred_prob, truth, ratio, 1, 0)
          cm <- caret::confusionMatrix(factor(pred_label, levels = c(0,1)), factor(truth, levels = c(0,1)), positive = "1")
          scad_res_df <- rbind(scad_res_df, data.frame(dataset = dataset_id, drug = drug, ratio = ratio, round = round, mode = mode, AUROC = auroc, AUPRC = auprc, NormAUPRC = norm_auprc, ACC = as.numeric(cm$overall["Accuracy"]), Precision = as.numeric(cm$byClass["Precision"]), Recall = as.numeric(cm$byClass["Recall"]), F1 = as.numeric(cm$byClass["F1"])))
        }
      }
    }
  }
}
scad_res_df <- scad_res_df %>%
  mutate(dataset = case_when(
    dataset == "GSE131984" & drug == "Pal" ~ "GSE131984_Pal_DMSO_rm",
    dataset == "GSE131984" & drug == "Pac" ~ "GSE131984_Pac_DMSO_rm",
    dataset == "GSE131984" & drug == "JQ1" ~ "GSE131984_JQ1_DMSO_rm",
    dataset == "GSE169246" & drug == "PacBlood" ~ "GSE169246_PacBlood",
    dataset == "GSE169246" & drug == "PacTissue" ~ "GSE169246_PacTissue",
    TRUE ~ dataset
  ))

scad_res_df_avg <- scad_res_df %>%
  group_by(dataset, ratio, round, mode) %>%
  summarise(AUROC = mean_na(AUROC), AUPRC = mean_na(AUPRC), NormAUPRC = mean_na(NormAUPRC), ACC = mean_na(ACC), Precision = mean_na(Precision), Recall = mean_na(Recall), F1 = mean_na(F1), .groups = "drop")

scad_best_df <- scad_res_df_avg %>%
  filter(!is.na(NormAUPRC)) %>%
  group_by(round, ratio, dataset) %>%
  slice_max(order_by = NormAUPRC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### scDEAL ####
scdeal_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path("/media/shenli/Expansion/scDrug_benchmark/scDEAL", paste0(round, "_", ratio, "_ratio_logs"))
    log_files <- list.files(res_folder_path, pattern = "\\.log$", full.names = TRUE)
    for (log_file in log_files) {
      log_name <- basename(log_file)
      dataset_id <- str_match(log_name, "^(.+?)_[^_]+_[^_]+_bottleneck\\d+_sampling[^_]+\\.log$")[,2]
      bottleneck <- str_match(log_name, "_bottleneck(\\d+)_")[,2]
      sampling <- str_match(log_name, "_sampling([^_\\.]+)")[,2]
      txt <- readChar(log_file, file.info(log_file)$size)
      auc <- str_match(txt, "auroc_score\\s+([0-9\\.eE+-]+)")[,2]
      ap <- str_match(txt, "ap_score\\s+([0-9\\.eE+-]+)")[,2]
      acc <- str_match(txt, "Accuracy:\\s*([0-9\\.eE+-]+)")[,2]
      recall <- str_match(txt, "recall\\s*([0-9\\.eE+-]+)")[,2]
      precision <- str_match(txt, "precision\\s*([0-9\\.eE+-]+)")[,2]
      f1 <- str_match(txt, "f1score\\s*([0-9\\.eE+-]+)")[,2]
      norm_ap <- str_match(txt, "Normalized\\s*AUPRC\\s*[:=]?\\s*([0-9eE+\\.-]+)")[,2]
      scdeal_res_df <- rbind(scdeal_res_df, data.frame(round = round, ratio = ratio, bottleneck = bottleneck, sampling_method = sampling, dataset = dataset_id, AUROC = as.numeric(auc), AUPRC = as.numeric(ap), NormAUPRC = as.numeric(norm_ap), ACC = as.numeric(acc), Recall = as.numeric(recall), Precision = as.numeric(precision), F1 = as.numeric(f1)))
    }
  }
}
scdeal_best_df <- scdeal_res_df %>%
  filter(!is.na(NormAUPRC)) %>%
  group_by(round, ratio, dataset) %>%
  slice_max(order_by = NormAUPRC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### scDr ####
scdr_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_unbalanced_results"), ratio, "scDr")
    res_file <- file.path(res_folder_path, "drmref_scDr_performance_unbalanced.csv")
    if (!file.exists(res_file)) next
    scdr_res <- read.csv(res_file, header = TRUE)
    scdr_res$ratio <- ratio
    scdr_res$round <- round
    scdr_res_df <- rbind(scdr_res_df, scdr_res)
  }
}

#### scIDUC ####
sciduc_res_df <- data.frame()
for (ratio in ratio_list) {
  for (round in round_list) {
    res_folder_path <- file.path(root_path, paste0(round, "_round_unbalanced_results"), ratio, "sciduc")
    csv_files <- list.files(res_folder_path, pattern = "\\.csv$", full.names = TRUE)
    pattern <- "^(.+?)_sciduc_res\\((.+)\\)_alpha[0-9.]+_ratio[0-9.]+\\.csv$"
    for (csv_file in csv_files) {
      csv_name <- basename(csv_file)
      m <- str_match(csv_name, pattern)
      if (!is.na(m[1,1])) {
        data_id <- m[1,2]
      } else {
        m2 <- str_match(csv_name, "^(.*?)_sciduc_res\\((\\w+)\\).*\\.csv$")
        if (is.na(m2[1,1])) next
        data_id <- m2[1,2]
      }
      pred_df <- read.csv(csv_file, header = TRUE, row.names = 1, check.names = FALSE)
      rownames(pred_df) <- gsub("\\.","-",pred_df$Cell_Name)
      tar_drug <- names(pred_df)[2]
      sc_data <- file.path(raw_data_path, paste0(round, "_round_unbalanced_sampling"), ratio, paste0(data_id, ".RData"))
      loaded_objects <- load(sc_data)
      meta_data <- get(loaded_objects[1])@meta.data
      rownames(meta_data) <- gsub("\\.","-",rownames(meta_data))
      common_cells <- intersect(rownames(pred_df), rownames(meta_data))
      pred_df <- pred_df[common_cells,,drop = FALSE]
      pred_df$condition <- meta_data[common_cells,"condition"]
      score <- pred_df[[tar_drug]]
      pred_df$pred_label <- ratio_pred_label(score, pred_df$condition, ratio, "resistant", "sensitive")
      true_label <- factor(pred_df$condition, levels = c("resistant","sensitive"))
      pred_label <- factor(pred_df$pred_label, levels = c("resistant","sensitive"))
      auc_value <- as.numeric(pROC::auc(pROC::roc(response = true_label, predictor = score, levels = c("sensitive","resistant"), quiet = TRUE)))
      pos_scores <- score[true_label == "resistant"]
      neg_scores <- score[true_label != "resistant"]
      auprc_value <- PRROC::pr.curve(scores.class0 = pos_scores, scores.class1 = neg_scores, curve = FALSE)$auc.integral
      baseline <- mean(true_label == "resistant")
      norm_auprc_value <- (auprc_value - baseline)/(1 - baseline)
      cm <- caret::confusionMatrix(pred_label, true_label, positive = "resistant")
      sciduc_res_df <- rbind(sciduc_res_df, data.frame(round = round, ratio = ratio, dataset = data_id, AUROC = auc_value, AUPRC = auprc_value, NormAUPRC = norm_auprc_value, ACC = as.numeric(cm$overall["Accuracy"]), Recall = as.numeric(cm$byClass["Recall"]), Precision = as.numeric(cm$byClass["Precision"]), F1 = as.numeric(cm$byClass["F1"])))
    }
  }
}
sciduc_best_df <- sciduc_res_df %>%
  group_by(round, ratio, dataset) %>%
  slice_max(order_by = NormAUPRC, n = 1, with_ties = FALSE) %>%
  ungroup()

#### Merge ####
all_res_list <- list(beyondcell_best_df, cadrres_sc_best_df, dreep_best_df, drugformer_res_df, precily_sc_res_df, scad_best_df, scdeal_best_df, scdr_res_df, sciduc_best_df)
names(all_res_list) <- c("Beyondcell","CaDRReS-sc","DREEP","DrugFormer","Precily","SCAD","scDEAL","scDr","scIDUC")

target_cols <- c("dataset","round","ratio","auroc","auprc","normauprc","recall","f1","precision","acc")
merged_df <- lapply(names(all_res_list), function(tool_name) {
  df <- all_res_list[[tool_name]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  colnames(df) <- tolower(colnames(df))
  keep_cols <- intersect(target_cols, colnames(df))
  df_sub <- df[,keep_cols,drop = FALSE]
  for (col in setdiff(target_cols, keep_cols)) df_sub[[col]] <- NA
  df_sub <- df_sub[,target_cols]
  df_sub$Tool <- tool_name
  df_sub
}) %>% bind_rows()

colnames(merged_df) <- c("Dataset","Round","Ratio","AUROC","AUPRC","NormAUPRC","Recall","F1","Precision","ACC","Tool")
merged_df$Ratio <- gsub("_",":",merged_df$Ratio)
merged_df[is.na(merged_df)] <- 0

save(merged_df, beyondcell_best_df, cadrres_sc_best_df, dreep_best_df, drugformer_res_df, precily_sc_res_df, scad_best_df, scdeal_best_df, scdr_res_df, sciduc_best_df, file = file.path(root_path, "unb_data_processing.RData"))
write.csv(merged_df, file.path(root_path, "unbalanced_final_df.csv"), row.names = FALSE)
