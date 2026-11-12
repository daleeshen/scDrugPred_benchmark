suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if(length(args) >= 1) args[1] else "./original_data"
output_root <- if(length(args) >= 2) args[2] else "./imbalanced_datasets"
if(!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

ratios <- c(3, 10, 30, 100)
max_total <- 10000
seed_list <- c(123, 456, 789)
dmso_remove_list <- c("GSE131984_JQ1", "GSE131984_Pac", "GSE131984_Pal")
exclude_datasets <- character(0)

normalize_condition <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("\\s+", "", x)
  x <- gsub("resisitant|resistence|resistantt", "resistant", x)
  x
}

get_dataset_id <- function(path) {
  sub("_seurat_afterAnno.*$", "", basename(path))
}

load_seurat_object <- function(path) {
  env <- new.env(parent = globalenv())
  loaded_objects <- load(path, envir = env)
  seurat_objects <- loaded_objects[vapply(loaded_objects, function(x) inherits(env[[x]], "Seurat"), logical(1))]
  if(length(seurat_objects) == 0) stop("No Seurat object found in: ", path)
  if(length(seurat_objects) > 1) warning("Multiple Seurat objects found in ", path, "; using ", seurat_objects[1])
  env[[seurat_objects[1]]]
}

get_counts <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

infer_classes <- function(obj, cond_col = "condition") {
  if(!cond_col %in% colnames(obj@meta.data)) stop("Missing condition column.")
  cond <- normalize_condition(obj@meta.data[[cond_col]])
  tab <- sort(table(cond), decreasing = FALSE)
  if(length(tab) != 2) stop("Condition must contain exactly two classes.")
  list(
    minority = names(tab)[1],
    majority = names(tab)[2],
    n_min = as.integer(tab[1]),
    n_maj = as.integer(tab[2])
  )
}

compute_counts <- function(n_min_available, n_maj_available, ratio, max_total = 10000) {
  n_min <- min(
    floor(max_total / (1 + ratio)),
    n_min_available,
    floor(n_maj_available / ratio)
  )
  if(n_min <= 0) return(list(n_min = 0L, n_maj = 0L))
  list(n_min = as.integer(n_min), n_maj = as.integer(ratio * n_min))
}

sample_imbalanced <- function(obj, cond_col, minority_label, majority_label, n_min, n_maj, seed) {
  cond <- normalize_condition(obj@meta.data[[cond_col]])
  min_cells <- Cells(obj)[which(cond == minority_label)]
  maj_cells <- Cells(obj)[which(cond == majority_label)]
  if(length(min_cells) < n_min || length(maj_cells) < n_maj) stop("Insufficient cells for sampling.")
  set.seed(seed)
  sampled_min <- sample(min_cells, n_min, replace = FALSE)
  sampled_maj <- sample(maj_cells, n_maj, replace = FALSE)
  subset(obj, cells = c(sampled_min, sampled_maj))
}

scdata_list <- list.files(input_dir, pattern = "_seurat_afterAnno.*\\.RData$", full.names = TRUE)
if(length(scdata_list) == 0) stop("No RData files found in: ", input_dir)

manifest <- data.frame()

for(scdata in scdata_list) {
  dataset_id <- get_dataset_id(scdata)
  if(dataset_id %in% exclude_datasets) next
  message("Loading ", dataset_id, "...")
  obj <- load_seurat_object(scdata)
  
  if(dataset_id %in% dmso_remove_list) {
    if(!"treatment" %in% colnames(obj@meta.data)) stop("Missing treatment column in ", dataset_id)
    keep_cells <- rownames(obj@meta.data)[tolower(trimws(as.character(obj@meta.data$treatment))) != "dmso"]
    obj <- subset(obj, cells = keep_cells)
  }
  
  if(!"condition" %in% colnames(obj@meta.data)) {
    warning("Missing condition column in ", dataset_id, "; skipping.")
    next
  }
  
  valid_cells <- rownames(obj@meta.data)[!is.na(obj@meta.data$condition)]
  obj <- subset(obj, cells = valid_cells)
  
  classes <- tryCatch(
    infer_classes(obj, "condition"),
    error = function(e) {
      warning(dataset_id, ": ", conditionMessage(e))
      NULL
    }
  )
  if(is.null(classes)) next
  
  for(ratio in ratios) {
    target <- compute_counts(classes$n_min, classes$n_maj, ratio, max_total)
    if(target$n_min == 0 || target$n_maj == 0) {
      warning(dataset_id, ": ratio 1:", ratio, " is not feasible; skipping.")
      next
    }
    
    for(round_idx in seq_along(seed_list)) {
      seed <- seed_list[round_idx]
      message(dataset_id, " | ratio=1:", ratio, " | round=", round_idx, " | seed=", seed)
      
      unbalanced_sampled_dataset <- sample_imbalanced(
        obj = obj,
        cond_col = "condition",
        minority_label = classes$minority,
        majority_label = classes$majority,
        n_min = target$n_min,
        n_maj = target$n_maj,
        seed = seed
      )
      
      cond_norm <- normalize_condition(unbalanced_sampled_dataset@meta.data$condition)
      final_min <- sum(cond_norm == classes$minority)
      final_maj <- sum(cond_norm == classes$majority)
      if(final_min != target$n_min || final_maj != target$n_maj) {
        stop("Unexpected sampling result for ", dataset_id, ", ratio 1:", ratio, ", round=", round_idx)
      }
      
      round_dir <- file.path(output_root, paste0(round_idx, "_round_unbalanced_sampling"))
      ratio_dir <- file.path(round_dir, paste0("1_", ratio))
      export_dir <- file.path(ratio_dir, "metadata_count")
      dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
      
      rdata_file <- file.path(ratio_dir, paste0(dataset_id, ".RData"))
      save(unbalanced_sampled_dataset, file = rdata_file)
      
      raw_counts <- get_counts(unbalanced_sampled_dataset)
      write.csv(raw_counts, file = file.path(export_dir, paste0(dataset_id, "_raw_count_matrix.csv")))
      write.csv(unbalanced_sampled_dataset@meta.data, file = file.path(export_dir, paste0(dataset_id, "_meta_data.csv")))
      
      manifest <- bind_rows(
        manifest,
        data.frame(
          Dataset = dataset_id,
          Ratio = paste0("1:", ratio),
          Round = round_idx,
          Seed = seed,
          Minority_label = classes$minority,
          Majority_label = classes$majority,
          Minority_cells = target$n_min,
          Majority_cells = target$n_maj,
          Total_cells = target$n_min + target$n_maj,
          Max_total = max_total,
          RData_file = rdata_file,
          stringsAsFactors = FALSE
        )
      )
      
      rm(unbalanced_sampled_dataset, raw_counts)
      gc()
    }
  }
  
  rm(obj)
  gc()
}

write.csv(manifest, file.path(output_root, "imbalanced_sampling_manifest.csv"), row.names = FALSE)
message("Imbalanced dataset generation completed.")