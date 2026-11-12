suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if(length(args) >= 1) args[1] else "./original_data"
output_root <- if(length(args) >= 2) args[2] else "./balanced_datasets"
if(!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

k_cell_list <- c(500, 1000, 2000, 3000, 5000)
seed_list <- c(123, 456, 789)
dmso_remove_list <- c("GSE131984_JQ1", "GSE131984_Pac", "GSE131984_Pal")
exclude_datasets <- character(0)

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
  obj@meta.data$ID <- rownames(obj@meta.data)
  condition_values <- as.character(obj@meta.data$condition)
  if(any(is.na(condition_values))) {
    keep_cells <- rownames(obj@meta.data)[!is.na(condition_values)]
    obj <- subset(obj, cells = keep_cells)
  }
  cell_counts <- as.data.frame(obj@meta.data) %>%
    count(condition, name = "n")
  if(nrow(cell_counts) != 2) {
    warning(dataset_id, " does not contain exactly two response conditions; skipping.")
    next
  }
  
  for(k_cell in k_cell_list) {
    target_n <- k_cell / 2
    if(any(cell_counts$n < target_n)) {
      message(dataset_id, ": insufficient cells for N=", k_cell, "; skipping.")
      next
    }
    min_cells_per_type <- min(cell_counts$n)
    
    for(seed_idx in seq_along(seed_list)) {
      seed <- seed_list[seed_idx]
      message(dataset_id, " | N=", k_cell, " | round=", seed_idx, " | seed=", seed)
      set.seed(seed)
      balanced_cells <- as.data.frame(obj@meta.data) %>%
        group_by(condition) %>%
        sample_n(size = min_cells_per_type) %>%
        pull(ID)
      balanced_seurat <- subset(obj, cells = balanced_cells)
      
      set.seed(seed)
      if(k_cell == 1000) {
        sampled_cells <- as.data.frame(balanced_seurat@meta.data) %>%
          group_by(condition) %>%
          sample_frac(size = 1000 / ncol(balanced_seurat)) %>%
          pull(ID)
      } else {
        sampled_cells <- as.data.frame(balanced_seurat@meta.data) %>%
          group_by(condition) %>%
          sample_n(size = target_n) %>%
          pull(ID)
      }
      
      balanced_sampled_datasets <- subset(balanced_seurat, cells = sampled_cells)
      final_counts <- table(balanced_sampled_datasets@meta.data$condition)
      if(ncol(balanced_sampled_datasets) != k_cell || length(final_counts) != 2 || any(final_counts != target_n)) {
        stop("Unexpected sampling result for ", dataset_id, ", N=", k_cell, ", round=", seed_idx)
      }
      
      save_dir <- file.path(output_root, paste0(seed_idx, "_round_raw_", k_cell, "cells"))
      export_dir <- file.path(save_dir, "metadata_count")
      dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)
      
      rdata_file <- file.path(save_dir, paste0("balanced_sampled_", dataset_id, "_seurat_afterAnno.RData"))
      save(balanced_sampled_datasets, file = rdata_file)
      
      raw_counts <- get_counts(balanced_sampled_datasets)
      write.csv(raw_counts, file = file.path(export_dir, paste0(dataset_id, "_raw_count_matrix.csv")))
      write.csv(balanced_sampled_datasets@meta.data, file = file.path(export_dir, paste0(dataset_id, "_meta_data.csv")))
      
      manifest <- bind_rows(
        manifest,
        data.frame(
          Dataset = dataset_id,
          Total_cells = k_cell,
          Round = seed_idx,
          Seed = seed,
          Condition_1 = names(final_counts)[1],
          Condition_1_cells = as.integer(final_counts[1]),
          Condition_2 = names(final_counts)[2],
          Condition_2_cells = as.integer(final_counts[2]),
          RData_file = rdata_file,
          stringsAsFactors = FALSE
        )
      )
      rm(balanced_cells, balanced_seurat, sampled_cells, balanced_sampled_datasets, raw_counts)
      gc()
    }
  }
  rm(obj)
  gc()
}

write.csv(manifest, file.path(output_root, "balanced_sampling_manifest.csv"), row.names = FALSE)
message("Balanced dataset generation completed.")