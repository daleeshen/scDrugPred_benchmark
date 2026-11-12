suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
example_dir <- if(length(args) >= 1) args[1] else "./example"

if(!dir.exists(example_dir)) stop("Example directory does not exist: ", example_dir)

load_seurat_object <- function(path) {
  if(grepl("\\.rds$", path, ignore.case = TRUE)) {
    obj <- readRDS(path)
    if(!inherits(obj, "Seurat")) stop("RDS file does not contain a Seurat object.")
    return(obj)
  }
  
  env <- new.env(parent = globalenv())
  loaded <- load(path, envir = env)
  seurat_names <- loaded[vapply(loaded, function(x) inherits(env[[x]], "Seurat"), logical(1))]
  
  if(length(seurat_names) == 0) stop("No Seurat object found in: ", path)
  if(length(seurat_names) > 1) warning("Multiple Seurat objects found; using ", seurat_names[1])
  
  env[[seurat_names[1]]]
}

get_counts <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

scenarios <- c("balanced_1000", "imbalanced_1_3", "imbalanced_1_30")

for(scenario in scenarios) {
  scenario_dir <- file.path(example_dir, scenario)
  
  if(!dir.exists(scenario_dir)) {
    warning("Scenario directory not found: ", scenario_dir)
    next
  }
  
  input_files <- list.files(
    scenario_dir,
    pattern = "\\.(RData|rds)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if(length(input_files) == 0) {
    warning("No Seurat files found in: ", scenario_dir)
    next
  }
  
  output_dir <- file.path(scenario_dir, "scdeal_input")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(input_file in input_files) {
    message("Processing: ", input_file)
    
    obj <- load_seurat_object(input_file)
    
    if(!"condition" %in% colnames(obj@meta.data)) {
      stop("Metadata column 'condition' not found in: ", input_file)
    }
    
    if(!any(c("percent.mt", "pct_counts_mt") %in% colnames(obj@meta.data))) {
      stop("Neither 'percent.mt' nor 'pct_counts_mt' found in metadata: ", input_file)
    }
    
    counts <- get_counts(obj)
    gene_names <- rownames(counts)
    cell_names <- colnames(counts)
    
    exp_counts <- as.data.frame(as.matrix(counts), check.names = FALSE)
    rownames(exp_counts) <- gene_names
    colnames(exp_counts) <- cell_names
    
    meta_data <- obj@meta.data[cell_names, , drop = FALSE]
    
    sample_id <- tools::file_path_sans_ext(basename(input_file))
    
    count_file <- file.path(
      output_dir,
      paste0(sample_id, "_raw_count_matrix.csv")
    )
    
    metadata_file <- file.path(
      output_dir,
      paste0(sample_id, "_meta_data.csv")
    )
    
    write.csv(
      exp_counts,
      count_file,
      row.names = TRUE,
      quote = FALSE
    )
    
    write.csv(
      meta_data,
      metadata_file,
      row.names = TRUE,
      quote = FALSE
    )
    
    message("Saved expression: ", count_file)
    message("Saved metadata: ", metadata_file)
    
    rm(obj, counts, exp_counts, meta_data)
    gc()
  }
}

message("scDEAL example input preparation completed.")