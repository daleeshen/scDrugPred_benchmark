suppressPackageStartupMessages({
  library(Seurat)
  library(GenomicFeatures)
  library(rtracklayer)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

example_dir <- if(length(args) >= 1) args[1] else "./example"
gtf_file <- if(length(args) >= 2) args[2] else "../methods/Precily/Precily/Homo_sapiens.GRCh38.110.gtf"
reference_file <- if(length(args) >= 3) args[3] else "./reference/precily_gene_lengths_GRCh38_110.rds"

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

calculate_tpm <- function(raw_counts, gene_lengths) {
  gene_lengths_kb <- gene_lengths / 1000
  rpk <- raw_counts / gene_lengths_kb
  tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6
  return(tpm)
}

build_gene_lengths <- function(gtf_file, reference_file) {
  if(!file.exists(gtf_file)) stop("GTF file does not exist: ", gtf_file)
  
  message("Building Precily gene-length reference from GRCh38.110 GTF...")
  
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  gene_ranges <- genes(txdb, columns = "gene_id")
  
  gtf_data <- rtracklayer::import(gtf_file)
  gene_symbols <- as.data.frame(gtf_data) %>%
    filter(type == "gene") %>%
    dplyr::select(gene_id, gene_name) %>%
    distinct()
  
  gene_lengths <- as.data.frame(gene_ranges)
  gene_lengths <- merge(
    gene_lengths,
    gene_symbols,
    by.x = "gene_id",
    by.y = "gene_id"
  )
  
  gene_lengths <- gene_lengths[, c("gene_name", "width")]
  colnames(gene_lengths) <- c("Gene", "Length")
  
  gene_lengths <- gene_lengths %>%
    group_by(Gene) %>%
    summarise(Length = max(Length), .groups = "drop")
  
  gene_lengths <- gene_lengths[!is.na(gene_lengths$Gene), ]
  gene_lengths <- as.data.frame(gene_lengths)
  rownames(gene_lengths) <- gene_lengths$Gene
  
  dir.create(dirname(reference_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(gene_lengths, reference_file)
  
  message("Gene-length reference saved to: ", reference_file)
  
  gene_lengths
}

if(file.exists(reference_file)) {
  message("Loading cached Precily gene-length reference...")
  gene_lengths <- readRDS(reference_file)
} else {
  gene_lengths <- build_gene_lengths(gtf_file, reference_file)
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
  
  output_dir <- file.path(scenario_dir, "precily_input")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(input_file in input_files) {
    message("Processing: ", input_file)
    
    obj <- load_seurat_object(input_file)
    counts <- get_counts(obj)
    
    common_genes <- intersect(rownames(counts), rownames(gene_lengths))
    
    if(length(common_genes) == 0) {
      stop("No genes overlap with the Precily GRCh38.110 annotation: ", input_file)
    }
    
    raw_counts <- as.matrix(counts[common_genes, , drop = FALSE])
    gene_lengths_filtered <- gene_lengths[common_genes, , drop = FALSE]
    
    tpm_matrix <- calculate_tpm(
      raw_counts,
      gene_lengths_filtered$Length
    )
    
    log_tpm_matrix <- log2(tpm_matrix + 1)
    
    sample_id <- tools::file_path_sans_ext(basename(input_file))
    
    logtpm_file <- file.path(
      output_dir,
      paste0(sample_id, "_logTPM.csv")
    )
    
    metadata_file <- file.path(
      output_dir,
      paste0(sample_id, "_meta_data.csv")
    )
    
    write.csv(
      log_tpm_matrix,
      logtpm_file,
      row.names = TRUE,
      quote = FALSE
    )
    
    write.csv(
      obj@meta.data[colnames(log_tpm_matrix), , drop = FALSE],
      metadata_file,
      row.names = TRUE,
      quote = FALSE
    )
    
    message("Saved logTPM: ", logtpm_file)
    message("Saved metadata: ", metadata_file)
    
    rm(obj, counts, raw_counts, tpm_matrix, log_tpm_matrix)
    gc()
  }
}

message("Precily example input preparation completed.")