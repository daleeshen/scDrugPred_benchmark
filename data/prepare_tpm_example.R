suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
example_dir <- if(length(args) >= 1) args[1] else "./example"
reference_file <- if(length(args) >= 2) args[2] else "./reference/gene_annotation_GRCh38_113.rds"
gtf_file <- if(length(args) >= 3) args[3] else "./reference/Homo_sapiens.GRCh38.113.gtf.gz"

if(!dir.exists(example_dir)) stop("Example directory does not exist: ", example_dir)

load_seurat_object <- function(path) {
  env <- new.env(parent = globalenv())
  loaded <- load(path, envir = env)
  seurat_names <- loaded[vapply(loaded, function(x) inherits(env[[x]], "Seurat"), logical(1))]
  if(length(seurat_names) == 0) stop("No Seurat object found in: ", path)
  env[[seurat_names[1]]]
}

get_counts <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

Count2TPM <- function(counts, effLen) {
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

build_reference <- function(gtf_file, reference_file) {
  if(!file.exists(gtf_file)) stop("GTF file does not exist: ", gtf_file)
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
  })
  
  message("Reference cache not found. Building annotation from local GTF...")
  gtf <- rtracklayer::import(gtf_file)
  
  gene_map <- as.data.frame(S4Vectors::mcols(gtf[gtf$type == "gene"])) %>%
    select(gene_id, gene_name) %>%
    filter(!is.na(gene_id), !is.na(gene_name), gene_name != "") %>%
    mutate(gene_id = sub("\\..*", "", gene_id)) %>%
    distinct(gene_id, gene_name) %>%
    rename(hgnc_symbol = gene_name)
  
  exons <- gtf[gtf$type == "exon"]
  exons$gene_id <- sub("\\..*", "", exons$gene_id)
  exons <- exons[!is.na(exons$gene_id)]
  exons_by_gene <- split(exons, exons$gene_id)
  
  gfe <- data.frame(
    gene_id = names(exons_by_gene),
    length = vapply(exons_by_gene, function(x) sum(GenomicRanges::width(GenomicRanges::reduce(x))), numeric(1)),
    stringsAsFactors = FALSE
  )
  
  dir.create(dirname(reference_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(gene_map = gene_map, gfe = gfe), reference_file)
  message("Reference cache saved to: ", reference_file)
  
  rm(gtf, exons, exons_by_gene)
  gc()
  list(gene_map = gene_map, gfe = gfe)
}

if(file.exists(reference_file)) {
  message("Loading cached gene annotation...")
  ref <- readRDS(reference_file)
} else {
  ref <- build_reference(gtf_file, reference_file)
}

gene_map <- ref$gene_map
gfe <- ref$gfe
rm(ref)
gc()

rdata_files <- list.files(example_dir, pattern = "\\.RData$", full.names = TRUE, recursive = TRUE)
if(length(rdata_files) == 0) stop("No RData files found in: ", example_dir)

for(rdata_file in rdata_files) {
  message("Processing: ", rdata_file)
  obj <- load_seurat_object(rdata_file)
  expr_matrix <- as.data.frame(get_counts(obj))
  
  expr_matrix$hgnc_symbol <- rownames(expr_matrix)
  expr_matrix <- inner_join(expr_matrix, gene_map, by = "hgnc_symbol") %>%
    group_by(gene_id) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
  
  expr_matrix <- as.data.frame(expr_matrix)
  rownames(expr_matrix) <- expr_matrix$gene_id
  expr_matrix <- expr_matrix[, -1, drop = FALSE]
  expr_matrix <- expr_matrix[rownames(expr_matrix) %in% gfe$gene_id, , drop = FALSE]
  
  gene_length_idx <- match(rownames(expr_matrix), gfe$gene_id)
  effLen <- gfe$length[gene_length_idx]
  
  expr_tpm <- apply(expr_matrix, 2, Count2TPM, effLen = effLen)
  rownames(expr_tpm) <- rownames(expr_matrix)
  
  expr_tpm2 <- as.data.frame(expr_tpm)
  expr_tpm2$gene_id <- rownames(expr_tpm2)
  expr_tpm2 <- inner_join(expr_tpm2, gene_map, by = "gene_id") %>%
    filter(hgnc_symbol != "") %>%
    group_by(hgnc_symbol) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  expr_tpm2 <- as.data.frame(expr_tpm2)
  rownames(expr_tpm2) <- expr_tpm2$hgnc_symbol
  expr_tpm2 <- expr_tpm2[, -1, drop = FALSE]
  
  meta_data <- obj@meta.data[colnames(expr_tpm2), , drop = FALSE]
  
  file_base <- sub("\\.RData$", "", basename(rdata_file))
  output_dir <- file.path(dirname(rdata_file), "TPM_metadata")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  write.csv(
    expr_tpm2,
    file.path(output_dir, paste0(file_base, "_tpm.csv")),
    row.names = TRUE
  )
  
  write.csv(
    meta_data,
    file.path(output_dir, paste0(file_base, "_meta_data.csv")),
    row.names = TRUE
  )
  
  message("Saved TPM and metadata to: ", output_dir)
  
  rm(obj, expr_matrix, expr_tpm, expr_tpm2, meta_data)
  gc()
}

message("CaDRReS-sc example input preparation completed.")