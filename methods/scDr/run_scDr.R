suppressPackageStartupMessages({
  library(Seurat)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4) {
  stop("Usage: Rscript run_scdr.R <input_file> <output_dir> <up_signature> <down_signature>")
}

input_file <- args[1]
output_dir <- args[2]
up_signature_file <- args[3]
down_signature_file <- args[4]

if(!file.exists(input_file)) stop("Input file does not exist: ", input_file)
if(!file.exists(up_signature_file)) stop("Up-regulated drug signature file does not exist: ", up_signature_file)
if(!file.exists(down_signature_file)) stop("Down-regulated drug signature file does not exist: ", down_signature_file)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

up_drug_signature <- readRDS(up_signature_file)
down_drug_signature <- readRDS(down_signature_file)

obj <- load_seurat_object(input_file)
sample_id <- tools::file_path_sans_ext(basename(input_file))

message("Input: ", input_file)
message("Cells: ", ncol(obj))
message("Running scDr...")

exp_counts <- get_counts(obj)

sceAnalyse_new <- CreateSeuratObject(
  counts = exp_counts,
  project = sample_id
)

sceAnalyse_new <- NormalizeData(
  sceAnalyse_new,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

exp_1 <- tryCatch(
  as.data.frame(GetAssayData(sceAnalyse_new, assay = "RNA", layer = "data")),
  error = function(e) as.data.frame(GetAssayData(sceAnalyse_new, assay = "RNA", slot = "data"))
)

score_last <- NULL

for(i in names(up_drug_signature)) {
  drug_signature_one_of <- c(
    up_drug_signature[[i]][1:200],
    down_drug_signature[[i]][1:200]
  )
  
  gene <- intersect(names(drug_signature_one_of), rownames(exp_1))
  exp <- as.matrix(exp_1[gene, , drop = FALSE])
  
  row_sd <- rowSds(exp, useNames = FALSE)
  valid_genes <- which(!is.na(row_sd) & row_sd != 0)
  
  exp <- exp[valid_genes, , drop = FALSE]
  drug_signature_one_of <- drug_signature_one_of[valid_genes]
  
  exp <- (exp - rowMeans(exp)) / row_sd[valid_genes]
  
  score <- drug_signature_one_of * exp
  score <- colSums(score) / length(valid_genes)
  
  score_last <- rbind(score_last, score)
}

rownames(score_last) <- names(up_drug_signature)

output_file <- file.path(
  output_dir,
  paste0(sample_id, ".scDr_score.tsv")
)

write.table(
  score_last,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = TRUE
)

message("Completed.")
message("Output: ", output_file)