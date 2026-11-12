suppressPackageStartupMessages({
  library(Seurat)
  library(beyondcell)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4) {
  stop("Usage: Rscript run_beyondcell.R <input_file> <drug> <output_dir> <PSc|SSc>")
}

input_file <- args[1]
drug <- toupper(args[2])
output_dir <- args[3]
collection_name <- toupper(args[4])

if(!file.exists(input_file)) stop("Input file does not exist: ", input_file)
if(!collection_name %in% c("PSC", "SSC")) stop("Collection must be PSc or SSc.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_input <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if(ext == "rds") {
    obj <- readRDS(path)
    if(!inherits(obj, "Seurat")) stop("RDS file does not contain a Seurat object.")
    return(obj)
  }
  
  if(ext == "rdata") {
    env <- new.env(parent = globalenv())
    loaded <- load(path, envir = env)
    seurat_names <- loaded[
      vapply(loaded, function(x) inherits(env[[x]], "Seurat"), logical(1))
    ]
    
    if(length(seurat_names) == 0) {
      stop("No Seurat object found in RData file.")
    }
    
    if(length(seurat_names) > 1) {
      warning("Multiple Seurat objects found; using ", seurat_names[1])
    }
    
    return(env[[seurat_names[1]]])
  }
  
  stop("Supported input formats: .RData and .rds")
}

get_counts <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

sample_id <- tools::file_path_sans_ext(basename(input_file))

message("Input: ", input_file)
message("Drug: ", drug)
message("Collection: ", collection_name)

obj <- load_input(input_file)
count_mat <- get_counts(obj)

sc <- CreateSeuratObject(
  counts = count_mat,
  project = sample_id
)

original_meta <- obj@meta.data[colnames(sc), , drop = FALSE]
for(col in colnames(original_meta)) {
  sc@meta.data[[col]] <- original_meta[[col]]
}

sc <- NormalizeData(
  sc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

sc <- ScaleData(sc, verbose = FALSE)
DefaultAssay(sc) <- "RNA"

if(collection_name == "PSC") {
  gs <- GetCollection(PSc)
} else {
  gs <- GetCollection(SSc)
}

if(drug == "TAE684") {
  drug <- "NVP-TAE684"
}

message("Cells: ", ncol(sc))
message("Running bcScore...")

bc <- bcScore(
  sc,
  gs,
  expr.thres = 0.1
)

bc.filtered <- bcSubset(
  bc,
  nan.cells = 0.95
)

bc.filtered@normalized[
  is.na(bc.filtered@normalized)
] <- 0

bc <- bcRecompute(
  bc.filtered,
  slot = "normalized"
)

message("Cells retained after BCS calculation: ", ncol(bc@normalized))
message("Regressing out nFeature_RNA...")

bc <- bcRegressOut(
  bc,
  vars.to.regress = "nFeature_RNA"
)

preferred_names <- toupper(
  as.character(drugInfo$IDs$preferred.drug.names)
)

drug_ids <- as.character(
  drugInfo$IDs[preferred_names == drug, 1]
)

drug_ids <- unique(
  drug_ids[!is.na(drug_ids)]
)

drug_ids <- drug_ids[
  drug_ids %in% rownames(bc@normalized)
]

if(length(drug_ids) == 0) {
  stop("No matching Beyondcell signature found for drug: ", drug)
}

score_mat <- t(
  as.matrix(
    bc@normalized[
      drug_ids,
      ,
      drop = FALSE
    ]
  )
)

score_df <- data.frame(
  Cell = rownames(score_mat),
  score_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

threshold_df <- data.frame(
  Signature = drug_ids,
  Switch_point = as.numeric(
    bc@switch.point[drug_ids]
  ),
  stringsAsFactors = FALSE
)

score_file <- file.path(
  output_dir,
  paste0(sample_id, "_Beyondcell_scores.csv")
)

threshold_file <- file.path(
  output_dir,
  paste0(sample_id, "_Beyondcell_switch_points.csv")
)

write.csv(
  score_df,
  score_file,
  row.names = FALSE
)

write.csv(
  threshold_df,
  threshold_file,
  row.names = FALSE
)

message("Completed.")
message("Matched signatures: ", paste(drug_ids, collapse = ", "))
message("Scores: ", score_file)
message("Switch points: ", threshold_file)