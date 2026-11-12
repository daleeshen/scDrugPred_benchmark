suppressPackageStartupMessages({
  library(Seurat)
  library(gficf)
  library(DREEP)
  library(dplyr)
  library(progress)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  stop("Usage: Rscript run_dreep.R <input_file> <drug> <output_dir> <reference_rds>")
}

input_file <- args[1]
target_drug <- toupper(args[2])
output_dir <- args[3]
reference_file <- args[4]

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
  env[[seurat_names[1]]]
}

get_counts <- function(obj) {
  tryCatch(
    GetAssayData(obj, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(obj, assay = "RNA", slot = "counts")
  )
}

tsmessage <- function(..., verbose = TRUE) {
  if(isTRUE(verbose)) cat(format(Sys.time(), "%T"), ..., "\n")
}

make.gene.sets <- function(xx, n.genes = 250, nc = 4) {
  sq <- seq(
    1,
    ncol(xx),
    ifelse(ncol(xx) <= 1000, floor(ncol(xx) / 2), 1000)
  )
  
  pb <- progress_bar$new(
    format = "  processing [:bar] :current/:total (:percent) [eta::eta]",
    total = length(sq),
    clear = FALSE,
    width = 60
  )
  pb$tick(0)
  
  cl <- snow::makeCluster(nc)
  l <- NULL
  
  for(i in 2:length(sq)) {
    x <- as.matrix(xx[, sq[i - 1]:(sq[i] - 1)])
    snow::clusterExport(cl, "x", envir = environment())
    snow::clusterExport(cl, "n.genes", envir = environment())
    
    l <- c(l, snow::parLapply(cl, seq_len(ncol(x)), function(idx) {
      names(x[order(x[, idx], decreasing = TRUE)[1:n.genes], idx])
    }))
    pb$tick()
  }
  
  if(sq[length(sq)] <= ncol(xx)) {
    x <- as.matrix(xx[, sq[i]:ncol(xx)])
    snow::clusterExport(cl, "x", envir = environment())
    snow::clusterExport(cl, "n.genes", envir = environment())
    
    l <- c(l, snow::parLapply(cl, seq_len(ncol(x)), function(idx) {
      names(x[order(x[, idx], decreasing = TRUE)[1:n.genes], idx])
    }))
  }
  
  snow::stopCluster(cl)
  pb$tick()
  
  names(l) <- colnames(xx)
  l
}

assignInNamespace("make.gene.sets", make.gene.sets, ns = "DREEP")

testDREEP <- function(
    M,
    n.markers = 250,
    cores = 0,
    gsea = "simple",
    gpds.signatures = c("CTRP2", "GDSC"),
    th.fdr = 0.1,
    verbose = TRUE,
    storeCellMarkers = FALSE,
    rescalePerc = FALSE,
    useUncorrectedPvals = FALSE,
    th.pval = 0.05
) {
  gpds.signatures <- toupper(gpds.signatures)
  gpds.signatures <- unique(base::match.arg(
    arg = gpds.signatures,
    choices = c("CTRP2", "GDSC", "PRISM"),
    several.ok = TRUE
  ))
  
  gsea <- base::match.arg(
    arg = gsea,
    choices = c("simple", "multilevel"),
    several.ok = FALSE
  )
  
  if(cores == 0) cores <- detectCores()
  
  tsmessage("Loading GPDS signatures..", verbose = verbose)
  
  meta.drug <- NULL
  M.drug.mkrs <- NULL
  
  pb <- progress_bar$new(
    format = "  processing [:bar] :current/:total (:percent) [eta::eta]",
    total = length(gpds.signatures),
    clear = FALSE,
    width = 60
  )
  pb$tick(0)
  
  for(i in gpds.signatures) {
    M.drug.mkrs <- cbind(
      M.drug.mkrs,
      readRDS(paste0(find.package("DREEP"), "/data/", i, ".gpds.rds"))
    )
    
    meta.drug <- rbind(
      meta.drug,
      readRDS(paste0(find.package("DREEP"), "/data/", i, ".drug.metadata.rds"))
    )
    pb$tick()
  }
  
  M <- M[rownames(M) %in% rownames(M.drug.mkrs), ]
  M.drug.mkrs <- M.drug.mkrs[rownames(M.drug.mkrs) %in% rownames(M), ]
  
  tsmessage("Extracting top scoring genes from each cell..", verbose = verbose)
  
  l.mkrs <- make.gene.sets(M, n.markers, 1)
  cells <- colnames(M)
  
  rm(M)
  gc(reset = TRUE)
  
  tsmessage("Running DREEP...", verbose = verbose)
  
  cl <- snow::makeCluster(cores)
  genes <- rownames(M.drug.mkrs)
  
  snow::clusterExport(
    cl,
    c("l.mkrs", "genes", "gsea"),
    envir = environment()
  )
  
  r <- snow::parApply(cl, M.drug.mkrs, 2, function(x) {
    names(x) <- genes
    
    if(gsea == "simple") {
      res <- fgsea::fgseaSimple(
        pathways = l.mkrs,
        stats = x,
        nproc = 1,
        gseaParam = 0,
        nperm = 1000
      )
    } else {
      res <- fgsea::fgseaMultilevel(
        pathways = l.mkrs,
        stats = x,
        nproc = 1,
        gseaParam = 0,
        eps = 0
      )
    }
    
    as.data.frame(res[, c("pval", "ES")])
  })
  
  snow::stopCluster(cl)
  tsmessage("DREEP Running finished!!", verbose = verbose)
  
  len_list <- sapply(r, nrow)
  min_len <- min(len_list)
  
  if(length(unique(len_list)) > 1) {
    message(
      "Different fgsea result lengths detected. Truncating to ",
      min_len,
      " cells."
    )
  }
  
  M.es <- sapply(r, function(x) x[1:min_len, "ES"])
  M.pval <- sapply(r, function(x) x[1:min_len, "pval"])
  
  rownames(M.es) <- cells[1:min_len]
  rownames(M.pval) <- cells[1:min_len]
  
  tsmessage("Adjusting pvalues..", verbose = verbose)
  
  M.fdr <- apply(M.pval, 2, p.adjust, method = "fdr")
  
  tmp <- reshape2::melt(M.es)
  
  if(!useUncorrectedPvals) {
    tmp$fdr <- reshape2::melt(M.fdr)$value
  } else {
    tmp$fdr <- reshape2::melt(M.pval)$value
    th.fdr <- th.pval
  }
  
  if(rescalePerc) {
    df <- tmp %>%
      group_by(Var2) %>%
      summarize(
        sens = sum(value < 0 & fdr < th.fdr) / sum(fdr < th.fdr),
        res = sum(value > 0 & fdr < th.fdr) / sum(fdr < th.fdr),
        med = median(value)
      )
  } else {
    df <- tmp %>%
      group_by(Var2) %>%
      summarize(
        sens = sum(value < 0 & fdr < th.fdr) / length(fdr),
        res = sum(value > 0 & fdr < th.fdr) / length(fdr),
        med = median(value)
      )
  }
  
  df <- as.data.frame(df)
  rm(tmp)
  
  rownames(df) <- df$Var2
  df <- df[, -1]
  
  df$conpound.id <- sapply(
    strsplit(names(r), "_", fixed = TRUE),
    function(x) x[[2]]
  )
  
  df$conpound <- sapply(
    strsplit(df$conpound.id, ":", fixed = TRUE),
    function(x) x[[1]]
  )
  
  df$dataset <- sapply(
    strsplit(names(r), "_", fixed = TRUE),
    function(x) x[[1]]
  )
  
  rownames(df) <- df$conpound.id
  df <- df[order(df$med, decreasing = FALSE), ]
  
  df$target <- meta.drug$target[
    fastmatch::fmatch(df$conpound.id, meta.drug$id)
  ]
  df$moa <- meta.drug$moa[
    fastmatch::fmatch(df$conpound.id, meta.drug$id)
  ]
  df$smiles <- meta.drug$smiles[
    fastmatch::fmatch(df$conpound.id, meta.drug$id)
  ]
  df$dataset <- meta.drug$dataset[
    fastmatch::fmatch(df$conpound.id, meta.drug$id)
  ]
  
  if(storeCellMarkers) {
    return(list(
      df = df,
      es.mtx = M.es,
      es.pval = M.pval,
      es.fdr = M.fdr,
      l.mkrs = l.mkrs
    ))
  }
  
  list(
    df = df,
    es.mtx = M.es,
    es.pval = M.pval,
    es.fdr = M.fdr
  )
}

message("Loading input: ", input_file)
obj <- load_seurat_object(input_file)

message("Loading gene annotation...")
ref <- readRDS(reference_file)
gene_info <- ref$gene_map
gene_map <- setNames(gene_info$gene_id, gene_info$hgnc_symbol)
rm(ref, gene_info)
gc()

message("Running GFICF...")
counts <- get_counts(obj)
gficf_res <- gficf::gficf(M = counts, verbose = TRUE)
ndata <- gficf_res$gficf

rm(gficf_res, counts)
gc()

message("Converting gene symbols to Ensembl IDs...")
mapped_ids <- gene_map[rownames(ndata)]
new_rownames <- ifelse(
  is.na(mapped_ids),
  rownames(ndata),
  mapped_ids
)
rownames(ndata) <- new_rownames

message("Running DREEP...")
dreep_data <- testDREEP(
  M = ndata,
  n.markers = 500,
  gsea = "simple",
  gpds.signatures = c("CTRP2", "GDSC")
)

matches <- grepl(
  target_drug,
  toupper(colnames(dreep_data$es.mtx)),
  fixed = TRUE
)

if(!any(matches)) {
  stop("No DREEP signature found for drug: ", target_drug)
}

matched_signatures <- colnames(dreep_data$es.mtx)[matches]

message("Matched signatures:")
print(matched_signatures)

pred_df <- data.frame(
  cells = rownames(dreep_data$es.mtx),
  dreep_data$es.mtx[, matches, drop = FALSE],
  check.names = FALSE
)

if("condition" %in% colnames(obj@meta.data)) {
  metadata_df <- data.frame(
    cells = colnames(obj),
    condition = as.character(obj$condition),
    stringsAsFactors = FALSE
  )
  
  pred_df <- merge(
    metadata_df,
    pred_df,
    by = "cells",
    all.x = TRUE,
    sort = FALSE
  )
}

sample_id <- tools::file_path_sans_ext(basename(input_file))
output_file <- file.path(
  output_dir,
  paste0(sample_id, "_DREEP_pred.csv")
)

write.csv(
  pred_df,
  output_file,
  row.names = FALSE
)

message("Saved: ", output_file)