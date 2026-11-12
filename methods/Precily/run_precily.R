suppressPackageStartupMessages({
  library(reticulate)
})

use_python("/opt/conda/envs/precily/bin/python", required = TRUE)

suppressPackageStartupMessages({
  library(keras)
  library(GSVA)
  library(impute)
  library(GSEABase)
  library(BiocParallel)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 5) {
  stop("Usage: Rscript run_precily.R <logTPM_file> <metadata_file> <drug> <cancer_type> <output_dir>")
}

logtpm_file <- args[1]
metadata_file <- args[2]
current_drug <- args[3]
cancerType <- args[4]
output_dir <- args[5]

if(!file.exists(logtpm_file)) stop("logTPM file does not exist: ", logtpm_file)
if(!file.exists(metadata_file)) stop("Metadata file does not exist: ", metadata_file)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

PRECILY_DIR <- "/media/shenli/Expansion/DRMREF/scDrugBenchmark/methods/Precily/Precily"
FEATURE_FILE <- file.path(PRECILY_DIR, "Features.csv")
MODEL_FILES <- file.path(PRECILY_DIR, paste0("Model", 1:5, ".hdf5"))
GMT_FILE <- file.path(PRECILY_DIR, "Fig2/Fig2b/c2.cp.v6.1.symbols.gmt")
GDSC_METADATA_FILE <- file.path(PRECILY_DIR, "Fig2/Fig2b/GDSC2_metadata.RData")

required_files <- c(FEATURE_FILE, MODEL_FILES, GMT_FILE, GDSC_METADATA_FILE)
missing_files <- required_files[!file.exists(required_files)]
if(length(missing_files) > 0) {
  stop("Missing Precily files:\n", paste(missing_files, collapse = "\n"))
}

load_env <- new.env()
load(GDSC_METADATA_FILE, envir = load_env)
if(!exists("metadata", envir = load_env)) {
  stop("Object 'metadata' was not found in: ", GDSC_METADATA_FILE)
}
gdsc_metadata <- load_env$metadata
rm(load_env)

geneSets <- getGmt(GMT_FILE)

drugPred_modified <- function(enrichment.scores, metadata, CancerType) {
  log_msg <- function(msg) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
    flush.console()
  }
  
  pred <- list()
  
  log_msg("Reading Features.csv ...")
  features <- read.table(FEATURE_FILE, sep = ",")[, 1]
  log_msg(sprintf("Number of features: %d", length(features)))
  
  log_msg("Loading 5 pretrained models ...")
  model1 <- load_model_hdf5(MODEL_FILES[1])
  model2 <- load_model_hdf5(MODEL_FILES[2])
  model3 <- load_model_hdf5(MODEL_FILES[3])
  model4 <- load_model_hdf5(MODEL_FILES[4])
  model5 <- load_model_hdf5(MODEL_FILES[5])
  log_msg("All models loaded successfully.")
  
  total_samples <- ncol(enrichment.scores)
  log_msg(sprintf("Start prediction for %d samples ...", total_samples))
  
  for(i in 1:ncol(enrichment.scores)) {
    expAve <- as.matrix(enrichment.scores[, i])
    common_features <- intersect(c(rownames(expAve), colnames(metadata)), features)
    Final_Index <- which(metadata$TCGA_DESC %in% CancerType)
    metainfo <- metadata[Final_Index, ]
    col <- dim(metainfo)[1]
    
    output <- sapply(features, function(x) {
      if(x %in% common_features) {
        if(x %in% colnames(metainfo)) {
          temp <- metainfo[, x]
        } else if(x %in% rownames(expAve)) {
          temp <- rep(expAve[x, ], col)
        }
      } else {
        temp <- rep(NA, col)
      }
      return(temp)
    })
    
    log_msg("  Feature matrix constructed.")
    
    impdata <- impute.knn(t(output[, 1:1329]))
    xtest <- cbind(t(impdata$data), output[, 1330:ncol(output)])
    
    log_msg(sprintf(
      "  Missing value imputation done. xtest dim: %d × %d",
      nrow(xtest), ncol(xtest)
    ))
    
    log_msg("  Running 5-model ensemble prediction ...")
    
    prediction1 <- model1 %>% predict(xtest)
    prediction2 <- model2 %>% predict(xtest)
    prediction3 <- model3 %>% predict(xtest)
    prediction4 <- model4 %>% predict(xtest)
    prediction5 <- model5 %>% predict(xtest)
    
    predictions <- apply(
      cbind.data.frame(prediction1, prediction2, prediction3, prediction4, prediction5),
      1,
      mean
    )
    predictions <- round(predictions, 4)
    
    data <- cbind.data.frame(metainfo$DRUG_NAME, predictions)
    colnames(data) <- c("DRUGS", "Predictedvalue")
    mat <- unique(data)
    
    colnames(mat)[2] <- colnames(enrichment.scores)[i]
    pred[[i]] <- mat[order(mat[, 2]), ]
  }
  
  log_msg("All samples processed successfully.")
  return(pred)
}

message("Loading logTPM matrix...")
log_tpm_matrix <- read.csv(
  logtpm_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

sample_metadata <- read.csv(
  metadata_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if(!"condition" %in% colnames(sample_metadata)) {
  stop("Metadata column 'condition' was not found.")
}

matched_drug <- unique(
  gdsc_metadata$DRUG_NAME[
    toupper(gdsc_metadata$DRUG_NAME) == toupper(current_drug)
  ]
)

if(length(matched_drug) == 0) {
  stop(current_drug, " not found in Precily GDSC2 metadata.")
}

if(!cancerType %in% gdsc_metadata$TCGA_DESC) {
  stop(cancerType, " not found in Precily TCGA_DESC.")
}

message("Drug: ", matched_drug[1])
message("Cancer type: ", cancerType)
message("Cells: ", ncol(log_tpm_matrix))

message("Running GSVA...")
param <- gsvaParam(
  exprData = as.matrix(log_tpm_matrix),
  geneSets = geneSets,
  kcdf = "Gaussian",
  maxDiff = TRUE
)

enrichment.scores <- gsva(
  param,
  verbose = TRUE,
  BPPARAM = SerialParam(progressbar = TRUE)
)

message("GSVA completed. Running Precily...")

df1 <- drugPred_modified(
  enrichment.scores,
  gdsc_metadata,
  cancerType
)

drug_score <- sapply(seq_along(df1), function(i) {
  x <- df1[[i]]
  idx <- toupper(as.character(x$DRUGS)) == toupper(current_drug)
  
  if(any(idx)) {
    as.numeric(x[idx, 2][1])
  } else {
    NA_real_
  }
})

names(drug_score) <- colnames(enrichment.scores)

pred_res <- data.frame(
  cells = names(drug_score),
  condition = sample_metadata[names(drug_score), "condition"],
  precily_score = as.numeric(drug_score),
  stringsAsFactors = FALSE
)

pred_res <- pred_res[!is.na(pred_res$precily_score), ]

sample_id <- sub("_logTPM\\.csv$", "", basename(logtpm_file), ignore.case = TRUE)
output_file <- file.path(
  output_dir,
  paste0(sample_id, "_Precily_pred.csv")
)

write.csv(
  pred_res,
  output_file,
  row.names = FALSE,
  quote = FALSE
)

message("Completed.")
message("Prediction file: ", output_file)