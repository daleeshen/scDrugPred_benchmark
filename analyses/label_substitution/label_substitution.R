library(tidyverse)

# SCAD architecture + scDEAL labels
# Prerequisite: SCAD must first be run using its original labels to generate the
# Source_exprs_resp_z.*.tsv files for each dataset-drug pair, geneset, scenario,
# and sampling round. This section reuses those SCAD source-expression files,
# retains the common source cell lines with available scDEAL labels, and replaces
# only the response column with the corresponding scDEAL labels.

pair_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"
scdeal_label_file <- "/media/shenli/Expansion/scDrug_benchmark/scDEAL/data/original_drug_data/ALL_label_binary_wf.csv"
scad_label_file <- "/media/shenli/Expansion/SCAD/data/response/GDSC_response.all_drugs.tsv"

out_scad <- "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source"
out_scad_scdeal_label <- "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source_scDEAL_label"
out_scdeal_scad_label <- "/media/shenli/Expansion/DRMREF/New_fig5/scDEAL_SCAD_label"

dir.create(out_scad, recursive = TRUE, showWarnings = FALSE)
dir.create(out_scad_scdeal_label, recursive = TRUE, showWarnings = FALSE)
dir.create(out_scdeal_scad_label, recursive = TRUE, showWarnings = FALSE)

normalize_cell_id <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^DATA\\.", "", x, ignore.case = TRUE)
  x <- sub("\\.0+$", "", x)
  x
}

normalize_drug <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("^X(?=[0-9])", "", x, perl = TRUE)
  gsub("[^A-Z0-9]", "", x)
}

normalize_label <- function(x, output = c("RS","word","numeric")) {
  output <- match.arg(output)
  x <- tolower(trimws(as.character(x)))
  z <- case_when(
    x %in% c("r","resistant","resistance","0") ~ "R",
    x %in% c("s","sensitive","sensitivity","1") ~ "S",
    TRUE ~ NA_character_
  )
  if(output == "word") return(ifelse(z == "R","resistant",ifelse(z == "S","sensitive",NA)))
  if(output == "numeric") return(ifelse(z == "R",0,ifelse(z == "S",1,NA)))
  z
}

convert_scad_drug <- function(drug) {
  case_when(
    drug == "PLX-4720" ~ "PLX4720",
    drug == "Palbociclib" ~ "PD-0332991",
    TRUE ~ drug
  )
}

pairs <- read.csv(pair_file, stringsAsFactors = FALSE, check.names = FALSE)
scdeal_raw <- read.csv(scdeal_label_file, stringsAsFactors = FALSE, check.names = FALSE)
scad_raw <- read.delim(scad_label_file, stringsAsFactors = FALSE, check.names = FALSE)

scdeal_first_col_name <- names(scdeal_raw)[1]
names(scdeal_raw)[1] <- "cell_id"
names(scad_raw)[1] <- "cell_id"

scdeal_raw$cell_id_raw <- scdeal_raw$cell_id
scdeal_raw$cell_id <- normalize_cell_id(scdeal_raw$cell_id)
scad_raw$cell_id <- normalize_cell_id(scad_raw$cell_id)

scdeal_drugs <- setdiff(names(scdeal_raw), c("cell_id","cell_id_raw"))
scad_drugs <- setdiff(names(scad_raw), "cell_id")
scdeal_keys <- normalize_drug(scdeal_drugs)
scad_keys <- normalize_drug(scad_drugs)

scenario_table <- tibble(
  scenario = c("balanced","balanced","balanced","imbalance_1_3","imbalance_1_3","imbalance_1_3","imbalance_1_30","imbalance_1_30","imbalance_1_30"),
  round = rep(1:3,3),
  base_path = c(
    "/media/shenli/Expansion/SCAD/data/split_norm/1st_round_history",
    "/media/shenli/Expansion/DRMREF/2_round_1000cells_res/1000cells/SCAD/data/split_norm",
    "/media/shenli/Expansion/DRMREF/3_round_1000cells_res/1000cells/SCAD/data/split_norm",
    "/media/shenli/Expansion/benchmark_data_unb/1_round_unbalanced_sampling/1_3/metadata_count/data_for_SCAD/data/split_norm",
    "/media/shenli/Expansion/benchmark_data_unb/2_round_unbalanced_sampling/1_3/metadata_count/data_for_SCAD/data/split_norm",
    "/media/shenli/Expansion/benchmark_data_unb/3_round_unbalanced_sampling/1_3/metadata_count/data_for_SCAD/data/split_norm",
    "/media/shenli/Expansion/benchmark_data_unb/1_round_unbalanced_sampling/1_30/metadata_count/data_for_SCAD/data/split_norm",
    "/media/shenli/Expansion/benchmark_data_unb/2_round_unbalanced_sampling/1_30/metadata_count/data_for_SCAD/data/split_norm",
    "/media/shenli/Expansion/benchmark_data_unb/3_round_unbalanced_sampling/1_30/metadata_count/data_for_SCAD/data/split_norm"
  )
)

genesets <- c("","_ppi","_tp4k")
scad_source_summary <- list()

for(s in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table$scenario[s]
  round_id <- scenario_table$round[s]
  base_path <- scenario_table$base_path[s]
  for(i in seq_len(nrow(pairs))) {
    dataset <- as.character(pairs[i,1])
    drug <- as.character(pairs[i,2])
    drug_scad <- convert_scad_drug(drug)
    key <- normalize_drug(drug)
    scdeal_col <- scdeal_drugs[match(key, scdeal_keys)]
    if(is.na(scdeal_col)) stop("Cannot match scDEAL drug: ", drug)
    scdeal_response <- normalize_label(scdeal_raw[[scdeal_col]], "numeric")
    names(scdeal_response) <- scdeal_raw$cell_id
    valid_scdeal_ids <- names(scdeal_response)[!is.na(scdeal_response)]
    for(gs in genesets) {
      file_name <- paste0("Source_exprs_resp_z.", drug_scad, gs, ".tsv")
      input_file <- file.path(base_path, paste0(dataset,"_",drug_scad), file_name)
      if(!file.exists(input_file)) stop("File does not exist: ", input_file)
      source <- read.delim(input_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
      source_ids <- normalize_cell_id(rownames(source))
      keep <- source_ids %in% valid_scdeal_ids
      source_common <- source[keep,,drop = FALSE]
      source_common_ids <- source_ids[keep]
      source_scdeal <- source_common
      source_scdeal$response <- as.numeric(scdeal_response[source_common_ids])
      output_folder_scad <- file.path(out_scad, scenario, paste0("round_",round_id), paste0(dataset,"_",drug_scad))
      output_folder_scdeal <- file.path(out_scad_scdeal_label, scenario, paste0("round_",round_id), paste0(dataset,"_",drug_scad))
      dir.create(output_folder_scad, recursive = TRUE, showWarnings = FALSE)
      dir.create(output_folder_scdeal, recursive = TRUE, showWarnings = FALSE)
      write.table(source_common, file.path(output_folder_scad,file_name), sep = "\t", quote = FALSE, col.names = NA)
      write.table(source_scdeal, file.path(output_folder_scdeal,file_name), sep = "\t", quote = FALSE, col.names = NA)
      scad_source_summary[[length(scad_source_summary)+1]] <- data.frame(
        scenario = scenario, round = round_id, dataset = dataset, drug = drug,
        geneset = ifelse(gs == "","all",sub("^_","",gs)),
        n_before = nrow(source), n_common = nrow(source_common),
        SCAD_R = sum(source_common$response == 0, na.rm = TRUE),
        SCAD_S = sum(source_common$response == 1, na.rm = TRUE),
        scDEAL_R = sum(source_scdeal$response == 0, na.rm = TRUE),
        scDEAL_S = sum(source_scdeal$response == 1, na.rm = TRUE),
        missing_scDEAL_label = sum(is.na(source_scdeal$response))
      )
    }
  }
}

scad_source_summary <- bind_rows(scad_source_summary)
write.csv(scad_source_summary, file.path(dirname(out_scad),"SCAD_source_label_substitution_summary.csv"), row.names = FALSE)

substitution_drugs <- c(
  "ERLOTINIB","PACLITAXEL","SORAFENIB","CRIZOTINIB","DASATINIB","BORTEZOMIB",
  "LAPATINIB","GEMCITABINE","VINORELBINE","X5.FLUOROURACIL","RUXOLITINIB",
  "CAMPTOTHECIN","VINBLASTINE","CISPLATIN","CYTARABINE","DOCETAXEL",
  "GEFITINIB","VORINOSTAT","NILOTINIB","OLAPARIB","AXITINIB","AFATINIB",
  "TAMOXIFEN","TRAMETINIB","DABRAFENIB","TEMOZOLOMIDE"
)

excluded_drugs <- c(
  "AZD6482","AZD7762","AZD8055","BI.2536","BMS.345541",
  "BMS.754807","MG.132","MK.2206","OSI.027"
)

substitution_keys <- normalize_drug(substitution_drugs)
scad_substitution_cols <- scad_drugs[match(substitution_keys, scad_keys)]
if(any(is.na(scad_substitution_cols))) stop("Some SCAD substitution drugs cannot be matched.")

common_cells <- intersect(scdeal_raw$cell_id, scad_raw$cell_id)
scad_substitution_df <- scad_raw[
  match(common_cells,scad_raw$cell_id),
  c("cell_id",scad_substitution_cols),
  drop = FALSE
]

for(j in 2:ncol(scad_substitution_df)) scad_substitution_df[[j]] <- normalize_label(scad_substitution_df[[j]], "word")
complete_idx <- complete.cases(scad_substitution_df[,-1])
selected_cells <- scad_substitution_df$cell_id[complete_idx]

if(length(selected_cells) != 224) stop("Unexpected number of selected cell lines: ", length(selected_cells))

scad_substitution_df <- scad_substitution_df[complete_idx,,drop = FALSE]
rownames(scad_substitution_df) <- scad_substitution_df$cell_id

scdeal_selected <- scdeal_raw[
  match(selected_cells,scdeal_raw$cell_id),
  ,
  drop = FALSE
]
rownames(scdeal_selected) <- scdeal_selected$cell_id

scdeal_scad_part <- data.frame(cell_id = selected_cells, stringsAsFactors = FALSE)
for(i in seq_along(substitution_drugs)) {
  scdeal_drug <- substitution_drugs[i]
  scad_drug <- scad_substitution_cols[i]
  scdeal_scad_part[[scdeal_drug]] <- scad_substitution_df[selected_cells,scad_drug]
}

remaining_drugs <- scdeal_drugs[
  !normalize_drug(scdeal_drugs) %in% substitution_keys &
    !normalize_drug(scdeal_drugs) %in% normalize_drug(excluded_drugs)
]

scdeal_remaining_part <- scdeal_selected[,c("cell_id",remaining_drugs),drop = FALSE]

scdeal_SCAD_label <- cbind(
  scdeal_scad_part,
  scdeal_remaining_part[,-1,drop = FALSE]
)

if(dim(scdeal_SCAD_label)[1] != 224 || dim(scdeal_SCAD_label)[2] != 75) {
  stop("Unexpected final scDEAL label matrix dimension: ", paste(dim(scdeal_SCAD_label), collapse = " x "))
}
if(anyNA(scdeal_SCAD_label)) stop("SCAD-labeled scDEAL matrix contains NA.")

raw_ids <- scdeal_raw$cell_id_raw[match(scdeal_SCAD_label$cell_id,scdeal_raw$cell_id)]
scdeal_SCAD_label[[1]] <- raw_ids
names(scdeal_SCAD_label)[1] <- scdeal_first_col_name

write.csv(
  scdeal_SCAD_label,
  file.path(out_scdeal_scad_label,"ALL_label_binary_wf.csv"),
  row.names = FALSE,
  quote = FALSE
)

substitution_check <- data.frame(
  drug = substitution_drugs,
  scad_drug = scad_substitution_cols,
  n_cells = length(selected_cells),
  n_R = sapply(scad_substitution_cols,function(x) sum(scad_substitution_df[[x]] == "resistant")),
  n_S = sapply(scad_substitution_cols,function(x) sum(scad_substitution_df[[x]] == "sensitive"))
)

write.csv(
  substitution_check,
  file.path(out_scdeal_scad_label,"SCAD_label_substitution_summary.csv"),
  row.names = FALSE
)

cat("SCAD architecture + SCAD labels:", out_scad, "\n")
cat("SCAD architecture + scDEAL labels:", out_scad_scdeal_label, "\n")
cat("scDEAL architecture + SCAD labels:", file.path(out_scdeal_scad_label,"ALL_label_binary_wf.csv"), "\n")
cat("scDEAL SCAD-label matrix dimension:", nrow(scdeal_SCAD_label), "x", ncol(scdeal_SCAD_label), "\n")