library(tidyverse)
library(scales)
library(ggrepel)

# ============================================================
# Paths
# ============================================================
scad_ic50_file <- "/media/shenli/Expansion/SCAD/data/response/GDSC_response.logIC50.all_drugs.tsv"
scad_label_file <- "/media/shenli/Expansion/SCAD/data/response/GDSC_response.all_drugs.tsv"
scdeal_label_file <- "/media/shenli/Expansion/scDrug_benchmark/scDEAL/data/original_drug_data/ALL_label_binary_wf.csv"
pair_file <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison/scDEAL_NewFig5C_dataset_drug.csv"

fig5b_dir <- "/media/shenli/Expansion/DRMREF/New_fig5/Figure5B_label_comparison"
out_scad <- "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source"
out_scdeal <- "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_common_source_scDEAL_label"
dir.create(fig5b_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_scad, recursive = TRUE, showWarnings = FALSE)
dir.create(out_scdeal, recursive = TRUE, showWarnings = FALSE)

min_common_n <- 30

# ============================================================
# Functions
# ============================================================
clean_missing <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN", "N/A", "na", "nan", "n/a", "NULL", "null", "None", "none", "-")] <- NA_character_
  x
}

normalize_cell_id <- function(x) {
  x <- clean_missing(x)
  x <- sub("^DATA\\.", "", x, ignore.case = TRUE)
  x <- sub("\\.0+$", "", x)
  trimws(x)
}

normalize_label <- function(x) {
  x <- tolower(clean_missing(x))
  case_when(
    x %in% c("r", "resistant", "resistance", "0") ~ "R",
    x %in% c("s", "sensitive", "sensitivity", "1") ~ "S",
    TRUE ~ NA_character_
  )
}

normalize_drug_name <- function(x) {
  x |> as.character() |> str_trim() |> str_to_lower() |> str_replace_all("[^a-z0-9]+", "")
}

convert_scad_drug <- function(x) {
  case_when(
    x == "PLX-4720" ~ "PLX4720",
    x == "Palbociclib" ~ "PD-0332991",
    TRUE ~ x
  )
}

read_wide <- function(file) {
  if (grepl("\\.tsv$|\\.txt$", file, ignore.case = TRUE))
    readr::read_tsv(file, col_types = cols(.default = col_character()), name_repair = "minimal", show_col_types = FALSE)
  else
    readr::read_csv(file, col_types = cols(.default = col_character()), name_repair = "minimal", show_col_types = FALSE)
}

standardize_matrix <- function(df) {
  names(df)[1] <- "cell_id"
  df$cell_id <- normalize_cell_id(df$cell_id)
  if (anyDuplicated(df$cell_id[!is.na(df$cell_id)])) stop("Duplicated cell IDs detected.")
  df
}

remove_duplicated_drugs <- function(dfs) {
  duplicated_drugs <- unique(unlist(lapply(dfs, function(df) {
    x <- trimws(names(df)[-1])
    unique(x[duplicated(x) | duplicated(x, fromLast = TRUE)])
  })))
  lapply(dfs, function(df) {
    drug_names <- trimws(names(df)[-1])
    df[, c(TRUE, !drug_names %in% duplicated_drugs), drop = FALSE]
  })
}

label_to_long <- function(df, label_name, drug_name) {
  out <- df |>
    pivot_longer(-cell_id, names_to = drug_name, values_to = label_name, values_drop_na = FALSE)
  out[[label_name]] <- normalize_label(out[[label_name]])
  out |> mutate(drug_key = normalize_drug_name(.data[[drug_name]]))
}

ic50_to_long <- function(df) {
  df |>
    pivot_longer(-cell_id, names_to = "drug_scad_ic50", values_to = "logIC50_raw", values_drop_na = FALSE) |>
    mutate(
      drug_key = normalize_drug_name(drug_scad_ic50),
      logIC50 = readr::parse_number(clean_missing(logIC50_raw))
    ) |>
    select(-logIC50_raw)
}

check_normalized_collision <- function(df, raw_col) {
  x <- df |>
    distinct(drug_key, drug_raw = .data[[raw_col]]) |>
    count(drug_key, name = "n") |>
    filter(n > 1)
  x$drug_key
}

check_duplicate_pairs <- function(df, name) {
  if (df |> count(cell_id, drug_key) |> filter(n > 1) |> nrow() > 0)
    stop(name, " contains duplicated drug-cell pairs.")
}

cohen_kappa <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  if (!length(x)) return(NA_real_)
  po <- mean(x == y)
  pe <- mean(x == "R") * mean(y == "R") + mean(x == "S") * mean(y == "S")
  if (pe == 1) NA_real_ else (po - pe) / (1 - pe)
}

resistant_jaccard <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]; y <- y[keep]
  u <- sum(x == "R" | y == "R")
  if (!u) NA_real_ else sum(x == "R" & y == "R") / u
}

mcnemar_exact <- function(scdeal, scad) {
  a <- sum(scdeal == "S" & scad == "R", na.rm = TRUE)
  b <- sum(scdeal == "R" & scad == "S", na.rm = TRUE)
  if (a + b == 0) 1 else binom.test(min(a, b), a + b, p = 0.5)$p.value
}

rank_auc <- function(score, label) {
  keep <- !is.na(score) & !is.na(label)
  score <- score[keep]; label <- label[keep]
  np <- sum(label == "R"); nn <- sum(label == "S")
  if (!np || !nn) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[label == "R"]) - np * (np + 1) / 2) / (np * nn)
}

# ============================================================
# 1. Read and normalize SCAD / scDEAL data
# ============================================================
raw <- list(
  ic50 = standardize_matrix(read_wide(scad_ic50_file)),
  scad = standardize_matrix(read_wide(scad_label_file)),
  scdeal = standardize_matrix(read_wide(scdeal_label_file))
)
raw <- remove_duplicated_drugs(raw)

scad_ic50_long <- ic50_to_long(raw$ic50)
scad_label_long <- label_to_long(raw$scad, "scad_label", "drug_scad")
scdeal_label_long <- label_to_long(raw$scdeal, "scdeal_label", "drug_scdeal")

collision_keys <- unique(c(
  check_normalized_collision(scad_ic50_long, "drug_scad_ic50"),
  check_normalized_collision(scad_label_long, "drug_scad"),
  check_normalized_collision(scdeal_label_long, "drug_scdeal")
))
collision_keys <- collision_keys[!is.na(collision_keys)]

if (length(collision_keys)) {
  scad_ic50_long <- filter(scad_ic50_long, !drug_key %in% collision_keys)
  scad_label_long <- filter(scad_label_long, !drug_key %in% collision_keys)
  scdeal_label_long <- filter(scdeal_label_long, !drug_key %in% collision_keys)
}

check_duplicate_pairs(scad_ic50_long, "SCAD logIC50")
check_duplicate_pairs(scad_label_long, "SCAD labels")
check_duplicate_pairs(scdeal_label_long, "scDEAL labels")

# ============================================================
# 2. Common drug-cell pairs and label comparison
# ============================================================
master_pairs <- scad_label_long |>
  inner_join(scdeal_label_long, by = c("cell_id", "drug_key")) |>
  left_join(scad_ic50_long |> select(cell_id, drug_key, logIC50), by = c("cell_id", "drug_key")) |>
  mutate(
    drug_display = coalesce(drug_scad, drug_scdeal),
    label_pair = case_when(
      scdeal_label == "S" & scad_label == "S" ~ "SS",
      scdeal_label == "S" & scad_label == "R" ~ "SR",
      scdeal_label == "R" & scad_label == "S" ~ "RS",
      scdeal_label == "R" & scad_label == "R" ~ "RR",
      TRUE ~ NA_character_
    )
  )

complete_pairs <- master_pairs |> filter(!is.na(scad_label), !is.na(scdeal_label))

comparison_by_drug <- complete_pairs |>
  group_by(drug_key, drug_display) |>
  summarise(
    n_complete = n(),
    n_logIC50_complete = sum(!is.na(logIC50)),
    resistant_fraction_scDEAL = mean(scdeal_label == "R"),
    resistant_fraction_SCAD = mean(scad_label == "R"),
    resistant_fraction_difference = resistant_fraction_scDEAL - resistant_fraction_SCAD,
    agreement = mean(scdeal_label == scad_label),
    cohen_kappa = cohen_kappa(scdeal_label, scad_label),
    resistant_jaccard = resistant_jaccard(scdeal_label, scad_label),
    n_SS = sum(label_pair == "SS"),
    n_SR = sum(label_pair == "SR"),
    n_RS = sum(label_pair == "RS"),
    n_RR = sum(label_pair == "RR"),
    mcnemar_p = mcnemar_exact(scdeal_label, scad_label),
    logIC50_AUROC_scDEAL = rank_auc(logIC50, scdeal_label),
    logIC50_AUROC_SCAD = rank_auc(logIC50, scad_label),
    .groups = "drop"
  ) |>
  mutate(mcnemar_FDR = p.adjust(mcnemar_p, method = "BH")) |>
  arrange(cohen_kappa, agreement)

write_csv(master_pairs, file.path(fig5b_dir, "common_drug_cell_line_master_table.csv"))
write_csv(comparison_by_drug, file.path(fig5b_dir, "label_comparison_metrics_by_drug.csv"))

eligible_metrics <- comparison_by_drug |> filter(n_complete >= min_common_n)
if (!nrow(eligible_metrics)) stop("No drug has >= ", min_common_n, " complete common cell lines.")

write_csv(eligible_metrics, file.path(fig5b_dir, "drugs_included_in_Figure5B.csv"))

# ============================================================
# 3. Figure 5B
# ============================================================
label_pair_colors <- c(SS = "#4DBBD5", SR = "#E64B35", RS = "#F39B7F", RR = "#3C5488")
label_pair_names <- c(SS = "S/S", SR = "S/R", RS = "R/S", RR = "R/R")

drug_labels <- eligible_metrics |>
  arrange(cohen_kappa, agreement) |>
  mutate(drug_plot_label = paste0(drug_display, " (n=", n_complete, ")"))

label_pair_data <- complete_pairs |>
  semi_join(eligible_metrics, by = "drug_key") |>
  count(drug_key, drug_display, label_pair) |>
  group_by(drug_key, drug_display) |>
  mutate(proportion = n / sum(n)) |>
  ungroup() |>
  left_join(drug_labels |> select(drug_key, drug_plot_label), by = "drug_key") |>
  mutate(
    drug_plot_label = factor(drug_plot_label, levels = rev(drug_labels$drug_plot_label)),
    label_pair = factor(label_pair, levels = c("SS", "SR", "RS", "RR"))
  )

p_label_pairs <- ggplot(label_pair_data, aes(proportion, drug_plot_label, fill = label_pair)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = label_pair_colors, labels = label_pair_names, drop = FALSE) +
  scale_x_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.005))) +
  labs(x = "Proportion of common cell lines", y = NULL, fill = "scDEAL/SCAD label") +
  theme_classic(base_size = 15) +
  theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11), legend.position = "bottom")

bubble_label_data <- eligible_metrics |>
  mutate(distance = abs(resistant_fraction_scDEAL - resistant_fraction_SCAD)) |>
  slice_max(distance, n = 10, with_ties = FALSE)

p_resistant_fraction <- ggplot(
  eligible_metrics,
  aes(resistant_fraction_SCAD, resistant_fraction_scDEAL, size = n_complete, fill = cohen_kappa)
) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.8) +
  geom_point(shape = 21, color = "black", stroke = 0.6, alpha = 0.95) +
  geom_text_repel(
    data = bubble_label_data,
    aes(label = drug_display),
    inherit.aes = FALSE,
    mapping = aes(x = resistant_fraction_SCAD, y = resistant_fraction_scDEAL, label = drug_display),
    size = 5, max.overlaps = Inf, seed = 123
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = percent_format(accuracy = 1)) +
  scale_fill_gradient2(low = "#3C5488", mid = "#F2F2F2", high = "#D7301F", midpoint = 0, limits = c(-1, 1), oob = squish) +
  scale_size_continuous(range = c(4, 11)) +
  labs(
    x = "Resistant fraction under SCAD labels",
    y = "Resistant fraction under scDEAL labels",
    fill = "Cohen's kappa",
    size = "Common cell lines"
  ) +
  theme_classic(base_size = 15) +
  theme(axis.title = element_text(size = 25), axis.text = element_text(size = 25),
        legend.title = element_text(size = 17), legend.text = element_text(size = 17))

representative_drugs <- eligible_metrics |>
  filter(n_logIC50_complete >= min_common_n, is.finite(cohen_kappa)) |>
  arrange(cohen_kappa)

if (!nrow(representative_drugs)) stop("No drug has sufficient complete logIC50 data.")

idx <- unique(c(1, ceiling(nrow(representative_drugs) / 2), nrow(representative_drugs)))
representative_drugs <- representative_drugs |>
  slice(idx) |>
  mutate(panel_label = paste0(drug_display, "\n\u03BA = ", sprintf("%.2f", cohen_kappa), "; n = ", n_complete))

representative_data <- complete_pairs |>
  semi_join(representative_drugs, by = "drug_key") |>
  left_join(representative_drugs |> select(drug_key, panel_label), by = "drug_key") |>
  filter(!is.na(logIC50)) |>
  mutate(label_pair = factor(label_pair, levels = c("SS", "SR", "RS", "RR")))

p_ic50 <- ggplot(representative_data, aes(label_pair, logIC50, fill = label_pair)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, alpha = 0.8, linewidth = 0.65) +
  geom_jitter(width = 0.13, alpha = 0.3, size = 0.9) +
  facet_wrap(~panel_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = label_pair_colors, labels = label_pair_names, drop = FALSE) +
  scale_x_discrete(labels = label_pair_names, drop = FALSE) +
  labs(x = "Label pair: scDEAL/SCAD", y = "GDSC logIC50") +
  theme_classic(base_size = 15) +
  theme(strip.text = element_text(size = 12.5, face = "bold"),
        axis.title = element_text(size = 15), axis.text.x = element_text(size = 11.5, angle = 35, hjust = 1),
        axis.text.y = element_text(size = 11.5), legend.position = "none")

ggsave(file.path(fig5b_dir, "Figure5B_left.png"), p_label_pairs, width = 8, height = 6, dpi = 300)
ggsave(file.path(fig5b_dir, "Figure5B_right_up.png"), p_resistant_fraction, width = 8, height = 7, dpi = 300)
ggsave(file.path(fig5b_dir, "Figure5B_right_down.png"), p_ic50, width = 8, height = 6, dpi = 300)

# ============================================================
# 4. Figure 5C: construct matched SCAD source datasets
# ============================================================
pairs <- read.csv(pair_file, stringsAsFactors = FALSE, check.names = FALSE)
names(pairs)[1:2] <- c("dataset", "drug")

scenario_table <- tribble(
  ~scenario, ~round, ~base_path,
  "balanced", 1, "/media/shenli/Expansion/SCAD/data/split_norm/1st_round_history",
  "balanced", 2, "/media/shenli/Expansion/DRMREF/2_round_1000cells_res/1000cells/SCAD/data/split_norm",
  "balanced", 3, "/media/shenli/Expansion/DRMREF/3_round_1000cells_res/1000cells/SCAD/data/split_norm",
  "imbalance_1_3", 1, "/media/shenli/Expansion/benchmark_data_unb/1_round_unbalanced_sampling/1_3/metadata_count/data_for_SCAD/data/split_norm",
  "imbalance_1_3", 2, "/media/shenli/Expansion/benchmark_data_unb/2_round_unbalanced_sampling/1_3/metadata_count/data_for_SCAD/data/split_norm",
  "imbalance_1_3", 3, "/media/shenli/Expansion/benchmark_data_unb/3_round_unbalanced_sampling/1_3/metadata_count/data_for_SCAD/data/split_norm",
  "imbalance_1_30", 1, "/media/shenli/Expansion/benchmark_data_unb/1_round_unbalanced_sampling/1_30/metadata_count/data_for_SCAD/data/split_norm",
  "imbalance_1_30", 2, "/media/shenli/Expansion/benchmark_data_unb/2_round_unbalanced_sampling/1_30/metadata_count/data_for_SCAD/data/split_norm",
  "imbalance_1_30", 3, "/media/shenli/Expansion/benchmark_data_unb/3_round_unbalanced_sampling/1_30/metadata_count/data_for_SCAD/data/split_norm"
)

genesets <- c(all = "", ppi = "_ppi", tp4k = "_tp4k")
summary_list <- list()
k <- 1

for (s in seq_len(nrow(scenario_table))) {
  scenario <- scenario_table$scenario[s]
  round_id <- scenario_table$round[s]
  base_path <- scenario_table$base_path[s]
  
  for (i in seq_len(nrow(pairs))) {
    dataset <- pairs$dataset[i]
    drug <- pairs$drug[i]
    drug_key <- normalize_drug_name(drug)
    drug_scad <- convert_scad_drug(drug)
    
    scdeal_response <- scdeal_label_long |>
      filter(drug_key == !!drug_key, !is.na(scdeal_label)) |>
      transmute(
        cell_id,
        response = if_else(scdeal_label == "R", 0, 1)
      )
    
    if (!nrow(scdeal_response))
      stop("No valid scDEAL labels found for drug: ", drug)
    
    if (anyDuplicated(scdeal_response$cell_id))
      stop("Duplicated scDEAL cell IDs for drug: ", drug)
    
    scdeal_vec <- setNames(scdeal_response$response, scdeal_response$cell_id)
    
    for (gs_name in names(genesets)) {
      gs <- genesets[[gs_name]]
      file_name <- paste0("Source_exprs_resp_z.", drug_scad, gs, ".tsv")
      input_file <- file.path(base_path, paste0(dataset, "_", drug_scad), file_name)
      
      if (!file.exists(input_file)) {
        warning("Missing source file: ", input_file)
        next
      }
      
      source <- read.delim(input_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
      source_ids <- normalize_cell_id(rownames(source))
      keep <- source_ids %in% names(scdeal_vec)
      source_common <- source[keep, , drop = FALSE]
      common_ids <- source_ids[keep]
      
      source_scdeal <- source_common
      source_scdeal$response <- as.numeric(scdeal_vec[common_ids])
      
      folder_scad <- file.path(out_scad, scenario, paste0("round_", round_id), paste0(dataset, "_", drug_scad))
      folder_scdeal <- file.path(out_scdeal, scenario, paste0("round_", round_id), paste0(dataset, "_", drug_scad))
      dir.create(folder_scad, recursive = TRUE, showWarnings = FALSE)
      dir.create(folder_scdeal, recursive = TRUE, showWarnings = FALSE)
      
      write.table(source_common, file.path(folder_scad, file_name), sep = "\t", quote = FALSE, col.names = NA)
      write.table(source_scdeal, file.path(folder_scdeal, file_name), sep = "\t", quote = FALSE, col.names = NA)
      
      summary_list[[k]] <- tibble(
        scenario = scenario,
        round = round_id,
        dataset = dataset,
        drug = drug,
        geneset = gs_name,
        n_before = nrow(source),
        n_common = nrow(source_common),
        SCAD_R = sum(source_common$response == 0, na.rm = TRUE),
        SCAD_S = sum(source_common$response == 1, na.rm = TRUE),
        scDEAL_R = sum(source_scdeal$response == 0, na.rm = TRUE),
        scDEAL_S = sum(source_scdeal$response == 1, na.rm = TRUE),
        missing_scDEAL_label = sum(is.na(source_scdeal$response))
      )
      k <- k + 1
    }
  }
}

summary_df <- bind_rows(summary_list)
write_csv(summary_df, "/media/shenli/Expansion/DRMREF/New_fig5/SCAD_source_label_comparison_summary_all_scenarios.csv")
print(summary_df)
