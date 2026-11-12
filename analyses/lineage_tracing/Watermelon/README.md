# Watermelon Lineage-Tracing Analysis

This folder contains the analysis scripts for the Watermelon lineage-tracing benchmark.

## Files

Analysis scripts:

    analyses/lineage_tracing/Watermelon/
    ├── performance_pres_vs_prer.R
    └── watermelon_analysis.R

Input data:

    data/lineage_tracing/Watermelon/
    ├── GSE150949_metaData_with_lineage.txt
    ├── GSE150949_pc9_count_matrix.csv
    ├── watermelon_D0.rds
    └── Watermelon_preS_vs_preR_performance.csv

## watermelon_analysis.R

This script reconstructs the Watermelon lineage-tracing analysis from the original count matrix and lineage metadata.

It performs:

- construction of the Seurat object
- lineage-based assignment of pre-sensitive, pre-resistant, resistant, and uncertain cells
- generation of D0, D3, D7, and D14 subsets
- normalization, cell-cycle scoring, UMAP, and Harmony integration where required
- differential expression analysis between pre-resistant and pre-sensitive cells
- Hallmark, KEGG, and GO-MF GSEA
- comparison between the intrinsic pre-resistant/pre-sensitive signal and treatment-associated changes at D3, D7, and D14

Before running, please update the working directory, input/output paths, and GMT file path in the script.

The required raw files are located in:

    data/lineage_tracing/Watermelon/

## performance_pres_vs_prer.R

This script summarizes model performance for distinguishing lineage-defined pre-sensitive and pre-resistant cells.

The benchmark uses:

    data/lineage_tracing/Watermelon/watermelon_D0.rds

as the reference Seurat object.

The script reads prediction results from the evaluated methods, calculates AUROC, ACC, F1, Precision, and Recall across three independent sampling rounds, and generates the performance and prediction-score UMAP visualizations.

Please update all method-specific prediction-result paths before running this script.

## Precomputed performance results

The final performance results from the three independent sampling rounds are provided in:

    data/lineage_tracing/Watermelon/Watermelon_preS_vs_preR_performance.csv

This table can be used directly to reproduce the performance visualization without rerunning all prediction methods.

## Re-running predictions

If independent prediction reproduction is required, the expression matrix and metadata can be exported from:

    data/lineage_tracing/Watermelon/watermelon_D0.rds

into the input format required by each method.

The corresponding methods can then be rerun using the scripts provided under:

    methods/

The newly generated prediction results can then be processed using:

    performance_pres_vs_prer.R

## Recommended workflow

    1. Update all paths in watermelon_analysis.R and performance_pres_vs_prer.R.
    2. Run watermelon_analysis.R to reproduce lineage labeling and biological analyses.
    3. Use Watermelon_preS_vs_preR_performance.csv directly for performance visualization, or rerun the individual methods if needed.
    4. Run performance_pres_vs_prer.R to recalculate performance and generate prediction-score UMAPs.
