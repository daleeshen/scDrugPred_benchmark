# Benchmark Dataset Construction

This directory contains scripts used to generate the balanced and imbalanced benchmark datasets.

## Scripts

- `generate_balanced_datasets.R`: generates balanced subsets with 500, 1000, 2000, 3000, and 5000 cells.
- `generate_imbalanced_datasets.R`: generates imbalanced subsets with ratios of 1:3, 1:10, 1:30, and 1:100, with a maximum of 10,000 cells per subset.

## Usage

Balanced datasets:
`Rscript generate_balanced_datasets.R <input_dir> <output_dir>`

Imbalanced datasets:
`Rscript generate_imbalanced_datasets.R <input_dir> <output_dir>`

Input files should be annotated Seurat objects containing a binary `condition` column in `meta.data`.
