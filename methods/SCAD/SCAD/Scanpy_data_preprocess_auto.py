import scanpy as sc
import pandas as pd
import numpy as np
import re
import os
import argparse
import sys

def data_preprocess(input_data, save_path, ID, drug):
    raw_counts = pd.read_csv(input_data, index_col=0)  # First column contains gene names and columns represent cells
    raw_counts = raw_counts.T  # Transpose to cells × genes

    if drug == 'KDM5-C70':
        drug = 'KDM5'
    else:
        drug = drug.replace("-", "")

    SAVE_DATA_TO = save_path + '/data/original_norm/' + ID + '_' + drug + '/'
    if not os.path.exists(SAVE_DATA_TO):
        os.makedirs(SAVE_DATA_TO)
        print("Directory ", SAVE_DATA_TO, " Created ")
    else:
        print("Directory ", SAVE_DATA_TO, " already exists")

    adata = sc.AnnData(raw_counts)
    adata.var_names_make_unique()

    '''
    # 2. Quality control
    # This step is omitted because the DRMREF datasets have already undergone quality control.
    adata.obs["n_genes"] = (adata.X > 0).sum(1)
    adata.obs["n_counts"] = adata.X.sum(1)
    adata.var["n_cells"] = (adata.X > 0).sum(0)

    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))  # Applicable to both human and mouse
    adata.obs["mito_frac"] = adata[:, adata.var["mt"]].X.sum(1) / adata.obs["n_counts"]

    min_genes = 200
    min_counts = 500
    max_counts = np.percentile(adata.obs["n_counts"], 99)  # Remove cells with abnormally high total counts
    max_mito = 0.2  # Maximum mitochondrial fraction: 20%

    adata = adata[adata.obs["n_genes"] > min_genes, :]
    adata = adata[adata.obs["n_counts"] > min_counts, :]
    adata = adata[adata.obs["n_counts"] < max_counts, :]
    adata = adata[adata.obs["mito_frac"] < max_mito, :]

    min_cells = 3
    adata = adata[:, adata.var["n_cells"] > min_cells]
    '''

    # 3. Normalization: total-count normalization followed by log transformation
    adata.layers["raw"] = adata.X.copy()  # Store raw counts
    sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize each cell to a total count of 1e4
    sc.pp.log1p(adata)  # Log transformation
    adata.layers["log_normalized"] = adata.X.copy()

    # 4. Z-score scaling
    adata.X = sc.pp.scale(adata.X, zero_center=True, copy=True)
    adata.layers["z_score"] = adata.X.copy()

    # 5. Highly variable gene selection
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, subset=False)
    hvgs = adata.var[adata.var["highly_variable"]].index

    # 6. Retain genes shared between the source and target domains
    ppi_gene_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "original", "2128_ppi_gene_list.txt")
    with open(ppi_gene_path, "r") as f:
        source_genes = set(line.strip() for line in f)  # Source-domain gene list, one gene per line

    common_genes = list(source_genes & set(adata.var_names))
    adata_common = adata[:, common_genes]

    # 7. Export processed expression matrices as TSV files
    pd.DataFrame(
        adata.layers["z_score"],
        index=adata.obs_names,
        columns=adata.var_names
    ).to_csv(
        SAVE_DATA_TO + 'Original_exprs_resp_z.' + drug + '_' + ID + '.tsv',
        sep="\t"
    )

    pd.DataFrame(
        adata[:, hvgs].layers["z_score"],
        index=adata.obs_names,
        columns=hvgs
    ).to_csv(
        SAVE_DATA_TO + 'Original_exprs_resp_z.' + drug + '_' + ID + '_tp4k.tsv',
        sep="\t"
    )

    pd.DataFrame(
        adata_common.layers["z_score"],
        index=adata_common.obs_names,
        columns=common_genes
    ).to_csv(
        SAVE_DATA_TO + 'Original_exprs_resp_z.' + drug + '_' + ID + '_ppi.tsv',
        sep="\t"
    )

def main():
    parser = argparse.ArgumentParser(description="Single-cell RNA-seq Data Preprocessing Script")

    # Input/output paths
    parser.add_argument('--input', type=str, required=True, help="Input CSV file or data folder path")
    parser.add_argument('--output_dir', type=str, default='./data/original_norm/', help="Output root directory")

    # Dataset-specific parameters
    parser.add_argument('--id', type=str, help="Dataset ID")
    parser.add_argument('--drug', type=str, help="Drug name")

    # Preprocessing parameters
    parser.add_argument('--hvg_num', type=int, default=4000, help="Number of highly variable genes")

    args = parser.parse_args()

    out_dir = args.output_dir
    input_csv = args.input
    tar_id = args.id
    tar_drug = args.drug

    data_preprocess(
        input_data=input_csv,
        save_path=out_dir,
        ID=tar_id,
        drug=tar_drug
    )

if __name__ == "__main__":
    main()
