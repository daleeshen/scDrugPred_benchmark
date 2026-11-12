import os
import argparse
import warnings

warnings.filterwarnings("ignore")

import pandas as pd

from cadrres_sc import pp, model, utility


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run CaDRReS-sc prediction from TPM input"
    )
    parser.add_argument(
        "--tpm", required=True,
        help="Path to TPM expression file (gene x cell)"
    )
    parser.add_argument(
        "--outdir", required=True,
        help="Directory to save prediction results"
    )
    parser.add_argument(
        "--model-dir", required=True,
        help="Directory containing pretrained CaDRReS-sc model files"
    )
    parser.add_argument(
        "--gdsc-exp", required=True,
        help="Path to GDSC expression file"
    )
    parser.add_argument(
        "--essential-genes", required=True,
        help="Path to essential gene list"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    tpm_file = args.tpm
    outdir = args.outdir
    model_dir = args.model_dir
    gdsc_exp_path = args.gdsc_exp
    ess_gene_path = args.essential_genes

    model_spec_name_list = [
        "cadrres-wo-sample-bias",
        "cadrres-wo-sample-bias-weight"
    ]

    os.makedirs(outdir, exist_ok=True)

    gdsc_exp_df = pd.read_csv(
        gdsc_exp_path,
        sep="\t",
        index_col=0
    ).groupby(level=0).mean()

    gdsc_log2_mean_fc_exp_df, _ = pp.gexp.normalize_log2_mean_fc(
        gdsc_exp_df
    )

    ess_gene_list = utility.get_gene_list(ess_gene_path)

    sample_name = os.path.basename(tpm_file)
    if sample_name.endswith("_TPM.csv"):
        sample_id = sample_name[:-8]
    else:
        sample_id = os.path.splitext(sample_name)[0]

    gene_exp_df = pd.read_csv(
        tpm_file,
        sep=",",
        index_col=0
    ).groupby(level=0).mean()

    cell_line_log2_mean_fc_exp_df, _ = pp.gexp.normalize_log2_mean_fc(
        gene_exp_df
    )

    test_kernel_df = pp.gexp.calculate_kernel_feature(
        cell_line_log2_mean_fc_exp_df,
        gdsc_log2_mean_fc_exp_df,
        gene_list=ess_gene_list
    )

    for model_spec_name in model_spec_name_list:
        print(f"Running CaDRReS-sc: {model_spec_name}")

        model_file = os.path.join(
            model_dir,
            f"{model_spec_name}_param_dict.pickle"
        )

        cadrres_model = model.load_model(model_file)

        pred_df, _ = model.predict_from_model(
            cadrres_model,
            test_kernel_df,
            model_spec_name
        )

        save_path = os.path.join(
            outdir,
            f"{sample_id}_{model_spec_name}_test_pred.csv"
        )

        pred_df.to_csv(save_path)
        print(f"Saved: {save_path}")


if __name__ == "__main__":
    main()
