import os
import sys
import pickle
import argparse

os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import pandas as pd
import numpy as np
np.set_printoptions(precision=2)
from collections import Counter

import importlib
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')


def parse_args():
    parser = argparse.ArgumentParser(description="Train CaDRReS-sc model")
    parser.add_argument("--cadrres-root", required=True, help="Path to the CaDRReS-sc source directory")
    parser.add_argument("--tissue-info", required=True, help="Path to GDSC_tissue_info.csv")
    parser.add_argument("--ic50", required=True, help="Path to gdsc_all_abs_ic50_bayesian_sigmoid_only9dosages.csv")
    parser.add_argument("--drug-stat", required=True, help="Path to drug_stat.csv")
    parser.add_argument("--gdsc-exp", required=True, help="Path to GDSC_exp.tsv")
    parser.add_argument("--essential-genes", required=True, help="Path to essential_genes.txt")
    parser.add_argument("--output-dir", required=True, help="Directory to save trained models")
    parser.add_argument(
        "--model-spec",
        required=True,
        choices=["cadrres-wo-sample-bias", "cadrres-wo-sample-bias-weight"],
        help="CaDRReS-sc model specification"
    )
    return parser.parse_args()


args = parse_args()

scriptpath = args.cadrres_root
sys.path.append(os.path.abspath(scriptpath))

from cadrres_sc import pp, model, evaluation, utility

tissue_sample_df = pd.read_csv(args.tissue_info, index_col=0)
tissue_sample_df.index = tissue_sample_df.index.astype(str)

cell_line_obs_df = pd.read_csv(args.ic50, index_col=0)
cell_line_obs_df.index = cell_line_obs_df.index.astype(str)

# cell lines list which your model will be trained on
cell_line_sample_list = cell_line_obs_df.index.astype(str)

# you can choose specific cancer types for downstream modeling
# cell_line_hn_sample_list = tissue_sample_df[tissue_sample_df['TCGA_CLASS']=='HNSC'].index
# cell_line_hn_obs_df = cell_line_obs_df.loc[cell_line_hn_sample_list]

dataset_drug_df = pd.read_csv(args.drug_stat, index_col=0)
dataset_drug_df.index = dataset_drug_df.index.astype(str)

dataset_drug_list = dataset_drug_df.index
dataset_drug_df.shape

print("Dataframe shape:", dataset_drug_df.shape, "\n")

gene_exp_df = pd.read_csv(args.gdsc_exp, sep='\t', index_col=0)
gene_exp_df = gene_exp_df.groupby(gene_exp_df.index).mean()
print("Dataframe shape:", gene_exp_df.shape, "\n")

cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(gene_exp_df)

ess_gene_list = utility.get_gene_list(args.essential_genes)

cell_line_sample_list = np.array([s for s in cell_line_sample_list if s in gene_exp_df.columns])

cell_line_log2_mean_fc_exp_df = cell_line_log2_mean_fc_exp_df[cell_line_sample_list]
cell_line_obs_df = cell_line_obs_df.loc[cell_line_sample_list, dataset_drug_list]
dataset_drug_df = dataset_drug_df.loc[dataset_drug_list]

print(cell_line_log2_mean_fc_exp_df.shape, cell_line_obs_df.shape, dataset_drug_df.shape)

kernel_feature_df = pp.gexp.calculate_kernel_feature(
    cell_line_log2_mean_fc_exp_df,
    cell_line_log2_mean_fc_exp_df,
    ess_gene_list
).loc[cell_line_sample_list]

# kernel feature based only on training samples
X_train = kernel_feature_df.loc[cell_line_sample_list, cell_line_sample_list]

# observed drug response
Y_train = cell_line_obs_df.loc[cell_line_sample_list]

print("Dataframe shape:", X_train.shape, "\n")
print("Dataframe shape:", Y_train.shape, "\n")

model_spec_name = args.model_spec
print('You have select:' + str(model_spec_name))

indication_specific_degree = 1

output_dir = args.output_dir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print('Results will be saved in ', output_dir)

sample_weights_logistic_x0_df = model.get_sample_weights_logistic_x0(
    dataset_drug_df,
    'log2_max_conc',
    X_train.index
)

indication_weight_df = pd.DataFrame(
    np.ones(Y_train.shape),
    index=Y_train.index,
    columns=Y_train.columns
)

cv_cell_line_sample_list = [
    cl for cl in cell_line_sample_list
    if cl in X_train.index
]

indication_weight_df.loc[cv_cell_line_sample_list, :] = (
    indication_weight_df.loc[cv_cell_line_sample_list, :] *
    indication_specific_degree
)

if model_spec_name in ['cadrres', 'cadrres-wo-sample-bias']:
    cadrres_model_dict, cadrres_output_dict = model.train_model(
        Y_train,
        X_train,
        Y_train,
        X_train,
        10,
        0.0,
        100000,
        0.01,
        model_spec_name=model_spec_name,
        save_interval=5000,
        output_dir=output_dir
    )
elif model_spec_name in ['cadrres-wo-sample-bias-weight']:
    cadrres_model_dict, cadrres_output_dict = model.train_model_logistic_weight(
        Y_train,
        X_train,
        Y_train,
        X_train,
        sample_weights_logistic_x0_df,
        indication_weight_df,
        10,
        0.0,
        100000,
        0.01,
        model_spec_name=model_spec_name,
        save_interval=5000,
        output_dir=output_dir
    )

print('Saving ' + os.path.join(output_dir, '{}_param_dict.pickle'.format(model_spec_name)))
pickle.dump(
    cadrres_model_dict,
    open(os.path.join(output_dir, '{}_param_dict.pickle'.format(model_spec_name)), 'wb')
)

print('Saving ' + os.path.join(output_dir, '{}_output_dict.pickle'.format(model_spec_name)))
pickle.dump(
    cadrres_output_dict,
    open(os.path.join(output_dir, '{}_output_dict.pickle'.format(model_spec_name)), 'wb')
)
