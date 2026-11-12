# Copyright 2023 BioMap (Beijing) Intelligence Technology Limited

import argparse
import random, os
import numpy as np
import pandas as pd
import argparse
import torch
from tqdm import tqdm
import scipy.sparse
from scipy.sparse import issparse
import scanpy as sc
from load import *
import pickle

####################################Settings#################################

parser = argparse.ArgumentParser(description='Drug_response_pre')

parser.add_argument(
    '--task_name',
    type=str,
    default='deepcdr',
    help='task name'
)

parser.add_argument(
    '--input_type',
    type=str,
    default='singlecell',
    choices=['singlecell', 'bulk'],
    help='input type; default: singlecell'
)

parser.add_argument(
    '--output_type',
    type=str,
    default='cell',
    choices=['cell', 'gene', 'gene_batch', 'gene_expression'],
    help='cell or gene embedding; default: cell the difference between gene and gene_batch is that in gene mode the gene embedding will be processed one by one. while in gene_batch mode, the gene embedding will be processed in batch. GEARS use gene_batch mode.'
)

parser.add_argument(
    '--pool_type',
    type=str,
    default='all',
    choices=['all', 'max'],
    help='pooling type of cell embedding; default: all only valid for output_type=cell'
)

parser.add_argument(
    '--tgthighres',
    type=str,
    default='t4',
    help='the targeted high resolution (start with t) or the fold change of the high resolution (start with f), or the addtion (start with a) of the high resoultion. only valid for input_type=singlecell'
)

# CHANGED:
# data_path now receives one TXT file instead of a directory.
parser.add_argument(
    '--data_path',
    type=str,
    required=True,
    help='path to one input TXT file'
)

parser.add_argument(
    '--save_path',
    type=str,
    default='./',
    help='save path'
)

parser.add_argument(
    '--pre_normalized',
    type=str,
    default='F',
    choices=['F', 'T', 'A'],
    help='if normalized before input; default: False (F). choice: True(T), Append(A) When input_type=bulk: pre_normalized=T means log10(sum of gene expression). pre_normalized=F means sum of gene expression without normalization. When input_type=singlecell: pre_normalized=T or F means gene expression is already normalized+log1p or not. pre_normalized=A means gene expression is normalized and log1p transformed. the total count is appended to the end of the gene expression matrix.'
)

parser.add_argument(
    '--demo',
    action='store_true',
    default=False,
    help='if demo, only infer 10 samples'
)

parser.add_argument(
    '--version',
    type=str,
    default='ce',
    help='only valid for output_type=cell. For read depth enhancemnet, version=rde For others, version=ce'
)

parser.add_argument(
    '--model_path',
    type=str,
    required=True,
    help='pre-trained model path'
)

parser.add_argument(
    '--ckpt_name',
    type=str,
    default='01B-resolution',
    help='checkpoint name'
)

parser.add_argument(
    '--gene_index',
    type=str,
    required=True,
    help='path to OS_scRNA_gene_index.19264.tsv'
)

parser.add_argument(
    '--sample_output',
    type=str,
    required=True,
    help='path to output *_samples.pkl file'
)

parser.add_argument(
    '--label_output',
    type=str,
    required=True,
    help='path to output *_labels.npy file'
)

args = parser.parse_args()


def main_gene_selection(X_df, gene_list):
    """
    Describe:
        rebuild the input adata to select target genes encode protein
    Parameters:
        adata->`~anndata.AnnData` object: adata with var index_name by gene symbol
        gene_list->list: wanted target gene
    Returns:
        adata_new->`~anndata.AnnData` object
        to_fill_columns->list: zero padding gene
    """

    to_fill_columns = list(set(gene_list) - set(X_df.columns))

    padding_df = pd.DataFrame(
        np.zeros((X_df.shape[0], len(to_fill_columns))),
        columns=to_fill_columns,
        index=X_df.index
    )

    X_df = pd.DataFrame(
        np.concatenate(
            [df.values for df in [X_df, padding_df]],
            axis=1
        ),
        index=X_df.index,
        columns=list(X_df.columns) + list(padding_df.columns)
    )

    X_df = X_df[gene_list]

    var = pd.DataFrame(index=X_df.columns)

    var['mask'] = [
        1 if i in to_fill_columns else 0
        for i in list(var.index)
    ]

    return X_df, to_fill_columns, var


gene_list_df = pd.read_csv(
    args.gene_index,
    header=0,
    delimiter='\t'
)

gene_list = list(gene_list_df['gene_name'])


def main():
    # Set random seed
    random.seed(0)
    np.random.seed(0)  # numpy random generator

    torch.manual_seed(0)
    torch.cuda.manual_seed_all(0)

    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    ###########################################################################
    # Load one input file
    ###########################################################################

    # CHANGED:
    # args.data_path is now the complete path to one TXT file.
    file_path = args.data_path
    file_name = os.path.basename(file_path)

    gexpr_feature = pd.read_csv(
        file_path,
        sep='\t',
        index_col=0,
        usecols=lambda column: column != 'Condition'
    )

    print(gexpr_feature.columns)

    lbl = list(
        pd.read_csv(
            file_path,
            sep='\t',
            usecols=['Condition']
        )['Condition']
    )

    label = [
        1 if element == 'sensitive' else 0
        for element in lbl
    ]

    print(gexpr_feature.shape)

    if gexpr_feature.shape[1] != 19264:
        print('covert gene feature into 19264')

        gexpr_feature, to_fill_columns, var = main_gene_selection(
            gexpr_feature,
            gene_list
        )

        assert gexpr_feature.shape[1] >= 19264

    print(gexpr_feature)

    if (
        args.pre_normalized == 'F'
        and args.input_type == 'bulk'
    ):
        adata = sc.AnnData(gexpr_feature)

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

        gexpr_feature = pd.DataFrame(
            adata.X,
            index=adata.obs_names,
            columns=adata.var_names
        )

    if args.demo:
        gexpr_feature = gexpr_feature.iloc[:10, :]

    print(gexpr_feature.shape)

    ###########################################################################
    # Load model
    ###########################################################################

    ckpt_path = args.model_path

    key = 'cell'

    pretrainmodel, pretrainconfig = load_model_frommmf(
        ckpt_path,
        key
    )

    # print(pretrainmodel)
    pretrainmodel.eval()

    ###########################################################################
    # Inference
    ###########################################################################

    data_dict = {}

    x_all = []
    x_padding_all = []
    pos_id_all = []

    for i in tqdm(range(gexpr_feature.shape[0])):
        with torch.no_grad():
            # pre-Normalization
            if args.pre_normalized == 'F':
                tmpdata = (
                    np.log1p(
                        gexpr_feature.iloc[i, :]
                        / gexpr_feature.iloc[i, :].sum()
                        * 1e4
                    )
                ).tolist()

            totalcount = gexpr_feature.iloc[i, :].sum()

            # select resolution
            pretrain_gene_x = torch.tensor(
                tmpdata
                + [
                    float(args.tgthighres[1:]),
                    np.log10(totalcount)
                ]
            ).unsqueeze(0).cuda(0)

            data_gene_ids = torch.arange(
                19266,
                device=pretrain_gene_x.device
            ).repeat(
                pretrain_gene_x.shape[0],
                1
            )

        value_labels = pretrain_gene_x > 0

        x, x_padding = gatherData(
            pretrain_gene_x,
            value_labels,
            pretrainconfig['pad_token_id']
        )

        position_gene_ids, _ = gatherData(
            data_gene_ids,
            value_labels,
            pretrainconfig['pad_token_id']
        )

        x_all.extend(x)
        x_padding_all.extend(x_padding)
        pos_id_all.extend(position_gene_ids)

        # Cell embedding

    data_dict['x'] = x_all
    data_dict['x_padding'] = x_padding_all
    data_dict['pos_id'] = pos_id_all

    ###########################################################################
    # Save results
    ###########################################################################

    with open(
        args.sample_output,
        'wb'
    ) as p:
        pickle.dump(data_dict, p)

    np.save(
        args.label_output,
        np.array(label)
    )

    print('saved' + file_name)


if __name__ == '__main__':
    main()