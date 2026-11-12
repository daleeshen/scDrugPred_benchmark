import scanpy as sc
import pandas as pd
import numpy as np
import re
import os
import argparse
import sys

def data_preprocess(input_data, save_path, ID, drug):
    raw_counts = pd.read_csv(input_data, index_col=0)  # 第一列是基因名，列名是细胞名
    raw_counts = raw_counts.T  # **转置，使得行是细胞，列是基因**
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
    # 2. 质控（这一步我省略了，因为DRMREF中的数据已经做过了质控）
    adata.obs["n_genes"] = (adata.X > 0).sum(1)
    adata.obs["n_counts"] = adata.X.sum(1)
    adata.var["n_cells"] = (adata.X > 0).sum(0)

    adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))  # 适用于人类和小鼠
    adata.obs["mito_frac"] = adata[:, adata.var["mt"]].X.sum(1) / adata.obs["n_counts"]

    min_genes = 200  # 适应具体数据集
    min_counts = 500
    max_counts = np.percentile(adata.obs["n_counts"], 99)  # 过滤异常高的细胞
    max_mito = 0.2  # 线粒体占比上限20%

    adata = adata[adata.obs["n_genes"] > min_genes, :]
    adata = adata[adata.obs["n_counts"] > min_counts, :]
    adata = adata[adata.obs["n_counts"] < max_counts, :]
    adata = adata[adata.obs["mito_frac"] < max_mito, :]

    min_cells = 3
    adata = adata[:, adata.var["n_cells"] > min_cells]

    '''

    # 3. 归一化：Total counts normalization & CPM
    adata.layers["raw"] = adata.X.copy()  # 备份原始count
    sc.pp.normalize_total(adata, target_sum=1e4)  # 归一化到1e4
    sc.pp.log1p(adata)  # 对数转换
    adata.layers["log_normalized"] = adata.X.copy()

    # 4. Z-score标准化
    adata.X = sc.pp.scale(adata.X, zero_center=True, copy=True)  # z-score标准化
    adata.layers["z_score"] = adata.X.copy()

    # 5. 高度可变基因筛选
    sc.pp.highly_variable_genes(adata, n_top_genes=4000, subset=False)
    hvgs = adata.var[adata.var["highly_variable"]].index

    # 6. 只保留Source和Target Domain共有基因
    with open("/media/shenli/Expansion/SCAD/data/original/2128_ppi_gene_list.txt", "r") as f:
        source_genes = set(line.strip() for line in f)  # 读取Source基因列表，每行一个基因

    common_genes = list(source_genes & set(adata.var_names))  # 取交集
    adata_common = adata[:, common_genes]  # 只保留Source和Target共有的基因

    # 7. 导出TSV文件
    pd.DataFrame(adata.layers["z_score"], index=adata.obs_names, columns=adata.var_names).to_csv(
        SAVE_DATA_TO + 'Original_exprs_resp_z.' + drug + '_'+ ID +'.tsv', sep="\t")
    pd.DataFrame(adata[:, hvgs].layers["z_score"], index=adata.obs_names, columns=hvgs).to_csv(
        SAVE_DATA_TO + 'Original_exprs_resp_z.' + drug + '_'+ ID +'_tp4k.tsv', sep="\t")
    pd.DataFrame(adata_common.layers["z_score"], index=adata_common.obs_names, columns=common_genes).to_csv(
        SAVE_DATA_TO + 'Original_exprs_resp_z.' + drug + '_'+ ID +'_ppi.tsv', sep="\t")


# 设置文件路径
data_folder_path_list = ['','']
data_folder = "E:/DRMREF/data_for_SCAD/"
info_file = 'C:/Users/Li Shen/Desktop/scDrugPredict_benchmark/dataset_drug_single.csv'  # 读取 dataset_info
data_info = pd.read_csv(info_file, index_col=0)  # 读取GSE和药物信息
dataset_info = dict(zip(data_info.index.astype(str), data_info.iloc[:, 0]))

# 遍历raw count数据
for file in os.listdir(data_folder):
    if file.endswith("_count_matrix.csv"):
        file_path = os.path.join(data_folder, file)
        match = re.search(r"sampled_(.+?)_raw", file)
        geo_id = match.group(1)
        drug = dataset_info.get(geo_id, "Unknown")
        data_preprocess(file_path, geo_id, drug)



