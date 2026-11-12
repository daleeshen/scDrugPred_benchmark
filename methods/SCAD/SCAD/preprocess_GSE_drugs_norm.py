from sklearn import preprocessing
from pathlib import Path
import os.path
import math
from sklearn.feature_selection import VarianceThreshold
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
#import GEOparse
from collections import Counter
import re
import pandas as pd
import numpy as np
import copy
import openpyxl
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)


def HGNCid2ENTERZid(genelist):
    dict_sym = {}
    HGNC = pd.read_csv("/scratch/project_2011055/SCAD/data/original/HGNC_symbol_all_genes.tsv", sep='\t', index_col=0, decimal='.')
    for index in range(len(HGNC['Approved symbol'])):
        if (not isinstance(HGNC['Approved symbol'][index], float)):  
            if HGNC['NCBI Gene ID'][index] == HGNC['NCBI Gene ID'][index]:
                dict_sym[HGNC['Approved symbol'][index]] = int(HGNC['NCBI Gene ID'][index])
    ncbilist = []
    indexkeep = []
    for index in range(len(genelist)):
        if genelist[index] in dict_sym:
            indexkeep.append(index)
            ncbilist.append(dict_sym[genelist[index]])
    return (ncbilist, indexkeep, dict_sym)

##creat a Dict for GSE number and drug for matching

folder_path = Path("/scratch/project_2011055/benchmark_data_unb")
data_path_list = [d for d in folder_path.iterdir() if d.is_dir()]
data_info = '/scratch/project_2011055/dataset_drug_single.csv'  # 读取 dataset_info
df = pd.read_csv("/scratch/project_2011055/SCAD/data/response/GDSC_response.all_drugs.tsv",sep='\t', index_col=0, decimal='.')
df_ic50 = pd.read_csv("/scratch/project_2011055/SCAD/data/response/GDSC_response.logIC50.all_drugs.tsv",sep='\t', index_col=0, decimal='.')
data_info_df = pd.read_csv(data_info, header=0)
key_col_name = data_info_df.columns[0]
value_col_name = data_info_df.columns[1:3]
value_col_name = data_info_df.columns[1:3]

ratio_list = ['1_3','1_10','1_30','1_100']

for data_path_general in data_path_list:
    for ratio in ratio_list:
        data_path = data_path_general / ratio
        print('目前处理路径为: ', str(data_path))
        matched_results = {}
        for file_name in os.listdir(data_path):
            match = re.search(r'(.+?)\.RData$', file_name)
            if match:
                extracted_text = match.group(1)
                matched_row = data_info_df[data_info_df[key_col_name] == extracted_text]
                if not matched_row.empty:
                    drug = matched_row[value_col_name].values[0][0]
                    if drug in df.columns:
                        matched_results[extracted_text] = matched_row[value_col_name].values[0]
                    else:
                        if drug == 'Palbociclib':
                            drug = 'PD-0332991'
                            matched_results[extracted_text] = np.array([drug, matched_row[value_col_name].values[0][1]], dtype=object)
                        elif drug == 'PLX-4720':
                            drug = 'PLX4720'
                            matched_results[extracted_text] = np.array([drug, matched_row[value_col_name].values[0][1]], dtype=object)
                        elif drug == 'TAE684':
                            drug = 'NVP-TAE684'
                            matched_results[extracted_text] = np.array([drug, matched_row[value_col_name].values[0][1]], dtype=object)
                        elif drug == 'JQ1':
                            drug = 'JQ1'
                            matched_results[extracted_text] = np.array([drug, matched_row[value_col_name].values[0][1]], dtype=object)
        for key, value in matched_results.items():
            DRUG = value[0]  # Set to the name of the drug being split
            GEO_id = key
            Disease = value[1]
            print('开始处理数据集: ', str(GEO_id))
            non_solid_tumor = ['Myeloma', 'leukemia']
            if any(word in Disease for word in non_solid_tumor):
                Solid="False"  ###use GDSC solid tumor only
            else:
                Solid="True"
            print(Disease, Solid)

            geneset = ['', '_ppi', '_tp4k']

            SAVE_ROOT_PATH = os.path.join(data_path, 'metadata_count/data_for_SCAD')
            ROOT_PATH = os.path.join(data_path, 'metadata_count/data_for_SCAD')

            SAVE_RESULT_TO = SAVE_ROOT_PATH + '/data/split_norm/' + GEO_id + '_' + DRUG + '/'
            if DRUG == 'NVP-TAE684':
                LOAD_DATA_FROM = ROOT_PATH + '/data/original_norm/' + GEO_id + '_TAE684/'
            elif DRUG == 'PD-0332991':
                LOAD_DATA_FROM = ROOT_PATH + '/data/original_norm/' + GEO_id + '_Palbociclib/'
            else:
                LOAD_DATA_FROM = ROOT_PATH + '/data/original_norm/' + GEO_id + '_' + DRUG + '/'
            dirName = SAVE_RESULT_TO

            if not os.path.exists(dirName):
                os.makedirs(dirName)
                print("Directory ", dirName, " Created ")
            else:
                print("Directory ", dirName, " already exists")

            DRUG_name = DRUG

            d1 = df_ic50.loc[:, [DRUG_name]]
            d1.columns = ["logIC50"]
            d2 = df.loc[:, [DRUG_name]]
            d2.columns = ["response"]

            print("original cell line count with logIC50 = {}".format(len(d1.index)))
            print("original cell line count with binarized label = {}".format(len(d2.index)))

            d1.dropna(inplace=True)
            d2.dropna(inplace=True)

            ####decide whether to use solid tumor cell lines only####
            COSMIC_ids = pd.read_excel("/scratch/project_2011055/SCAD/data/response/TableS1E.xlsx")
            COSMIC_ids = COSMIC_ids.iloc[3:, [1, 2, 8]]
            COSMIC_ids = COSMIC_ids.iloc[:-1, ]
            COSMIC_ids.columns = ["name", 'COSMIC', 'Tissue_descriptor1']
            COSMIC_ids.set_index("name", inplace=True, drop=True)  ## 1001 * 2

            # 过滤出实体瘤（solid tumors）细胞系
            solid_tumors = COSMIC_ids['COSMIC'][
                (COSMIC_ids['Tissue_descriptor1'] != "leukemia") &
                (COSMIC_ids['Tissue_descriptor1'] != "lymphoma") &
                (COSMIC_ids['Tissue_descriptor1'] != "myeloma")
                ]

            if Solid == "True":
                d1 = d1[d1.index.isin(solid_tumors)]
                d2 = d2[d2.index.isin(solid_tumors)]
                print("dataset: {}, disease: {}, type: Solid".format(GEO_id, Disease))
            else:
                print("dataset: {}, disease: {}, type: non-Solid".format(GEO_id, Disease))

            print("intersect cell line count with logIC50 = {}".format(len(d1.index)))
            print("intersect cell line count with binarized label = {}".format(len(d2.index)))

            ##check ids
            if DRUG == 'JQ1':
                d_drug = pd.concat([d2, d1], axis=1, join='inner')
                d_drug.loc[:, "drug"] = DRUG
                d_drug.index.name = "sample_name"
            else:
                d_drug = pd.concat([d2, d1], axis=1)
                d_drug.loc[:, "drug"] = DRUG
                d_drug.index.name = "sample_name"

            gdsc_rma_all = pd.read_csv("/scratch/project_2011055/SCAD/data/original/gdsc-rma_gene-expression.csv", sep=',', index_col=0, decimal='.')

            etz = HGNCid2ENTERZid(gdsc_rma_all.columns)

            ncbilist = etz[0]
            indexkeep = etz[1]
            dict_ = etz[2]

            hgnckeep = gdsc_rma_all.columns[indexkeep]

            gdsc_rma_all.shape
            gdsc_rma_all = gdsc_rma_all.iloc[:, indexkeep]
            gdsc_rma_all.columns = ncbilist

            ########### Get intersect DRUG cell lines with both rma and drug DRUG IC50 value
            drug_cells = []
            for i in (gdsc_rma_all.index):
                if i in d_drug.index:
                    drug_cells.append(i)

            drug_cells.sort()  ## sort by cell line name

            GDSC_exprs_drug = gdsc_rma_all.loc[drug_cells]
            print("GDSC Sample Size: {}".format(GDSC_exprs_drug.shape))

            d_drug_intersect = d_drug.loc[drug_cells]
            Counter(GDSC_exprs_drug.index == d_drug_intersect.index) #
            print("GDSC cell line index checking: {}".format(Counter(GDSC_exprs_drug.index == d_drug_intersect.index)))

            d_drug_intersect['logIC50'].describe()

            print("GDSC: Num of resistant and Num of Sensitive = {}".format(Counter(d_drug_intersect['response'].values)))
            d_drug_intersect['response'].values[d_drug_intersect['response'].values == 'R'] = 0  ## 0 means resistant
            d_drug_intersect['response'].values[d_drug_intersect['response'].values == 'S'] = 1  ## 1 means sensitive
            print("GDSC: Num of resistant and Num of Sensitive = {}".format(Counter(d_drug_intersect['response'].values)))

            scaler1 = StandardScaler()
            GDSC_exprs_drug[GDSC_exprs_drug.columns] = scaler1.fit_transform(GDSC_exprs_drug[GDSC_exprs_drug.columns])
            GDSC_exprs_drug_z = GDSC_exprs_drug

            #######################################################################
            ##      data preprocessing for GSE108383 expr.z.resp tsv table       ##
            #######################################################################
            # Loading Normalized data:
            metadata_ROOT_PATH = os.path.join(data_path, 'metadata_count')
            GSE_metadata = pd.read_csv(metadata_ROOT_PATH + '/' + GEO_id + '_meta_data.csv')
            GSE_metadata.iloc[:, 0] = GSE_metadata.iloc[:, 0].str.replace('-', '.', regex=False)
            GSE_response_stat = GSE_metadata['condition'].tolist()
            GSE_response = []
            for i in GSE_response_stat:
                if i == 'resistant':
                    response = 0
                    GSE_response.append(response)
                else:
                    response = 1
                    GSE_response.append(response)
            print("{}: Num of resistant(0) and Num of Sensitive(1) = {}".format(GEO_id, Counter(GSE_response)))
            GDSC_exprs_drug_z = GDSC_exprs_drug

            for gs in geneset:
                if DRUG_name == 'NVP-TAE684':
                    GSE_exprs_z = pd.read_csv(
                        LOAD_DATA_FROM + 'Original_exprs_resp_z.' + 'TAE684_' + GEO_id + gs + '.tsv',
                        sep='\t', index_col=0, decimal='.')
                elif DRUG_name == 'PD-0332991':
                    GSE_exprs_z = pd.read_csv(
                        LOAD_DATA_FROM + 'Original_exprs_resp_z.' + 'Palbociclib_' + GEO_id + gs + '.tsv',
                        sep='\t', index_col=0, decimal='.')
                else:
                    GSE_exprs_z = pd.read_csv(
                        LOAD_DATA_FROM + 'Original_exprs_resp_z.' + DRUG_name + '_' + GEO_id + gs + '.tsv',
                        sep='\t', index_col=0, decimal='.')

                GSE_etz = HGNCid2ENTERZid(GSE_exprs_z.columns)
                GSE_ncbilist = GSE_etz[0]
                GSE_indexkeep = GSE_etz[1]

                GSE_hgnckeep = GSE_exprs_z.columns[GSE_indexkeep]

                GSE_exprs_z.shape
                GSE_exprs_z = GSE_exprs_z.iloc[:, GSE_indexkeep]
                GSE_exprs_z.columns = GSE_ncbilist
                ls = GSE_exprs_z.columns.intersection(GDSC_exprs_drug.columns)
                valid_genes = [gene for gene in ls if gene in GDSC_exprs_drug_z.columns]

                GDSC_exprs_drug_z = GDSC_exprs_drug_z.loc[:, valid_genes]
                GSE_exprs_z = GSE_exprs_z.loc[:, valid_genes]
                print(ls)
                print("GDSC sample size and number of genes = {}".format(GDSC_exprs_drug_z.shape))
                print("{} sample size and number of genes = {}".format(GEO_id, GSE_exprs_z.shape))

                GDSC_exprs_drug_z_genes = list(GDSC_exprs_drug_z.columns.values)
                GSE_exprs_z_genes = list(GSE_exprs_z.columns.values)

                if not (GDSC_exprs_drug_z_genes == GSE_exprs_z_genes):
                    print("\nWARNING: genes do not have the same order\n")

                GDSC_exprs_drug_z_resp = GDSC_exprs_drug_z.copy(deep=True)
                GDSC_exprs_drug_z_resp.insert(0, "response", list(d_drug_intersect['response']), allow_duplicates=False)
                GDSC_exprs_drug_z_resp.insert(1, "logIC50", list(d_drug_intersect['logIC50']), allow_duplicates=False)

                GDSC_exprs_drug_z_resp = shuffle(GDSC_exprs_drug_z_resp, random_state=42)

                gdsc_save_path = os.path.join(SAVE_RESULT_TO, 'Source_exprs_resp_z.' + DRUG + gs + '.tsv')
                GDSC_exprs_drug_z_resp.to_csv(path_or_buf=gdsc_save_path,
                                              sep='\t', index=True)


                print(" - - successfully saved GDSC {} Source data - - \n".format(DRUG))

                GSE_exprs_z_resp = GSE_exprs_z.copy(deep=True)
                GSE_exprs_z_resp.insert(0, "response", GSE_response, allow_duplicates=False)
                GSE_exprs_z_resp = shuffle(GSE_exprs_z_resp, random_state=42)
                GSE_exprs_z_resp.to_csv(
                    path_or_buf=os.path.join(SAVE_RESULT_TO, 'Target_expr_resp_z.' + DRUG + gs + '.tsv'),
                    sep='\t', index=True)

                print(" - - successfully saved preprocessed {}{} {} Target data - - \n".format(GEO_id, gs, DRUG))


