import csv
import os
import pickle
import pandas as pd
from datasets import Dataset
csv.field_size_limit(10000000)

# 读取基因token映射字典
with open('./token_gene_dictionary.pickle', 'rb') as file:
    dictionary_genename_token_pair = pickle.load(file)

# 遍历数据集
root_path = '/media/shenli/Expansion/scDrug_benchmark/DrugFormer'
k_ratio_list = ['1_3','1_10','1_30','1_100']
dataset_info = pd.read_csv('./dataset_drug_single.csv')

for round_idx in range(1, 4):
    for k_ratio in k_ratio_list:
        data_path = f'{root_path}/{round_idx}_round_{k_ratio}_countMatrix/'
        for ID in dataset_info.iloc[:, 0]:
            project_name = ID
            print(f'Start processing {project_name} | round = {round_idx}, k_ratio = {k_ratio}')

            csv_file_path = os.path.join(data_path, f'{project_name}.countMatrix.txt')
            labels_cell = []
            samples = []

            try:
                with open(csv_file_path, 'r', encoding='utf-8') as file:
                    csv_reader = csv.reader(file)
                    pattern = [-99999]
                    headr = True

                    for row in csv_reader:
                        row_data = row[0].split('\t')
                        if headr:
                            for gene in row_data[:-1]:
                                gene = gene.replace('"', '')
                                token = dictionary_genename_token_pair.get(gene, -99999)
                                pattern.append(token)
                            pattern.append(-99999)
                            headr = False
                        else:
                            if len(pattern) != len(row_data):

                                print(f"[Warning] Pattern length mismatch in {project_name}, skipping row.")
                                continue

                            seq_pattern_order_id_EXPscore = []
                            for j in range(len(row_data)):
                                if j == 0:
                                    continue
                                elif j == len(row_data) - 1:
                                    cond = row_data[j]
                                    if 'sensitive' in row_data[j]:
                                        labels_cell.append(1)
                                    elif 'resistant' in row_data[j]:
                                        labels_cell.append(0)
                                else:
                                    if row_data[j] != '0' and pattern[j] != -99999:
                                        seq_pattern_order_id_EXPscore.append((pattern[j], row_data[j]))

                            sorted_seq = sorted(seq_pattern_order_id_EXPscore, key=lambda x: x[1], reverse=True)
                            sample = [item[0] for item in sorted_seq][:2048]
                            sample += [0] * (2048 - len(sample))  # pad
                            samples.append(sample)

                if len(samples) == 0:
                    print(f"[Warning] No valid samples found for {project_name}, skipping save.")
                    continue

                if len(samples) != len(labels_cell):
                    print(f"[Error] Sample-label mismatch in {project_name}, skipping save.")
                    continue

                dt = {'input_ids': samples, 'cell_label': labels_cell}
                save_path = os.path.join(root_path, f'{round_idx}_{k_ratio}_cell_dt', project_name)
                os.makedirs(save_path, exist_ok=True)
                print(f'Saving dataset for {project_name} to {save_path}...')

                my_dataset = Dataset.from_dict(dt)
                my_dataset.save_to_disk(save_path)
                print(f'Dataset for {project_name} saved successfully.')

            except FileNotFoundError:
                print(f'[Error] File not found: {csv_file_path}')





