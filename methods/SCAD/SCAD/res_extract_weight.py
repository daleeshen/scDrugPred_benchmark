import os
import re
import pandas as pd
import glob
from sklearn.metrics import accuracy_score, precision_recall_curve, average_precision_score

GSE_folder = ['GSE108397_PLX4720', 'GSE115251_NVP-TAE684', 'GSE131984_JQ1_DMSO_rm_JQ1', 'GSE131984_Pac_DMSO_rm_Paclitaxel', 'GSE131984_Pal_DMSO_rm_PD-0332991', 'GSE140440_Docetaxel', 'GSE149214_Erlotinib', 'GSE158457_Dasatinib', 'GSE163836_Paclitaxel', 'GSE169246_PacBlood_Paclitaxel', 'GSE169246_PacTissue_Paclitaxel', 'GSE175716_Sorafenib', 'GSE186960_Gemcitabine']

def extract_AUC_val(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        content = f.read()

    match = re.search(r'Average Test AUC = ([\d\.]+)', content)
    if match:
        return float(match.group(1))
    else:
        return "AUC not Found"



geneset = ['_norm','_ppi','_tp4k']
df_list = []
for folder_name in GSE_folder:
    for gs in geneset:
        LOAD_PATH = 'D:/SCAD/best_param_res_weight/best_param_results' + gs + '/' + folder_name + '_' + 'train_test_summary_5folds_5seeds.txt'
        AUC = extract_AUC_val(LOAD_PATH)
        GEO_id = folder_name.rsplit('_', 1)[0]
        Drug = folder_name.rsplit('_', 1)[1]

        pred_val_path = 'D:/SCAD/best_param_res_weight/best_param_results' + gs + '/' + folder_name + '/bimodal/best_models/FX/test/'
        for tar_path in glob.glob(os.path.join(pred_val_path, "*")):
            if os.path.isdir(tar_path):
                pred_val_file = glob.glob(os.path.join(tar_path, "*_y.tsv"))
                with open(pred_val_file[0], 'r', encoding='utf-8') as f:
                    pred_df = pd.read_csv(f, sep='\t')

                true_lab = pred_df.iloc[:, 2]
                pred_prob = pred_df.iloc[:, 3]

                predicted_labels = (pred_prob >= 0.5).astype(int)
                acc = accuracy_score(true_lab, predicted_labels)

        df_list.append([GEO_id, Drug, AUC, acc, 'SCAD_weight'+gs])





df = pd.DataFrame(df_list, columns=['Dataset','Drug','AUC', 'ACC','Mode'])
#print(df)
#df.to_csv("D:/SCAD/best_param_res_weight/SCAD_weight_best_param_AUC.csv", index=False)
#print(df)
#print('------------------------------ ↓MAX↓ ------------------------------')
new_df = df.loc[df.groupby("Dataset")['AUC'].idxmax()]
#print(new_df)
#new_df.to_csv("D:/SCAD/best_param_res_weight/SCAD_weight_best_param_AUC_max.csv", index=False)


