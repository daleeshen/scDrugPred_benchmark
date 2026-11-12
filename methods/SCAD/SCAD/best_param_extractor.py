import re


def extract_best_auc_params(log_text):
    # Regular expression to extract necessary information
    pattern = re.compile(
        r"Drug = ([\w-]+), extractor = (\w+), source average validation AUC of iter(\d+) = ([\d.]+), "
        r"hdim(\d+)_zdim(\d+)_lr([\d.]+)_epoch(\d+)_lambda1([\d.]+)_"
        r"dropouts([\d.]+)_([\d.]+)_([\d.]+)_mbS(\d+)_mbT(\d+)"
    )

    best_iter = None
    best_auc = float('-inf')
    best_params = {}

    # Iterate through matches and find the best AUC
    for match in pattern.finditer(log_text):
        drug, extractor, iter_num, auc, hdim, zdim, lr, epoch, lambda1, d1, d2, d3, mbS, mbT = match.groups()
        auc = float(auc)
        if auc > best_auc:
            best_auc = auc
            best_iter = iter_num
            best_params = {
                "Drug": drug,
                "Extractor": extractor,
                "AUC": auc,
                "hdim": int(hdim),
                "zdim": int(zdim),
                "lr": float(lr),
                "epoch": int(epoch),
                "lambda1": float(lambda1),
                "dropouts": (float(d1), float(d2), float(d3)),
                "mbS": int(mbS),
                "mbT": int(mbT),
                "iter": int(iter_num)
            }

    return best_params


def extract_best_auc_params_smote(log_text, drug):
    # Regular expression to extract necessary information
    pattern = re.compile(
        rf"Drug = {drug}, extractor = (\w+), source average validation AUC of iter(\d+) = ([\d.]+), "
        r"hdim(\d+)_zdim(\d+)_lr([\d.]+)_epoch(\d+)_lambda([\d.]+)_"
        r"dropouts([\d.]+)_([\d.]+)_([\d.]+)_mbS(\d+)_mbT(\d+)"
    )

    best_iter = None
    best_auc = float('-inf')
    best_params = {}

    # Iterate through matches and find the best AUC
    for match in pattern.finditer(log_text):
        extractor, iter_num, auc, hdim, zdim, lr, epoch, lambda1, d1, d2, d3, mbS, mbT = match.groups()
        auc = float(auc)
        if auc > best_auc:
            best_auc = auc
            best_iter = iter_num
            best_params = {
                "Drug": drug,
                "Extractor": extractor,
                "AUC": auc,
                "hdim": int(hdim),
                "zdim": int(zdim),
                "lr": float(lr),
                "epoch": int(epoch),
                "lambda1": float(lambda1),
                "dropouts": (float(d1), float(d2), float(d3)),
                "mbS": int(mbS),
                "mbT": int(mbT),
                "iter": int(iter_num)
            }

    return best_params


GSE_folder = ['GSE108397_PLX4720', 'GSE115251_NVP-TAE684', 'GSE131984_JQ1_DMSO_rm_JQ1', 'GSE131984_Pac_DMSO_rm_Paclitaxel', 'GSE131984_Pal_DMSO_rm_PD-0332991', 'GSE140440_Docetaxel', 'GSE149214_Erlotinib', 'GSE158457_Dasatinib', 'GSE163836_Paclitaxel', 'GSE169246_PacBlood_Paclitaxel', 'GSE169246_PacTissue_Paclitaxel', 'GSE175716_Sorafenib', 'GSE186960_Gemcitabine']

path = 'D:/SCAD/results_smote_norm/train_eval_on_source_hypers_summary_5fds_smote_fixlr.txt'
with open(path, 'r') as f:
    file = f.read()

GEO_id = GSE_folder[0].rsplit('_', 1)[0]
Drug = GSE_folder[0].rsplit('_', 1)[1]
print(Drug)
out = extract_best_auc_params_smote(file, Drug)
print(out)


