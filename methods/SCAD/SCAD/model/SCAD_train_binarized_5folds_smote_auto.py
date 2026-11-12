import argparse
import torch
import numpy as np
import pandas as pd
import os
import random
import itertools
import re
from itertools import cycle
from sklearn.metrics import *
import sklearn.preprocessing as sk
from imblearn.over_sampling import SMOTE
from imblearn.pipeline import Pipeline
from imblearn.under_sampling import RandomUnderSampler
from SCADmodules import *

def get_args():
    parser = argparse.ArgumentParser(description="SCAD Best Parameter Evaluation with SMOTE")
    parser.add_argument('--id', type=str, required=True, help='Dataset ID (e.g., GSE140440)')
    parser.add_argument('--drug', type=str, required=True, help='Drug name (e.g., Docetaxel)')
    parser.add_argument('--solid_tumor', type=int, choices=[0, 1], required=True, help='1 for solid tumor, 0 for non-solid tumor')
    parser.add_argument('--base_path', type=str, required=True, help='Base directory containing processed data and results')
    parser.add_argument('--gpu', type=int, default=1, choices=[0, 1], help='1 for CUDA, 0 for CPU')
    return parser.parse_args()

def reset_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True

def extract_best_auc_params(log_content):
    lines = log_content.strip().split('\n')
    best_auc = -1.0
    best_params_str = ""

    for line in lines:
        if not line.strip():
            continue
        auc_match = re.search(r"Avg_Val_AUC=([\d.]+)", line)
        if auc_match:
            current_auc = float(auc_match.group(1))
            if current_auc > best_auc:
                params_match = re.search(r"Params=(.+)", line)
                if params_match:
                    best_auc = current_auc
                    best_params_str = params_match.group(1).strip()

    if not best_params_str:
        raise ValueError("No valid parameter information could be extracted from the log.")

    print(f"--- Best validation AUC: {best_auc} ---")
    print(f"--- Corresponding parameter string: {best_params_str} ---")

    param_dict = {}
    keys = ['hdim', 'zdim', 'lr', 'epoch', 'lambda1', 'mbS', 'mbT']

    for key in keys:
        match = re.search(rf"{key}([\d.]+)", best_params_str)
        if not match:
            raise ValueError(f"Parameter {key} was not found in: {best_params_str}")
        val_str = match.group(1)
        param_dict[key] = float(val_str) if '.' in val_str else int(val_str)

    return param_dict

@torch.no_grad()
def predict_label(XTestCells, gen_model, map_model, device, extractor):
    gen_model.eval()
    gen_model.to(device)
    map_model.eval()
    map_model.to(device)

    if extractor in ['FX', 'FX_SFMX', 'FX_MLP']:
        F_xt_test = gen_model(XTestCells.to(device))
    else:
        raise NotImplementedError(f'{extractor} is not a valid extractor')

    return map_model(F_xt_test)

def evaluate_model(XTestCells, YTestCells, gen_model, map_model, device, extractor):
    XTestCells = XTestCells.to(device)
    YTestCells = YTestCells.to(device)

    y_predicted = predict_label(
        XTestCells,
        gen_model,
        map_model,
        device,
        extractor
    )

    C_loss_eval = torch.nn.BCELoss()
    closs_test = C_loss_eval(y_predicted, YTestCells)

    yt_true_test = YTestCells.view(-1, 1).cpu()
    y_predicted_cpu = y_predicted.cpu()

    AUC_test = roc_auc_score(
        yt_true_test.detach().numpy(),
        y_predicted_cpu.detach().numpy()
    )
    APR_test = average_precision_score(
        yt_true_test.detach().numpy(),
        y_predicted_cpu.detach().numpy()
    )

    return closs_test, AUC_test, APR_test, y_predicted_cpu

def main():
    args = get_args()

    device = "cuda" if args.gpu == 1 and torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")

    folder_name = f"{args.id}_{args.drug}"
    Drug = args.drug
    extractor = "FX"
    genesets = ['_norm', '_ppi', '_tp4k']
    init_seed = 42

    for gs in genesets:
        # GitHub runs each geneset in a separate process.
        reset_seed(init_seed)

        print("================================================================================")
        print(f"Final SCAD SMOTE model: {folder_name} | geneset={gs}")
        print("================================================================================")

        log_file_path = os.path.join(
            args.base_path,
            f"results_smote{gs}",
            'train_eval_on_source_hypers_summary_5fds_smote_fixlr.txt'
        )

        if not os.path.exists(log_file_path):
            raise FileNotFoundError(
                f"Hyperparameter search log not found: {log_file_path}"
            )

        with open(log_file_path, 'r') as f:
            log_content = f.read()

        best_params = extract_best_auc_params(log_content)

        h_dim = best_params['hdim']
        z_dim = best_params['zdim']
        epoch = best_params['epoch']
        lam1 = best_params['lambda1']
        mbS = best_params['mbS']
        mbT = best_params['mbT']
        lr = best_params['lr']

        # Step 4 uses only 8/8 and 32/32, matching the GitHub final script.
        if mbS != mbT:
            raise ValueError(
                f"GitHub SCAD SMOTE final model assumes mbS == mbT, "
                f"but received mbS={mbS}, mbT={mbT}"
            )

        dropout_gen = 0.5
        dropout_mtl = dropout_gen
        dropout_dg = dropout_gen

        print(
            f"h_dim={h_dim}, z_dim={z_dim}, lr={lr}, "
            f"epoch={epoch}, lambda1={lam1}, mbS={mbS}, mbT={mbT}"
        )

        results_save_base = os.path.join(
            args.base_path,
            f"best_param_results_smote{gs}",
            folder_name
        )

        mode_path = os.path.join(
            results_save_base,
            "bimodal",
            "best_models",
            extractor
        )

        os.makedirs(
            os.path.join(mode_path, "test"),
            exist_ok=True
        )

        tr_ts_summary_path = os.path.join(
            args.base_path,
            f"best_param_results_smote{gs}",
            f"{args.id}_{Drug}_train_test_summary_5folds_5seeds.txt"
        )

        load_data_from = os.path.join(
            args.base_path,
            'data',
            'split_norm',
            folder_name,
            f'stratified{gs}'
        )

        if gs == '_norm':
            source_dir = 'source_5_folds_norm'
            target_dir = 'target_5_folds_norm'
        else:
            source_dir = f'source_5_folds_{gs}'
            target_dir = f'target_5_folds_{gs}'

        ls_splits = ['split1', 'split2', 'split3', 'split4', 'split5']

        param_str = (
            f"hdim{h_dim}_zdim{z_dim}_lr{lr}_epoch{epoch}_"
            f"lambda1{lam1}_mbS{mbS}_mbT{mbT}"
        )

        test_res_dir = os.path.join(
            mode_path,
            "test",
            f"mbS{mbS}_mbT{mbT}"
        )

        os.makedirs(
            test_res_dir,
            exist_ok=True
        )

        test_results_file = os.path.join(
            test_res_dir,
            f"{param_str}_smote_seed{init_seed}.tsv"
        )

        test_pred_y_file = os.path.join(
            test_res_dir,
            f"{param_str}_seed{init_seed}_predicted_y.tsv"
        )

        if os.path.isfile(test_results_file):
            os.remove(test_results_file)

        if os.path.isfile(test_pred_y_file):
            os.remove(test_pred_y_file)

        with open(test_results_file, 'a') as f:
            f.write("-- Parameters --\n\n")
            f.write(
                f"h_dim: {h_dim}\n"
                f"z_dim: {z_dim}\n"
                f"lr: {lr}\n"
                f"epoch: {epoch}\n"
                f"lambda1: {lam1}\n"
                f"mbS: {mbS}\n"
                f"mbT: {mbT}\n"
            )

        # Same fixed SMOTE pipeline as the GitHub implementation.
        over = SMOTE(sampling_strategy=0.5)
        under = RandomUnderSampler(sampling_strategy=0.5)
        pipeline = Pipeline(
            steps=[
                ('o', over),
                ('u', under)
            ]
        )

        auc_test_total = []
        apr_test_total = []

        for split in ls_splits:
            print(f"Processing {folder_name} | {gs} | {split}...")

            XTrainGDSC = pd.read_csv(
                os.path.join(
                    load_data_from,
                    source_dir,
                    split,
                    "X_train_source.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            YTrainGDSC = pd.read_csv(
                os.path.join(
                    load_data_from,
                    source_dir,
                    split,
                    "Y_train_source.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            # Loaded by the original GitHub implementation.
            XValGDSC = pd.read_csv(
                os.path.join(
                    load_data_from,
                    source_dir,
                    split,
                    "X_val_source.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            YValGDSC = pd.read_csv(
                os.path.join(
                    load_data_from,
                    source_dir,
                    split,
                    "Y_val_source.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            # Fixed SMOTE(0.5) + undersampling(0.5), no fallback.
            XTrainGDSC, YTrainGDSC = pipeline.fit_resample(
                XTrainGDSC,
                YTrainGDSC
            )

            XTrainCells = pd.read_csv(
                os.path.join(
                    load_data_from,
                    target_dir,
                    split,
                    "X_train_target.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            YTrainCells = pd.read_csv(
                os.path.join(
                    load_data_from,
                    target_dir,
                    split,
                    "Y_train_target.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            XTestCells = pd.read_csv(
                os.path.join(
                    load_data_from,
                    target_dir,
                    split,
                    "X_test_target.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            YTestCells = pd.read_csv(
                os.path.join(
                    load_data_from,
                    target_dir,
                    split,
                    "Y_test_target.tsv"
                ),
                sep='\t',
                index_col=0,
                decimal='.'
            )

            # Fit scaler on SMOTE-processed source train + target train.
            XTrainCombined = pd.concat([
                XTrainGDSC,
                XTrainCells
            ])

            scalerTrain = sk.StandardScaler()
            scalerTrain.fit(
                XTrainCombined.values
            )

            XTrainGDSC_N = scalerTrain.transform(
                XTrainGDSC.values
            )

            XTrainCells_N = scalerTrain.transform(
                XTrainCells.values
            )

            XTestCells_N = scalerTrain.transform(
                XTestCells.values
            )

            TXTestCells_N = torch.FloatTensor(
                XTestCells_N
            ).to(device)

            TYTestCells = torch.FloatTensor(
                YTestCells.values.astype(int)
            ).to(device)

            CDataset = torch.utils.data.TensorDataset(
                torch.FloatTensor(XTrainGDSC_N),
                torch.FloatTensor(YTrainGDSC.values)
            )

            PDataset = torch.utils.data.TensorDataset(
                torch.FloatTensor(XTrainCells_N),
                torch.FloatTensor(
                    YTrainCells.values.astype(int)
                )
            )

            # GitHub SMOTE version uses shuffle=True for both domains.
            CLoader = torch.utils.data.DataLoader(
                dataset=CDataset,
                batch_size=mbS,
                shuffle=True
            )

            PLoader = torch.utils.data.DataLoader(
                dataset=PDataset,
                batch_size=mbT,
                shuffle=True
            )

            _, IE_dim = XTrainGDSC_N.shape

            if extractor == 'FX':
                Gen = FX(
                    dropout_gen,
                    IE_dim,
                    h_dim,
                    z_dim
                )
            else:
                raise NotImplementedError(
                    f'{extractor} is not a valid extractor'
                )

            Map = MTLP(
                dropout_mtl,
                h_dim,
                z_dim
            )

            DG = Discriminator(
                dropout_dg,
                h_dim,
                z_dim
            )

            Gen.to(device)
            Map.to(device)
            DG.to(device)

            optimizer = torch.optim.Adagrad(
                itertools.chain(
                    Gen.parameters(),
                    Map.parameters(),
                    DG.parameters()
                ),
                lr=lr
            )

            C_loss = torch.nn.BCELoss()

            for it in range(epoch):
                for data in zip(
                    cycle(CLoader),
                    PLoader
                ):
                    DataS = data[0]
                    DataT = data[1]

                    xs = DataS[0].to(device)
                    ys = DataS[1].view(-1, 1).to(device)

                    xt = DataT[0].to(device)
                    yt = DataT[1].view(-1, 1).to(device)

                    if xs.size()[0] < 5 or xt.size()[0] < 5:
                        continue

                    Gen.train()
                    Map.train()
                    DG.train()

                    F_xs = Gen(xs)
                    F_xt = Gen(xt)

                    yhat_xs = Map(F_xs)

                    closs = C_loss(
                        yhat_xs,
                        ys
                    )

                    Labels = torch.ones(
                        F_xs.size(0),
                        1
                    )

                    Labelt = torch.zeros(
                        F_xt.size(0),
                        1
                    )

                    Lst = torch.cat(
                        [Labels, Labelt],
                        0
                    ).to(device)

                    Xst = torch.cat(
                        [F_xs, F_xt],
                        0
                    ).to(device)

                    yhat_DG = DG(Xst)

                    DG_loss = C_loss(
                        yhat_DG,
                        Lst
                    )

                    loss1 = closs
                    loss2 = lam1 * DG_loss
                    Loss = loss1 + loss2

                    optimizer.zero_grad()
                    Loss.backward()
                    optimizer.step()

            test_loss, test_auc, test_apr, test_predicted_y = evaluate_model(
                TXTestCells_N,
                TYTestCells,
                Gen,
                Map,
                device,
                extractor
            )

            print(
                f"{split}: Test AUC={test_auc:.6f}, "
                f"Test APR={test_apr:.6f}"
            )

            auc_test_total.append(
                test_auc
            )

            apr_test_total.append(
                test_apr
            )

            with open(test_results_file, 'a') as f:
                f.write(
                    f"-- Split {split} --\n"
                )
                f.write(
                    f"Test loss: {test_loss}\t"
                    f"Test AUC: {test_auc}\t"
                    f"Test APR: {test_apr}\n\n"
                )

            # Additional benchmark output:
            # split, cell_id, true_response, predicted_score
            pred_np = test_predicted_y.detach().cpu().numpy()

            with open(test_pred_y_file, 'a') as pf:
                for i in range(len(YTestCells)):
                    pf.write(
                        f"{split}\t"
                        f"{YTestCells.index[i]}\t"
                        f"{YTestCells['response'].values[i]}\t"
                        f"{pred_np[i][0]}\n"
                    )

        auc_test_total = np.array(
            auc_test_total
        )
        apr_test_total = np.array(
            apr_test_total
        )

        avg_auc = np.mean(
            auc_test_total
        )
        std_auc = np.std(
            auc_test_total
        )
        avg_apr = np.mean(
            apr_test_total
        )
        std_apr = np.std(
            apr_test_total
        )

        with open(test_results_file, 'a') as f:
            f.write("\n-- Average Test AUC --\n")
            f.write(
                f"Mean: {avg_auc}\t"
                f"Standard Deviation: {std_auc}\n"
            )

            f.write("\n-- Average Test Precision Recall --\n")
            f.write(
                f"Mean: {avg_apr}\t"
                f"Standard Deviation: {std_apr}\n"
            )

        with open(tr_ts_summary_path, 'a') as f:
            f.write(
                f"Average Test AUC = {avg_auc}; "
                f"Average Test Precision Recall = {avg_apr}; "
                f"path = {test_results_file}\n"
            )

        print(
            f"{gs}: Avg Test AUC={avg_auc:.6f}, "
            f"Avg Test APR={avg_apr:.6f}"
        )

if __name__ == "__main__":
    main()